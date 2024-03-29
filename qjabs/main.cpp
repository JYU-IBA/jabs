#include <QApplication>
#include <QLocale>
#include <QFileOpenEvent>
#include <QtDebug>
#include <cstdarg>
#include <clocale>
#include "mainwindow.h"
extern "C" {
#include <stdio.h>
#include "message.h"
#include "win_compat.h"
int fit_iter_callback(fit_stats stats);
}

jabs_msg_level jabs_message_verbosity;


class QJaBSApplication : public QApplication
{
public:
    QJaBSApplication(int &argc, char **argv)
        : QApplication(argc, argv)
    {
    }

    bool event(QEvent *event) override
    {
        if (event->type() == QEvent::FileOpen) {
           QFileOpenEvent *openEvent = static_cast<QFileOpenEvent *>(event);
           mainWindow->openFile(openEvent->file());
        }
        return QApplication::event(event);
    }
    void setMainWindow(MainWindow *mw) {mainWindow = mw;}
private:
    MainWindow *mainWindow;
};


MainWindow *mainWindow; /* Dirty trick! jabs_message needs to know the address of MainWindow, but we don't want to mess with the function signature (to pass a pointer directly), since jabs_message() has an alternate implementation in the CLI version. */

int main(int argc, char *argv[])
{
    jabs_message_verbosity = JABS_DEFAULT_VERBOSITY;
    setlocale(LC_NUMERIC,"C");
    QJaBSApplication a(argc, argv);
    setlocale(LC_NUMERIC,"C");  /* This should force C part of JaBS to use "." as a decimal separator even when locale says it is ",". */
    QCoreApplication::setOrganizationName("University of Jyväskylä");
    QCoreApplication::setOrganizationDomain("jyu.fi");
    QCoreApplication::setApplicationName("JaBS");
    QCoreApplication::setApplicationVersion(jabs_version());
    MainWindow w;
    mainWindow = &w;
    a.setMainWindow(&w);
    w.show();
    if(argc == 2) {
        w.openFile(argv[1]);
    }
    return a.exec();
}

static const char *jabs_msg_levels[MSG_ERROR+1] = {"Debug", "Verbose", "Default", "Important", "Warning", "Error"};

const char *jabs_message_level_str(jabs_msg_level level) {
    if(level <= MSG_ERROR) {
        return jabs_msg_levels[level];
    }
    return NULL;
}

void jabs_message(jabs_msg_level level, const char * format, ...) {
    if(level < jabs_message_verbosity) {
        return;
    }
    va_list argp;
    va_start(argp, format);
    char *str_out;
    vasprintf(&str_out, format, argp);
    mainWindow->addMessage(level, str_out);
    fputs(str_out, stderr);
    free(str_out);
    va_end(argp);
}

void jabs_message_printf(jabs_msg_level level, FILE *f, const char * format, ...) {
    if(level < jabs_message_verbosity) {
        return;
    }
    va_list argp;
    va_start(argp, format);
    if(f == stderr)     {
        char *str_out;
        vasprintf(&str_out, format, argp);
        mainWindow->addMessage(level, str_out);
        fputs(str_out, stderr);
        free(str_out);
    } else {
        vfprintf(f, format, argp);
    }
    va_end(argp);
}

int fit_iter_callback(fit_stats stats) {
    return mainWindow->fitCallback(stats);
}
