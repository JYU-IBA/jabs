#include <QApplication>
#include <QLocale>
#include <cstdarg>
#include <clocale>
#include "mainwindow.h"
extern "C" {
#include <stdio.h>
#include "message.h"
#include "win_compat.h"
}




MainWindow *mainWindow; /* Dirty trick! jabs_message needs to know the address of MainWindow, but we don't want to mess with the function signature (to pass a pointer directly), since jabs_message() has an alternate implementation in the CLI version. */

int main(int argc, char *argv[])
{
    setlocale(LC_NUMERIC,"C");
    QApplication a(argc, argv);
    setlocale(LC_NUMERIC,"C");  /* This should force C part of JaBS to use "." as a decimal separator even when locale says it is ",". */
    MainWindow w;
    mainWindow = &w;
    w.show();
    return a.exec();
}

void jabs_message(jabs_msg_level level, FILE *f, const char * format, ...) {
    va_list argp;
    va_start(argp, format);
    if(f == stderr)     {
        char *str_out;
        vasprintf(&str_out, format, argp);
        mainWindow->addMessage(str_out);
        fputs(str_out, stderr);
        free(str_out);
    } else {
        vfprintf(f, format, argp);
    }
    va_end(argp);
}
