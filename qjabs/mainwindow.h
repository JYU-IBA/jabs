#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "highlighter.h"

extern "C" {
#include <jibal_units.h>
#include "jabs.h"
#include "sample.h"
#include "script.h"
#include "generic.h"
#include "message.h"
#include "defaults.h"
#include <jibal_defaults.h>
#ifdef WIN32
#include "win_compat.h"
#endif
}


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    void addMessage(const char *msg);
    ~MainWindow();

private slots:
    void on_action_Run_triggered();

    void on_plotSpinBox_valueChanged(int arg1);

    void plotSpectrum(size_t i_det);

    void on_action_Open_File_triggered();

    void on_action_Quit_triggered();

    void on_action_New_File_triggered();

    void on_plainTextEdit_textChanged();

    void on_action_Save_File_as_triggered();

    bool saveScriptToFile(const QString &filename);

    void on_action_Save_File_triggered();

    void clearPlotAndOutput();

    void on_actionAbout_triggered();

private:
    void updateWindowTitle();
    void setFilename(const QString &filename);
    Ui::MainWindow *ui;
    struct jibal *jibal;
    script_session *session;
    Highlighter *highlighter;
    QString filename;
    QString filebasename;
    QString originalPath;
    bool needsSaving;
};
#endif // MAINWINDOW_H
