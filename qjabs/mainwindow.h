#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>


extern "C" {
#include <jibal_units.h>
#include "jabs.h"
#include "sample.h"
#include "script.h"
#include "generic.h"
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
    void on_actionRun_triggered();

    void on_plotSpinBox_valueChanged(int arg1);

    void plotSpectrum(size_t i_det);

private:
    Ui::MainWindow *ui;
    struct jibal *jibal;
    script_session *session;

};
#endif // MAINWINDOW_H
