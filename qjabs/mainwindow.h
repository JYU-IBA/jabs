#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSettings>
#include "highlighter.h"
#include "fitdialog.h"
#include "plotdialog.h"

extern "C" {
#include <jibal_units.h>
#include "jabs.h"
#include "sample.h"
#include "script_generic.h"
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
    void addMessage(jabs_msg_level level, const char *msg);
    int fitCallback(fit_stats stats);
    ~MainWindow();

public slots:
    void openFile(const QString &filename);

private slots:
    void on_actionPreferences_triggered();

private slots:
    void plotDialogClosed();

    void readSettings(); /* does not include plot settings, see readPlotSettings() */

    void readPlotSettings();

    void on_action_Plot_triggered();

    void on_action_Run_triggered();

    void on_plotSpinBox_valueChanged(int arg1);

    void plotSpectrum(size_t i_det);

    void on_action_Open_File_triggered();

    void on_action_Quit_triggered();

    void on_action_New_File_triggered();

    void on_editor_textChanged();

    void on_action_Save_File_as_triggered();

    bool saveScriptToFile(const QString &filename);

    void on_action_Save_File_triggered();

    void resetAll();

    void on_actionAbout_triggered();

    void on_commandLineEdit_returnPressed();

    void on_msgTextBrowser_textChanged();

    void updateRecentFileActions();

    void openRecentFile();

    void openLink(const QUrl &link);

    void postInit();

private:
    void updateWindowTitle();
    void setFilename(const QString &filename);
    int runLine(const QString &line);
    void plotSession(bool error = false);
    bool askToSave();
    void closeEvent(QCloseEvent *event);
    void enableRun(bool enabled);
    void closeFitDialog();
    int initSession();
    int closeSession();
    static QString makeFileLink(const QString &filename);
    void setNeedsSaving(bool value);
    void updateListOfVisibleGraphs();
    Ui::MainWindow *ui;
    FitDialog *fitDialog;
    PlotDialog *plotDialog;
    struct jibal *jibal;
    script_session *session;
    Highlighter *highlighter;
    QString filename;
    QString filebasename;
    QString originalPath;
    bool needsSaving;
    bool firstRun;
    QString aboutString;
    int maxRecentFiles;
    QAction *recentFileActs; /* array */
    bool showIsotopes;
    int plotIsotopesZ;
    QSettings settings;
    QStringList visibleGraphs;
};
#endif // MAINWINDOW_H
