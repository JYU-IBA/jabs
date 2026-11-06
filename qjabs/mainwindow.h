#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSettings>
#include <QElapsedTimer>
#include <QProgressBar>
#include <QFrame>
#include <QLabel>
#include <QCheckBox>
#include <QPushButton>
#include <QThread>
#include <QSharedPointer>
#include <atomic>
#include <vector>
#include "highlighter.h"
#include "plotdialog.h"
#include "jabsworker.h"
#include "colorprofile.h"

extern "C" {
#include <jibal_units.h>
#include "jabs.h"
#include "fit.h"
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

struct PlotSnapshot { /* Copy of spectra (of all detectors), made by the worker thread during a fit, for plotting in the GUI thread */
    std::vector<result_spectra> spectra;
    jabs_histogram *ref = nullptr;
    PlotSnapshot() = default;
    PlotSnapshot(const PlotSnapshot &) = delete;
    PlotSnapshot &operator=(const PlotSnapshot &) = delete;
    ~PlotSnapshot() {
        for(result_spectra &s : spectra) {
            result_spectra_free(&s);
        }
        jabs_histogram_free(ref);
    }
};

Q_DECLARE_METATYPE(fit_stats)
Q_DECLARE_METATYPE(jabs_msg_level)
Q_DECLARE_METATYPE(QSharedPointer<PlotSnapshot>)

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    void addMessage(jabs_msg_level level, const char *msg);
    int fitCallback(fit_stats stats); /* Called by the worker thread during a fit, must not touch widgets */
    void reportProgress(double fraction, const char *msg); /* Can be called from any thread */
    void reportProgressDone(); /* Can be called from any thread */
    ~MainWindow();

signals:
    void runFinished();
    void sendMessage(jabs_msg_level level, QString msg);
    void startRun(const QStringList &lines); /* Handled by the worker thread */
    void fitProgress(fit_stats stats, QSharedPointer<PlotSnapshot> snapshot); /* Emitted by the worker thread */
    void progressChanged(double fraction, QString msg, bool done); /* Emitted by any thread */

public slots:
    void openFile(const QString &filename);
    void runRoi(const QString &roi);
    void updateListOfVisibleGraphs();
    void scrollMsgBoxToBottom();

private slots:
    void receiveMessage(jabs_msg_level level, QString msg);

    void onRunFinished(bool error, qsizetype lineno);

    void onFitProgress(fit_stats stats, QSharedPointer<PlotSnapshot> snapshot);

    void stopRun();

    void onProgressChanged(double fraction, QString msg, bool done);

    void applyColors(); /* Applies the active color profile (settings) to widgets, plot and icons */

    void scheduleApplyColors();

    void on_actionPreferences_triggered();

    void plotDialogClosed();

    void readSettings(); /* does not include plot settings, see readPlotSettings() */

    void readPlotSettings();

    void on_action_Plot_triggered();

    void on_action_Run_triggered();

    void plotSpectrum(size_t i_det);

    void on_action_Open_File_triggered();

    void on_action_Quit_triggered();

    void on_action_New_File_triggered();

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

    void on_actionIDF_triggered();

    void onEnergyAxisSet(bool value); /* SpectrumPlot x axis change from channels to energy or vice versa. value true for energy. */

    void onSpectrumLegendMoved(bool outside);

    void updateDetectorList();

    void on_comboBox_currentIndexChanged(int index);

    void on_actionNext_detector_triggered();

    void on_actionPrevious_detector_triggered();

private:
    enum RunMode {RUN_SCRIPT, RUN_COMMAND_LINE};
    void updateWindowTitle();
    void setFilename(const QString &filename);
    int runLine(const QString &line);
    void startWorker(const QStringList &lines, RunMode mode);
    void stopWorkerThread();
    void plotSession(bool error = false);
    void plotSpectra(const result_spectra *spectra, const jabs_histogram *ref_histo);
    bool askToSave();
    void closeEvent(QCloseEvent *event);
    void changeEvent(QEvent *event) override;
    void enableRun(bool enabled);
    void updateRunAction(); /* Run action is a stop action while running */
    void setupRunPanel();
    void showRunPanel(bool show);
    int initSession();
    int closeSession();
    static QString makeFileLink(const QString &filename, const QString &linkname = "");
    void setNeedsSaving(bool value);
    void showInitialMessages();

    Ui::MainWindow *ui;
    PlotDialog *plotDialog;
    struct jibal *jibal;
    script_session *session;
    Highlighter *highlighter;
    QString filename;
    QString filebasename;
    QString originalPath;
    QString filedirpath;
    bool needsSaving;
    bool firstRun;
    QString aboutString;
    int maxRecentFiles;
    QAction *recentFileActs; /* array */
    bool showIsotopes;
    int plotIsotopesZ;
    bool confidenceLimits;
    QSettings settings;
    QStringList visibleGraphs;
    int warningCounter;
    int defaultVerbosity;
    ColorProfile colorProfile; /* Active profile */
    bool applyingColors;
    bool applyColorsPending;
    Qt::ColorScheme forcedColorScheme; /* Color scheme set by us (Unknown = follow OS) */
    QThread *workerThread;
    JabsWorker *worker;
    bool running; /* Worker is running. GUI must not access session while this is true. */
    RunMode runMode;
    QElapsedTimer runTimer;
    QAction *actionStop; /* In run panel (and shortcut) */
    QIcon runIcon; /* Icons of run action (themed) */
    QIcon stopIcon;
    QFrame *runPanel; /* Shown during a run: progress, fit statistics, stop */
    QLabel *runStatusLabel;
    QProgressBar *progressBar;
    QWidget *fitRow; /* Part of run panel shown during fits */
    QLabel *fitPhaseLabel;
    QLabel *fitIterLabel;
    QLabel *fitChisqLabel;
    QCheckBox *plotWhileFittingCheckBox;
    QPushButton *endFitButton;
    std::atomic<bool> fitAbort; /* Graceful abort of fit (after current iteration) */
    std::atomic<bool> plotWhileFitting;
    QSharedPointer<PlotSnapshot> snapshot; /* Latest spectra received during a run */
    bool snapshotPlotted; /* Snapshot has been plotted during this run */
};
#endif // MAINWINDOW_H
