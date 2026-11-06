#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFontDatabase>
#include <QStyleHints>
#include <QToolButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QDebug>
#include <QTextStream>
#include "preferencesdialog.h"

#include "script.h"
#include "script_session.h"
#include "script_command.h"
#include "idf2jbs.h"
extern "C" {
int fit_iter_callback(fit_stats stats);
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow),
      maxRecentFiles(5), jibal(NULL), session(NULL), highlighter(NULL), plotDialog(NULL),
      showIsotopes(false), plotIsotopesZ(JIBAL_ANY_Z), workerThread(nullptr), worker(nullptr), running(false), runMode(RUN_SCRIPT),
      actionStop(nullptr), runPanel(nullptr), runStatusLabel(nullptr), progressBar(nullptr), fitRow(nullptr), fitPhaseLabel(nullptr), fitIterLabel(nullptr),
      fitChisqLabel(nullptr), plotWhileFittingCheckBox(nullptr), endFitButton(nullptr), fitAbort(false), plotWhileFitting(false), applyingColors(false), applyColorsPending(false), snapshotPlotted(false), forcedColorScheme(Qt::ColorScheme::Unknown)
{
    qRegisterMetaType<fit_stats>("fit_stats");
    qRegisterMetaType<jabs_msg_level>("jabs_msg_level"); /* Messages can come from other threads (worker, OpenMP), queued connection needs this */
    qRegisterMetaType<QSharedPointer<PlotSnapshot>>("QSharedPointer<PlotSnapshot>");
    aboutString = QString("JaBS version ") + jabs_version() + "\n\n"
                       + "Using JIBAL version "+ jibal_version() + ", compiled using version " + JIBAL_VERSION + "\n\n"
                       + "Using Qt version " + qVersion() + ", compiled using version " + QT_VERSION_STR + "\n\n"
                       + "Copyright 2021 - 2026 Jaakko Julin <jaakko.julin@jyu.fi>\n";
    ui->setupUi(this);
    QIcon icon(":/icons/jabs.svg");
    QApplication::setWindowIcon(icon);
    highlighter = new Highlighter(ui->editor->document());
    ui->widget->setVisible(false);
    ui->detectorFrame->setVisible(false); /* Will be made visible if necessary */
    setWindowIcon(icon);
    ui->msgTextBrowser->setOpenLinks(false); /* We use the connection below to handle links */
    connect(this, &MainWindow::sendMessage, this, &MainWindow::receiveMessage);
    connect(ui->msgTextBrowser, &QTextBrowser::anchorClicked, this, &MainWindow::openLink);
    //connect(this, &MainWindow::runFinished, this, &MainWindow::scrollMsgBoxToBottom);
    originalPath = QDir::currentPath();
    ui->splitter->setSizes(QList<int>() << 1 << 3 << 1);
    ui->splitter_2->setSizes(QList<int>() << 1 << 3);
    showInitialMessages();
    ui->action_Run->setShortcutContext(Qt::ApplicationShortcut);
    ui->action_New_File->setShortcut(QKeySequence::New);
    ui->action_Open_File->setShortcut(QKeySequence::Open);
    ui->action_Save_File->setShortcut(QKeySequence::Save);
    ui->action_Save_File_as->setShortcut(QKeySequence::SaveAs);
    ui->action_Run->setShortcut(QKeySequence::Refresh);
    ui->action_Run->setShortcutVisibleInContextMenu(true);
    ui->action_Quit->setShortcut(QKeySequence::Quit);
    ui->actionPrevious_detector->setShortcut(QKeySequence::Back);
    ui->actionNext_detector->setShortcut(QKeySequence::Forward);
    setNeedsSaving(false);
    firstRun = true;
    statusBar()->showMessage(QString("JaBS ") + jabs_version() + ", cwd: " +  QDir::currentPath(), 2000);
    recentFileActs = new QAction[maxRecentFiles];
    updateRecentFileActions();
    for(int i = 0; i < maxRecentFiles; i++) {
        QMenu *m = ui->menuRecent_Files;
        QAction *a = &recentFileActs[i];
        if(m)
            m->addAction(a);
        connect(a, &QAction::triggered, this, &MainWindow::openRecentFile);
    }
    ui->menuRecent_Files->setToolTipsVisible(true);
    actionStop = new QAction(QIcon(":/icons/stop-box.svg"), QString("Stop"), this);
    actionStop->setToolTip(QString("Stop running simulation or fit"));
    actionStop->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Period));
    actionStop->setEnabled(false);
    connect(actionStop, &QAction::triggered, this, &MainWindow::stopRun);
    addAction(actionStop); /* Shortcut works in the whole window. Toolbar and menu have the run action, which stops a run. */
    connect(this, &MainWindow::fitProgress, this, &MainWindow::onFitProgress); /* Emitted from worker thread, queued */
    connect(this, &MainWindow::progressChanged, this, &MainWindow::onProgressChanged); /* Emitted from worker thread, queued */
    setupRunPanel();
    ColorProfile::migrateOldSettings(settings);
    readSettings(); /* Also applies colors */
    readPlotSettings();
    connect(QGuiApplication::styleHints(), &QStyleHints::colorSchemeChanged, this, &MainWindow::scheduleApplyColors); /* OS light / dark mode */
    int status = initSession();
    if(status != 0) {
        enableRun(false);
        ui->action_Plot->setEnabled(false);
    } else {
        workerThread = new QThread(this);
        workerThread->setStackSize(8 * 1024 * 1024); /* Default for secondary threads is small on some platforms (512 kB on macOS) */
        worker = new JabsWorker(session);
        worker->moveToThread(workerThread);
        connect(workerThread, &QThread::finished, worker, &QObject::deleteLater);
        connect(this, &MainWindow::startRun, worker, &JabsWorker::run);
        connect(worker, &JabsWorker::finished, this, &MainWindow::onRunFinished);
        connect(worker, &JabsWorker::statusMessage, this, [this](const QString &msg) {runStatusLabel->setText(msg); statusBar()->showMessage(msg);});
        workerThread->start();
    }
    connect(ui->widget, &SpectrumPlot::rangeSelected, this, &MainWindow::runRoi);
    connect(ui->widget, &SpectrumPlot::energyAxisSet, this, &MainWindow::onEnergyAxisSet);
    connect(ui->widget, &SpectrumPlot::legendMoved, this, &MainWindow::onSpectrumLegendMoved);
    connect(ui->widget, &SpectrumPlot::graphVisibilityChanged, this, &MainWindow::updateListOfVisibleGraphs);

    updateDetectorList();

    ui->msgTextBrowser->append("\n");
    ui->msgTextBrowser->ensureCursorVisible();
    /* Modification flag of the document, not textChanged(), which is emitted also when only formats change (e.g. syntax highlighting after colors are applied) */
    connect(ui->editor->document(), &QTextDocument::modificationChanged, this, &MainWindow::setNeedsSaving);
}

void MainWindow::addMessage(jabs_msg_level level, const char *msg)
{
    emit sendMessage(level, QString(msg));
}

void MainWindow::reportProgress(double fraction, const char *msg)
{
    emit progressChanged(fraction, msg ? QString(msg) : QString(), false);
}

void MainWindow::reportProgressDone()
{
    emit progressChanged(0.0, QString(), true);
}

void MainWindow::setupRunPanel()
{ /* Panel below the main area, shown during runs. The central widget (from .ui) is put in a container together with the panel. */
    QWidget *main = takeCentralWidget();
    QWidget *container = new QWidget(this);
    QVBoxLayout *containerLayout = new QVBoxLayout(container);
    containerLayout->setContentsMargins(0, 0, 0, 0);
    containerLayout->setSpacing(0);
    containerLayout->addWidget(main, 1);
    runPanel = new QFrame(container);
    runPanel->setFrameShape(QFrame::StyledPanel);
    QVBoxLayout *panelLayout = new QVBoxLayout(runPanel);
    QHBoxLayout *progressRow = new QHBoxLayout();
    runStatusLabel = new QLabel(runPanel);
    runStatusLabel->setMinimumWidth(100);
    runStatusLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred); /* Long text does not widen the window */
    progressRow->addWidget(runStatusLabel, 3);
    progressBar = new QProgressBar(runPanel);
    progressBar->setRange(0, 0);
    progressBar->setTextVisible(true);
    progressBar->setFormat("%p%");
    progressRow->addWidget(progressBar, 2);
    QToolButton *stopButton = new QToolButton(runPanel);
    stopButton->setDefaultAction(actionStop);
    stopButton->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    progressRow->addWidget(stopButton);
    panelLayout->addLayout(progressRow);
    fitRow = new QWidget(runPanel);
    QHBoxLayout *fitLayout = new QHBoxLayout(fitRow);
    fitLayout->setContentsMargins(0, 0, 0, 0);
    fitPhaseLabel = new QLabel(fitRow);
    fitIterLabel = new QLabel(fitRow);
    fitChisqLabel = new QLabel(fitRow);
    fitLayout->addWidget(fitPhaseLabel);
    fitLayout->addSpacing(20);
    fitLayout->addWidget(fitIterLabel);
    fitLayout->addSpacing(20);
    fitLayout->addWidget(fitChisqLabel);
    fitLayout->addStretch(1);
    plotWhileFittingCheckBox = new QCheckBox(QString("Plot while fitting"), fitRow);
    connect(plotWhileFittingCheckBox, &QCheckBox::toggled, this, [this](bool checked) {
        plotWhileFitting = checked;
        settings.setValue("plotWhileFitting", checked);
        if(!checked && running) {
            ui->widget->setVisible(false);
        }
    });
    fitLayout->addWidget(plotWhileFittingCheckBox);
    endFitButton = new QPushButton(QString("End fit after this iteration"), fitRow);
    endFitButton->setToolTip(QString("Fit stops after the current iteration, result so far is kept"));
    connect(endFitButton, &QPushButton::clicked, this, [this]() {
        fitAbort = true;
        endFitButton->setEnabled(false);
        endFitButton->setText(QString("Ending after this iteration..."));
    });
    fitLayout->addWidget(endFitButton);
    panelLayout->addWidget(fitRow);
    containerLayout->addWidget(runPanel);
    setCentralWidget(container);
    showRunPanel(false);
}

void MainWindow::showRunPanel(bool show)
{
    runPanel->setVisible(show);
    if(show) {
        runStatusLabel->setText(QString("Running."));
        progressBar->setRange(0, 0); /* Busy indicator */
        fitRow->setVisible(false);
        plotWhileFittingCheckBox->blockSignals(true);
        plotWhileFittingCheckBox->setChecked(plotWhileFitting);
        plotWhileFittingCheckBox->blockSignals(false);
        endFitButton->setEnabled(true);
        endFitButton->setText(QString("End fit after this iteration"));
    }
}

void MainWindow::onProgressChanged(double fraction, QString msg, bool done)
{
    if(!running) {
        return;
    }
    if(done) {
        progressBar->setRange(0, 0); /* Something is still running, but no progress information */
        runStatusLabel->setText(QString("Running.")); /* Previous progress text (e.g. fit phase) may no longer apply */
        runStatusLabel->setToolTip(QString());
        return;
    }
    if(fraction < 0.0) {
        progressBar->setRange(0, 0); /* Busy indicator */
    } else {
        progressBar->setRange(0, 1000);
        progressBar->setValue((int)(1000.0 * qBound(0.0, fraction, 1.0)));
    }
    if(!msg.isEmpty()) {
        runStatusLabel->setText(msg);
        runStatusLabel->setToolTip(msg);
    }
}

MainWindow::~MainWindow()
{
    stopWorkerThread();
    closeSession();
    delete [] recentFileActs;
    delete highlighter;
    delete ui;
}

void MainWindow::openFile(const QString &filename)
{
    if(running) {
        QMessageBox::warning(this, "Run in progress", QString("Can not open file %1 while a simulation or fit is running.").arg(filename));
        return;
    }
    QFile file(filename);
    if(file.open(QFile::ReadOnly | QFile::Text)) {
        QByteArray data = file.readAll();
        if(data.endsWith('\n')) {
            ui->editor->setPlainText(data.sliced(0, data.size() -1));
        } else {
            ui->editor->setPlainText(data);
        }
        file.close();
        setFilename(filename);
        setNeedsSaving(false);
        resetAll();
        plotSession();  /* Will hide plot etc. */
    }
}

void MainWindow::runRoi(const QString &roi) {
    if(running) {
        return;
    }
    runLine(QString("roi %1 %2").arg(ui->comboBox->currentIndex() + 1).arg(roi));
}

void MainWindow::plotDialogClosed()
{
    plotDialog = NULL;
}

void MainWindow::readPlotSettings()
{
    showIsotopes = settings.value("showIsotopes", QVariant(false)).toBool();
    plotIsotopesZ = settings.value("plotIsotopesZ", QVariant(JIBAL_ANY_Z)).toInt();
    ui->widget->setLegendVisible(settings.value("showLegend", QVariant(true)).toBool());
    ui->widget->setLegendOutside(settings.value("legendOutside", QVariant(false)).toBool());
    ui->widget->setEnergyAxis(settings.value("energyAxis", QVariant(false)).toBool());
    confidenceLimits = settings.value("confidenceLimits", QVariant(true)).toBool();
    plotSession();
}


int MainWindow::runLine(const QString &line) { /* Runs line in GUI thread, only for quick commands (e.g. roi) */
    if(running) {
        return SCRIPT_COMMAND_FAILURE;
    }
    return JabsWorker::runLine(session, line);
}

void MainWindow::startWorker(const QStringList &lines, RunMode mode)
{
    running = true;
    runMode = mode;
    fitAbort = false;
    plotWhileFitting = settings.value("plotWhileFitting", QVariant(false)).toBool();
    snapshot.clear();
    snapshotPlotted = false;
    jabs_abort_set(FALSE);
    enableRun(false);
    showRunPanel(true);
    runTimer.start();
    statusBar()->showMessage(QString("Running."));
    emit startRun(lines);
}

void MainWindow::stopWorkerThread()
{
    if(!workerThread) {
        return;
    }
    jabs_abort_set(TRUE);
    fitAbort = true;
    workerThread->quit(); /* Event loop of thread quits after current run (if any) has been aborted */
    workerThread->wait();
    workerThread = nullptr;
    worker = nullptr;
}

void MainWindow::stopRun()
{
    if(!running) {
        return;
    }
    statusBar()->showMessage(QString("Stopping..."));
    fitAbort = true;
    jabs_abort_set(TRUE);
}

void MainWindow::plotSession(bool error)
{
    if(running) { /* Session can not be accessed, but latest snapshot (if any) can be plotted */
        if(snapshot && !error) {
            plotSpectrum(ui->comboBox->currentIndex());
        }
        return;
    }
    updateDetectorList();
    if(session && session->fit && session->fit->sim && !error) {
        size_t i_det = ui->comboBox->currentIndex();
        plotSpectrum(i_det);
        if(firstRun) {
            ui->widget->resetZoom();
        }
    } else {
        ui->detectorFrame->setVisible(false);
        ui->widget->setVisible(false);
    }
}

bool MainWindow::askToSave()
{
    if(needsSaving) {
        QMessageBox msgBox;
        msgBox.setText("The file has been modified.");
        if(filename.isEmpty()) {
            msgBox.setInformativeText("You haven't saved your changes. Cancel and try again, or discard changes.");
            msgBox.setStandardButtons(QMessageBox::Discard | QMessageBox::Cancel);
            msgBox.setDefaultButton(QMessageBox::Cancel);
        } else {
            msgBox.setInformativeText("Do you want to save your changes?");
            msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
            msgBox.setDefaultButton(QMessageBox::Save);
        }
        int ret = msgBox.exec();
        if(ret == QMessageBox::Cancel) {
            return false;
        }
        if(ret == QMessageBox::Save) {
            on_action_Save_File_triggered();
        }
    }
    return true;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    if(running) {
        QMessageBox::StandardButton reply = QMessageBox::question(this, "Run in progress", "A simulation or fit is running. Stop it and quit?", QMessageBox::Yes | QMessageBox::Cancel);
        if(reply != QMessageBox::Yes) {
            event->ignore();
            return;
        }
        stopRun();
    }
    if(askToSave()) {
        event->accept();
    } else {
        event->ignore();
    }
}

void MainWindow::enableRun(bool enabled)
{
    ui->action_Run->setEnabled(enabled || running); /* Stops a run */
    updateRunAction();
    ui->commandLineEdit->setEnabled(enabled);
    ui->action_New_File->setEnabled(enabled);
    ui->action_Open_File->setEnabled(enabled);
    ui->menuRecent_Files->setEnabled(enabled);
    ui->action_Save_File->setEnabled(enabled && needsSaving);
    ui->action_Save_File_as->setEnabled(enabled);
    ui->actionIDF->setEnabled(enabled);
    ui->menuDetector->setEnabled(enabled);
    ui->comboBox->setEnabled(enabled);
    ui->editor->setReadOnly(!enabled); /* Highlighter uses session (commands), which a running script can modify */
    if(actionStop) {
        actionStop->setEnabled(!enabled && running);
    }
}

void MainWindow::updateRunAction()
{
    if(running) {
        ui->action_Run->setIcon(stopIcon);
        ui->action_Run->setText(QString("&Stop"));
        ui->action_Run->setToolTip(QString("Stop running simulation or fit"));
        ui->action_Run->setShortcut(QKeySequence()); /* Pressing run shortcut twice does not stop a run */
    } else {
        ui->action_Run->setIcon(runIcon);
        ui->action_Run->setText(QString("&Run script"));
        ui->action_Run->setToolTip(QString("Reset and run script"));
        ui->action_Run->setShortcut(QKeySequence::Refresh);
    }
}


int MainWindow::initSession()
{
    QString license_filename_str;
#if defined(Q_OS_MACOS)
    license_filename_str = QApplication::applicationDirPath() + "/../Resources/LICENSE.txt";
#elif defined(Q_OS_WIN)
    license_filename_str = QApplication::applicationDirPath() + "\\LICENSE.txt";
#elif defined(Q_OS_LINUX)
    license_filename_str = "/usr/share/common-licenses/GPL-2";
#endif
    if(!license_filename_str.isEmpty() && QFile::exists(license_filename_str)) {
            ui->msgTextBrowser->insertHtml(QString("<p>%1</p>\n").arg(MainWindow::makeFileLink(license_filename_str, "License file")));
            ui->msgTextBrowser->insertPlainText("\n\n");
    }

    QString config_filename_str = settings.value("jibalConfigurationFile").toString();
    if(config_filename_str.isEmpty()) {
#if defined(Q_OS_MACOS)
    config_filename_str = QApplication::applicationDirPath() + "/../Resources/jibal.conf";
#elif defined(Q_OS_WIN)
    config_filename_str = QApplication::applicationDirPath() + "\\jibal.conf";
#else
    config_filename_str = QApplication::applicationDirPath() + "/jibal.conf";
#endif
    } else {
        if(!QFile::exists(config_filename_str)) {
            QMessageBox::information(this, "Warning", QString("JIBAL configuration file is supposed to be \"%1\", but this file does not exist. Ignoring it. Change the value in preferences.").arg(config_filename_str));
        }
    }
    if(config_filename_str.isEmpty() || !QFile::exists(config_filename_str) ) {
         jibal = jibal_init(NULL);
         ui->msgTextBrowser->insertHtml(QString("<p>JIBAL configuration file as determined by JIBAL: %1</p>\n\n").arg(MainWindow::makeFileLink(jibal_config_filename(jibal))));
    } else {
        ui->msgTextBrowser->insertHtml(QString("<p>JIBAL configuration file as determined by JaBS: %1</p>\n").arg(MainWindow::makeFileLink(config_filename_str)));
        jibal = jibal_init(qPrintable(config_filename_str));
    }
    ui->msgTextBrowser->insertPlainText("\n");
    if(!jibal) {
        QMessageBox::critical(this, "Error", "Could not initialize JIBAL (a NULL pointer was returned).\n");
        return -1;
    }
    ui->msgTextBrowser->insertPlainText("\n");
    ui->msgTextBrowser->insertPlainText(jibal_status_string(jibal));
    if(jibal->error) {
        QMessageBox::critical(this, "Error", QString("Could not initialize JIBAL:\n\n%1\n").arg(jibal_status_string(jibal)));
        return -1;
    }
    ui->msgTextBrowser->insertHtml(QString("<p>GSTO files defined in: %1</p>\n").arg(MainWindow::makeFileLink(jibal->config->files_file)));
    if(jibal->gsto->n_files == 0) {
        QString warning = "No GSTO files (stopping, straggling) have been defined. Please be adviced that JaBS is unable to do any real work.";
        warning.append(QString("You could try adding %1 to %2").arg(JIBAL_FILES_FILE, jibal->config->userdatadir));
        if(jibal->config->files_file) {
            warning.append(QString(" or you should check contents of file %1.").arg(jibal->config->files_file));
        } else {
            warning.append(".");
        }
        QMessageBox::warning(this, "Warning", warning);
    }
    session = script_session_init(jibal, NULL);
    if(!session)  {
        QMessageBox::critical(this, "Error", "Could not initialize JaBS session.\n");
        return -1;
    }
    session->fit_iter_callback = &fit_iter_callback;
    highlighter->setSession(session);
    return 0;
}

int MainWindow::closeSession()
{
    script_session_free(session);
    jibal_free(jibal);
    return 0;
}

QString MainWindow::makeFileLink(const QString &filename, const QString &linkname)
{
    if(filename.isEmpty() && linkname.isEmpty()) {
        return QString("(no filename)");
    }
    if(filename.isEmpty()) {
        return QString(linkname.toHtmlEscaped() + " (no filename)"); // This should not hopefully never be seen by the user
    }
    return QString("<a href=\"%1\">%2</a>").arg(QUrl::fromLocalFile(filename).toString(), linkname.isEmpty() ? filename.toHtmlEscaped() : linkname.toHtmlEscaped());
}

void MainWindow::readSettings()
{
    static const QString defaultFontFamily = QFontDatabase::systemFont(QFontDatabase::FixedFont).family();
    QString editorFontFamily = settings.value("editorFontFamily").toString();
    if(editorFontFamily.isEmpty()) {
        editorFontFamily = defaultFontFamily;
    }
    QFont editorFont;
    editorFont.setFamily(editorFontFamily);
    editorFont.setPointSize(settings.value("editorFontSize", 11).toInt());
    ui->editor->setFont(editorFont);

    QString messageFontFamily = settings.value("messageFontFamily").toString();
    if(messageFontFamily.isEmpty()) {
        messageFontFamily = defaultFontFamily;
    }
    QFont messageFont;
    messageFont.setFamily(messageFontFamily);
    messageFont.setPointSize(settings.value("messageFontSize", 10).toInt());
    ui->msgTextBrowser->setFont(messageFont);
    defaultVerbosity = settings.value("defaultVerbosity", JABS_DEFAULT_VERBOSITY).toInt();
    applyColors();
}

void MainWindow::applyColors()
{
    if(applyingColors) {
        return;
    }
    applyingColors = true;
    colorProfile = ColorProfile::active(settings);
#if QT_VERSION >= QT_VERSION_CHECK(6, 8, 0)
    const Qt::ColorScheme scheme = colorProfile.colorScheme(); /* Light or dark appearance of the whole application (window, toolbar, dialogs), or follow OS */
    if(scheme != forcedColorScheme) {
        forcedColorScheme = scheme;
        QGuiApplication::styleHints()->setColorScheme(scheme); /* Palette change follows, colors will be applied again */
    }
#endif
    const QPalette appPalette = QApplication::palette();

    /* Only roles set in profile are set in widget palettes, others ("System") are inherited from application palette and follow its changes */
    QPalette editorPalette;
    if(colorProfile.color("editorBackground").isValid()) {
        const QColor bg = colorProfile.color("editorBackground");
        editorPalette.setColor(QPalette::Base, bg);
        editorPalette.setColor(QPalette::Button, bg.lightness() < 128 ? bg.lighter(130) : bg.darker(108)); /* Line number area */
    }
    if(colorProfile.color("editorText").isValid()) {
        const QColor text = colorProfile.color("editorText");
        editorPalette.setColor(QPalette::Text, text);
        QColor lineNumbers = text;
        lineNumbers.setAlpha(150);
        editorPalette.setColor(QPalette::ButtonText, lineNumbers);
    }
    ui->editor->setPalette(editorPalette);
    QTextCharFormat format;
    format.setForeground(colorProfile.color("editorComment"));
    highlighter->setCommentFormat(format);
    format.setForeground(colorProfile.color("editorCommand"));
    format.setFontWeight(QFont::Bold);
    highlighter->setCommandFormat(format);
    format.setForeground(colorProfile.color("editorVariable"));
    format.setFontWeight(QFont::StyleItalic);
    highlighter->setVariableFormat(format);
    highlighter->rehighlight();

    QPalette messagePalette;
    if(colorProfile.color("messageBackground").isValid()) {
        messagePalette.setColor(QPalette::Base, colorProfile.color("messageBackground"));
    }
    if(colorProfile.color("messageText").isValid()) {
        messagePalette.setColor(QPalette::Text, colorProfile.color("messageText"));
    }
    ui->msgTextBrowser->setPalette(messagePalette);

    ui->widget->setColors(colorProfile.color("plotBackground"), colorProfile.color("plotForeground"), colorProfile.color("plotGrid"), colorProfile.color("plotLegendBackground"));

    const QColor iconColor = colorProfile.color("iconColor", appPalette.color(QPalette::ButtonText));
    runIcon = ColorProfile::themedIcon(":/icons/play-box.svg", iconColor);
    stopIcon = ColorProfile::themedIcon(":/icons/stop-box.svg", iconColor);
    updateRunAction();
    const QList<QPair<QAction *, QString>> icons = {
        {ui->action_New_File, ":/icons/file-document.svg"},
        {ui->action_Open_File, ":/icons/folder-open.svg"},
        {ui->action_Save_File, ":/icons/content-save.svg"},
        {ui->actionNext_detector, ":/icons/chevron-right.svg"},
        {ui->actionPrevious_detector, ":/icons/chevron-left.svg"},
        {actionStop, ":/icons/stop-box.svg"}
    };
    for(const auto &icon : icons) {
        if(icon.first) {
            icon.first->setIcon(ColorProfile::themedIcon(icon.second, iconColor));
        }
    }
    plotSession();
    applyingColors = false;
}

void MainWindow::scheduleApplyColors()
{
    if(applyColorsPending) {
        return;
    }
    applyColorsPending = true;
    QTimer::singleShot(100, this, [this]() {
        applyColorsPending = false;
        applyColors();
    });
}

void MainWindow::changeEvent(QEvent *event)
{
    if(event->type() == QEvent::ApplicationPaletteChange || event->type() == QEvent::PaletteChange || event->type() == QEvent::ThemeChange) { /* System colors may have changed */
        scheduleApplyColors();
    }
    QMainWindow::changeEvent(event);
}

void MainWindow::setNeedsSaving(bool value)
{
    needsSaving = value;
    if(!value) { /* Text as it is now is saved (or new) */
        ui->editor->document()->setModified(false);
    }
    updateWindowTitle();
}

void MainWindow::showInitialMessages()
{
    ui->msgTextBrowser->insertPlainText(aboutString);
    ui->msgTextBrowser->insertPlainText(QString(COPYRIGHT_STRING));
    ui->msgTextBrowser->insertHtml("<p>Visit the <a href=\"https://github.com/JYU-IBA/jabs\">GitHub page</a> for latest information.</p>");
    ui->msgTextBrowser->insertPlainText("\n\n");
}

void MainWindow::updateListOfVisibleGraphs()
{
    if(ui->widget->isVisible()) {
        visibleGraphs = ui->widget->visibleGraphs();
    }
}

void MainWindow::scrollMsgBoxToBottom()
{
    ui->msgTextBrowser->ensureCursorVisible();
    //ui->msgTextBrowser->verticalScrollBar()->setValue(ui->msgTextBrowser->verticalScrollBar()->maximum());
}

void MainWindow::receiveMessage(jabs_msg_level level, QString msg)
{
    ui->msgTextBrowser->moveCursor(QTextCursor::End);
    switch(level) {
    case MSG_ERROR:
        ui->msgTextBrowser->setTextColor(colorProfile.color("messageError"));
        break;
    case MSG_WARNING:
        ui->msgTextBrowser->setTextColor(colorProfile.color("messageWarning"));
        warningCounter++;
        break;
    case MSG_DEBUG:
        ui->msgTextBrowser->setTextColor(colorProfile.color("messageDebug"));
        break;
    default:
        ui->msgTextBrowser->setTextColor(ui->msgTextBrowser->palette().color(QPalette::Text));
        break;
    }
    ui->msgTextBrowser->insertPlainText(msg);
}

int MainWindow::fitCallback(fit_stats stats)
{ /* Called by the worker thread. Session can be accessed here (the worker is not using it), widgets can not. */
    QSharedPointer<PlotSnapshot> snap;
    if(plotWhileFitting && session && session->fit && session->fit->spectra) {
        snap = QSharedPointer<PlotSnapshot>::create();
        snap->spectra.resize(session->fit->n_det_spectra);
        for(size_t i = 0; i < session->fit->n_det_spectra; i++) {
            result_spectra_copy(&snap->spectra[i], &session->fit->spectra[i]);
        }
        snap->ref = session->fit->ref ? jabs_histogram_clone(session->fit->ref) : nullptr;
    }
    emit fitProgress(stats, snap);
    return fitAbort ? 1 : 0;
}

void MainWindow::onFitProgress(fit_stats stats, QSharedPointer<PlotSnapshot> snap)
{
    if(!running) { /* Late signal, run has already finished */
        return;
    }
    if(stats.phase == FIT_PHASE_NONE) { /* Fit is over, script may continue (e.g. simulate) */
        fitRow->setVisible(false);
        return;
    }
    fitRow->setVisible(true);
    fitPhaseLabel->setText(QString("Phase: %1").arg(fit_phase_to_str((fit_phase) stats.phase)));
    fitIterLabel->setText(QString("Iteration: %1").arg(stats.iter));
    fitChisqLabel->setText(stats.chisq_dof > 0.0 ? QString("χ²/dof: %1").arg(stats.chisq_dof, 0, 'g', 6) : QString("χ²/dof: -"));
    if(plotWhileFittingCheckBox->isChecked()) {
        if(snap) {
            snapshot = snap;
            ui->widget->setVisible(true);
            plotSpectrum(ui->comboBox->currentIndex());
            if(!snapshotPlotted) { /* First plot during this run: axes are set to data (like plotSession() does on first run), later the user's zoom is kept */
                snapshotPlotted = true;
                if(firstRun) {
                    ui->widget->resetZoom();
                }
            }
        }
    } else {
        ui->widget->setVisible(false);
    }
}

void MainWindow::on_action_Run_triggered()
{
    if(!session || !worker) {
        QMessageBox::critical(this, "Error", "Session not initialized.\n");
        return;
    }
    if(running) { /* Run action is a stop action while running */
        stopRun();
        return;
    }
    if(firstRun) {
        resetAll();
    } else {
        updateListOfVisibleGraphs();
        ui->msgTextBrowser->clear();
        script_reset(session, 0, NULL);
    }
    if(!filedirpath.isEmpty()) {
        QDir::setCurrent(filedirpath);
    }
    QString text = ui->editor->toPlainText();
    QTextStream stream = QTextStream(&text, QIODevice::ReadOnly);
    QStringList lines;
    while(!stream.atEnd()) {
        lines.append(stream.readLine());
    }
    warningCounter = 0;
    jabs_message_verbosity = (jabs_msg_level) defaultVerbosity;
    startWorker(lines, RUN_SCRIPT);
}

void MainWindow::onRunFinished(bool error, qsizetype lineno)
{
    running = false;
    showRunPanel(false);
    snapshot.clear();
    bool aborted = jabs_abort_requested() || fitAbort;
    jabs_abort_set(FALSE);
    enableRun(true);
    QString message;
    if(runMode == RUN_COMMAND_LINE) {
        if(!error) {
            ui->commandLineEdit->clear();
        }
        plotSession();
        message = error ? QString("Command failed.") : QString("Command successful in %1 seconds.").arg(runTimer.elapsed()/1000.0);
    } else {
        plotSession();
        if(error) {
            message = aborted ? QString("Run stopped on line %1.").arg(lineno) : QString("Error on line %1. Run aborted.").arg(lineno);
        } else {
            message = QString("Run successful in %1 seconds, %2 lines processed.").arg(runTimer.elapsed()/1000.0).arg(lineno);
            firstRun = false;
        }
    }
    if(warningCounter == 1) {
        message.append(QString(" 1 warning."));
    } else if(warningCounter > 1) {
        message.append(QString(" %1 warnings.").arg(warningCounter));
    }
    statusBar()->showMessage(message, 5000);
    QTimer::singleShot(0, this, &MainWindow::scrollMsgBoxToBottom); /* Messages come from queue somewhere, this makes scrolling happen after that */
    emit runFinished();
}


void MainWindow::plotSpectrum(size_t i_det)
{
    if(running) {
        if(snapshot && i_det < snapshot->spectra.size()) {
            plotSpectra(&snapshot->spectra[i_det], snapshot->ref);
        }
        return;
    }
    if(!session || !session->fit || !session->fit->spectra || i_det >= session->fit->n_det_spectra) {
        ui->widget->clearAll();
        ui->widget->setVisible(false);
        return;
    }
    plotSpectra(&session->fit->spectra[i_det], session->fit->ref);
}

void MainWindow::plotSpectra(const result_spectra *spectra, const jabs_histogram *ref_histo)
{
    ui->widget->clearAll();
    if(ref_histo) {
        QCPGraph *g = ui->widget->addGraphHisto("Reference", colorProfile.color("plotReference"));
        ui->widget->drawDataToGraph(g, ref_histo);
    }
    if(spectra) {
        if(spectra->n_spectra > RESULT_SPECTRA_SIMULATED) {
            QCPGraph *g = ui->widget->addGraphHisto("Simulated", colorProfile.color("plotSimulated"));
            ui->widget->drawDataToGraph(g, result_spectra_simulated_histo(spectra));
            ui->widget->setGraphVisibility(g, false);
        }
        if(spectra->n_spectra > RESULT_SPECTRA_EXPERIMENTAL) {
            QCPGraph *g = ui->widget->addGraphHisto("Experimental", colorProfile.color("plotExperimental"));
            ui->widget->drawDataToGraph(g, result_spectra_experimental_histo(spectra));
            ui->widget->setGraphVisibility(g, false);
        }
        if(confidenceLimits && spectra->n_spectra > RESULT_SPECTRA_UNCERTAINTY_POSITIVE) {
            ui->widget->drawFilledDataToGraphs("Confidence limit", result_spectra_uncertainty_histo_negative(spectra), result_spectra_uncertainty_histo_positive(spectra), colorProfile.color("plotConfidence"));
        }
        size_t n_spectra = spectra->n_spectra;
        jabs_histogram *histo = NULL;
        int colorindex = 0;
        for(int i = RESULT_SPECTRA_N_FIXED; i < n_spectra; ++i) {
            const result_spectrum *s = &spectra->s[i];
            if(!s || !s->histo || !s->histo->n) {
                continue;
            }
            int Z_this = s->target_isotope ? s->target_isotope->Z : JIBAL_ANY_Z;
            /* TODO: spectrum can be empty (e.g. if reaction is not possible), detect and continue without plotting it */
            const result_spectrum *s_next;
            int Z_next;
            if(i+1 < n_spectra) {
                s_next = &spectra->s[i + 1];
                Z_next = s_next->target_isotope ? s_next->target_isotope->Z : JIBAL_ANY_Z;
            } else {
                s_next = NULL;
                Z_next = JIBAL_ANY_Z;
            }

            if(histo) {
                for(size_t i = 0; i < histo->n && i < s->histo->n; i++) {
                    /* ignores different ranges in GSL histograms so energy histograms could, in theory, be different */
                    histo->bin[i] += s->histo->bin[i];
                }
            } else {
                histo = jabs_histogram_clone(s->histo);
            }
            bool nextMatches = s_next && s->type == s_next->type && s->ds == s_next->ds && Z_this == Z_next && Z_this != JIBAL_ANY_Z;
            bool showThisIsotope = showIsotopes &&  (plotIsotopesZ == JIBAL_ANY_Z || plotIsotopesZ == Z_this);
            if(nextMatches && !showThisIsotope) {
                continue;
            }
            if(showThisIsotope || !nextMatches) {
            /* Plot histo (accumulated histogram) if:
                 * 1. If show isotopes is set and Z is JIBAL_ANY_Z or we have a match with plotIsotopesZ
                 * 2. If this is the last reaction (there is no s_next), therefore the last histogram
                 * 3. Next reaction has different type (e.g. RBS vs ERD) or different Z
                 * */
                QString name;
                if(s->target_isotope == NULL || s->type == 0) {
                    name = s->name;
                } else {
                    name = QString("") + reaction_type_to_string(s->type) + " ";
                    if(showThisIsotope) {
                        name.append(s->target_isotope->name);
                    } else {
                        name.append(jibal_element_name(jibal->elements, Z_this));
                    }
                    if(s->ds) {
                        name.append(" (DS)");
                    }
                }
                QColor color = colorProfile.reactionColor(colorindex);
                colorindex++;
                QCPGraph *g = ui->widget->addGraphHisto(name, color);
                ui->widget->drawDataToGraph(g, histo);
                ui->widget->setGraphVisibility(g, false);
                jabs_histogram_free(histo);
                histo = NULL;
            }
        }
    }


    if(ui->widget->graphCount() == 0) {
        ui->widget->setVisible(false);
        return;
    }
    ui->widget->setVisible(true);
    if(ui->widget->visibleGraphs().isEmpty()) {
        visibleGraphs.append("Simulated");
        visibleGraphs.append("Experimental");
        if(confidenceLimits) {
            visibleGraphs.append("Confidence limit");
        }
    }
    for(int i = 0; i < ui->widget->graphCount(); ++i) {
        QString name = ui->widget->graph(i)->name();
        ui->widget->setGraphVisibility(ui->widget->graph(i), visibleGraphs.contains(name));

    }
    ui->widget->updateVerticalRange();
    ui->widget->replot(QCustomPlot::rpQueuedReplot); /* Widget may have just been made visible (or resized), replot after layout has been updated */
}


void MainWindow::on_action_Open_File_triggered()
{
    if(!MainWindow::askToSave())
        return;
    QString filename = QFileDialog::getOpenFileName(this, "Open script file", "", tr("JaBS script files (*.jbs)"));
    if(filename.isEmpty()) {
        return;
    }
    openFile(filename);
}


void MainWindow::on_action_Quit_triggered()
{
    close();
}

void MainWindow::updateWindowTitle()
{
    bool newfile;
    QString title = QString("JaBS ") + jabs_version_simple() + " - ";
    if(filename.isEmpty()) {
        title.append("New file");
        newfile = true;
    } else {
        title.append(filebasename);
        newfile = false;
    }
    if(needsSaving) {
        title.append("*");
    } else {
        ui->action_Save_File->setEnabled(false);
    }
    ui->action_Save_File->setEnabled(needsSaving && !newfile);
    setWindowTitle(title);
}

void MainWindow::setFilename(const QString &filename)
{
    MainWindow::filename = filename;
    QFile file(filename);
    QFileInfo fi(file);
    MainWindow::filebasename = fi.baseName();
    filedirpath = fi.absoluteDir().path();
    setWindowFilePath(fi.absoluteFilePath());
    if(!QDir::setCurrent(filedirpath)) {
        qDebug() << "Can't set the current working directory!";
    }
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(filename);
    files.prepend(filename);
    while (files.size() > maxRecentFiles)
        files.removeLast();

    settings.setValue("recentFileList", files);
    updateRecentFileActions();
}


void MainWindow::on_action_New_File_triggered()
{
    if(!askToSave()) {
        return;
    }
    ui->editor->clear();
    resetAll();
    plotSession(); /* Will hide plot */
    filename.clear();
    filedirpath.clear();
    filebasename.clear();
    QDir::setCurrent(originalPath);
    setNeedsSaving(false);
    showInitialMessages();
}



void MainWindow::on_action_Save_File_as_triggered()
{
    QString filename = QFileDialog::getSaveFileName(this, "Save script to file", MainWindow::filename, tr("JaBS script files (*.jbs)"));
    if(!filename.isEmpty()) {
        if(saveScriptToFile(filename)) {
            setFilename(filename);
            setNeedsSaving(false);
        }
    }
}

bool MainWindow::saveScriptToFile(const QString &filename)
{
    QFile file(filename);
    if(file.open(QFile::WriteOnly | QFile::Text)) {
        QTextStream stream(&file);
        stream << ui->editor->toPlainText();
        stream << '\n';
        file.close();
    } else {
         QMessageBox::critical(this, tr("Error"), tr("Can not save to file."));
         return false;
    }
    return true;
}


void MainWindow::on_action_Save_File_triggered()
{
    if(filename.isEmpty())
        return;
    if(saveScriptToFile(filename)) {
        setNeedsSaving(false);
    }
}

void MainWindow::resetAll()
{
    ui->widget->clearAll();
    ui->widget->replot();
    ui->msgTextBrowser->clear();
    firstRun = true;
    script_reset(session, 0, NULL);
    visibleGraphs.clear();
    visibleGraphs.append("Simulated");
    visibleGraphs.append("Experimental");
    warningCounter = 0;
}


void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("QJaBS"), aboutString + "\n" + QString(COPYRIGHT_STRING).simplified());
}

void MainWindow::on_commandLineEdit_returnPressed()
{
    if(running || !worker) {
        return;
    }
    warningCounter = 0;
    updateListOfVisibleGraphs();
    startWorker(QStringList(ui->commandLineEdit->text()), RUN_COMMAND_LINE);
}

void MainWindow::on_msgTextBrowser_textChanged()
{
    ui->msgTextBrowser->moveCursor(QTextCursor::End);
    ui->msgTextBrowser->ensureCursorVisible();
}

void MainWindow::updateRecentFileActions()
{
    QSettings settings;
    QStringList files = settings.value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), maxRecentFiles);

    for(int i = 0; i < numRecentFiles; ++i) {
        QAction *a = &recentFileActs[i];
        QString text = QFileInfo(files[i]).fileName();
        a->setText(text);
        a->setData(files[i]);
        a->setVisible(true);
        a->setToolTip(files[i]);
    }
    for(int j = numRecentFiles; j < maxRecentFiles; ++j)
        recentFileActs[j].setVisible(false);
}

void MainWindow::openRecentFile()
{
    if(!askToSave()) {
        return;
    }
    QAction *action = qobject_cast<QAction *>(sender());
    if(action)
        openFile(action->data().toString());
}

void MainWindow::openLink(const QUrl &link)
{
    QDesktopServices::openUrl(link);
}

void MainWindow::on_action_Plot_triggered()
{
    if(plotDialog) {
        plotDialog->show();
        plotDialog->raise();
        return;
    }
    if(!jibal)
        return;
    plotDialog = new PlotDialog(this, jibal);
    connect(plotDialog, &PlotDialog::closed, this, &MainWindow::plotDialogClosed);
    connect(plotDialog, &PlotDialog::settingsSaved, this, &MainWindow::readPlotSettings);
    plotDialog->show();
}


void MainWindow::on_actionPreferences_triggered()
{
    PreferencesDialog *dialog = new PreferencesDialog(this);
    dialog->setModal(true);
    connect(dialog, &PreferencesDialog::settingsSaved, this, &MainWindow::readSettings);
    dialog->open();
}


void MainWindow::on_actionIDF_triggered()
{
    if(!MainWindow::askToSave())
        return;
    QMessageBox::StandardButton reply;
    QString filename = QFileDialog::getOpenFileName(this, "Import IDF file", "", tr("IDF files (*.idf;*.xml;*.xnra)"));
    if(filename.isEmpty()) {
        return;
    }
    QString question = QString("The program will now convert the IDF file %1 to a JaBS script (file extension .jbs) and a spectra (if applicable; extension(s) .dat).\n\nThese files, if they already exist, may be overwritten.\n\nContinue?").arg(filename);
    reply = QMessageBox::question(this, "Possible overwrite", question, QMessageBox::Ok|QMessageBox::Cancel);
    if(reply != QMessageBox::Ok) {
        return;
    }
    char *filename_out = nullptr;
    idf_error idferr = idf_parse(qPrintable(filename), &filename_out);
    if(idferr == IDF2JBS_SUCCESS) {
        openFile(filename_out);
    } else {
        QMessageBox::critical(this, "Error", QString("Could not import file %1.\n\nidf2jbs failed with error: %2.\n").arg(filename).arg(idf_error_code_to_str(idferr)));
    }
}

void MainWindow::onEnergyAxisSet(bool value)
{
    size_t i_det = ui->comboBox->currentIndex();
    plotSpectrum(i_det);
    ui->widget->resetZoom();
    settings.setValue("energyAxis", value);
}

void MainWindow::onSpectrumLegendMoved(bool outside)
{
    settings.setValue("legendOutside", outside);
}

void MainWindow::updateDetectorList()
{
    if(!session || !session->fit || !session->fit->sim) {
        ui->detectorFrame->setVisible(false);
        ui->actionNext_detector->setVisible(false);
        ui->actionPrevious_detector->setVisible(false);
        return;
    }
    bool multiple_detectors = (session->fit->sim->n_det > 1);
    ui->detectorFrame->setVisible(multiple_detectors);
    ui->actionPrevious_detector->setVisible(multiple_detectors);
    ui->actionNext_detector->setVisible(multiple_detectors);
    ui->menuDetector->setEnabled(multiple_detectors);
    int new_index = ui->comboBox->currentIndex() < session->fit->sim->n_det ? ui->comboBox->currentIndex() : 0;
    ui->comboBox->blockSignals(true);
    ui->comboBox->clear();
    for(size_t i_det = 0; i_det < session->fit->sim->n_det; i_det++) {
        ui->comboBox->addItem(QString(session->fit->sim->det[i_det]->name));
    }
    ui->comboBox->update();
    ui->comboBox->setCurrentIndex(new_index);
    ui->comboBox->blockSignals(false);
}


void MainWindow::on_comboBox_currentIndexChanged(int index)
{
    if(running) {
        return;
    }
    ui->actionPrevious_detector->setEnabled(index > 0);
    ui->actionNext_detector->setEnabled(index < session->fit->sim->n_det - 1);
    plotSpectrum(index);
}


void MainWindow::on_actionNext_detector_triggered()
{
    ui->comboBox->setCurrentIndex(ui->comboBox->currentIndex() + 1);
}


void MainWindow::on_actionPrevious_detector_triggered()
{
    ui->comboBox->setCurrentIndex(ui->comboBox->currentIndex() - 1);
}

