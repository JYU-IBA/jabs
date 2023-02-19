#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFontDatabase>
#include <QDebug>
#include <QTextStream>
#include "preferencesdialog.h"


extern "C" {
#include "script.h"
#include "script_session.h"
#include "script_command.h"
#include "idf2jbs.h"
int fit_iter_callback(fit_stats stats);
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow),
      maxRecentFiles(5), fitDialog(NULL), jibal(NULL), session(NULL), highlighter(NULL), plotDialog(NULL),
      showIsotopes(false), plotIsotopesZ(JIBAL_ANY_Z)
{
    aboutString = QString("JaBS version ") + jabs_version() + "\n\n"
                       + "Using JIBAL version "+ jibal_version() + ", compiled using version " + JIBAL_VERSION + "\n\n"
                       + "Using Qt version " + qVersion() + ", compiled using version " + QT_VERSION_STR + "\n\n"
                       + "Copyright 2021 - 2023 Jaakko Julin <jaakko.julin@jyu.fi>\n";
    ui->setupUi(this);
    QIcon icon(":/icons/jabs.svg");
    QApplication::setWindowIcon(icon);
    ui->widget->setVisible(false);
    ui->detectorFrame->setVisible(false); /* Will be made visible if necessary */
    setWindowIcon(icon);
    ui->msgTextBrowser->setOpenLinks(false); /* We use the connection below to handle links */
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
    readSettings();
    readPlotSettings();
    int status = initSession();
    if(status != 0) {
        enableRun(false);
        ui->action_Plot->setEnabled(false);
    }
    connect(ui->widget, &SpectrumPlot::rangeSelected, this, &MainWindow::runRoi);
    connect(ui->widget, &SpectrumPlot::energyAxisSet, this, &MainWindow::onEnergyAxisSet);
    connect(ui->widget, &SpectrumPlot::legendMoved, this, &MainWindow::onSpectrumLegendMoved);
    connect(ui->widget, &SpectrumPlot::graphVisibilityChanged, this, &MainWindow::updateListOfVisibleGraphs);

    ui->msgTextBrowser->append("\n");
    ui->msgTextBrowser->ensureCursorVisible();
    ui->editor->blockSignals(true); /* extremely dirty hack:
    on Linux (for whatever reason) textChanged() gets emitted although the text does not change.
    Maybe it is related to highlighter or clipboard, who knows.

    The workaround is to wait 100 ms (0 was not enough, 50 ms was not enough) and then allow signals again. */
    QTimer::singleShot(100, this, &MainWindow::postInit);
}

void MainWindow::addMessage(jabs_msg_level level, const char *msg)
{
    ui->msgTextBrowser->moveCursor(QTextCursor::End);
    switch(level) {
    case MSG_ERROR:
        ui->msgTextBrowser->setTextColor(Qt::red);
        break;
    case MSG_WARNING:
        ui->msgTextBrowser->setTextColor(Qt::darkYellow);
        warningCounter++;
        break;
    case MSG_DEBUG:
        ui->msgTextBrowser->setTextColor(Qt::gray);
        break;
    default:
        ui->msgTextBrowser->setTextColor(Qt::black); /* TODO: defaults from theme? */
        break;
    }
    ui->msgTextBrowser->insertPlainText(msg);
    //repaint();
}

MainWindow::~MainWindow()
{
    closeSession();
    delete [] recentFileActs;
    delete highlighter;
    delete ui;
}

void MainWindow::openFile(const QString &filename)
{
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
    runLine(QString("roi %1 %2").arg(ui->comboBox->currentIndex()).arg(roi));
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
    plotSession();
}


int MainWindow::runLine(const QString &line) {
    if(line.isEmpty()) {
        return SCRIPT_COMMAND_SUCCESS;
    }
    jabs_message(MSG_INFO, stderr, "%s%s\n", PROMPT, qPrintable(line));
    int status = script_execute_command(session, qPrintable(line));
    return status;
}

void MainWindow::plotSession(bool error)
{
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
    if(askToSave()) {
        event->accept();
    } else {
        event->ignore();
    }
}

void MainWindow::enableRun(bool enabled)
{
    ui->action_Run->setEnabled(enabled);
    ui->commandLineEdit->setEnabled(enabled);
    ui->action_New_File->setEnabled(enabled);
    ui->action_Open_File->setEnabled(enabled);
    ui->menuRecent_Files->setEnabled(enabled);
    ui->action_Save_File->setEnabled(enabled && needsSaving);
    ui->action_Save_File_as->setEnabled(enabled);
}

void MainWindow::closeFitDialog()
{
    if(fitDialog) {
        fitDialog->close();
        delete fitDialog;
        fitDialog = NULL;
    }
}

int MainWindow::initSession()
{
    QString config_filename_str = settings.value("jibalConfigurationFile").toString();
    if(config_filename_str.isEmpty()) {
#if defined(Q_OS_OSX)
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
         ui->msgTextBrowser->insertHtml(QString("JIBAL configuration file (determined by JIBAL): %1\n\n").arg(MainWindow::makeFileLink(jibal_config_filename(jibal))));
    } else {
        ui->msgTextBrowser->insertHtml(QString("<p>JIBAL configuration file (determined by JaBS): %1</p>\n").arg(MainWindow::makeFileLink(config_filename_str)));
        ui->msgTextBrowser->insertPlainText("\n");
        jibal = jibal_init(qPrintable(config_filename_str));
    }
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
    highlighter = new Highlighter(ui->editor->document());
    highlighter->setSession(session);
    return 0;
}

int MainWindow::closeSession()
{
    script_session_free(session);
    jibal_free(jibal);
    return 0;
}

QString MainWindow::makeFileLink(const QString &filename)
{
    if(filename.isEmpty()) {
        return QString("(no filename)");
    }
    return QString("<a href=\"%1\">%2</a>").arg(QUrl::fromLocalFile(filename).toString(), filename.toHtmlEscaped());
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
}

void MainWindow::setNeedsSaving(bool value)
{
    needsSaving = value;
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

int MainWindow::fitCallback(fit_stats stats)
{
    if(!fitDialog) {
        fitDialog = new FitDialog();
        fitDialog->show();
    }
    fitDialog->updateStats(stats);
    if(fitDialog->isPlotWhileFitting()) {
        plotSession(/*stats.error < 0*/);
    } else {
        ui->widget->setVisible(false);
    }
    QCoreApplication::processEvents();
    return fitDialog->isAborted();
}

void MainWindow::on_action_Run_triggered()
{
    if(!session) {
        QMessageBox::critical(this, "Error", "Session not initialized.\n");
        return;
    }
    enableRun(false);
    if(firstRun) {
        resetAll();
    } else {
        updateListOfVisibleGraphs();
        ui->msgTextBrowser->clear();
        script_reset(session, 0, NULL);
    }
    QString text = ui->editor->toPlainText();
    QTextStream stream = QTextStream(&text, QIODevice::ReadOnly);
    size_t lineno = 0;
    bool error = false;
    warningCounter = 0;
    jabs_message_verbosity = (jabs_msg_level) defaultVerbosity;
    QElapsedTimer refreshTimer;
    refreshTimer.start();
    QElapsedTimer runTimer;
    runTimer.start();
    while(!stream.atEnd() && !error) {
        QString line = stream.readLine();
        lineno++;
        QString command_str = line.split(" ").first();
        const script_command *c = script_command_find(session->commands, qPrintable(command_str));
        if(refreshTimer.elapsed() > 200 || c && (c->f == script_fit || c->f == script_simulate)) { /* Make sure stuff gets shown before sim or fit starts, or occasionally with a timer */
            refreshTimer.restart();
            if(c && c->f == script_fit) {
                statusBar()->showMessage(QString("Running a fit."));
            } else if(c && c->f == script_simulate) {
                statusBar()->showMessage(QString("Running a simulation."));
            } else {
                statusBar()->showMessage(QString("Running."));
            }
            repaint();
            QCoreApplication::processEvents();
        }
        if(runLine(line) < 0) {
            error = true;
            break;
        }
        if(session->file_depth > 0) {
            if(script_process(session) < 0) {
                error = true;
                break;
            }
        }
        //closeFitDialog(); /* if fit dialog was opened (by callback resulting in running this line), close it immediately*/
    }
    closeFitDialog();
    enableRun(true);
    plotSession(/*error */);
    QString message;
    if(error) {
        message = QString("Error on line %1. Run aborted.").arg(lineno);
    } else {
        message = QString("Run successful in %1 seconds, %2 lines processed.").arg(runTimer.elapsed()/1000.0).arg(lineno);
        firstRun = false;
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
    ui->widget->clearAll();
    if(!session || !session->fit || !session->fit->spectra || i_det >= session->fit->spectra->n_spectra) {
        ui->widget->setVisible(false);
        return;
    }
    gsl_histogram *ref_histo = session->fit->ref;
    if(ref_histo) {
        ui->widget->drawDataToChart("Reference", ref_histo->range, ref_histo->bin, ref_histo->n, QColor("Gray"));
    }
    result_spectra *spectra = &session->fit->spectra[i_det];
    if(spectra) {
        size_t n_spectra = spectra->n_spectra;
        gsl_histogram *histo = NULL;
        int colorindex = 0;
        for(int i = 0; i < n_spectra; ++i) {
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
                histo = gsl_histogram_clone(s->histo);
            }
            bool nextMatches = s_next && s->type == s_next->type && Z_this == Z_next && Z_this != JIBAL_ANY_Z;
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
                }
                QColor color;
                switch(i) {
                case RESULT_SPECTRA_SIMULATED:
                    color = QColor("Blue");
                    break;
                case RESULT_SPECTRA_EXPERIMENTAL:
                    color = QColor("Black");
                    break;
                default:
                    color = SpectrumPlot::getColor(colorindex);
                    colorindex++;
                    break;
                }
                ui->widget->drawDataToChart(name, histo->range, histo->bin, histo->n, color);
                ui->widget->setGraphVisibility(ui->widget->graph(), false);
                gsl_histogram_free(histo);
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
    }
    for(int i = 0; i < ui->widget->graphCount(); ++i) {
        QString name = ui->widget->graph(i)->name();
        ui->widget->setGraphVisibility(ui->widget->graph(i), visibleGraphs.contains(name));

    }
    ui->widget->updateVerticalRange();
    ui->widget->replot();
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
    setWindowFilePath(fi.absoluteFilePath());
    if(!QDir::setCurrent(fi.absolutePath())) {
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
    filebasename.clear();
    QDir::setCurrent(originalPath);
    setNeedsSaving(false);
    showInitialMessages();
}


void MainWindow::on_editor_textChanged()
{
    if(!needsSaving) {
        setNeedsSaving(true);
    }
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
    warningCounter = 0;
    updateListOfVisibleGraphs();
    int status = runLine(ui->commandLineEdit->text());
    if(status == SCRIPT_COMMAND_SUCCESS) {
            ui->commandLineEdit->clear();
            plotSession();
    }
    closeFitDialog();
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

void MainWindow::postInit()
{
    ui->editor->blockSignals(false);
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
        return;
    }
    ui->detectorFrame->setVisible(session->fit->sim->n_det > 1);
    int old_i_det = ui->comboBox->currentIndex();
    ui->comboBox->clear();
    for(size_t i_det = 0; i_det < session->fit->sim->n_det; i_det++) {
        ui->comboBox->addItem(QString(session->fit->sim->det[i_det]->name));
    }
    ui->comboBox->setCurrentIndex(old_i_det < session->fit->sim->n_det ? old_i_det : 0);
}


void MainWindow::on_comboBox_currentIndexChanged(int index)
{
    plotSpectrum(index);
}

