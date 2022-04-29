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
                       + "Copyright 2021 - 2022 Jaakko Julin <jaakko.julin@jyu.fi>\n";
    ui->setupUi(this);
    ui->widget->setVisible(false);
    ui->plotSettingsGroupBox->setVisible(false); /* Will be made visible if necessary */
    QIcon icon(":/icons/jabs.svg");
    QApplication::setWindowIcon(icon);
    setWindowIcon(icon);
    ui->msgTextBrowser->setOpenLinks(false); /* We use the connection below to handle links */
    connect(ui->msgTextBrowser, &QTextBrowser::anchorClicked, this, &MainWindow::openLink);
    originalPath = QDir::currentPath();
    ui->splitter->setSizes(QList<int>() << 1 << 3 << 1);
    ui->splitter_2->setSizes(QList<int>() << 1 << 3);
    ui->msgTextBrowser->insertPlainText(aboutString);
    ui->msgTextBrowser->insertPlainText(QString(COPYRIGHT_STRING));
    ui->msgTextBrowser->insertHtml("<p>Visit the <a href=\"https://github.com/JYU-IBA/jabs\">GitHub page</a> for latest information.</p>");
    ui->msgTextBrowser->insertPlainText("\n\n");
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
        ui->editor->setPlainText(file.readAll());
        file.close();
        setFilename(filename);
        setNeedsSaving(false);
        resetAll();
        plotSession();  /* Will hide plot etc. */
    }
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
    plotSession();
}


int MainWindow::runLine(const QString &line) {
    jabs_message(MSG_INFO, stderr, "%s%s\n", PROMPT, qPrintable(line));
    int status = script_execute_command(session, qPrintable(line));
    return status;
}

void MainWindow::plotSession(bool error)
{
    if(session && session->fit && session->fit->sim && session->fit->n_ws && !error) {
        ui->plotSpinBox->setMaximum(session->fit->n_ws);
        ui->plotSettingsGroupBox->setVisible(session->fit->n_ws > 1);
        plotSpectrum(ui->plotSpinBox->value() - 1);
        if(firstRun) {
            ui->widget->resetZoom();
            firstRun = false;
        }
        ui->widget->setVisible(true);
    } else {
        ui->plotSettingsGroupBox->setVisible(false);
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
    } else {
        ui->msgTextBrowser->insertHtml(QString("<p>JIBAL configuration file: %1</p>\n").arg(MainWindow::makeFileLink(config_filename_str)));
        ui->msgTextBrowser->insertPlainText("\n");
        jibal = jibal_init(qPrintable(config_filename_str));
    }
    if(!jibal) {
        QMessageBox::critical(this, "Error", "Could not initialize JIBAL (a NULL pointer was returned).\n");
        return -1;
    }
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
#ifdef WIN32
    static const QString defaultFontFamily = "Courier New";
#else
    static const QString defaultFontFamily = QFontDatabase::systemFont(QFontDatabase::FixedFont).family();
#endif
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
}

void MainWindow::setNeedsSaving(bool value)
{
    needsSaving = value;
    updateWindowTitle();
}

int MainWindow::fitCallback(fit_stats stats)
{
    QCoreApplication::processEvents();
    if(!fitDialog) {
        fitDialog = new FitDialog();
        fitDialog->show();
    }
    fitDialog->updateStats(stats);
    if(fitDialog->isPlotWhileFitting()) {
        plotSession(stats.error < 0);
    }
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
        ui->widget->clearAll();
        ui->msgTextBrowser->clear();
        script_reset(session, 0, NULL);
    }
    QString text = ui->editor->toPlainText();
    QTextStream stream = QTextStream(&text, QIODevice::ReadOnly);
    size_t lineno = 0;
    bool error = false;
    while(!stream.atEnd() && !error) { /* TODO: this needs a loop to process script files. Loading script files has currently no effect (other than files getting opened) since the execution of session->files is not implemented! */
        QString line = stream.readLine();
        lineno++;
#ifdef DEBUG
        qDebug() << "Processing line " << lineno << line;
#endif
        if(line.isEmpty())
                continue;
        if(line.at(0) == '#')
            continue;
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
        QCoreApplication::processEvents();
        closeFitDialog(); /* if fit dialog was opened (by callback resulting in running this line), close it immediately*/
    }
    closeFitDialog();
    enableRun(true);
    ui->msgTextBrowser->ensureCursorVisible();
    plotSession(error);
}


void MainWindow::on_plotSpinBox_valueChanged(int arg1)
{
    plotSpectrum(arg1 - 1);
}

void MainWindow::plotSpectrum(size_t i_det)
{
    ui->widget->clearAll();
    if(!session || !session->fit || !session->fit->sim) {
        ui->widget->replot();
        return;
    }
    if(i_det >= session->fit->sim->n_det)
        return;
    gsl_histogram *sim_histo = fit_data_sim(session->fit, i_det);
    gsl_histogram *exp_histo = fit_data_exp(session->fit, i_det);
    sim_workspace *ws = fit_data_ws(session->fit, i_det);
    if(exp_histo) {
        ui->widget->drawDataToChart("Experimental", exp_histo->bin, exp_histo->n, QColor("Black"));
    }
    if(sim_histo) {
        ui->widget->drawDataToChart("Simulated", sim_histo->bin, sim_histo->n, QColor("Blue"));
    }
    if(ws) {
        gsl_histogram *histo = NULL;
        int colorindex = 0;
        for(int i = 0; i < session->fit->sim->n_reactions; ++i) {
            sim_reaction *r = &ws->reactions[i];
            sim_reaction *r_next = NULL;
            if(i+1 < ws->n_reactions) {
                r_next = &ws->reactions[i+1];
            }
            if(!r)
                continue;
            if(r->last_brick == 0)
                continue;
            if(r->n_bricks > 0 && r->histo->n > 0) {
                if(histo) {
                    for(size_t i = 0; i < histo->n && i < r->histo->n; i++) {
                        /* ignores different ranges in GSL histograms */
                        histo->bin[i] += r->histo->bin[i];
                    }
                } else {
                    histo = gsl_histogram_clone(r->histo);
                }
                bool showThisIsotope = showIsotopes &&  (plotIsotopesZ == JIBAL_ANY_Z || plotIsotopesZ == r->r->target->Z);
                if(showThisIsotope || !r_next || (r->r->type != r_next->r->type || r->r->target->Z != r_next->r->target->Z)) {
                    /* Plot histo if:
                     * 1. If show isotopes is set and Z is JIBAL_ANY_Z or matches with target->Z
                     * 2. If this is the last reaction (there is no r_next), therefore the last histogram
                     * 3. Next reaction has different type (e.g. RBS vs ERD) or different Z
                     * */
                    if(histo) {
                        QString name = QString("") + reaction_name(r->r) + " ";
                        if(showThisIsotope) {
                            name.append(r->r->target->name);
                        } else {
                            name.append(jibal_element_name(jibal->elements, r->r->target->Z));
                        }
                        ui->widget->drawDataToChart(name, histo->bin, histo->n, SpectrumPlot::getColor(colorindex));
                        colorindex++;
                        ui->widget->setGraphVisibility(ui->widget->graph(), false);
                        gsl_histogram_free(histo);
                        histo = NULL;
                    }
                }
            }
        }
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
}


void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("QJaBS"), aboutString + "\n" + QString(COPYRIGHT_STRING).simplified());
}

void MainWindow::on_commandLineEdit_returnPressed()
{
    int status = runLine(ui->commandLineEdit->text());
    if(status == SCRIPT_COMMAND_SUCCESS) {
            ui->commandLineEdit->clear();
    }
    plotSession(status < 0);
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

void MainWindow::preferencesDialogFinished()
{
    qDebug() << "Finished, I suppose.";
}


void MainWindow::on_action_Plot_triggered()
{
    if(plotDialog) {
        qDebug() << "Plot dialog already exists.";
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

