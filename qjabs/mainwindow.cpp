#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFontDatabase>
#include <QDebug>
#include <QTextStream>

extern "C" {
#include "script.h"
#include "script_session.h"
#include "script_command.h"
void fit_iter_callback(fit_stats stats);
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    aboutString = QString("JaBS version ") + jabs_version() + "\n\n"
                       + "Using JIBAL version "+ jibal_version() + ", compiled using version " + JIBAL_VERSION + "\n\n"
                       + "Using Qt version " + qVersion() + ", compiled using version " + QT_VERSION_STR + "\n\n"
                       + "Copyright 2021 - 2022 Jaakko Julin <jaakko.julin@jyu.fi>\n";
    ui->setupUi(this);
    ui->plotSettingsGroupBox->setVisible(false); /* Will be made visible if necessary */
    QIcon icon(":/icons/icon.svg");
    QApplication::setWindowIcon(icon);
    QApplication::setApplicationName("QJaBS");
    setWindowIcon(icon);
    originalPath = QDir::currentPath();
    ui->splitter->setSizes(QList<int>() << 1 << 3 << 1);
    ui->splitter_2->setSizes(QList<int>() << 1 << 2);
#ifdef WIN32
    QFont fixedFont = QFont("Courier New");
#else
    QFont fixedFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
#endif
    fixedFont.setPointSize(12);
    ui->plainTextEdit->setFont(fixedFont);
    fixedFont.setPointSize(11);
    ui->msgTextEdit->setFont(fixedFont);
    ui->msgTextEdit->insertPlainText(aboutString);
    ui->msgTextEdit->insertPlainText(QString(COPYRIGHT_STRING));
    QString config_filename_str;
#if defined(Q_OS_OSX)
    config_filename_str = QApplication::applicationDirPath() + "/../Resources/jibal.conf";
#endif
    if(config_filename_str.isEmpty() || !QFile::exists(config_filename_str) ) {
         jibal = jibal_init(NULL);
    } else {
         jibal = jibal_init(qPrintable(config_filename_str));
    }
    ui->msgTextEdit->insertPlainText(jibal_status_string(jibal));
    session = script_session_init(jibal, NULL);
    session->fit->fit_iter_callback = &fit_iter_callback;
    ui->action_Run->setShortcutContext(Qt::ApplicationShortcut);
    highlighter = new Highlighter(ui->plainTextEdit->document());
    highlighter->setSession(session);
    ui->action_New_File->setShortcut(QKeySequence::New);
    ui->action_Open_File->setShortcut(QKeySequence::Open);
    ui->action_Save_File->setShortcut(QKeySequence::Save);
    ui->action_Save_File_as->setShortcut(QKeySequence::SaveAs);
    ui->action_Run->setShortcut(QKeySequence::Refresh);
    ui->action_Run->setShortcutVisibleInContextMenu(true);
    ui->action_Quit->setShortcut(QKeySequence::Quit);
    needsSaving = false;
    updateWindowTitle();
    firstRun = true;
    statusBar()->showMessage(QString("JaBS ") + jabs_version() + ", cwd: " +  QDir::currentPath(), 2000);
}

void MainWindow::addMessage(jabs_msg_level level, const char *msg)
{
    ui->msgTextEdit->moveCursor(QTextCursor::End);
    switch(level) {
    case MSG_ERROR:
        ui->msgTextEdit->setTextColor(Qt::red);
        break;
    case MSG_WARNING:
        ui->msgTextEdit->setTextColor(Qt::darkYellow);
        break;
    case MSG_DEBUG:
        ui->msgTextEdit->setTextColor(Qt::gray);
        break;
    default:
        ui->msgTextEdit->setTextColor(Qt::black); /* TODO: defaults from theme? */
        break;
    }
    ui->msgTextEdit->insertPlainText(msg);
    //repaint();
}

MainWindow::~MainWindow()
{
    script_session_free(session);
    jibal_free(jibal);
    delete highlighter;
    delete ui;
}

void MainWindow::openFile(const QString &filename)
{
    QFile file(filename);
    if(file.open(QFile::ReadOnly | QFile::Text)) {
        ui->plainTextEdit->setPlainText(file.readAll());
        file.close();
        resetAll();
        setFilename(filename);
        needsSaving = false;
        updateWindowTitle();
    }
}


int MainWindow::runLine(const QString &line) {
    jabs_message(MSG_INFO, stderr, "%s%s\n", PROMPT, qPrintable(line));
    int status = script_execute_command(session, qPrintable(line));
    return status;
}

void MainWindow::plotSession()
{
    if(session && session->fit && session->fit->sim) {
        ui->plotSpinBox->setMaximum(session->fit->sim->n_det);
        ui->plotSettingsGroupBox->setVisible(session->fit->sim->n_det > 1);
        plotSpectrum(ui->plotSpinBox->value() - 1);
        if(firstRun) {
            ui->widget->resetZoom();
            firstRun = false;
        }
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

void MainWindow::fitCallback(fit_stats stats)
{
    QCoreApplication::processEvents();
}

void MainWindow::on_action_Run_triggered()
{
    if(!session) {
        return;
    }
    if(firstRun) {
        resetAll();
    } else {
        ui->widget->clearAll();
        ui->msgTextEdit->clear();
        script_reset(session, 0, NULL);
    }
    QString text = ui->plainTextEdit->toPlainText();
    QTextStream stream = QTextStream(&text, QIODevice::ReadOnly);
    size_t lineno = 0;
    while(!stream.atEnd()) { /* TODO: this needs a loop to process script files. Loading script files has currently no effect (other than files getting opened) since the execution of session->files is not implemented! */
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
            return; /* or break? */
        }
        if(session->file_depth > 0) {
            if(script_process(session) < 0) {
                return;
            }
        }
        QCoreApplication::processEvents();
    }
    plotSession();
}


void MainWindow::on_plotSpinBox_valueChanged(int arg1)
{
    plotSpectrum(arg1 - 1);
}

void MainWindow::plotSpectrum(size_t i_det)
{
    ui->widget->clearAll();
    if(!session || !session->fit || !session->fit->sim) {
        qDebug() << "Nothing to plot.";
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
                if(!r_next || (r->r->type != r_next->r->type || r->r->target->Z != r_next->r->target->Z)) {
                    if(histo) {
                      ui->widget->drawDataToChart(QString("") + reaction_name(r->r) + " " + jibal_element_name(jibal->elements, r->r->target->Z), histo->bin, histo->n, QColor("Red"));
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
}


void MainWindow::on_action_New_File_triggered()
{
    ui->plainTextEdit->clear();
    resetAll();
    filename.clear();
    filebasename.clear();
    QDir::setCurrent(originalPath);
    needsSaving = false;
    updateWindowTitle();
}


void MainWindow::on_plainTextEdit_textChanged()
{
    if(!needsSaving) {
        needsSaving = true;
        updateWindowTitle();
    }
}


void MainWindow::on_action_Save_File_as_triggered()
{
    QString filename = QFileDialog::getSaveFileName(this, "Save script to file", MainWindow::filename, tr("JaBS script files (*.jbs)"));
    if(!filename.isEmpty()) {
        if(saveScriptToFile(filename)) {
            setFilename(filename);
            needsSaving = false;
            updateWindowTitle();
        }
    }
}

bool MainWindow::saveScriptToFile(const QString &filename)
{
    QFile file(filename);
    if(file.open(QFile::WriteOnly | QFile::Text)) {
        QTextStream stream(&file);
        stream << ui->plainTextEdit->toPlainText();
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
        needsSaving = false;
        updateWindowTitle();
    }
}

void MainWindow::resetAll()
{
    ui->widget->clearAll();
    ui->widget->replot();
    ui->msgTextEdit->clear();
    firstRun = true;
    script_reset(session, 0, NULL);
}


void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("QJaBS"), aboutString + "\n" + QString(COPYRIGHT_STRING).simplified());
}

void MainWindow::on_commandLineEdit_returnPressed()
{
    if(runLine(ui->commandLineEdit->text()) == SCRIPT_COMMAND_SUCCESS) {
            ui->commandLineEdit->clear();
    }
    plotSession();
}

void MainWindow::on_msgTextEdit_textChanged()
{
    ui->msgTextEdit->moveCursor(QTextCursor::End);
    ui->msgTextEdit->ensureCursorVisible();
}

