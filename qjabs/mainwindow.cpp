#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFontDatabase>
#include <QDebug>
#include <QTextStream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->plotSettingsGroupBox->setVisible(false); /* Will be made visible if necessary */
    QIcon icon(":/icons/icon.svg");
    QApplication::setWindowIcon(icon);
    QApplication::setApplicationName("QJaBS");
    setWindowIcon(icon);
    originalPath = QDir::currentPath();
    ui->splitter->setSizes(QList<int>() << 1 << 3);
    ui->splitter_2->setSizes(QList<int>() << 1 << 2);
#ifdef WIN32
    QFont fixedFont = QFont("Courier New");
#else
    QFont fixedFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
#endif
    fixedFont.setPointSize(12);
    ui->plainTextEdit->setFont(fixedFont);
    fixedFont.setPointSize(11);
    ui->msgPlainTextEdit->setFont(fixedFont);
    jibal = jibal_init(nullptr);
    jibal_status_print(stderr, jibal);
    session = script_session_init(jibal, NULL);
    ui->action_Run->setShortcutContext(Qt::ApplicationShortcut);
    highlighter = new Highlighter(ui->plainTextEdit->document());
    ui->action_New_File->setShortcut(QKeySequence::New);
    ui->action_Open_File->setShortcut(QKeySequence::Open);
    ui->action_Save_File->setShortcut(QKeySequence::Save);
    ui->action_Save_File_as->setShortcut(QKeySequence::SaveAs);
    ui->action_Run->setShortcut(QKeySequence::Refresh);
    ui->action_Run->setShortcutVisibleInContextMenu(true);
    ui->action_Quit->setShortcut(QKeySequence::Quit);
    needsSaving = false;
    updateWindowTitle();
}

void MainWindow::addMessage(const char *msg)
{
    ui->msgPlainTextEdit->moveCursor(QTextCursor::End);
    ui->msgPlainTextEdit->insertPlainText(msg);
    repaint();
    //ui->msgTextEdit->append(msg);
}

MainWindow::~MainWindow()
{
    script_session_free(session);
    jibal_free(jibal);
    delete highlighter;
    delete ui;
}


void MainWindow::on_action_Run_triggered()
{
    if(!session) {
        return;
    }
    script_reset(session, 0, NULL);
    QString text = ui->plainTextEdit->toPlainText();
    QTextStream stream = QTextStream(&text, QIODevice::ReadOnly);
    size_t lineno = 0;
    int status = 0;
    int exit_session = FALSE;
    while(!stream.atEnd()) {
        QString line = stream.readLine();
        qDebug() << "Processing line:" << line;
        lineno++;
        if(line.isEmpty())
                continue;
        if(line.at(0) == '#')
            continue;

        char **argv = string_to_argv(qPrintable(line));
        if(!argv) {
            fprintf(stderr, "Something went wrong in parsing arguments.\n");
            continue;
        }
        char **a = argv;
        int argc = 0;
        while(*a != NULL) {
            a++;
            argc++;
        }
#ifdef DEBUG
        for(int i = 0; i < argc; i++) {
            fprintf(stderr, "args: %i: \"%s\"\n", i, argv[i]);
        }
#endif
        if(argc) {
            int found = FALSE;
            for(const struct script_command *c = script_commands; c->name != NULL; c++) {
                if(strcmp(c->name, argv[0]) == 0) {
                    found = TRUE;
                    if(c->f == NULL) {
                        exit_session = TRUE;
                        break;
                    }
                    status = c->f(session, argc - 1, argv + 1);
                    if(c->f == script_load || c->f == script_reset) {
                        free(session->cf->vars);
                        session->cf->vars = NULL;
                        jibal_config_file_set_vars(session->cf, script_make_vars(session)); /* Loading and resetting things can reset some pointers (like fit->det, so we need to update those to the vars */
                    }
                    break;
                }
            }
            if(!found) {
                    jabs_message(MSG_ERROR, stderr, "Command \"%s\" not recognized on line %zu\n", argv[0], lineno);
                    exit_session = TRUE;
            }
            free(argv[0]);
            free(argv);
        }
        if(status) {
            qDebug() << "Error!!!! Status code" << status;
            exit_session = TRUE;
        }
        if(exit_session)
            break;
    }
    if(session && session->fit && session->fit->sim) {
        ui->plotSpinBox->setMaximum(session->fit->sim->n_det - 1);
        ui->plotSettingsGroupBox->setVisible(session->fit->sim->n_det > 1);
        plotSpectrum(ui->plotSpinBox->value());
    }
}


void MainWindow::on_plotSpinBox_valueChanged(int arg1)
{
    plotSpectrum(arg1);
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
            if(i+1 < session->fit->sim->n_reactions) {
                r_next = &ws->reactions[i+1];
            }
            if(!r)
                continue;
            if(r->last_brick == 0)
                continue;
            if(r->n_bricks > 0 && r->histo->n > 0) {
                if(histo) {
                    // TODO: check binning
                    gsl_histogram_add(histo, r->histo);
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
    ui->widget->replot();
}


void MainWindow::on_action_Open_File_triggered()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open script file", "", tr("JaBS script files (*.jbs)"));
    if(filename.isEmpty()) {
        return;
    }
    QFile file(filename);
    if(file.open(QFile::ReadOnly | QFile::Text)) {
        ui->plainTextEdit->setPlainText(file.readAll());
        file.close();
        clearPlotAndOutput();
        setFilename(filename);
        needsSaving = false;
        updateWindowTitle();
    }
}


void MainWindow::on_action_Quit_triggered()
{
    close();
}

void MainWindow::updateWindowTitle()
{
    bool newfile;
    QString title = QString("JaBS ") + jabs_version() + " - ";
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
    setWindowFilePath(filename);
    QDir::setCurrent(fi.absolutePath());
}


void MainWindow::on_action_New_File_triggered()
{
    ui->plainTextEdit->clear();
    clearPlotAndOutput();
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

void MainWindow::clearPlotAndOutput()
{
    ui->widget->clearAll();
    ui->widget->replot();
    ui->msgPlainTextEdit->clear();
}


void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("QJaBS"), QString("JaBS version ") + jabs_version() + "\n\n"
                       + "Using JIBAL version "+ jibal_version() + ", compiled using version " + JIBAL_VERSION + "\n\n"
                       + "Using Qt version " + qVersion() + ", compiled using version " + QT_VERSION_STR + "\n\n"
                       + "Copyright 2021 Jaakko Julin <jaakko.julin@jyu.fi>\n\n"
                       + QString(COPYRIGHT_STRING).simplified());
}


void MainWindow::on_autoRangeCheckBox_stateChanged(int arg1)
{
    ui->widget->setAutoRange(arg1 == Qt::Checked);
    ui->widget->replot();
}


void MainWindow::on_logScaleCheckBox_stateChanged(int arg1)
{
    ui->widget->setLogScale(arg1 == Qt::Checked);
    ui->widget->replot();
}

