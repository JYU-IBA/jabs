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
    ui->splitter->setSizes(QList<int>() << 1 << 3);
    ui->splitter_2->setSizes(QList<int>() << 1 << 2);
    setWindowTitle(QString("JaBS ") + jabs_version());
#ifdef WIN32
    QFont fixedFont = QFont("Courier New");
#else
    QFont fixedFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
#endif
    fixedFont.setPointSize(14);
    ui->plainTextEdit->setFont(fixedFont);
    jibal = jibal_init(nullptr);
    jibal_status_print(stderr, jibal);
    session = script_session_init(jibal, NULL);
    ui->action_Run->setShortcutContext(Qt::ApplicationShortcut);
    highlighter = new Highlighter(ui->plainTextEdit->document());
    ui->action_New_file->setShortcut(QKeySequence::New);
    ui->action_Open->setShortcut(QKeySequence::Open);
    ui->action_Run->setShortcut(QKeySequence::Refresh);
    ui->action_Run->setShortcutVisibleInContextMenu(true);
    ui->action_Quit->setShortcut(QKeySequence::Quit);
}

void MainWindow::addMessage(const char *msg)
{
    ui->msgTextEdit->append(msg);
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
                    jabs_message(MSG_ERROR, "Command \"%s\" not recognized on line %zu\n", argv[0], lineno);
                    exit_session = TRUE;
            }
            free(argv[0]);
            free(argv);
        }
        if(status) {
            qDebug() << "Error!!!! Status code" << status;
        }
        if(exit_session)
            break;
    }
    if(session && session->fit && session->fit->sim) {
        ui->plotSpinBox->setMaximum(session->fit->sim->n_det - 1);
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
        ui->widget->drawDataToChart("Experimental", exp_histo->bin, exp_histo->n, QColor("Black"), sim_histo?false:true);
    }
    if(sim_histo) {
        ui->widget->drawDataToChart("Simulated", sim_histo->bin, sim_histo->n, QColor("Blue"), true);
    }
    if(ws) {
        for(int i = 0; i < session->fit->sim->n_reactions; ++i) {
            sim_reaction *r = &ws->reactions[i];
            if(!r)
                continue;
            if(r->last_brick == 0)
                continue;
            if(r->n_bricks > 0 && r->histo->n > 0) {
                ui->widget->drawDataToChart(QString("") + reaction_name(&session->fit->sim->reactions[i]) + " " + session->fit->sim->reactions[i].target->name, r->histo->bin, r->histo->n, QColor("Red"), false);
                ui->widget->setGraphVisibility(ui->widget->graph(), false);
            }
        }
    }
    ui->widget->replot();
}


void MainWindow::on_action_Open_triggered()
{
    QString filename = QFileDialog::getOpenFileName(nullptr, "Import measurements from file", "", tr("JaBS script files (*.jbs)"));
    if(filename.isEmpty()) {
        return;
    }
    QFile file(filename);
    if (file.open(QFile::ReadOnly | QFile::Text)) {
        ui->plainTextEdit->setPlainText(file.readAll());
    }
}


void MainWindow::on_action_Quit_triggered()
{
    close();
}

