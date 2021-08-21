#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFontDatabase>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->splitter->setSizes(QList<int>() << 1 << 3);
    ui->splitter_2->setSizes(QList<int>() << 1 << 2);
    setWindowTitle(QString("JaBS ") + jabs_version());
    QFont fixedFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    fixedFont.setPointSize(14);
    ui->textEdit->setFont(fixedFont);
    jibal = jibal_init(nullptr);
    jibal_status_print(stderr, jibal);
    session = script_session_init(jibal, NULL);
}

void MainWindow::addMessage(const char *msg)
{
    ui->msgTextEdit->append(msg);
}

MainWindow::~MainWindow()
{
    script_session_free(session);
    jibal_free(jibal);
    delete ui;
}


void MainWindow::on_actionRun_triggered()
{
    if(!session) {
        return;
    }
    script_reset(session, 0, NULL);
    QString text = ui->textEdit->toPlainText();
    QTextStream stream = QTextStream(&text , QIODevice::ReadOnly);
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
                        jibal_config_file_set_vars(session->cf, script_make_vars(session)); /* Loading and resetting things can reset some pointers (like fit->det, so we need to update those to the vars */
                    }
                    break;
                }
            }
            if(!found) {
                    fprintf(stderr, "Command \"%s\" not recognized on line %zu\n", argv[0], lineno);
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
    if(!session || !session->fit || !session->fit->sim) {
        qDebug() << "Nothing to plot.";
        return;
    }
    if(i_det >= session->fit->sim->n_det)
        return;
    sim_workspace *ws = fit_data_ws(session->fit, i_det);
    if(ws) {
        ui->widget->drawDataToChart(session->fit->ws[i_det]->histo_sum->bin, session->fit->ws[i_det]->histo_sum->n);
    }
}

