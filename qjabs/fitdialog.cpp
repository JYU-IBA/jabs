#include "fitdialog.h"
#include "ui_fitdialog.h"

FitDialog::FitDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FitDialog), abort(false), plot(false)
{
    ui->setupUi(this);
}

bool FitDialog::isAborted()
{
    return abort;
}

bool FitDialog::isPlotWhileFitting()
{
    return plot;
}

void FitDialog::updateStats(fit_stats stats)
{
    ui->phaseSpinBox->setValue(stats.phase);
    ui->iterSpinBox->setValue(stats.iter);
    ui->chisqDoubleSpinBox->setValue(stats.chisq_dof);
}

FitDialog::~FitDialog()
{
    delete ui;
}

void FitDialog::on_buttonBox_rejected()
{
    abort = true;
}


void FitDialog::on_plotCheckBox_stateChanged(int arg1)
{
    plot = (arg1 == Qt::Checked);
}

