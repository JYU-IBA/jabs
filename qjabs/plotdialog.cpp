#include "plotdialog.h"
#include "ui_plotdialog.h"

PlotDialog::PlotDialog(QWidget *parent, struct jibal *jibal) :
    QDialog(parent),
    ui(new Ui::PlotDialog),
    jibal(jibal)
{
    setAttribute(Qt::WA_DeleteOnClose);
    ui->setupUi(this);
    connect(ui->buttonBox, &QDialogButtonBox::clicked, this, &PlotDialog::buttonClicked);
    connect(this, &PlotDialog::accepted, this, &PlotDialog::saveSettings);
    connect(this, &PlotDialog::finished, this, &PlotDialog::finish);
    for(int Z = JIBAL_ANY_Z; Z <= jibal->config->Z_max; Z++) {
        ui->comboBox->addItem(jibal_element_name(jibal->elements, Z), QVariant(Z));
    }
    readSettings();
}

PlotDialog::~PlotDialog()
{
    delete ui;
}

void PlotDialog::buttonClicked(const QAbstractButton *button)
{
    const QPushButton *b = qobject_cast<const QPushButton *>(button);
    if(!b)
        return;
    if (b == ui->buttonBox->button(QDialogButtonBox::Apply)) {
        saveSettings();
    }
}

void PlotDialog::finish(int result)
{
    emit closed();
}

void PlotDialog::readSettings()
{
    ui->plotWhileFittingCheckBox->setChecked(settings.value("plotWhileFitting", QVariant(false)).toBool());
    ui->isotopesCheckBox->setChecked(settings.value("showIsotopes", QVariant(false)).toBool());
    int index = ui->comboBox->findData(settings.value("plotIsotopesZ").toInt());
    if(index != -1) {
        ui->comboBox->setCurrentIndex(index);
    }
    ui->showlegendCheckBox->setChecked(settings.value("showLegend", QVariant(true)).toBool());
    ui->legendOutsideCheckBox->setChecked(settings.value("legendOutside", QVariant(false)).toBool());
    ui->confidencelimitsCheckBox->setChecked(settings.value("confidenceLimits", QVariant(true)).toBool());
}

void PlotDialog::saveSettings()
{
    settings.setValue("plotWhileFitting", ui->plotWhileFittingCheckBox->isChecked());
    settings.setValue("showIsotopes", ui->isotopesCheckBox->isChecked());
    settings.setValue("plotIsotopesZ", ui->comboBox->currentData().toInt());
    settings.setValue("showLegend", ui->showlegendCheckBox->isChecked());
    settings.setValue("legendOutside", ui->legendOutsideCheckBox->isChecked());
    settings.setValue("confidenceLimits", ui->confidencelimitsCheckBox->isChecked());
    emit settingsSaved();
}

void PlotDialog::on_isotopesCheckBox_stateChanged(int arg1)
{
    ui->comboBox->setEnabled(arg1 == Qt::Checked);
}
