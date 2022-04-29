#include "preferencesdialog.h"
#include "ui_preferencesdialog.h"
#include <QFileDialog>
#include <QMessageBox>

PreferencesDialog::PreferencesDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PreferencesDialog)
{
    setAttribute(Qt::WA_DeleteOnClose);
    ui->setupUi(this);
    readSettings();
    connect(this, &QDialog::accepted, this, &PreferencesDialog::saveSettings);
    ui->editorFontComboBox->setMinimumWidth(200);
    ui->messageFontComboBox->setMinimumWidth(200);
}

PreferencesDialog::~PreferencesDialog()
{
    delete ui;
}

void PreferencesDialog::on_jibalconfigurationfilechangePushButton_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open JIBAL configuration file", "", tr("JIBAL configuration files (*.conf)"));
    if(filename.isEmpty()) {
        return;
    }
    if(!QFile::exists(filename)) {
        return;
    }
    ui->jibalconfigurationfileLineEdit->setText(filename);
}

void PreferencesDialog::saveSettings()
{
    settings.setValue("editorFontFamily", ui->editorFontComboBox->currentFont().family());
    settings.setValue("editorFontSize", ui->editorfontsizeSpinBox->value());
    settings.setValue("messageFontFamily", ui->messageFontComboBox->currentFont().family());
    settings.setValue("messageFontSize", ui->messagefontsizeSpinBox->value());
    QString newJibalConfFilename = ui->jibalconfigurationfileLineEdit->text();
    if(jibalConfFilename != newJibalConfFilename) {
        QString msg = "Please restart JaBS for the change in JIBAL configuration to take effect.";
        if(!newJibalConfFilename.isEmpty() && !QFile::exists(newJibalConfFilename)) {
            msg.append(" Note that the file you have given does not exist!");
        }
        QMessageBox::information(this, "JIBAL configuration file changed", msg);
        settings.setValue("jibalConfigurationFile", newJibalConfFilename);
    }
    emit settingsSaved();
}

void PreferencesDialog::readSettings()
{
    static QString defaultFontFamily = QFontDatabase::systemFont(QFontDatabase::FixedFont).family();
    jibalConfFilename = settings.value("jibalConfigurationFile").toString();
    ui->jibalconfigurationfileLineEdit->setText(jibalConfFilename);
    ui->editorFontComboBox->setCurrentFont(QFont(settings.value("editorFontFamily", defaultFontFamily).toString()));
    ui->editorfontsizeSpinBox->setValue(settings.value("editorFontSize", 11).toInt());
    ui->messageFontComboBox->setCurrentFont(QFont(settings.value("messageFontFamily", defaultFontFamily).toString()));
    ui->messagefontsizeSpinBox->setValue(settings.value("messageFontSize", 10).toInt());
}


