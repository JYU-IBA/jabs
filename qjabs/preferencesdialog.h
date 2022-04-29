#ifndef PREFERENCESDIALOG_H
#define PREFERENCESDIALOG_H

#include <QDialog>
#include <QSettings>

namespace Ui {
class PreferencesDialog;
}

class PreferencesDialog : public QDialog
{
    Q_OBJECT

public:
    explicit PreferencesDialog(QWidget *parent = nullptr);
    ~PreferencesDialog();

private slots:
    void on_jibalconfigurationfilechangePushButton_clicked();
    void saveSettings();
    void readSettings();

signals:
    void settingsSaved();

private:
    Ui::PreferencesDialog *ui;
    QSettings settings;
};

#endif // PREFERENCESDIALOG_H
