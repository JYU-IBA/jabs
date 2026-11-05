#ifndef PREFERENCESDIALOG_H
#define PREFERENCESDIALOG_H

#include <QDialog>
#include <QSettings>
#include <QMap>
#include "colorprofile.h"

class QComboBox;
class QPushButton;
class QToolButton;
class QLabel;

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
    void onProfileSelected(int index);
    void newProfile();
    void renameProfile();
    void deleteProfile();

signals:
    void settingsSaved();

private:
    Ui::PreferencesDialog *ui;
    QSettings settings;
    QString jibalConfFilename;
    QList<ColorProfile> customProfiles; /* Working copies, saved when dialog is accepted */
    QComboBox *profileComboBox;
    QPushButton *renameButton;
    QPushButton *deleteButton;
    QLabel *profileNote;
    QMap<QString, QToolButton *> colorButtons;
    QMap<QString, QPushButton *> systemButtons;
    void setupColors();
    void fillProfileComboBox(const QString &select);
    QString selectedProfileName() const;
    ColorProfile displayedProfile() const;
    ColorProfile *editableProfile(); /* Custom profile being shown, nullptr if built-in */
    void updateColorButtons();
    bool validProfileName(const QString &name);
};

#endif // PREFERENCESDIALOG_H
