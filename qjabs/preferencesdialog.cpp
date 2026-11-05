#include "preferencesdialog.h"
#include "ui_preferencesdialog.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QColorDialog>
#include <QComboBox>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QPainter>
#include <QPushButton>
#include <QScrollArea>
#include <QToolButton>
#include <QGridLayout>
extern "C" {
#include "message.h"
}

PreferencesDialog::PreferencesDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PreferencesDialog)
{
    setAttribute(Qt::WA_DeleteOnClose);
    ui->setupUi(this);
    ui->editorFontComboBox->setMinimumWidth(200);
    ui->messageFontComboBox->setMinimumWidth(200);
    for(int i = 0; i <= MSG_ERROR; i++) {
        ui->verbosityComboBox->addItem(QString(jabs_message_level_str((jabs_msg_level) i)));
    }
    setupColors();
    readSettings();
    connect(this, &QDialog::accepted, this, &PreferencesDialog::saveSettings);
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
    settings.setValue("defaultVerbosity", ui->verbosityComboBox->currentIndex());
    ColorProfile::saveCustomProfiles(settings, customProfiles);
    settings.setValue("colorProfile", selectedProfileName());
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
    ui->verbosityComboBox->setCurrentIndex(settings.value("defaultVerbosity", JABS_DEFAULT_VERBOSITY).toInt());
}

void PreferencesDialog::setupColors()
{
    for(const QString &name : ColorProfile::customProfileNames(settings)) {
        customProfiles.append(ColorProfile::load(settings, name));
    }
    QGroupBox *box = new QGroupBox(QString("Colors"), this);
    QVBoxLayout *boxLayout = new QVBoxLayout(box);
    QHBoxLayout *profileLayout = new QHBoxLayout();
    profileLayout->addWidget(new QLabel(QString("Profile"), box));
    profileComboBox = new QComboBox(box);
    profileComboBox->setMinimumWidth(200);
    profileLayout->addWidget(profileComboBox, 1);
    QPushButton *newButton = new QPushButton(QString("New..."), box);
    newButton->setToolTip(QString("New profile, a copy of the profile currently shown"));
    renameButton = new QPushButton(QString("Rename..."), box);
    deleteButton = new QPushButton(QString("Delete"), box);
    profileLayout->addWidget(newButton);
    profileLayout->addWidget(renameButton);
    profileLayout->addWidget(deleteButton);
    boxLayout->addLayout(profileLayout);
    profileNote = new QLabel(box);
    profileNote->setWordWrap(true);
    boxLayout->addWidget(profileNote);

    QWidget *grid = new QWidget();
    QGridLayout *gridLayout = new QGridLayout(grid);
    int row = 0;
    QString group;
    for(const ColorRole &role : ColorProfile::roles()) {
        if(role.group != group) {
            group = role.group;
            QLabel *groupLabel = new QLabel(QString("<b>%1</b>").arg(group), grid);
            gridLayout->addWidget(groupLabel, row++, 0, 1, 3);
        }
        gridLayout->addWidget(new QLabel(role.label, grid), row, 0);
        QToolButton *button = new QToolButton(grid);
        button->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
        button->setMinimumWidth(120);
        connect(button, &QToolButton::clicked, this, [this, role]() {
            ColorProfile *p = editableProfile();
            if(!p) {
                return;
            }
            QColor initial = p->color(role.key);
            if(!initial.isValid()) {
                initial = palette().color(QPalette::Base);
            }
            QColor c = QColorDialog::getColor(initial, this, QString("%1: %2").arg(role.group, role.label), QColorDialog::ShowAlphaChannel);
            if(c.isValid()) {
                p->setColor(role.key, c);
                updateColorButtons();
            }
        });
        gridLayout->addWidget(button, row, 1);
        colorButtons[role.key] = button;
        if(role.systemAllowed) {
            QPushButton *systemButton = new QPushButton(QString("System"), grid);
            systemButton->setToolTip(QString("Use color of the operating system"));
            connect(systemButton, &QPushButton::clicked, this, [this, role]() {
                ColorProfile *p = editableProfile();
                if(p) {
                    p->setColor(role.key, QColor());
                    updateColorButtons();
                }
            });
            gridLayout->addWidget(systemButton, row, 2);
            systemButtons[role.key] = systemButton;
        }
        row++;
    }
    QScrollArea *scroll = new QScrollArea(box);
    scroll->setWidget(grid);
    scroll->setWidgetResizable(true);
    scroll->setMinimumHeight(250);
    boxLayout->addWidget(scroll, 1);

    ui->gridLayout_2->removeWidget(ui->buttonBox);
    ui->gridLayout_2->addWidget(box, 1, 0);
    ui->gridLayout_2->addWidget(ui->buttonBox, 2, 0);
    ui->gridLayout_2->setRowStretch(1, 1);
    resize(width(), 650);

    connect(profileComboBox, qOverload<int>(&QComboBox::currentIndexChanged), this, &PreferencesDialog::onProfileSelected);
    connect(newButton, &QPushButton::clicked, this, &PreferencesDialog::newProfile);
    connect(renameButton, &QPushButton::clicked, this, &PreferencesDialog::renameProfile);
    connect(deleteButton, &QPushButton::clicked, this, &PreferencesDialog::deleteProfile);
    fillProfileComboBox(settings.value("colorProfile", ColorProfile::automaticName).toString());
}

void PreferencesDialog::fillProfileComboBox(const QString &select)
{
    profileComboBox->blockSignals(true);
    profileComboBox->clear();
    profileComboBox->addItem(QString("Automatic (follows system light/dark mode)"), ColorProfile::automaticName);
    profileComboBox->addItem(ColorProfile::lightName, ColorProfile::lightName);
    profileComboBox->addItem(ColorProfile::darkName, ColorProfile::darkName);
    for(const ColorProfile &p : customProfiles) {
        profileComboBox->addItem(p.name, p.name);
    }
    int index = profileComboBox->findData(select);
    profileComboBox->setCurrentIndex(index >= 0 ? index : 0);
    profileComboBox->blockSignals(false);
    onProfileSelected(profileComboBox->currentIndex());
}

QString PreferencesDialog::selectedProfileName() const
{
    return profileComboBox->currentData().toString();
}

ColorProfile PreferencesDialog::displayedProfile() const
{
    QString name = selectedProfileName();
    for(const ColorProfile &p : customProfiles) {
        if(p.name == name) {
            return p;
        }
    }
    if(name == ColorProfile::lightName) {
        return ColorProfile::light();
    }
    if(name == ColorProfile::darkName) {
        return ColorProfile::dark();
    }
    return ColorProfile::automatic();
}

ColorProfile *PreferencesDialog::editableProfile()
{
    QString name = selectedProfileName();
    for(ColorProfile &p : customProfiles) {
        if(p.name == name) {
            return &p;
        }
    }
    return nullptr;
}

void PreferencesDialog::onProfileSelected(int index)
{
    Q_UNUSED(index);
    bool editable = (editableProfile() != nullptr);
    renameButton->setEnabled(editable);
    deleteButton->setEnabled(editable);
    if(editable) {
        profileNote->setText(QString("Click a color to change it. \"System\" uses the color of the operating system."));
    } else if(selectedProfileName() == ColorProfile::automaticName) {
        profileNote->setText(QString("Built-in profile, currently %1 colors with backgrounds and text of the operating system. Create a new profile to edit colors.").arg(ColorProfile::systemIsDark() ? ColorProfile::darkName : ColorProfile::lightName));
    } else {
        profileNote->setText(QString("Built-in profile, the whole application is %1 regardless of the operating system. Create a new profile to edit colors.").arg(selectedProfileName() == ColorProfile::darkName ? "dark" : "light"));
    }
    updateColorButtons();
}

void PreferencesDialog::updateColorButtons()
{
    const ColorProfile p = displayedProfile();
    const bool editable = (editableProfile() != nullptr);
    for(const ColorRole &role : ColorProfile::roles()) {
        QToolButton *button = colorButtons.value(role.key);
        if(!button) {
            continue;
        }
        QColor c = p.color(role.key);
        QPixmap swatch(24, 16);
        if(c.isValid()) {
            swatch.fill(Qt::white);
            QPainter painter(&swatch);
            painter.fillRect(0, 0, 12, 8, Qt::lightGray); /* Checkerboard shows transparency */
            painter.fillRect(12, 8, 12, 8, Qt::lightGray);
            painter.fillRect(swatch.rect(), c);
            painter.setPen(palette().color(QPalette::Mid));
            painter.drawRect(swatch.rect().adjusted(0, 0, -1, -1));
            painter.end();
            button->setText(c.alpha() == 255 ? c.name(QColor::HexRgb) : c.name(QColor::HexArgb));
        } else {
            swatch.fill(Qt::transparent);
            button->setText(QString("System"));
        }
        button->setIcon(QIcon(swatch));
        button->setEnabled(editable);
        QPushButton *systemButton = systemButtons.value(role.key);
        if(systemButton) {
            systemButton->setEnabled(editable && c.isValid());
        }
    }
}

bool PreferencesDialog::validProfileName(const QString &name)
{
    if(name.isEmpty() || name.contains('/') || name.contains('\\')) {
        QMessageBox::warning(this, "Invalid name", "Profile name can not be empty or contain slashes.");
        return false;
    }
    if(ColorProfile::isBuiltinName(name)) {
        QMessageBox::warning(this, "Invalid name", QString("\"%1\" is a built-in profile.").arg(name));
        return false;
    }
    for(const ColorProfile &p : customProfiles) {
        if(p.name == name) {
            QMessageBox::warning(this, "Invalid name", QString("Profile \"%1\" already exists.").arg(name));
            return false;
        }
    }
    return true;
}

void PreferencesDialog::newProfile()
{
    ColorProfile p = displayedProfile();
    bool ok;
    QString name = QInputDialog::getText(this, "New color profile", QString("Name of the new profile (a copy of %1):").arg(p.name), QLineEdit::Normal, QString("%1 (custom)").arg(p.name), &ok).trimmed();
    if(!ok || !validProfileName(name)) {
        return;
    }
    p.name = name;
    p.builtin = false;
    customProfiles.append(p);
    fillProfileComboBox(name);
}

void PreferencesDialog::renameProfile()
{
    ColorProfile *p = editableProfile();
    if(!p) {
        return;
    }
    bool ok;
    QString name = QInputDialog::getText(this, "Rename color profile", "New name:", QLineEdit::Normal, p->name, &ok).trimmed();
    if(!ok || name == p->name || !validProfileName(name)) {
        return;
    }
    p->name = name;
    fillProfileComboBox(name);
}

void PreferencesDialog::deleteProfile()
{
    ColorProfile *p = editableProfile();
    if(!p) {
        return;
    }
    if(QMessageBox::question(this, "Delete color profile", QString("Delete profile \"%1\"? The change is saved when you press OK.").arg(p->name)) != QMessageBox::Yes) {
        return;
    }
    QString name = p->name;
    customProfiles.removeIf([&name](const ColorProfile &cp) {return cp.name == name;});
    fillProfileComboBox(ColorProfile::automaticName);
}
