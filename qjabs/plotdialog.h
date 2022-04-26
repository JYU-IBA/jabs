 #ifndef PLOTDIALOG_H
#define PLOTDIALOG_H

#include <QDialog>
#include <QSettings>
#include <QCloseEvent>
#include <QAbstractButton>
#include <QPushButton>

extern "C" {
#include <jibal.h>
#include <jibal_units.h>
#include <jibal_masses.h>
}

namespace Ui {
class PlotDialog;
}

class PlotDialog : public QDialog
{
    Q_OBJECT

public:
    explicit PlotDialog(QWidget *parent, struct jibal *jibal);
    ~PlotDialog();

signals:
    void closed();
    void settingsSaved();

private slots:
    void on_isotopesCheckBox_stateChanged(int arg1);

private slots:
    void buttonClicked(const QAbstractButton *button);
    void finish(int result);
    void readSettings();
    void saveSettings();

private:
    Ui::PlotDialog *ui;
    QSettings settings;
    jibal *jibal;
};

#endif // PLOTDIALOG_H
