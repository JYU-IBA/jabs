#ifndef FITDIALOG_H
#define FITDIALOG_H

#include <QDialog>
extern "C" {
#include "fit.h"
}

namespace Ui {
class FitDialog;
}

class FitDialog : public QDialog
{
    Q_OBJECT

public:
    explicit FitDialog(QWidget *parent = nullptr);
    bool isAborted();
    bool isPlotWhileFitting();
    void updateStats(fit_stats stats);
    ~FitDialog();

private slots:
    void on_plotCheckBox_stateChanged(int arg1);

private slots:
    void on_buttonBox_rejected();

private slots:

private:
    Ui::FitDialog *ui;
    bool abort;
    bool plot;
};

#endif // FITDIALOG_H
