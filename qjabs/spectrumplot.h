#ifndef SPECTRUMPLOT_H
#define SPECTRUMPLOT_H

#include <QObject>
#include <QVector>
#include "qcustomplot.h"

class SpectrumPlot : public QCustomPlot
{
    Q_OBJECT
public:
    explicit SpectrumPlot(QWidget *parent=nullptr);
    void drawDataToChart(const QString &name, double *data, int n, const QColor &color, bool rescale);
    void clearAll();
    ~SpectrumPlot();

private slots:
    void plotxRangeChanged(const QCPRange &range);
    void plotyRangeChanged(const QCPRange &range);
private:
    double xmin, xmax, ymin, ymax;
};

#endif // SPECTRUMPLOT_H
