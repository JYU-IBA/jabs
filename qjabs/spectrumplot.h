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
    void setGraphVisibility(QCPGraph *graph, bool visible);
    void clearAll();
    ~SpectrumPlot();

private slots:
    void plotxRangeChanged(const QCPRange &range);
    void plotyRangeChanged(const QCPRange &range);
    void legendClicked(QCPLegend *legend, QCPAbstractLegendItem *item, QMouseEvent *event);

private:
    double xmin, xmax, ymin, ymax;
    QFont legendFont;
};

#endif // SPECTRUMPLOT_H
