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
    void drawDataToChart(const QString &name, double *data, int n, const QColor &color);
    void setGraphVisibility(QCPGraph *graph, bool visible);
    void clearAll();
    void updateVerticalRange(bool force = false);
    ~SpectrumPlot();

protected:
    bool eventFilter(QObject *obj, QEvent *event) override;

public slots:
    void setLogScale(bool value);
    void setAutoRange(bool value);
    void updateMaxima();
    void resetZoom();

private slots:
    void plotxRangeChanged(const QCPRange &range);
    void plotyRangeChanged(const QCPRange &range);
    void legendClicked(QCPLegend *legend, QCPAbstractLegendItem *item, QMouseEvent *event);
    void contextMenuRequest(const QPoint &pos);
    void hideSelectedGraph();

private:
    double getVerticalMaximum();
    double xmin, xmax;
    double data_ymax;
    bool logscale;
    bool autorange;
    QFont legendFont;
    QAction *logAction;
    QAction *autoRangeAction;
};

#endif // SPECTRUMPLOT_H
