#include "spectrumplot.h"

SpectrumPlot::SpectrumPlot(QWidget *parent) : QCustomPlot(parent) {
    setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    xAxis->setLabel("Channel");
    yAxis->setLabel("Counts");
    QCPLayoutGrid *subLayout = new QCPLayoutGrid;
    plotLayout()->addElement(1, 0, subLayout);
    subLayout->setMargins(QMargins(0, 0, 5, 5));
    subLayout->addElement(0, 0, legend);
    legend->setFillOrder(QCPLegend::foColumnsFirst);
    plotLayout()->setRowStretchFactor(1, 0.001);
    legend->setSelectableParts(QCPLegend::spItems);
    legend->setVisible(true);
    clearAll();
    connect(xAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this, &SpectrumPlot::plotxRangeChanged);
    connect(yAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this, &SpectrumPlot::plotyRangeChanged);
}

void SpectrumPlot::drawDataToChart(const QString &name, double *data, int n, const QColor &color, bool rescale)
{
    addGraph();
    graph()->setLineStyle(QCPGraph::lsStepLeft); /*  GSL histograms store the left edge of bin */
    graph()->setPen(color);
    graph()->setName(name);

    double maxy = 0.0;
    for(int i = 0; i < n; ++i) {
        if(data[i] > maxy)
            maxy = data[i];
        graph()->addData(i, data[i]);
    }
    if(maxy > ymax)
        ymax = maxy;
    if(n > xmax)
        xmax = n;
    if(rescale) {
        yAxis->setRange(0.0, maxy);
        xAxis->setRange(0.0, n);
    }
    setVisible(true);
}

void SpectrumPlot::clearAll()
{
    clearPlottables();
    xmin = 0.0;
    ymin = 0.0;
    ymin = 0.0;
    ymax = 0.0;
    setVisible(false);
}

SpectrumPlot::~SpectrumPlot() {

}

void SpectrumPlot::plotxRangeChanged(const QCPRange &range)
{
    QCPRange newrange(range);
    double w = range.upper - range.lower;
    if(w > xmax - xmin) {
        w = xmax - xmin;
    }
    if(range.lower < xmin) {
        newrange.lower = xmin;
        newrange.upper = xmin + w;
    }
    if(range.upper > xmax) {
        newrange.upper = xmax;
        newrange.lower = xmax - w;
    }
    xAxis->setRange(newrange);
}

void SpectrumPlot::plotyRangeChanged(const QCPRange &range)
{
    QCPRange newrange(range);
    double w = range.upper - range.lower;
    if(w > 2*(ymax - ymin)) {
        w = 2*(ymax - ymin);
    }
    if(range.lower < ymin) {
        newrange.lower = ymin;
        newrange.upper = ymin+w;
    }
    yAxis->setRange(newrange);
}
