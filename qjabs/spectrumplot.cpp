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
    legend->setWrap(5);
    plotLayout()->setRowStretchFactor(1, 0.001);
    legend->setSelectableParts(QCPLegend::spItems);
    legend->setVisible(true);
    legendFont = QFont(font());
    clearAll();
    connect(xAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this, &SpectrumPlot::plotxRangeChanged);
    connect(yAxis, qOverload<const QCPRange &>(&QCPAxis::rangeChanged), this, &SpectrumPlot::plotyRangeChanged);
    connect(this, &SpectrumPlot::legendClick, this, &SpectrumPlot::legendClicked);
    setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this, &QCustomPlot::customContextMenuRequested, this, &SpectrumPlot::contextMenuRequest);
    xAxis->setRange(0.0, 8192);
    logAction = new QAction("Log", this);
    logAction->setCheckable(true);
    connect(logAction, &QAction::triggered, this, &SpectrumPlot::setLogScale);
    setLogScale(false);
    autoRangeAction = new QAction("Autorange", this);
    autoRangeAction->setCheckable(true);
    connect(autoRangeAction, &QAction::triggered, this, &SpectrumPlot::setAutoRange);
    setAutoRange(true);
    setFocusPolicy(Qt::ClickFocus);
    installEventFilter(this);
}

void SpectrumPlot::drawDataToChart(const QString &name, double *data, int n, const QColor &color)
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
    if(autorange) {
        recalculateVerticalRange();
    }
    setVisible(true);
}

void SpectrumPlot::clearAll()
{
    clearPlottables();
    xmin = 0.0;
    xmax = 0.0;
    ymin = 0.0;
    ymax = 0.0;
    setVisible(false);
}

bool SpectrumPlot::eventFilter(QObject *obj, QEvent *event)
{
    if(event->type() == QEvent::KeyPress) {
        QKeyEvent *ke = static_cast<QKeyEvent *>(event);
        //qDebug() << "You pushed key " << ke->key();
        if(ke->key() == Qt::Key_Right) {
            xAxis->moveRange((xAxis->range().upper-xAxis->range().lower)*0.1);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Left) {
            xAxis->moveRange((xAxis->range().upper-xAxis->range().lower)*(-0.1));
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Up) {
            yAxis->moveRange((yAxis->range().upper-yAxis->range().lower)*0.1);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Down) {
            yAxis->moveRange((yAxis->range().upper-yAxis->range().lower)*(-0.1));
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Plus) {
            yAxis->scaleRange(0.5);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Minus) {
            yAxis->scaleRange(2.0);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_0) {
            xAxis->setRange(xmin, xmax);
            recalculateVerticalRange();
            yAxis->setRange(ymin, ymax);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_L) {
            setLogScale(!logscale);
            return true;
        }
        if(ke->key() == Qt::Key_A) {
            setAutoRange(!autorange);
            return true;
        }
        return false;
    } else {
        return QCustomPlot::eventFilter(obj, event);
    }
}

void SpectrumPlot::setLogScale(bool value)
{

    logscale = value;
    if(logscale) {
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        yAxis->setTicker(logTicker);
        yAxis->setScaleType(QCPAxis::stLogarithmic);
    } else {
        QSharedPointer<QCPAxisTickerFixed> linTicker(new QCPAxisTickerFixed);
        yAxis->setTicker(linTicker);
        yAxis->setScaleType(QCPAxis::stLinear);
        yAxis->ticker()->setTickStepStrategy(QCPAxisTicker::tssReadability);
        linTicker->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);
        linTicker->setTickCount(15);
    }
    recalculateVerticalRange();
    logAction->setChecked(logscale);
    replot();
}

void SpectrumPlot::setAutoRange(bool value)
{
    autorange = value;
    axisRect()->setRangeDragAxes(xAxis, autorange?nullptr:yAxis);
    if(autorange) {
        axisRect()->setRangeZoom(Qt::Horizontal);
        recalculateVerticalRange();
    } else {
        axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    }
    axisRect()->setRangeZoomAxes(xAxis, autorange?nullptr:yAxis);
    autoRangeAction->setChecked(autorange);
    replot();
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
    if(autorange) {
        recalculateVerticalRange();
    }
}

void SpectrumPlot::plotyRangeChanged(const QCPRange &range)
{
    if(autorange) {
        recalculateVerticalRange();
        return;
    }
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

void SpectrumPlot::legendClicked(QCPLegend *legend, QCPAbstractLegendItem *item, QMouseEvent *event)
{
    for (int i=0; i< graphCount(); ++i) {
        if(item == legend->itemWithPlottable(graph(i))) {
            setGraphVisibility(graph(i), !graph(i)->visible());
            replot();
            return;
        }
    }
}

void SpectrumPlot::recalculateVerticalRange()
{
    QCPRange r=xAxis->range();
    ymax=0.0;
    for(int i = 0; i < graphCount(); ++i) {
        QCPGraph *gr = graph(i);
        if(!gr->visible())
            continue;
        for(QCPGraphDataContainer::const_iterator it = gr->data()->constBegin(); it != gr->data()->constEnd(); ++it) {
            if(it->key < r.lower || it->key > r.upper)
                continue;
            if(it->value > ymax)
                ymax = it->value;
        }
    }
    ymin=logscale?0.5:0.0;
    if(logscale) {
        ymax = pow(10, ceil(log10(ymax))); /* lowest power of ten larger than actual maximum, .e.g     100000 for 12345 */
    } else {
        ymax = ceil(ymax/10)*10.0; /* lowest multiple of ten larger than actual maximum */
    }
    if(autorange) {
        yAxis->setRange(ymin, ymax);
    }
}

void SpectrumPlot::setGraphVisibility(QCPGraph *g, bool visible) {

    g->setVisible(visible);
    QCPAbstractLegendItem *item = legend->itemWithPlottable(g);
    if(item) {
        legendFont.setStrikeOut(!visible);
        item->setFont(legendFont);
    }
}

void SpectrumPlot::contextMenuRequest(const QPoint &pos)
{
    QMenu *menu = new QMenu(this);
    menu->setAttribute(Qt::WA_DeleteOnClose);
    if(legend->selectTest(pos, false) >= 0) { // legend

    } if(yAxis->selectTest(pos, false) >=0) {
        menu->addAction(logAction);
        menu->addAction(autoRangeAction);
    } else {
      if(selectedGraphs().size() > 0)
        menu->addAction("Hide selected graph", this, &SpectrumPlot::hideSelectedGraph);
    }
    menu->popup(mapToGlobal(pos));
}

void SpectrumPlot::hideSelectedGraph()
{
    if (selectedGraphs().size() > 0) {
        QCPGraph *graph = selectedGraphs().at(0);
        setGraphVisibility(graph, false);
        graph->setSelection(QCPDataSelection());
        replot();
    }
}
