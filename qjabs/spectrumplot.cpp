#include "spectrumplot.h"
extern "C" {
#include <jibal_units.h>
}

SpectrumPlot::SpectrumPlot(QWidget *parent) : QCustomPlot(parent) {
    subLayout = NULL;
    setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    setSelectRect(false);
    yAxis->setLabel("Counts");
    moveLegend(false);
    legend->setBorderPen(QPen(Qt::NoPen));
    legend->setBrush(QBrush(QColor(255,255,255,127)));
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
    legendOutsideAction = new QAction("Legend outside", this);
    legendOutsideAction->setCheckable(true);
    connect(legendOutsideAction, &QAction::triggered, this, &SpectrumPlot::setLegendOutside);
    zoomAction = new QAction("Zoom", this);
    zoomAction->setCheckable(true);
    connect(zoomAction, &QAction::triggered, this, &SpectrumPlot::setZoom);
    setZoom(false);
    energyScaleAction = new QAction("Energy axis", this);
    energyScaleAction->setCheckable(true);
    connect(energyScaleAction, &QAction::triggered, this, &SpectrumPlot::setEnergyAxis);
    setEnergyAxis(false);
    coordinatesText = new QCPItemText(this);
    coordinatesText->setText("");
    coordinatesText->setPositionAlignment(Qt::AlignLeft | Qt::AlignBottom);
    coordinatesText->setFont(QFont(font().family(), 10));
    connect(this, &QCustomPlot::mouseMove, this, &SpectrumPlot::onMouseMove);
    setAutoRange(true);
    setFocusPolicy(Qt::ClickFocus);
    installEventFilter(this);
    connect(this->selectionRect(), &QCPSelectionRect::accepted, this, &SpectrumPlot::selectionAccepted);
}

void SpectrumPlot::drawDataToChart(const QString &name, double *range, double *bin, int n, const QColor &color)
{
    addGraph();
    //graph()->setSelectable(QCP::SelectionType::stDataRange);
    graph()->setLineStyle(QCPGraph::lsStepLeft); /*  GSL histograms store the left edge of bin */
    QPen graphPen;
    graphPen.setColor(color);
    graph()->setPen(graphPen);
    graphPen.setWidth(2);
    graph()->selectionDecorator()->setPen(graphPen);
    graph()->setName(name);
    double maxy = 0.0;
    int max_i = 0;
    for(int i = 0; i < n; ++i) {
        if(bin[i] > maxy) {
            maxy = bin[i];
        }
        if(bin[i] > 0.0) {
            max_i = i;
        }
        if(energyAxis) {
            graph()->addData(range[i] / C_KEV, bin[i]);
        } else {
            graph()->addData(i, bin[i]);
        }
    }
    double xmax_new;
    max_i++; /* We want to see the right-hand edge of the highest bin, too.*/
    if(energyAxis) {
        xmax_new = range[max_i] / C_KEV; /* There are n+1 ranges, so this is always safe, even with max_i++. */
    } else {
        xmax_new = max_i;
    }
    graph()->addData(xmax_new, 0.0); /* Add an extra point (zero) */

    if(xmax_new > xmax) { /* Only increase, since multiple spectra are plotted. clearAll() has reset xmax. */
        xmax = xmax_new;
    }
    xmax += 1.0;  /* Add one channel or one keV, just so that the maximum is a little bit better visible */
    if(maxy > data_ymax) {
        data_ymax = maxy;
    }
    setVisible(true);
}

void SpectrumPlot::clearAll()
{
    clearPlottables();
    xmin = 0.0;
    xmax = 0.0;
    data_ymax = 0.0;
}

bool SpectrumPlot::eventFilter(QObject *obj, QEvent *event)
{
    if(event->type() == QEvent::KeyPress) {
        QKeyEvent *ke = static_cast<QKeyEvent *>(event);
        if(ke->key() == Qt::Key_Escape) {
            selectionRect()->cancel();
        }
        if(ke->key() == Qt::Key_Right|| ke->key() == Qt::Key_D) {
            xAxis->moveRange((xAxis->range().upper-xAxis->range().lower)*0.1);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Left|| ke->key() == Qt::Key_A) {
            xAxis->moveRange((xAxis->range().upper-xAxis->range().lower)*(-0.1));
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Up|| ke->key() == Qt::Key_W) {
            yAxis->moveRange((yAxis->range().upper-yAxis->range().lower)*0.1);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Down || ke->key() == Qt::Key_S) {
            yAxis->moveRange((yAxis->range().upper-yAxis->range().lower)*(-0.1));
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Plus || ke->key() == Qt::Key_E)  {
            yAxis->setRange(yAxis->range().lower, (yAxis->range().upper-yAxis->range().lower)*0.5 + yAxis->range().lower);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_Minus || ke->key() == Qt::Key_Q) {
            yAxis->setRange(yAxis->range().lower, (yAxis->range().upper-yAxis->range().lower)*2.0 + yAxis->range().lower);
            replot();
            return true;
        }
        if(ke->key() == Qt::Key_0 || ke->key() == Qt::Key_X) {
            resetZoom();
            return true;
        }
        if(ke->key() == Qt::Key_L) {
            setLogScale(!logscale);
            return true;
        }
        if(ke->key() == Qt::Key_1) {
            setAutoRange(!autorange);
            return true;
        }
        if(ke->key() == Qt::Key_Z) {
            setZoom(!zoom);
            return true;
        }
        if(ke->key() == Qt::Key_R) {
            startSelection();
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
    updateVerticalRange();
    logAction->setChecked(logscale);
    replot(QCustomPlot::rpQueuedReplot);
}

void SpectrumPlot::setAutoRange(bool value)
{
    autorange = value;
    axisRect()->setRangeDragAxes(xAxis, autorange?nullptr:yAxis);
    if(autorange) {
        axisRect()->setRangeZoom(Qt::Horizontal);
        updateVerticalRange();
    } else {
        axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    }
    axisRect()->setRangeZoomAxes(xAxis, autorange?nullptr:yAxis);
    autoRangeAction->setChecked(autorange);
    replot();
}

void SpectrumPlot::setZoom(bool value)
{
    zoom = value;
    setSelectRect(zoom);
    zoomAction->setChecked(zoom);
}

void SpectrumPlot::setEnergyAxis(bool value)
{
    energyScaleAction->setChecked(value);
    if(value) {
        xAxis->setLabel("Energy (keV)");
    } else {
        xAxis->setLabel("Channel");
    }
    setZoom(false); /* Abort zooming when axis changes */
    setSelectRect(false); /* Abort selection when axis changes */
    if(energyAxis != value) {
        energyAxis = value;
        emit energyAxisSet(value);
    }
}

void SpectrumPlot::updateMaxima()
{
    free(subLayout);
}

void SpectrumPlot::resetZoom()
{
    xAxis->setRange(xmin, xmax);
    updateVerticalRange(true);
    replot();
}

void SpectrumPlot::saveAsFile()
{
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    dialog.setNameFilter(tr("PDF files (*.pdf)"));
    dialog.setDefaultSuffix(".pdf");
    if (!dialog.exec())
        return;
    if(dialog.selectedFiles().isEmpty())
        return;
    QString fileName = dialog.selectedFiles().at(0);
    int w = QInputDialog::getInt(this, "Width", "Width", this->width(), 100, 10000);
    int h = QInputDialog::getInt(this, "Height", "Height", this->height(), 100, 10000);
    this->savePdf(fileName, w, h);
}

SpectrumPlot::~SpectrumPlot() {
    delete logAction;
    delete autoRangeAction;
    delete legendOutsideAction;
    delete zoomAction;
    delete energyScaleAction;
}

const QColor SpectrumPlot::getColor(int index)
{
    return colors[index % n_colors];
}

void SpectrumPlot::setSelectRect(bool value)
{
    selectRect = value;
    setSelectionRectMode(selectRect ? QCP::srmCustom : QCP::srmNone);
}

void SpectrumPlot::startSelection()
{
    setZoom(false);
    setSelectRect(true);
}

void SpectrumPlot::startZoom()
{
    zoom = true;
    setSelectRect(true);
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
        updateVerticalRange();
    }
}

void SpectrumPlot::plotyRangeChanged(const QCPRange &range)
{
    if(autorange) {
        updateVerticalRange();
        return;
    }
    QCPRange newrange(range);
    double w = range.upper - range.lower;

    double ymin = logscale?0.5:0.0;
    if(w > 2*(data_ymax - ymin)) {
        w = 2*(data_ymax - ymin);
    }
    if(range.lower < ymin) {
        newrange.lower = ymin;
        newrange.upper = ymin+w;
    }
    if(newrange.upper > data_ymax*20.0)
        newrange.upper = data_ymax*20.0;
    if(newrange.upper < 10.0)
        newrange.upper = 10.0;
    yAxis->setRange(newrange);
}

void SpectrumPlot::legendClicked(QCPLegend *legend, QCPAbstractLegendItem *item, QMouseEvent *event)
{
    QCPGraph *g = graphWithLegendItem(item);
    if(!g)
        return;
    setGraphVisibility(g, !g->visible());
    if(autorange) {
        updateVerticalRange();
    }
    replot(QCustomPlot::rpQueuedReplot);
}

void SpectrumPlot::updateVerticalRange(bool force)
{
    if(!force && !autorange)
        return;
    double ymin;
    double ymax = getVerticalMaximum();
    ymax *= (legendOutside || !legend->visible()) ? 1.01 : 1.1; /* leave more room at the top when legend is inside */
    if(logscale) {
        ymin = 0.5;
        ymax = pow(10, ceil(log10(ymax))); /* lowest power of ten larger than actual maximum, .e.g     100000 for 12345 */
    } else {
        ymin = 0.0;
        ymax = ceil(ymax/10)*10.0; /* lowest multiple of ten larger than actual maximum */
    }
    if(ymax < 10.0)
        ymax = 10.0;
    yAxis->setRange(ymin, ymax);
}

double SpectrumPlot::getVerticalMaximum()
{
    QCPRange r=xAxis->range();
    double ymax = 0.0;
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
    return ymax;
}

QStringList SpectrumPlot::visibleGraphs()
{
    QStringList result;
    for(int i = 0; i < graphCount(); ++i) {
        if(graph(i)->visible()) {
            result.append(graph(i)->name());
        }
    }
    return result;
}

void SpectrumPlot::moveLegend(bool outside)
{
    if(outside) {
        if(subLayout) {
            return;
        }
        subLayout = new QCPLayoutGrid;
        plotLayout()->addElement(0, 1, subLayout);
        subLayout->setMargins(QMargins(5, 0, 5, 5));
        subLayout->addElement(0, 0, legend);
        plotLayout()->setColumnStretchFactor(1, 0.001);
        legend->setWrap(20);
        legend->setFillOrder(QCPLegend::foRowsFirst);
    } else { /* Inside */
        axisRect()->insetLayout()->addElement(legend, Qt::AlignRight|Qt::AlignTop);
        legend->setWrap(5);
        legend->setFillOrder(QCPLegend::foColumnsFirst);
        delete subLayout; /* This is no longer needed */
        subLayout = NULL;
    }
    if(legendOutside != outside) { /* Actual change */
        legendOutside = outside;
        emit legendMoved(outside);
    }
}

void SpectrumPlot::setGraphVisibility(QCPGraph *g, bool visible) {

    updateVerticalRange();
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
        menu->addAction(legendOutsideAction);
        for(int i = 0; i < legend->itemCount(); ++i) {
            const QCPAbstractLegendItem *item = legend->item(i);
            if(item && item->selectTest(pos, false) >= 0) { // legend item

            }
        }
    } if(yAxis->selectTest(pos, false) >= 0) {
        menu->addAction(logAction);
        menu->addAction(autoRangeAction);
    } else if(xAxis->selectTest(pos, false) >= 0) {
        menu->addAction(energyScaleAction);
    } else {
      if(selectedGraphs().size() > 0) {
        QString name = selectedGraphs().at(0)->name();
        menu->addAction(QString("Hide selected graph (%1)").arg(name), this, &SpectrumPlot::hideSelectedGraph);
      }
    }
    menu->addAction("Reset zoom", this, &SpectrumPlot::resetZoom);
    menu->addAction(zoomAction);
    if(!energyAxis) { /* Selection of ranges is only implemented when X axis is channels */
        menu->addAction("Select range and copy to clipboard", this, &SpectrumPlot::startSelection);
    }
#ifdef SAVING_QCUSTOMPLOT_PDF_DOES_NOT_CRASH
    menu->addAction("Save as file...", this, &SpectrumPlot::saveAsFile);
#endif
    menu->popup(mapToGlobal(pos));

}

void SpectrumPlot::hideSelectedGraph()
{
    if (selectedGraphs().size() > 0) {
        QCPGraph *graph = selectedGraphs().at(0);
        setGraphVisibility(graph, false);
        graph->setSelection(QCPDataSelection());
        replot(QCustomPlot::rpQueuedReplot);
    }
}

void SpectrumPlot::onMouseMove(QMouseEvent *event)
{
    int x = qFloor(xAxis->pixelToCoord(event->pos().x()));
    double y = yAxis->pixelToCoord(event->pos().y());
    if(zoom) {
        coordinatesText->setText(QString("Zoom (%1, %2)").arg(x).arg(y, 6, 'f', 1));
    } else {
        coordinatesText->setText(QString("(%1, %2)").arg(x).arg(y, 6, 'f', 1));
    }
    coordinatesText->position->setCoords(QPointF(x, y));
    replot(QCustomPlot::rpQueuedReplot);
}

void SpectrumPlot::replotAll() /* Wrapper for replot() */
{
    replot(QCustomPlot::rpQueuedReplot);
}

void SpectrumPlot::selectionAccepted(const QRect &rect, QMouseEvent *event)
{
    Q_UNUSED(event)

    int left = qFloor(xAxis->pixelToCoord(qMin(rect.left(), rect.right())));
    int right = qCeil(xAxis->pixelToCoord(qMax(rect.left(), rect.right())));
    if(zoom) {
        xAxis->setRange(xAxis->pixelToCoord(rect.left()), xAxis->pixelToCoord(rect.right()));
        if(autorange) {
            updateVerticalRange();
        } else {
            yAxis->setRange(yAxis->pixelToCoord(rect.bottom()), yAxis->pixelToCoord(rect.top()));
        }
        replot(QCustomPlot::rpQueuedReplot);
    } else {
        QClipboard *clipboard = QGuiApplication::clipboard();
        QString rangeText = QString("[%1:%2]").arg(left).arg(right);
        clipboard->setText(rangeText);
        if(selectRect) {
            setSelectRect(false);
        }
        emit rangeSelected(rangeText);
    }
}

QCPGraph *SpectrumPlot::graphWithLegendItem(const QCPAbstractLegendItem *item)
{
    for (int i=0; i< graphCount(); ++i) {
        if(item == legend->itemWithPlottable(graph(i))) {
            return graph(i);
        }
    }
    return nullptr;
}

void SpectrumPlot::setLegendOutside(bool value)
{
    if(!legend->visible())
        value = false; /* When legend is not visible, override. This prevents an empty box outside the plot. */
    legendOutsideAction->blockSignals(true);
    legendOutsideAction->setChecked(value);
    legendOutsideAction->blockSignals(false);
    if(value == legendOutside) /* Legend is already where it is supposed to be */
        return;
    moveLegend(value);
    if(autorange) { /* We have different y-range scaling depending on the location of the legend */
        updateVerticalRange();
    }
    replot(QCustomPlot::rpQueuedReplot);
}

void SpectrumPlot::setLegendVisible(bool value)
{
    legend->setVisible(value);
    if(!value && legendOutside) { /* When legend is not visible, move it "inside" */
        moveLegend(false);
    }
    legendOutsideAction->setEnabled(value);
}
