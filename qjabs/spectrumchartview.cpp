#include "spectrumchartview.h"

SpectrumChartView::SpectrumChartView(QWidget *parent) : QChartView(parent),
    spectrumChart(new QChart()),
    series(new QLineSeries()),
    xAxis(new QValueAxis()),
    yAxis(new QValueAxis())
{
    setChart(spectrumChart);
    chart()->legend()->hide();
    chart()->addSeries(series);
    chart()->addAxis(xAxis, Qt::AlignBottom);
    chart()->addAxis(yAxis, Qt::AlignLeft);
    series->attachAxis(xAxis);
    series->attachAxis(yAxis);
    xAxis->setRange(0.0, 8191);
    xAxis->setLabelFormat("%.0f");
    yAxis->setRange(0.0, 1000.0);
    xAxis->setTitleText("Channel");
    yAxis->setTitleText("Counts");
    setRenderHint(QPainter::Antialiasing);
    setRubberBand(QChartView::RectangleRubberBand);
}

SpectrumChartView::~SpectrumChartView()
{
    delete series;
    delete spectrumChart;
}

void SpectrumChartView::drawDataToChart(double *data, int n)
{
    dataList.clear();
    dataList.reserve(n);
    double maxy = 0.0;
    for(int i = 0; i < n; ++i) {
        if(data[i] > maxy)
            maxy = data[i];
        dataList.append((QPointF(i, data[i])));
    }
    series->replace(dataList);
    yAxis->setRange(0.0, maxy);
    xAxis->setRange(0.0, n);
    xAxis->setTickCount(10);
    yAxis->setTickCount(10);
    xAxis->applyNiceNumbers();
    yAxis->applyNiceNumbers();
}

void SpectrumChartView::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Plus:
        chart()->zoomIn();
        break;
    case Qt::Key_Minus:
        chart()->zoomOut();
        break;
    case Qt::Key_Left:
        chart()->scroll(-10, 0);
        break;
    case Qt::Key_Right:
        chart()->scroll(10, 0);
        break;
    case Qt::Key_Up:
        chart()->scroll(0, 10);
        break;
    case Qt::Key_Down:
        chart()->scroll(0, -10);
        break;
    case Qt::Key_0:
        chart()->zoomReset();
        break;
    default:
        QGraphicsView::keyPressEvent(event);
        break;
    }
}

void SpectrumChartView::mouseMoveEvent(QMouseEvent *event)
{
    QPointF p = chart()->mapToValue(event->position());
    emit coordsChanged(p);
    QChartView::mouseMoveEvent(event);
}

void SpectrumChartView::wheelEvent(QWheelEvent *event)
{
    chart()->zoomReset();
    mFactor *= event->angleDelta().y() > 0 ? (3.0/4.0) : (4.0/3.0);
    QRectF rect = chart()->plotArea();
    //QPointF c = chart()->plotArea().center();
    QPointF c = event->position();
    rect.setWidth(mFactor*rect.width());
    rect.setHeight(mFactor*rect.height());
    rect.moveCenter(c);
    chart()->zoomIn(rect);
    QChartView::wheelEvent(event);
}
