#ifndef BPMCHARTVIEW_H
#define BPMCHARTVIEW_H

#include <QChartView>
#include <QChart>
#include <QLineSeries>
#include <QValueAxis>

class SpectrumChartView : public QChartView
{
    Q_OBJECT
public:
    SpectrumChartView(QWidget *parent = 0);
    ~SpectrumChartView();
    void drawDataToChart(double *data, int n);

protected:
    void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
    void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
    void wheelEvent(QWheelEvent *event) Q_DECL_OVERRIDE;

private:
    QChart *spectrumChart;
    QLineSeries *series;
    QValueAxis *xAxis;
    QValueAxis *yAxis;
    double mFactor = 1.0;
    QList<QPointF> dataList;
signals:
    void coordsChanged(QPointF p);
};

#endif // BPMCHARTVIEW_H
