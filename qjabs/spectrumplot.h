#ifndef SPECTRUMPLOT_H
#define SPECTRUMPLOT_H

#include <QObject>
#include <QVector>
#include <QClipboard>
#include "qcustomplot.h"
#include "histogram.h"


class SpectrumPlot : public QCustomPlot
{
    Q_OBJECT
public:
    explicit SpectrumPlot(QWidget *parent=nullptr);
    QCPGraph *addGraphHisto(const QString &name, const QColor &color);
    void drawDataToGraph(QCPGraph *g, const jabs_histogram *h, bool endzero = true, bool left = true);
    void drawFilledDataToGraphs(const QString &name, const jabs_histogram *h_low, const jabs_histogram *h_high, const QColor &color);
    void setGraphVisibility(QCPGraph *graph, bool visible);
    void clearAll();
    void updateVerticalRange(bool force = false);
    ~SpectrumPlot();
    static const int n_colors = 10;
    constexpr static const QColor colors[n_colors] = {
        QColor(0xfb, 0x9a, 0x99),
        QColor(0xe3, 0x1a, 0x1c),
        QColor(0xfd, 0xbf, 0x6f),
        QColor(0xff, 0x7f, 0x00),
        QColor(0xca, 0xb2, 0xd6),
        QColor(0x6a, 0x3d, 0x9a),
        QColor(0xa6, 0xce, 0xe3),
        QColor(0x1f, 0x78, 0xb4),
        QColor(0xb2, 0xdf, 0x8a),
        QColor(0x33, 0xa0, 0x2c),
    };
    static const QColor getColor(int index);
    QStringList visibleGraphs();

protected:
    bool eventFilter(QObject *obj, QEvent *event) override;

signals:
    void rangeSelected(QString);
    void energyAxisSet(bool);
    void legendMoved(bool outside);
    void graphVisibilityChanged(); /* One of the graphs was enabled / disabled (by clicking on the legend) */

public slots:
    void setLogScale(bool value);
    void setAutoRange(bool value);
    void setZoom(bool value);
    void setLegendOutside(bool value);
    void setLegendVisible(bool value);
    void setSelectRect(bool value);
    void startSelection();
    void startZoom();
    void setEnergyAxis(bool value);

    void updateMaxima();
    void resetZoom();
    void saveAsFile();

private slots:
    void plotxRangeChanged(const QCPRange &range);
    void plotyRangeChanged(const QCPRange &range);
    void legendClicked(QCPLegend *legend, QCPAbstractLegendItem *item, QMouseEvent *event);
    void contextMenuRequest(const QPoint &pos);
    void hideSelectedGraph();
    void onMouseMove(QMouseEvent *event);
    void replotAll();
    void selectionAccepted(const QRect &rect, QMouseEvent *event);

private:
    QCPGraph *graphWithLegendItem(const QCPAbstractLegendItem *item);
    double getVerticalMaximum();
    double xmin, xmax;
    double data_ymax;
    bool logscale;
    bool autorange;
    QFont legendFont;
    QAction *logAction;
    QAction *autoRangeAction;
    QAction *legendOutsideAction;
    QAction *zoomAction;
    QAction *energyScaleAction;
    QCPLayoutGrid *subLayout;
    QCPItemText *coordinatesText;
    void moveLegend(bool outside);
    bool legendOutside;
    bool selectRect;
    int zoom;
    bool energyAxis;
};

#endif // SPECTRUMPLOT_H
