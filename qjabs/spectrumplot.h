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

protected:
    bool eventFilter(QObject *obj, QEvent *event) override;

public slots:
    void setLogScale(bool value);
    void setAutoRange(bool value);
    void setLegendOutside(bool value);
    void setLegendVisible(bool value);

    void updateMaxima();
    void resetZoom();

private slots:
    void plotxRangeChanged(const QCPRange &range);
    void plotyRangeChanged(const QCPRange &range);
    void legendClicked(QCPLegend *legend, QCPAbstractLegendItem *item, QMouseEvent *event);
    void contextMenuRequest(const QPoint &pos);
    void hideSelectedGraph();
    void replotAll();

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
    QCPLayoutGrid *subLayout;
    void moveLegendOutside();
    void moveLegendInside();
    bool legendOutside;
};

#endif // SPECTRUMPLOT_H
