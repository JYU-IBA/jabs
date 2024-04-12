#!/bin/bash
url="https://www.qcustomplot.com/release/2.1.1/QCustomPlot.tar.gz";
echo "Downloading qcustomplot using curl from $url"
curl "https://www.qcustomplot.com/release/2.1.1/QCustomPlot.tar.gz" -o QCustomPlot.tar.gz
if echo "$(cat QCustomPlot.tar.gz.sha256)  QCustomPlot.tar.gz" | shasum -a 256 --check --status; then
    echo "SHA sum check OK, copying qcustomplot.cpp and .h files to qjabs directory."
    cp qcustomplot/qcustomplot.h qcustomplot/qcustomplot.cpp qjabs/
else 
    echo "QCustomPlot.tar.gz fails SHA sum check (see QCustomPlot.tar.gz.sha256), install shasum if not installed or verify some other way that the downloaded file is legitimate."
fi
