#!/bin/bash
qcustomplot_version="2.1.1"
release_url="https://www.qcustomplot.com/release/${qcustomplot_version}/";
for file in QCustomPlot.tar.gz; do
    shafile="${file}.sha256"
    url="${release_url}${file}"
    echo "Downloading $file using curl from $url"
    curl "$url" -o "$file"
    if echo "$(cat $shafile)  QCustomPlot.tar.gz" | shasum -a 256 --check --status; then
        echo "$file passes SHA sum check."
        if [[ "$file" == "QCustomPlot.tar.gz" ]]; then
            tar -xf "$file"
            echo "Copying qcustomplot.cpp and .h files to qjabs directory."
            cp qcustomplot/qcustomplot.h qcustomplot/qcustomplot.cpp qjabs/
            rm -rf qcustomplot
        fi
    else 
        echo "$file fails SHA sum check (see $shafile), install shasum if not installed or verify some other way that the downloaded file is legitimate."
    fi
done
