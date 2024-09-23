#!/bin/bash
for file in QCustomPlot.tar.gz; do
    shafile="${file}.sha256"
    url="$(head -n 1 "$shafile")"
    if [ -f "${file}" ]; then
        echo "File $file already exists, not downloading."
    else
        echo "Downloading $file using curl from $url"
        curl "$url" -o "$file"
    fi
    if echo "$(tail -n 1 $shafile)  $file" | shasum -a 256 --check --status; then
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
