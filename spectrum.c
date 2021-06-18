/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include "spectrum.h"
#include <string.h>


gsl_histogram *spectrum_read(const char *filename, const detector *det) {
    char *line=NULL;
    size_t line_size=0;
    FILE *in;
    if(!filename)
        return NULL;
    if(strcmp(filename, "-") == 0) {
        in = stdin;
    } else {
        in = fopen(filename, "r");
    }
    if(!in)
        return NULL;
    gsl_histogram *h = gsl_histogram_alloc(det->channels);
    gsl_histogram_reset(h);
    h->n = 0; /* We will calculate the real number of channels based on input. */
    size_t lineno = 0;
    size_t n_columns = 0; /* Number of columns (largest in file) */
    char **columns = NULL; /* Will be (re)allocated later */
#ifdef DEBUG
    fprintf(stderr, "Reading experimental spectrum from file %s. Detector number is %lu and it can have up to %lu channels.\n", filename, det->number, det->channels);
#endif
    while(getline(&line, &line_size, in) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strips all kinds of newlines! */
        if(strlen(line) >= 1 && *line == '#') /* Comment */
            continue;
        char *line_split = line;
        char *col;
        size_t n = 0; /* Number of columns on this row */
        while ((col = strsep(&line_split, " \t")) != NULL) {
            if(*col == '\0') /* Multiple separators are treated as one */
                continue;
            if(n == n_columns) {
                n_columns++;
#ifdef DEBUG
                fprintf(stderr, "(Re)allocating columns. New number %lu.\n", n_columns);
#endif
                columns = realloc(columns, n_columns*sizeof(char *));
                if(!columns) {
                    break;
                }
            }
            columns[n] = col;
            n++;
        }
        size_t ch;
        size_t column;
        char *end;
        if(det->number == 0) {
            ch = lineno-1;
            column = 0;
        } else {
            ch = strtoul(columns[0], &end, 10);
            if(end == columns[0]) {
                fprintf(stderr, "Error converting %s to channel number. Issue on line %lu of file %s.\n", columns[0], lineno, filename);
                break;
            }
            column = det->number;
        }
        if(n < column) {
            fprintf(stderr, "Not enough columns in experimental spectra on line %lu. Expected %lu, got %lu.\n", lineno, column, n);
        }

        if(ch >= det->channels) {
            fprintf(stderr, "Channel %lu is too large for detector (%lu channels). Issue on line %lu of file %s.\n", ch, det->channels, lineno, filename);
            break;
        }
        ch /= det->compress;
        double y = strtod(columns[column], &end);
        if(end == columns[column]) {
            fprintf(stderr, "Error converting column %lu \"%s\" to histogram value. Issue on line %lu of file %s.\n", column, columns[column], lineno, filename);
            break;
        }
        h->bin[ch] += y;
        if(ch > h->n)
            h->n = ch;
    }
    if(h->n == 0) {
        fprintf(stderr, "Experimental spectrum could be read from file \"%s\". Read %lu lines before stopping.\n", filename, lineno);
        gsl_histogram_free(h);
        h = NULL;
    } else {
        h->n++;
#ifdef DEBUG
        fprintf(stderr, "Read %lu lines from \"%s\", probably %lu channels. Allocation of %lu channels.\n",
                lineno, filename, h->n, det->channels);
#endif
    }
    free(line);
    free(columns);
    if(in != stdin) {
        fclose(in);
    }
    spectrum_set_calibration(h, det);
    return h;
}

void spectrum_set_calibration(gsl_histogram *h, const detector *det) {
    if(!h || !det)
        return;
    for(size_t i = 0; i < h->n + 1; i++) {
        h->range[i] = detector_calibrated(det, i);
    }
}

double spectrum_roi(gsl_histogram *h, size_t low, size_t high) {
    if(!h || h->n == 0)
        return 0.0;
    if(low >= h->n)
        return 0.0;
    if(high >= h->n)
        high = h->n - 1;
    double sum = 0.0;
    for(size_t i = low; i <= high; i++) {
        sum += h->bin[i];
    }
    return sum;
}

size_t spectrum_channels_in_range(gsl_histogram *h, size_t low, size_t high) {
    if(!h || h->n == 0)
        return 0;
    if(low >= h->n)
        return 0;
    if(high >= h->n)
        high = h->n - 1;
    return high - low + 1;
}
