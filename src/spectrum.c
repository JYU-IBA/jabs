/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <string.h>

#include "generic.h"
#include "spectrum.h"
#include "message.h"
#include "win_compat.h"

gsl_histogram *spectrum_read(const char *filename, size_t skip, size_t channels_max, size_t column, size_t compress) {
    char *line=NULL;
    size_t line_size=0;
    FILE *in = fopen_file_or_stream(filename, "r");
    int error = FALSE;
    if(!in)
        return NULL;
    size_t l = strlen(filename);
    char *delim;
    if(l > 4 && strncmp(filename+l-4, ".csv", 4) == 0) { /* File ending is .csv. TODO: CSV support is not complete (e.g. quotes are not supported) */
        if(skip == 0)
            skip = 1; /* TODO: this assumes CSV files always have headers. In worst case we skip the first line unintentionally (is it a big deal?) */
        delim = ",";
    } else {
        delim = " \t";
    }
    gsl_histogram *h = gsl_histogram_alloc(channels_max);
    gsl_histogram_reset(h);
    h->n = 0; /* We will calculate the real number of channels based on input. */
    size_t lineno = 0;
    size_t n_columns = 0; /* Number of columns (largest in file) */
    char **columns = NULL; /* Will be (re)allocated later */
#ifdef DEBUG
fprintf(stderr, "Reading experimental spectrum from file %s. Detector column is %zu and it can have up to %zu channels.\n", filename, column, channels_max);
#endif
    while(getline(&line, &line_size, in) > 0) {
        lineno++;
        if(skip) {
            skip--;
            continue;
        }
        if(jabs_line_is_comment(line)) {
            continue;
        }
        jabs_strip_newline(line);
        char *line_split = line;
        char *col_str;
        size_t n = 0; /* Number of columns on this row */
        while ((col_str = strsep(&line_split, delim)) != NULL) {
            if(*col_str == '\0') {/* Multiple separators are treated as one */
                continue;
            }
            if(n == n_columns) {
                n_columns++;
#ifdef DEBUG
                fprintf(stderr, "(Re)allocating columns. New number %zu.\n", n_columns);
#endif
                columns = realloc(columns, n_columns*sizeof(char *));
                if(!columns) {
                    error = TRUE;
                    break;
                }
            }
            columns[n] = col_str;
            n++;
        }
        size_t ch;
        char *end;
        if(column == 0) {
            ch = lineno-1;
        } else {
            ch = strtoul(columns[0], &end, 10);
            if(*end != '\0') {
                jabs_message(MSG_ERROR, stderr, "Error converting %s to channel number (must be unsigned integer). Issue on line %zu of file %s.\n", columns[0], lineno, filename);
                error = TRUE;
                break;
            }
        }
        if(column >= n) {
            jabs_message(MSG_ERROR, stderr, "Not enough columns in experimental spectra on line %zu. Expected %zu, got %zu.\n", lineno, column, n);
            error = TRUE;
            break;
        }

        if(ch >= channels_max) {
            jabs_message(MSG_ERROR, stderr, "Channel %zu is too large (max %zu channels). Issue on line %zu of file %s.\n", ch, channels_max, lineno, filename);
            error = TRUE;
            break;
        }
        ch /= compress;
        double y = strtod(columns[column], &end);
        if(*end != '\0') {
            jabs_message(MSG_ERROR, stderr, "Error converting col %zu \"%s\" to histogram value (floating point). Issue on line %zu of file %s.\n", column, columns[column], lineno, filename);
            error = TRUE;
            break;
        }
        h->bin[ch] += y;
        if(ch > h->n)
            h->n = ch;
    }
    if(h->n == 0 || error) {
        jabs_message(MSG_ERROR, stderr, "Experimental spectrum could not be read from file \"%s\". Read %zu lines before stopping.\n", filename, lineno);
        gsl_histogram_free(h);
        h = NULL;
    } else {
        h->n++;
#ifdef DEBUG
        fprintf(stderr, "Read %zu lines from \"%s\", probably %zu channels. Allocation of %zu channels.\n", lineno, filename, h->n, channels_max);
#endif
    }
    free(line);
    free(columns);
    fclose_file_or_stream(in);
    return h;
}

gsl_histogram *spectrum_read_detector(const char *filename, const detector *det) {
    if(!det)
        return NULL;
    return spectrum_read(filename, 0, det->channels, det->column, det->compress);
}

void spectrum_set_calibration(gsl_histogram *h, const detector *det, int Z) {
    if(!h || !det)
        return;
#ifdef DEBUG
    fprintf(stderr, "Updating spectrum %p, detector %p calibration. Z = %i\n", (void *) h, (void *) det, Z);
#endif
    for(size_t i = 0; i < h->n + 1; i++) {
        h->range[i] = detector_calibrated(det, Z, i);
    }
}

double spectrum_roi(const gsl_histogram *h, size_t low, size_t high) {
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

size_t spectrum_channels_in_range(const gsl_histogram *h, size_t low, size_t high) {
    if(!h || h->n == 0)
        return 0;
    if(low >= h->n)
        return 0;
    if(high >= h->n)
        high = h->n - 1;
    return high - low + 1;
}

int spectrum_compare(const gsl_histogram *h1, const gsl_histogram *h2, size_t low, size_t high, double *out) {
    if(!h1 || !h2)
        return EXIT_FAILURE;
    if(h1->n == 0 || h2->n == 0)
        return EXIT_FAILURE;
    if(high >= h1->n || high >= h2->n)
        return EXIT_FAILURE;
    if(low >= high)
        return EXIT_FAILURE;
    double sum = 0.0;
    size_t n = 0;
    for(size_t i = low; i <= high; i++) {
        if(h2->bin[i] < 1e-4) /* TODO: arbitrary cutoff to prevent div by zero */
            continue;
        sum += pow2((h1->bin[i] - h2->bin[i]))/h2->bin[i];
        n++;
    }
    *out = sqrt(sum)/(1.0*n);
    return EXIT_SUCCESS;
}
