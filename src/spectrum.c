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
#include <assert.h>
#include "jabs_debug.h"
#include "generic.h"
#include "spectrum.h"
#include "message.h"
#include "histogram.h"
#include "win_compat.h"

result_spectra *result_spectra_alloc(size_t n) {
    result_spectra *spectra = malloc(sizeof(result_spectra));
    spectra->n_spectra = n;
    spectra->s = calloc(n, sizeof(result_spectrum));
    return spectra;
}

void result_spectra_free(result_spectra *spectra) {
    if(!spectra) {
        return;
    }
    for(size_t i = 0; i < spectra->n_spectra; i++) {
        result_spectrum *s = &spectra->s[i];
        if(!s) {
            continue;
        }
        jabs_histogram_free(s->histo);
        free(s->name);
    }
    free(spectra->s);
    spectra->s = NULL;
    spectra->n_spectra = 0;
}

int result_spectra_copy(result_spectra *dest, const result_spectra *src) {
    dest->n_spectra = src->n_spectra;
    dest->s = calloc(dest->n_spectra, sizeof(result_spectrum));
    for(size_t i = 0; i < src->n_spectra; i++) {
        result_spectrum_copy(&dest->s[i], &src->s[i]);
    }
    return EXIT_SUCCESS;
}

int result_spectrum_copy(result_spectrum *dest, const result_spectrum *src) {
    result_spectrum_set(dest, src->histo, src->name, src->target_isotope, src->type);
    return EXIT_SUCCESS;
}

int result_spectrum_set(result_spectrum *dest, const jabs_histogram *h, const char *name, const jibal_isotope *target_isotope, reaction_type type) {
    dest->histo = h ? jabs_histogram_clone(h) : NULL;
    dest->name = name ? strdup(name) : NULL;
    dest->target_isotope = target_isotope;
    dest->type = type;
    return EXIT_SUCCESS;
}

size_t result_spectra_n_ch(const result_spectra *spectra) {
    size_t n_ch = 0;
    for(size_t i = 0; i < spectra->n_spectra; i++) {
        if(!spectra->s[i].histo) {
            continue;
        }
        n_ch = GSL_MAX(n_ch, spectra->s[i].histo->n);
    }
    return n_ch;
}

jabs_histogram *result_spectrum_histo(const result_spectra *spectra, size_t i_spectrum) {
    if(i_spectrum < spectra->n_spectra) {
        return spectra->s[i_spectrum].histo;
    }
    return NULL;
}
jabs_histogram *result_spectra_reaction_histo(const result_spectra *spectra, size_t i_reaction) {
    return result_spectrum_histo(spectra, RESULT_SPECTRA_N_FIXED + i_reaction);
}

jabs_histogram *result_spectra_simulated_histo(const result_spectra *spectra) {
    return result_spectrum_histo(spectra, RESULT_SPECTRA_SIMULATED);
}

jabs_histogram *result_spectra_experimental_histo(const result_spectra *spectra) {
    return result_spectrum_histo(spectra, RESULT_SPECTRA_EXPERIMENTAL);
}

jabs_histogram *spectrum_read(const char *filename, size_t skip, size_t channels_max, size_t column, size_t compress) {
    char *line=NULL;
    size_t line_size=0;
    FILE *in = fopen_file_or_stream(filename, "r");
    int error = FALSE;
    if(!in)
        return NULL;
    char *delim;
    if(strncmp(jabs_file_extension_const(filename), ".csv", 4) == 0) { /* File ending is .csv, assume CSV */
        if(skip == 0)
            skip = 1; /* TODO: this assumes CSV files always have headers. In worst case we skip the first line unintentionally (is it a big deal?) */
        delim = ",";
    } else {
        delim = " \t";
    }
    jabs_histogram *h = jabs_histogram_alloc(channels_max);
    jabs_histogram_reset(h);
    h->n = 0; /* We will calculate the real number of channels based on input. */
    size_t lineno = 0;
    size_t n_columns = 0; /* Number of columns (largest in file) */
    char **columns = NULL; /* Will be (re)allocated later */
    DEBUGMSG("Reading experimental spectrum from file %s. Detector column is %zu and it can have up to %zu channels.", filename, column, channels_max);
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
                DEBUGMSG("(Re)allocating columns. New number %zu.", n_columns);
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
        if(column >= n) {
            jabs_message(MSG_ERROR, "Not enough columns in experimental spectra on line %zu. Expected %zu, got %zu.\n", lineno, column, n);
            error = TRUE;
            break;
        }
        if(column == 0) {
            ch = lineno-1;
        } else {
            ch = strtoul(columns[0], &end, 10);
            if(*end != '\0') {
                jabs_message(MSG_ERROR, "Error converting %s to channel number (must be unsigned integer). Issue on line %zu of file %s.\n", columns[0], lineno, filename);
                error = TRUE;
                break;
            }
        }
        if(ch >= channels_max) {
            jabs_message(MSG_ERROR, "Channel %zu is too large (max %zu channels). Issue on line %zu of file %s.\n", ch, channels_max, lineno, filename);
            error = TRUE;
            break;
        }
        ch /= compress;
        double y = strtod(columns[column], &end);
        if(*end != '\0') {
            jabs_message(MSG_ERROR, "Error converting col %zu \"%s\" to histogram value (floating point). Issue on line %zu of file %s.\n", column, columns[column], lineno, filename);
            error = TRUE;
            break;
        }
        h->bin[ch] += y;
        if(ch > h->n)
            h->n = ch;
    }
    if(h->n == 0 || error) {
        jabs_message(MSG_ERROR, "Experimental spectrum could not be read from file \"%s\". Read %zu lines before stopping.\n", filename, lineno);
        jabs_histogram_free(h);
        h = NULL;
    } else {
        h->n++;
        DEBUGMSG("Read %zu lines from \"%s\", probably %zu channels. Allocation of %zu channels.", lineno, filename, h->n, channels_max);
    }
    free(line);
    free(columns);
    fclose_file_or_stream(in);
    return h;
}

jabs_histogram *spectrum_read_detector(const char *filename, const detector *det) {
    if(!det)
        return NULL;
    return spectrum_read(filename, 0, det->channels, det->column, det->compress);
}
