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


gsl_histogram *read_experimental_spectrum(const char *filename, size_t n) {
    char *line=NULL;
    size_t line_size=0;
    size_t ch=0, i=0;
    double y;
    FILE *in;
    if(strcmp(filename, "-") == 0) {
        in = stdin;
    } else {
        in = fopen(filename, "r");
    }
    if(!in)
        return NULL;
    gsl_histogram *h = gsl_histogram_alloc(n);
    size_t lineno = 0;
    while(getline(&line, &line_size, in) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strips all kinds of newlines! */
        if(line_size < 2)
            continue;
        int cols  = sscanf(line, "%lu %lf", &ch, &y);
        if(cols == 2) {
            if (ch >= n)
                continue; /* Silently? */
            h->bin[ch] = y;
            i = ch;
        } else if(cols == 1 && i < n) {
            h->bin[i] = ch;
        } else {
            fprintf(stderr, "Error while reading file \"%s\": line %lu garbled: \"%s\"", filename, lineno, line);
            break;
        }
        i++;
    }
#ifdef DEBUG
    fprintf(stderr, "Read %lu lines from \"%s\", probably %lu channels. Allocation of %lu channels.\n", lineno, filename, i, n);
#endif
    h->n = i;
    free(line);
    if(in != stdin) {
        fclose(in);
    }
    return h;
}

void set_experimental_spectrum_calibration(gsl_histogram *h, const simulation *sim) {
    size_t i;
    for(i = 0; i < h->n + 1; i++) {
        h->range[i] = sim->energy_offset + i*sim->energy_slope;
    }
}
