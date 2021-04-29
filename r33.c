/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "r33.h"

r33_file *r33_alloc() {
    r33_file *rfile = malloc(sizeof(r33_file));
    memset(rfile, 0, sizeof(r33_file));
    return rfile;
}

r33_file *r33_read(const char *filename) {
    FILE *f = fopen(filename, "r");
    if(!f)
        return NULL;
    r33_file *rfile = r33_alloc();
    char *line = NULL;
    size_t line_size = 0;
    size_t lineno = 0;
    while(getline(&line, &line_size, f) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strips all kinds of newlines! */
        char *line_split = line;
#if 0
        while((col = strsep(&line_split, " ,;:")) != NULL) {

        }
#endif
    }
    return rfile;
}

reaction *r33_file_to_reaction(const jibal_isotope *isotopes, const r33_file *rfile) {
    reaction *r = malloc(sizeof(reaction));
    r->cs = JIBAL_CS_NONE;
    r->incident = jibal_isotope_find(isotopes, NULL, rfile->zeds[0], (int)(floor(rfile->masses[0])));
    r->target = jibal_isotope_find(isotopes, NULL, rfile->zeds[1], (int)(floor(rfile->masses[1])));
    r->product = jibal_isotope_find(isotopes, NULL, rfile->zeds[2], (int)(floor(rfile->masses[2])));
    r->product_nucleus = jibal_isotope_find(isotopes, NULL, rfile->zeds[3], (int)(floor(rfile->masses[3])));
    if(rfile->filename) {
        r->filename = strdup(rfile->filename);
    } else {
        r->filename = NULL;
    }
    /* TODO: actual data */
    return r;
};
