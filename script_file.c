/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#include "script_file.h"
#include "generic.h"

script_file *script_file_open(const char *filename) {
    script_file *sfile = malloc(sizeof(script_file));
    sfile->filename = strdup_non_null(filename);
    sfile->line_size = 0;
    sfile->lineno = 0;
    sfile->line = NULL;
    sfile->f = fopen_file_or_stream(filename, "r");
    if(!sfile->f) {
        script_file_close(sfile);
        return NULL;
    }
    sfile->interactive = (sfile->f == stdin && isatty(fileno(stdin)));
    return sfile;
}

void script_file_close(script_file *sfile) {
    if(!sfile)
        return;
    if(sfile->f) {
        fclose_file_or_stream(sfile->f);
    }
#ifdef DEBUG
    fprintf(stderr, "Closed file %s\n", sfile->filename);
#endif
    free(sfile->filename);
    free(sfile->line);
    free(sfile);
}

ssize_t script_file_getline(script_file *sfile) {
    while(1) {
        ssize_t n = getline(&sfile->line, &sfile->line_size, sfile->f);
        if(n <= 0) {
            return n;
        }
        sfile->lineno++;
        sfile->line[strcspn(sfile->line, "\r\n")] = 0; /* Strip newlines */
#ifdef DEBUG
        fprintf(stderr, "File %s: line %zu: %s\n", sfile->filename, sfile->lineno, sfile->line);
#endif
        if(*sfile->line == '#') {/* Comment */
            continue;
        }
        return n;
    }
}
