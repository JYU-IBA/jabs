/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

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
#include "win_compat.h"
#else
#include <unistd.h>
#endif
#include "jabs_debug.h"
#include "script_command.h"
#include "script_file.h"
#include "generic.h"

#ifdef _READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

script_file *script_file_open(const char *filename) {
    script_file *sfile = malloc(sizeof(script_file));
    sfile->filename = strdup_non_null(filename);
    sfile->line_size = 0;
    sfile->lineno = 0;
    sfile->line = NULL;
    sfile->f = fopen_file_or_stream(filename, "r");
    if(sfile->f == stdin) {
        sfile->filename = strdup("(standard input)");
    }
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
    DEBUGMSG("Closed file %s", sfile->filename);
    free(sfile->filename);
    free(sfile->line);
    free(sfile);
}

ssize_t script_file_getline(script_file *sfile) {
#ifdef _READLINE
    rl_completer_quote_characters = "\"";
    char *line = NULL;
    size_t line_size = 0;
    while(line_size == 0) {
        line = readline(PROMPT);
        if(line) {
            line_size = strlen(line);
        } else {
            line_size = 0;
        }
    }
    add_history(line);
    free(sfile->line);
    sfile->line = line;
    sfile->line_size = line_size;
    return sfile->line_size;
#else
    while(1) {
        if(sfile->interactive) {
            fputs(PROMPT, stderr);
        }
        ssize_t n = getline(&sfile->line, &sfile->line_size, sfile->f);
        if(n <= 0) {
            return n;
        }
        sfile->lineno++;
        if(jabs_line_is_comment(sfile->line)) {
            continue;
        }
        jabs_strip_newline(sfile->line);
        if(!sfile->interactive) {
            jabs_message(MSG_INFO, "%s:%zu jabs> %s\n", sfile->filename, sfile->lineno, sfile->line);
        }
        return n;
    }
#endif
}
