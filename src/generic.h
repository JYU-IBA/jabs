/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_GENERIC_H
#define JABS_GENERIC_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>

char *strsep_with_quotes(char **stringp, const char *delim);
char **string_to_argv(const char *str, int *argc);
void argv_free(char **argv, int argc);
char *argv_to_string(int argc, char * const *argv);

FILE *fopen_file_or_stream(const char *filename, const char *mode); /* opens file and returns file pointer, returns NULL if fails, stderr if filename is NULL, stdout if filename is "-" */
void fclose_file_or_stream(FILE *f); /* fclose() if f is not stdout or stderr */

char *strdup_non_null(const char *s); /* if s is NULL, returns NULL, else strdup(s) */
int asprintf_append(char **ret, const char *format, ...);
int is_match(const char *candidate, const char *pattern); /* stolen from https://stackoverflow.com/questions/23457305/compare-strings-with-wildcard because of laziness. License unknown. */
char *jabs_strip_newline(char *str);
int jabs_line_is_comment(const char *line);
#endif //JABS_GENERIC_H
