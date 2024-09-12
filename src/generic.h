/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

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
#include <jibal_units.h>

#ifdef __cplusplus
extern "C" {
#endif

#define JABS_MAX(a,b) ((a) > (b) ? (a) : (b))
#define JABS_MIN(a,b) ((a) < (b) ? (a) : (b))

char **string_to_argv(const char *str, int *argc, char **s_out); /* Turns string str into an argument vector (output array allocated or NULL), number of arguments stored in argc. All argument strings point to one string, allocated and location stored to s_out */
void argv_free(char **argv, char *s_out);
char *argv_to_string(int argc, char * const *argv);

FILE *fopen_file_or_stream(const char *filename, const char *mode); /* opens file and returns file pointer, returns NULL if fails, stderr if filename is NULL, stdout if filename is "-" */
void fclose_file_or_stream(FILE *f); /* fclose() if f is not stdout or stderr */

char *strdup_non_null(const char *s); /* if s is NULL, returns NULL, else strdup(s) */
int asprintf_append(char **ret, const char *format, ...);
int is_match(const char *candidate, const char *pattern); /* stolen from https://stackoverflow.com/questions/23457305/compare-strings-with-wildcard because of laziness. License unknown. */
char *jabs_strip_newline(char *str);
int jabs_line_is_comment(const char *line);
char *jabs_file_extension(char *filename); /* Returns pointer in filename to the LAST dot, '.' in string filename or beginning of filename if not found. */
const char *jabs_file_extension_const(const char *filename);
int jabs_unit_convert(const jibal_units *units, char type, const char *str, double *out); /* Wrapper for jibal_unit_convert with error reporting using jabs_message and sanity checker (jabs_unit_sanity_check()). Negative values are errors that should not be ignored. */
int jabs_unit_sanity_check(double value, int type); /* returns 1 if no issue was found, returns 0 and prints a warning via jabs_message() if value is suspicious, and returns -1 for really crazy stuff */
int jabs_str_to_size_t(const char *str, size_t *out);
double jabs_clock(void);
#ifdef __cplusplus
}
#endif
#endif //JABS_GENERIC_H
