/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

    Some parts of this source file under different license, see below!

 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE /* Needed by vasprintf() on Linux, since it is a GNU extension */
#endif
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif
#include <stdint.h>
#include <jibal_generic.h>
#include "jabs_debug.h"
#include "win_compat.h"
#include "generic.h"
#include "defaults.h"
#include "message.h"

char **string_to_argv(const char *str, int *argc, char **s_out) { /* Returns allocated array of allocated strings, needs to be free'd (see argv_free()) */
    char *s = strdup(str);
    char *s_split = s;
    jabs_strip_newline(s);
    size_t len = strlen(s);
    size_t n = 0;
    char *col;
    while((col = jibal_strsep_with_quotes(&s_split, " \t")) != NULL) {
        if(*col == '\0') {
            continue;
        }
        n++;
    }
    char **out = malloc(sizeof(char *) * (n + 1));
    out[0] = s;
    size_t i = 0;
    for(size_t pos = 0; pos < len && i < n; pos++) {
        if(s[pos] == '\0') {
            continue;
        }
        out[i] = s + pos;
        size_t l = strlen(out[i]);
        pos += l;
        i++;
    }
    if(n == 0) {
        free(s); /* Passing the string that needs to be freed later as the first element and NULL terminating are not compatible if n == 0, unless we free the string already here. */
        s = NULL;
    }
    out[n] = NULL;
    if(argc) {
        *argc = n;
    }
    *s_out = s;
    return out;
}

void argv_free(char **argv, char *s_out) {
    free(s_out);
    free(argv);
}

int argc_from_argv(const char * const *argv) {
    if(!argv)
        return 0;
    const char * const *a = argv;
    int argc = 0;
    while(*a != NULL) {
        DEBUGMSG("got \"%s\" from string_to_argv", *a);
        a++;
        argc++;
    }
    return argc;
}


char *argv_to_string(int argc, char * const *argv) {
    if(argc < 1)
        return NULL;
    char *s = NULL;
    size_t len = 0;
    for(int i = 0; i < argc; i++) {
        len += strlen(argv[i]);
        len++;
    }
    s = malloc(sizeof(char) * len);
    *s = '\0';
    char *sp = s;
    for(int i = 0; i < argc; i++) {
        strcat(sp, argv[i]);
        sp += strlen(argv[i]);
        *sp = ' ';
        sp++;
    }
    *sp = '\0';
    return s;
}

FILE *fopen_file_or_stream(const char *filename, const char *mode) {
    FILE *f = NULL;
    if(!filename) {
        if(*mode == 'r')
            f = stdin;
        else
            f = stderr;
    } else if(strlen(filename) == 1 && *filename == '-')
        f = stdout;
    else {
        f = fopen(filename, mode);
    }
    if(!f && filename) {
        fprintf(stderr, "Can not open file \"%s\" (mode %s)\n", filename, mode);
    }
    return f;
}

void fclose_file_or_stream(FILE *f) {
    if(f == stdout || f == stderr || f == stdin)
        return;
    fclose(f);
}

char *strdup_non_null(const char *s) {
    if(s)
        return strdup(s);
    else
        return NULL;
}

int asprintf_append(char **ret, const char * restrict format, ...) {
    if(!ret)
        return -1;
    va_list argp;
    char *s = NULL;
    va_start(argp, format);
    int len = vasprintf(&s, format, argp);
    va_end(argp);
    if(len < 0)
        return len;
    size_t len_input = *ret ? strlen(*ret) : 0;
    len += len_input;
    *ret = realloc(*ret, sizeof(char) * (len + 1));
    if(*ret) {
        if(len_input == 0) {
            strncpy(*ret, s, len);
        } else {
            strncat(*ret, s, len);
        }
        (*ret)[len] = '\0'; /* Guarantee termination :) */
    } else {
        len = -1; /* Failure code */
    }
    free(s); /* Allocated by vasprintf() */
    return len;
}

int is_match(const char* candidate, const char* pattern) {
    int wildcard = 0;
    const char *last_pattern_start = 0;
    const char *last_line_start = 0;
    do
    {
        if (*pattern == *candidate)
        {
            if(wildcard == 1)
                last_line_start = candidate + 1;
            candidate++;
            pattern++;
            wildcard = 0;
        } else if (*pattern == '?') {
            if(*(candidate) == '\0') // the line is ended but char was expected
                return 0;
            if(wildcard == 1)
                last_line_start = candidate + 1;
            candidate++;
            pattern++;
            wildcard = 0;
        } else if (*pattern == '*') {
            if (*(pattern+1) == '\0') {
                return 1;
            }
            last_pattern_start = pattern;
            wildcard = 1;
            pattern++;
        } else if (wildcard) {
            if (*candidate == *pattern) {
                wildcard = 0;
                candidate++;
                pattern++;
                last_line_start = candidate + 1 ;
            } else {
                candidate++;
            }
        } else {
            if ((*pattern) == '\0' && (*candidate) == '\0')  // end of mask
                return 1; // if the line also ends here then the pattern match
            else
            {
                if (last_pattern_start != 0) {// try to restart the mask on the rest
                    pattern = last_pattern_start;
                    candidate = last_line_start;
                    last_line_start = 0;
                } else {
                    return 0;
                }
            }
        }
    } while (*candidate);

    if (*pattern == '\0') {
        return 1;
    }
    else {
        return 0;
    }
}

char *jabs_strip_newline(char *str) { /* Replaces newline with null, so everything after newline is discarded in practice */
    str[strcspn(str, "\r\n")] = 0;
    return str;
}

int jabs_line_is_comment(const char *line) {
    if(!line) {
        return FALSE;
    }
    if(*line == '#') {
        return TRUE;
    }
    if(strncmp(line, "REM ", 4) == 0) {
        return TRUE;
    }
    return FALSE;
}

char *jabs_file_extension(char *filename) {
    return (char *)jabs_file_extension_const(filename);
}

const char *jabs_file_extension_const(const char *filename) {
    const char *ext;
    for(ext = filename + strlen(filename); ext > filename && *ext != '.'; ext--) {}
    return ext;
}

int jabs_unit_convert(const jibal_units *units, char type, const char *str, double *out) {
    int error;
    DEBUGMSG("JaBS unit conversion, type '%c' (%i), string %s", type, type, str);
    if((error = jibal_unit_convert(units, type, str, out)) < 0) {
        jabs_message(MSG_ERROR, stderr, "Conversion of value \"%s\" failed: %s\n", str, jibal_unit_conversion_error_string(error));
        return error;
    }
    DEBUGMSG("Converted output to %p should have been set to %g successfully, running sanity check.", (void *)out, *out);
    return jabs_unit_sanity_check(*out, type);
}

int jabs_unit_sanity_check(double value, int type) {
    switch(type) {
        case JIBAL_UNIT_TYPE_LAYER_THICKNESS:
            if(value < 0.1 * C_TFU) {
                jabs_message(MSG_WARNING, stderr, "Layer thickness %g tfu (1e15 at./cm2) sounds quite small. Check your units.\n", value / C_TFU);
                return 0;
            }
            break;
        case JIBAL_UNIT_TYPE_ANGLE:
            if(value < -C_2PI || value> C_2PI) {
                jabs_message(MSG_WARNING, stderr, "Angle of %g degrees given. Check your units.\n", value / C_DEG);
                return 0;
            }
            break;
        case JIBAL_UNIT_TYPE_ENERGY:
            if(value> E_MAX) {
                jabs_message(MSG_ERROR, stderr, "Energy of %g MeV given. Check your units.\n", value / C_MEV);
                return -1;
            }
            break;
        case JIBAL_UNIT_TYPE_DENSITY:
            if(value< 0.025 * C_G_CM3 || value> 25 * C_G_CM3) {
                jabs_message(MSG_WARNING, stderr, "Density %g g/cm3 given. Check your units.\n", value/ C_G_CM3);
                return 0;
            }
            break;
        case JIBAL_UNIT_TYPE_SOLID_ANGLE:
            if(value < 0.001 * C_MSR || value > 4.0 * C_PI) {
                jabs_message(MSG_WARNING, stderr, "Solid angle %g msr given. Check your units.\n", value / C_MSR);
                return 0;
            }
            break;
        default:
            break;
    }
    return 1;
}


double jabs_clock() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (1.0 * clock() / CLOCKS_PER_SEC); /* TODO: replace with something that returns something comparable to omp_get_wtime() */
#endif
}

int jabs_str_to_size_t(const char *str, size_t *out) {
    if(!str) {
        jabs_message(MSG_ERROR, stderr, "Conversion error, null pointer.\n");
        return -1;
    }
    if(*str == '\0') {
        jabs_message(MSG_ERROR, stderr, "Conversion error, empty string.\n");
        return -1;
    }
    char *end;
    long val = strtol(str, &end, 10);
    DEBUGMSG("strtol reports %li\n", val);
    if(*end != '\0') { /* Not entirely valid */
        jabs_message(MSG_ERROR, stderr, "Conversion error, \"%s\" is not a valid unsigned number.\n", str);
        return -1;
    }
    if(val < 0) {
        jabs_message(MSG_ERROR, stderr, "Positive value expected, you gave \"%s\"\n", str);
        return -1;
    }
    if(val == LONG_MAX) {
        jabs_message(MSG_ERROR, stderr, "Overflow: \"%s\"\n", str);
        return -1;
    }
    if((unsigned long long)val < SIZE_MAX) { /* SIZE_MAX could theoretically be quite small (65535) */
        *out = val;
    } else {
        *out = SIZE_MAX;
    }
    return 0;
}
