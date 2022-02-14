#ifndef _GNU_SOURCE
#define _GNU_SOURCE /* Needed by vasprintf() on Linux, since it is a GNU extension */
#endif
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "win_compat.h"
#include "generic.h"

/* strsep from NetBSD, modified by Jaakko Julin. Original copyright note below */

/*      $NetBSD: strsep.c,v 1.14 2003/08/07 16:43:52 agc Exp $  */

/*-
 * Copyright (c) 1990, 1993
 *      The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

char *strsep_with_quotes(char **stringp, const char *delim) {
    char *s;
    const char *spanp;
    int c, sc;
    char *tok;

    if ((s = *stringp) == NULL)
        return (NULL);

    tok = s;
    if(*s == '"') {
        *s = '\0';
        s++;
        size_t end = strcspn(s, "\""); /* Find ending quote (or end of string) */
        if(s[end] == '"') {
            s[end] = 0;
            *stringp = s+end+1;
            return tok+1;
        }
        else {
            return tok; /* We return the position of the starting quote (now a '\0') because no terminating quote was found.*/
        }
    }
    for (tok = s;;) {
        c = *s++;
        spanp = delim;
        do {
            if ((sc = *spanp++) == c) {
                if (c == 0)
                    s = NULL;
                else
                    s[-1] = 0;
                *stringp = s;
                return (tok);
            } else {
            }
        } while (sc != 0);
    }
}


/* End of code from NetBSD */

char **string_to_argv(const char *str, int *argc) { /* Returns allocated array of allocated strings, needs to be free'd (see argv_free()) */
    char *s = strdup(str);
    char *s_split = s;
    s[strcspn(s, "\r\n")] = 0;
    size_t len = strlen(s);
    size_t n = 0;
    char *col;
    while((col = strsep_with_quotes(&s_split, " \t")) != NULL) { /* TODO: parse quotation marks so that 'foo "bar baz"' becomes out[0] == "foo", out[1] == "bar baz"; */
        if(*col == '\0') {
            continue;
        }
        n++;
    }
    char **out = malloc(sizeof(char *) * (n + 1));
    out[0] = s;
    size_t pos;
    size_t i = 1;
    for(pos = 0; pos < len && i < n; pos++) {
        if(s[pos] == '\0' /*&& s[pos-1] != '\0'*/) {
            out[i] = s+pos+1;
            if(*out[i] != '\0') { /* Consecutive delimeters (turned to '\0' by strsep above) are ignored */
                i++;
            }
        }
    }
    for(i = 0; i < n; i++) { /* We should have our final strings now, let's make deep copies */
#ifdef ARGV_DEEP_COPY
        out[i] = strdup(out[i]);
#endif
#ifdef DEBUG
        fprintf(stderr, "argv[%zu] = %s\n", i, out[i]);
#endif
    }
    out[n] = NULL;
    if(argc) {
        *argc = n;
    }
#ifdef ARGV_DEEP_COPY
    free(s);
#endif
    return out;
}

void argv_free(char **argv, int argc) {
#ifdef ARGV_DEEP_COPY
    for(int i = 0; i < argc; i++) {
        free(argv[i]);
    }
#else
    (void) argc;
    free(argv[0]);
#endif
    free(argv);
}

int argc_from_argv(const char * const *argv) {
    if(!argv)
        return 0;
    const char * const *a = argv;
    int argc = 0;
    while(*a != NULL) {
#ifdef DEBUG
        fprintf(stderr, "got \"%s\" from string_to_argv\n", *a);
#endif
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
    if(ret) {
        if(len_input == 0) {
            strncpy(*ret, s, len);
        } else {
            strncat(*ret, s, len);
        }
    }
    free(s);
    return len;
}
