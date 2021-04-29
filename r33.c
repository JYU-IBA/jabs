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

typedef enum {
    R33_PARSE_STATE_INIT = 0,
    R33_PARSE_STATE_COMMENT = 1,
    R33_PARSE_STATE_HEADERS = 2,
    R33_PARSE_STATE_DATA = 3,
    R33_PARSE_STATE_END = 4
} r33_parser_state;

r33_file *r33_file_alloc() {
    r33_file *rfile = malloc(sizeof(r33_file));
    memset(rfile, 0, sizeof(r33_file));
    return rfile;
}

void r33_file_free(r33_file *rfile) {
    free(rfile->filename);
    free(rfile->comment);
    free(rfile->source);
    free(rfile->name);
    for(size_t i = 0; i < R33_N_ADDRESS_FIELDS; i++) {
        free(rfile->address[i]);
    }
    free(rfile->reaction);
    free(rfile->composition);
    free(rfile->data);
    free(rfile);
}

int r33_file_data_realloc(r33_file *rfile, size_t n) {
    if(n == 0) {
        n = R33_INITIAL_ALLOC;
    }
    if(n < rfile->n_data) {
        fprintf(stderr, "Warning, possible data loss in R33 file.\n");
        rfile->n_data = n;
    }
    rfile->data = realloc(rfile->data, sizeof(r33_data) * n);
    if(!rfile->data)
        return -1;
    rfile->n_data_alloc = n;
    return 0;
}


r33_file *r33_file_read(const char *filename) {
    FILE *f = fopen(filename, "r");
    if(!f) {
        fprintf(stderr, "Could open file \"%s\".", filename);
        return NULL;
    }
    r33_file *rfile = r33_file_alloc();
    char *line = NULL;
    size_t line_size = 0;
    size_t lineno = 0;
    r33_parser_state state = R33_PARSE_STATE_INIT;
    int valid = TRUE; /* Unless proven otherwise. */
    while(getline(&line, &line_size, f) > 0 && valid) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strips all kinds of newlines! R33 requires CR LF sequences, but we happily accept just LF. */
        if(state == R33_PARSE_STATE_INIT || state == R33_PARSE_STATE_HEADERS) {
            char *line_split = line;
            strsep(&line_split, ":");
            if(line_split == NULL) {
                fprintf(stderr, "Line %zu is not valid. Stopping.\n", lineno);
                valid = FALSE;
                break;
            }
            if(line_split[0] == ' ') {
                line_split++; /* Skip space. This parser doesn't care if there is no space! */
            }
            r33_header_type type = r33_header_type_find(line);
            if(state == R33_PARSE_STATE_INIT) {
                if(type == R33_HEADER_COMMENT) {
                    r33_string_append(&rfile->comment, line_split);
                    state = R33_PARSE_STATE_COMMENT;
                    continue;
                } else {
                    fprintf(stderr, "Expected comment field, got something else instead.\n");
                    valid = FALSE;
                    break;
                }
            }
            /* This point is reached only when state is R33_PARSE_STATE_HEADERS */
            if(type == R33_HEADER_NONE) {
                fprintf(stderr, "Line %zu: unhandled header type \"%s\"\n", lineno, line);
                continue;
            }
            if(type == R33_HEADER_VERSION) {
                char *s = r33_string_upper(line_split);
                if(strcmp(line_split, "R33A") == 0) {
                    rfile->version = R33_VERSION_R33a;
                } else if(strcmp(line_split, "R33") == 0) {
                    rfile->version = R33_VERSION_R33;
                } else {
                    fprintf(stderr, "\"%s\" is not a valid or supported R33 file version. I'll assume you meant \"R33\".\n", line_split);
                    rfile->version = R33_VERSION_R33; /* Grudgingly */
                }
                free(s);
                continue;
            }
            if(type == R33_HEADER_SOURCE) {
                r33_string_overwrite(&rfile->source, line_split);
                continue;
            }
            if(type == R33_HEADER_NAME) {
                r33_string_overwrite(&rfile->name, line_split);
                continue;
            }
            if(type == R33_HEADER_REACTION) {
                r33_string_overwrite(&rfile->reaction, line_split);
                continue;
            }
            if(type == R33_HEADER_DATA) {
                state = R33_PARSE_STATE_DATA;
                continue;
            }
            if(type >= R33_HEADER_ADDRESS1 && type <= R33_HEADER_ADDRESS9) {
                r33_string_overwrite(&rfile->address[type-R33_HEADER_ADDRESS1], line_split);
                continue;
            }
            fprintf(stderr, "Line %zu: unhandled valid header type %i (%s)\n", lineno, type, r33_header_string(type));
        } else if(state == R33_PARSE_STATE_COMMENT) {
            if(strlen(line) == 0) {
                state = R33_PARSE_STATE_HEADERS;
                continue;
            }
            r33_string_append(&rfile->comment, line);
        } else if (state == R33_PARSE_STATE_DATA) {
            if(strncmp(line, "End", 3) == 0) { /* "EndData:" or "Enddata" or "eNdDAta", we take the hint... */
                state = R33_PARSE_STATE_END;
                continue;
            }
            if(rfile->n_data == rfile->n_data_alloc) { /* Full, allocate more memory */
                if(r33_file_data_realloc(rfile, rfile->n_data * 2)) {
                    valid = FALSE;
                    fprintf(stderr, "Issues with R33 file reallocation.\n");
                }
            }
            if(r33_values_read(line, rfile->data[rfile->n_data], R33_N_DATA_COLUMNS) != R33_N_DATA_COLUMNS) {
                fprintf(stderr, "Not enough data on line %zu.\n", lineno);
                valid = FALSE;
                break;
            }
            rfile->n_data++;
        } else if (state == R33_PARSE_STATE_END) {
            fprintf(stderr, "End reached, %zu lines read.\n", lineno); /* This will NOT be printed out if data ends with Enddata */
            break;
        } else {
            fprintf(stderr, "R33 parser state machine broken.\n");
            break;
        }
    }
    /* TODO: check that all required headers are given.. or don't, because they probably aren't. This is not a validator. */
    /* TODO: check that no conflicts exist (mutually exclusive or contradictory), handle some conflicts gracefully */
    if(!valid) {
        r33_file_free(rfile);
        return NULL;
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
    /* TODO: copy and convert data (check units etc) */
    return r;
}

const char *r33_header_string(r33_header_type type) {
    for(const r33_header *h = r33_headers; h->type != 0; h++) {
        if(type == h->type) {
            return h->s;
        }
    }
    return NULL;
}

r33_header_type r33_header_type_find(const char *s) {
    if(!s)
        return R33_HEADER_NONE;
    for(const r33_header *h = r33_headers; h->s != 0; h++) {
        if(strcmp(s, h->s) == 0) {
            return h->type;
        }
    }
    return R33_HEADER_NONE;
}

void r33_string_append(char **dest, const char *src) {
    if(!src || src[0] == '\0') /* Appending NULL or null-terminated string is NOP */
        return;
    if(!*dest) {
        *dest = strdup(src);
        return;
    }
    size_t len = strlen(*dest) + 1 + strlen(src) + 1;
    *dest = realloc(*dest, sizeof(char)*len+1);
    strcat(*dest, "\n");
    strcat(*dest, src);
}

void r33_string_overwrite(char **dest, const char *src) {
    if(!src || src[0] == '\0')
        return;
    if(*dest) {
        free(*dest);
    }
    *dest = strdup(src);
    return;
}

size_t r33_values_read(const char *str, double *dest, size_t n) {
    char *str_orig = strdup(str);
    char *split = strdup(str);
    char *col;
    size_t i = 0;
    while((col = strsep(&split, " ,;:")) != NULL) {
        if(*col == '\0') {/* Multiple separators are treated as one */
            continue;
        }
        if(i==n)
            break;
        dest[i] = strtod(split, NULL);
        i++;
    }
    free(str_orig);
    return i;
}

char *r33_string_upper(const char *str) {
    char *s_orig = strdup(str);
    for(char *s = s_orig; *(s) != '\0'; s++) {
        if(*s >= 'a' && *s < 'z') {
            *s += 'A' - 'a';
        }
    }
    return s_orig;
}
