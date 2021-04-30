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
    rfile->energy = 1.0; /* Weird default*/
    rfile->enfactors[0] = 1.0;
    rfile->sigfactors[0] = 1.0;
    rfile->distribution = R33_DIST_ENERGY;
    rfile->unit = R33_UNIT_MB;
    for(size_t i = 0; i < R33_N_NUCLEI; i++) {
        rfile->masses[i] = 1.0;
        rfile->zeds[i] = 1.0;
    }
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
    for(size_t i = 0; i < R33_N_NUCLEI; i++) {
        free(rfile->reaction_nuclei[i]);
    }
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

int r33_parse_header_content(r33_file *rfile, r33_header_type type, const char *line_split) {
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
        return 0;
    }
    if(type == R33_HEADER_SOURCE) {
        r33_string_overwrite(&rfile->source, line_split);
        return 0;
    }
    if(type == R33_HEADER_NAME) {
        r33_string_overwrite(&rfile->name, line_split);
        return 0;
    }

    if(type >= R33_HEADER_ADDRESS1 && type <= R33_HEADER_ADDRESS9) {
        r33_string_overwrite(&rfile->address[type-R33_HEADER_ADDRESS1], line_split);
        return 0;
    }
    if(type == R33_HEADER_SERIAL) {
        rfile->serial = (int)strtol(line_split, NULL, 10);
        return 0;
    }
    if(type == R33_HEADER_REACTION) {
        r33_string_overwrite(&rfile->reaction, line_split);
        return 0;
    }
    if(type == R33_HEADER_MASSES) {
        r33_values_read(line_split, rfile->masses, R33_N_NUCLEI);
        return 0;
    }
    if(type == R33_HEADER_ZEDS) {
        r33_values_read(line_split, rfile->zeds, R33_N_NUCLEI);
        return 0;
    }
    if(type == R33_HEADER_COMPOSITION) {
        r33_string_overwrite(&rfile->composition, line_split); /* TODO: parse? */
        return 0;
    }
    if(type == R33_HEADER_QVALUE) {
        r33_values_read(line_split, rfile->Qvalues, R33_N_QVALUES);
        return 0;
    }
    if(type == R33_HEADER_DISTRIBUTION) {
        if(strcmp(line_split, "Angle") == 0) {
            rfile->distribution = R33_DIST_ANGLE;
        } else if(strcmp(line_split, "Energy") == 0) {
            rfile->distribution = R33_DIST_ENERGY;
        } else {
            fprintf(stderr, "Unknown distribution \"%s\", assuming defaults.\n", line_split);
        }
        return 0;
    }
    if(type == R33_HEADER_THETA) {
        rfile->theta = strtod(line_split, NULL); /* Note: no unit conversion! */
        return 0;
    }
    if(type == R33_HEADER_ENERGY) {
        rfile->energy = strtod(line_split, NULL); /* Note: no unit conversion! */
        return 0;
    }
    if(type == R33_HEADER_SIGFACTORS) {
        r33_values_read(line_split, rfile->sigfactors, R33_N_SIGFACTORS);
        return 0;
    }
    if(type == R33_HEADER_UNITS) {
        if(strcmp(line_split, "tot") == 0) {
            rfile->unit = R33_UNIT_TOT;
        } else if(strcmp(line_split, "rr") == 0) {
            rfile->unit = R33_UNIT_RR;
        } else if(strcmp(line_split, "mb") == 0) {
            rfile->unit = R33_UNIT_MB;
        } else {
            fprintf(stderr, "Unit \"%s\" not recognized, assuming defaults.\n", line_split);
        }
        return 0;
    }
    if(type == R33_HEADER_ENFACTORS) {
        r33_values_read(line_split, rfile->enfactors, R33_N_ENFACTORS);
        return 0;
    }
    return -1;
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
                fprintf(stderr, "Line %zu: unhandled header type \"%s\"\n", lineno, line); /* This is not fatal */
                continue;
            }
            if(type == R33_HEADER_DATA) {
                state = R33_PARSE_STATE_DATA;
                continue;
            }
            if(type == R33_HEADER_NVALUES) {
                rfile->nvalues = strtol(line_split, NULL, 10);
                if(rfile->nvalues < 0) {
                    fprintf(stderr, "Negative Nvalues are not allowed.\n");
                    valid = FALSE;
                    break;
                }
                if(rfile->version == R33_VERSION_R33a && rfile->nvalues > 0) {
                    fprintf(stderr, "R33a version does not allow non-zero Nvalues.\n");
                    valid = FALSE;
                    break;
                }
                state = R33_PARSE_STATE_DATA;
                continue;
            }
            if(r33_parse_header_content(rfile, type, line_split)) { /* Headers that don't cause state changes are parsed here */
                fprintf(stderr, "Line %zu: unhandled valid header type %i (%s)\n", lineno, type, r33_header_string(type));
            }
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
            if(rfile->nvalues > 0 && rfile->n_data == (unsigned long)rfile->nvalues) { /* Reached predetermined number of values, stop. */
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
    if(rfile->unit == R33_UNIT_TOT && rfile->distribution == R33_DIST_ANGLE) {
        fprintf(stderr, "R33 file is not valid, since unit is 'tot', but distribution is angle.\n");
        valid = FALSE;
    }
    /* TODO: check that no conflicts exist (mutually exclusive or contradictory), handle some conflicts gracefully */
    r33_parse_reaction_string(rfile);
    if(!valid) {
        r33_file_free(rfile);
        return NULL;
    }
    return rfile;
}

reaction *r33_file_to_reaction(const jibal_isotope *isotopes, const r33_file *rfile) {
    const jibal_isotope *nuclei[R33_N_NUCLEI];
    for(size_t i = 0; i < R33_N_NUCLEI; i++) {
#ifdef R33_IGNORE_REACTION_STRING
        nuclei[i] = jibal_isotope_find(isotopes, NULL, r33_double_to_int(rfile->zeds[i]), r33_double_to_int(rfile->masses[i]));
#else
        nuclei[i] = jibal_isotope_find(isotopes, rfile->reaction_nuclei[i], 0, 0);
#endif
        if(!nuclei[i]) {
#ifdef R33_IGNORE_REACTION_STRING
            fprintf(stderr, "Could not parse an isotope from Z=%g, mass=%g.", rfile->zeds[i], rfile->masses[i]);
#else
            fprintf(stderr, "Could not parse an isotope from \"%s\".\n", rfile->reaction_nuclei[i]);
#endif
            return NULL;
        }
    }
    if(rfile->composition) {
        fprintf(stderr, "This program does not currently support \"Composition\" in R33 files.\n");
        return NULL;
    }
    reaction *r = malloc(sizeof(reaction));
    r->cs = JIBAL_CS_NONE;
    r->incident = nuclei[0];
    r->target = nuclei[1];
    r->product = nuclei[2];
    r->product_nucleus = nuclei[3];
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
        dest[i] = strtod(col, NULL);
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

void r33_parse_reaction_string(r33_file *rfile) {
    if(!rfile->reaction)
        return;
    char *p_str[R33_N_NUCLEI]; /* Should be 4. */
    char *s;
    p_str[0] = strdup(rfile->reaction);
    s = p_str[0];
    strsep(&s, "(");
    p_str[1] = s;
    strsep(&s, ",");
    p_str[2] = s;
    strsep(&s, ")");
    p_str[3] = s;
#ifdef DEBUG
    fprintf(stderr, "Reaction parser got %s(%s,%s)%s\n", p_str[0], p_str[1], p_str[2], p_str[3]);
#endif
    for(size_t i = 0; i < R33_N_NUCLEI; i++) {
        rfile->reaction_nuclei[i] = strdup(p_str[i]);
    }
    free(p_str[0]);
}
int r33_double_to_int(double d) {
    d = round(d);
    d += 0.5;
    return ((int) d);
}
