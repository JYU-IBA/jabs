/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_R33_H
#define JABS_R33_H

#include <jibal_masses.h>
#include "reaction.h"

#define R33_INITIAL_ALLOC 128 /* Allocate memory this many data points, more memory will be allocated as is necessary */

#define R33_N_DATA_COLUMNS 4
#define R33_N_ADDRESS_FIELDS 9 /* Note that if this is ever changed, changes are required elsewhere too. */
#define R33_N_QVALUES 5 /* Up to 5 Q-values are allowed */
#define R33_N_NUCLEI 4 /* Four nuclei are involved in a reaction. No exceptions! */
#define R33_N_SIGFACTORS 2
#define R33_N_ENFACTORS 3
#define R33_UNKNOWN "Unknown"

typedef enum {
    R33_VERSION_NONE = 0,
    R33_VERSION_R33 = 1,
    R33_VERSION_R33a = 2
} r33_version;

typedef enum {
    R33_HEADER_NONE = 0,
    R33_HEADER_COMMENT = 1,
    R33_HEADER_VERSION = 2,
    R33_HEADER_SOURCE = 3,
    R33_HEADER_NAME = 4,
    R33_HEADER_ADDRESS1 = 5, /* It is imperative for the code that the numbering of these address fields is the way it is (consecutive). */
    R33_HEADER_ADDRESS2 = 6,
    R33_HEADER_ADDRESS3 = 7,
    R33_HEADER_ADDRESS4 = 8,
    R33_HEADER_ADDRESS5 = 9,
    R33_HEADER_ADDRESS6 = 10,
    R33_HEADER_ADDRESS7 = 11,
    R33_HEADER_ADDRESS8 = 12,
    R33_HEADER_ADDRESS9 = 13,
    R33_HEADER_SERIAL = 14,
    R33_HEADER_REACTION = 15,
    R33_HEADER_MASSES = 16,
    R33_HEADER_ZEDS = 17,
    R33_HEADER_COMPOSITION = 18,
    R33_HEADER_QVALUE = 19,
    R33_HEADER_DISTRIBUTION = 20,
    R33_HEADER_THETA = 21,
    R33_HEADER_ENERGY = 22,
    R33_HEADER_SIGFACTORS = 23,
    R33_HEADER_UNITS = 24,
    R33_HEADER_ENFACTORS = 25,
    R33_HEADER_NVALUES = 26,
    R33_HEADER_DATA = 27,
    R33_HEADER_ENDDATA = 28
} r33_header_type;

typedef enum {
    R33_HEADER_OPTIONAL = 0,
    R33_HEADER_REQUIRED = 1,
    R33_HEADER_MUTEX1   = 2, /* Mutually exclusive: "Theta" and "Energy". One is required, if both are present we determine which one to use (based on "Distribution") */
    R33_HEADER_MUTEX2   = 3 /* Mutually exclusive: "Nvalues" and "Data". */
} r33_optional;

typedef enum {
    R33_UNIT_MB = 0, /* mb/sr */
    R33_UNIT_RR = 1, /* ratio to Rutherford */
    R33_UNIT_TOT = 2 /* mb */
} r33_unit;

typedef enum {
    R33_DIST_ENERGY = 0,
    R33_DIST_ANGLE = 1
} r33_distribution;

typedef struct {
    r33_header_type type;
    const char *s;
    r33_optional opt;
} r33_header;

static const r33_header r33_headers[] = {
         {R33_HEADER_COMMENT,   "Comment", R33_HEADER_REQUIRED},
         {R33_HEADER_VERSION,   "Version", R33_HEADER_OPTIONAL}, /* Required for R33a */
         {R33_HEADER_SOURCE,    "Source", R33_HEADER_REQUIRED},
         {R33_HEADER_NAME,      "Name", R33_HEADER_REQUIRED}, /* Author or institution. */
         {R33_HEADER_ADDRESS1,         "Address1", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS2,         "Address2", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS3,         "Address3", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS4,         "Address4", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS5,         "Address5", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS6,         "Address6", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS7,         "Address7", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS8,         "Address8", R33_HEADER_OPTIONAL},
         {R33_HEADER_ADDRESS9,         "Address9", R33_HEADER_OPTIONAL},
         {R33_HEADER_SERIAL,           "Serial Number", R33_HEADER_REQUIRED},
         {R33_HEADER_REACTION,         "Reaction", R33_HEADER_REQUIRED},
         {R33_HEADER_MASSES,           "Masses", R33_HEADER_REQUIRED},
         {R33_HEADER_ZEDS,             "Zeds", R33_HEADER_REQUIRED},
         {R33_HEADER_COMPOSITION,      "Composition", R33_HEADER_OPTIONAL},
         {R33_HEADER_QVALUE,           "Qvalue", R33_HEADER_REQUIRED},
         {R33_HEADER_DISTRIBUTION,     "Distribution", R33_HEADER_REQUIRED},
         {R33_HEADER_THETA,            "Theta", R33_HEADER_MUTEX1},
         {R33_HEADER_ENERGY,           "Energy", R33_HEADER_MUTEX1},
         {R33_HEADER_SIGFACTORS,       "Sigfactors", R33_HEADER_OPTIONAL},
         {R33_HEADER_UNITS,            "Units", R33_HEADER_OPTIONAL},
         {R33_HEADER_ENFACTORS,        "Enfactors", R33_HEADER_OPTIONAL},
         {R33_HEADER_NVALUES,          "Nvalues", R33_HEADER_MUTEX2},
         {R33_HEADER_DATA,             "Data", R33_HEADER_MUTEX2},
         {R33_HEADER_ENDDATA,          "Enddata", R33_HEADER_OPTIONAL},
         {R33_HEADER_NONE, 0, R33_HEADER_OPTIONAL}
};

typedef double r33_data[R33_N_DATA_COLUMNS];

typedef struct r33_file {
    char *filename;
    char *comment;
    r33_version version;
    char *source;
    char *name;
    char *address[R33_N_ADDRESS_FIELDS];
    long int serial;
    char *reaction;
    char *reaction_nuclei[R33_N_NUCLEI];
    double masses[R33_N_NUCLEI]; /* Reaction masses m2(m1,m3)m4 */
    double zeds[R33_N_NUCLEI]; /* Same order as masses. Stored as doubles so we don't have to write extra code. Storing (small-ish) ints as doubles works well enough :) */
    char *composition;
    double Qvalues[R33_N_QVALUES];
    r33_distribution distribution;
    double theta;
    double energy;
    double sigfactors[R33_N_SIGFACTORS];
    r33_unit unit;
    double enfactors[R33_N_ENFACTORS];
    r33_data *data;
    size_t n_data;
    size_t n_data_alloc;
    long int nvalues;
} r33_file;

r33_file *r33_file_alloc();
void r33_file_free(r33_file *rfile);
int r33_file_data_realloc(r33_file *rfile, size_t n);
int r33_parse_header_content(r33_file *rfile, r33_header_type type, const char *line_split); /* Internal */
r33_file *r33_file_read(const char *filename);
reaction *r33_file_to_reaction(const jibal_isotope *isotopes, const r33_file *rfile);
const char *r33_header_string(r33_header_type type);
r33_header_type r33_header_type_find(const char *s);
void r33_string_append(char **dest, const char *src);
void r33_string_overwrite(char **dest, const char *src);
size_t r33_values_read(const char *str, double *dest, size_t n);
char *r33_string_upper(const char *str);
void r33_parse_reaction_string(r33_file *rfile);
int r33_double_to_int(double d);
#endif // JABS_R33_H
