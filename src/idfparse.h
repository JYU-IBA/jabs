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
#ifndef IDFPARSER_H
#define IDFPARSER_H
#include <libxml/tree.h>
#include <jibal_units.h>

#define IDF_BUF_MSG_MAX (1024) /* Maximum length of input string to idf_output_printf. Longer strings will be truncated. */
#define IDF_BUF_SIZE_INITIAL (8*IDF_BUF_MSG_MAX)
#define JABS_FILE_SUFFIX ".jbs"
#define EXP_SPECTRUM_FILE_SUFFIX ".dat"
#define SAVE_SPECTRUM_FILE_SUFFIX ".csv"

typedef enum idf_error {
    IDF2JBS_SUCCESS = 0,
    IDF2JBS_FAILURE = -1,
    IDF2JBS_FAILURE_COULD_NOT_READ = -2,
    IDF2JBS_FAILURE_NOT_IDF_FILE = -3,
    IDF2JBS_FAILURE_WRONG_EXTENSION = -4,
    IDF2JBS_FAILURE_COULD_NOT_WRITE_OUTPUT = -5,
    IDF2JBS_FAILURE_NO_SAMPLES_DEFINED = -6
} idf_error;

typedef struct idf_unit {
    const char *unit;
    double factor;
} idf_unit;

static const idf_unit idf_units[] = {
        {"fraction", 1.0},
        {"degree", C_DEG},
        {"radian", 1.0},
        {"1e15at/cm2", C_TFU},
        {"#particles", 1.0},
        {"amu", C_U},
        {"keV", C_KEV},
        {"keV/channel", C_KEV},
        {"keV/channel^2", C_KEV},
        {"MeV", C_MEV},
        {"m",  1.0},
        {"cm", C_CM},
        {"mm", C_MM},
        {"us", C_US},
        {"ns", C_NS},
        {"ps", C_PS},
        {"msr", C_MSR},
        {"sr", 1.0},
        {0, 0}};

typedef struct idf_parser {
    char *filename;
    char *basename; /* Same as filename, but without the extension */
    xmlDoc *doc;
    xmlNode *root_element;
    char *buf;
    size_t buf_size;
    size_t pos_write;
    idf_error error;
    size_t n_samples;
    size_t i_sample; /* Keep track on how many samples are defined in the file */
    char *sample_basename;
    size_t n_spectra; /* Number of spectra in a sample (the one that is currently being parsed) */
    size_t i_spectrum;
} idf_parser;

xmlNode *idf_findnode_deeper(xmlNode *root, const char *path, const char **path_next); /* Called by idf_findnode(), goes one step deeper in path */
xmlNode *idf_findnode(xmlNode *root, const char *path); /* Finds the first node with given path. Path should look "like/this/thing" */
size_t idf_foreach(idf_parser *idf, xmlNode *node, const char *name, idf_error (*f)(idf_parser *idf, xmlNode *node)); /* Runs f(parser,child_node) for each child element of node element name "name" */
char *idf_node_content_to_str(const xmlNode *node); /* will always return a char * which should be freed by free() */
const xmlChar *idf_xmlstr(const char *s);
double idf_node_content_to_double(const xmlNode *node); /* performs unit conversion */
int idf_node_content_to_boolean(const xmlNode *node);
double idf_unit_string_to_SI(xmlChar *unit);
double idf_unit_mode(xmlChar *mode); /* FWHM etc */
int idf_stringeq(const void *a, const void *b);
int idf_stringneq(const void *a, const void *b, size_t n);
idf_error idf_write_simple_data_to_file(const char *filename, const char *x, const char *y);
idf_error idf_output_printf(idf_parser *idf, const char *format, ...);
idf_error idf_buffer_realloc(idf_parser *idf);
idf_parser *idf_file_read(const char *filename);
void idf_file_free(idf_parser *idf);
idf_error idf_write_buf_to_file(const idf_parser *idf, char **filename_out); /* Name is generated automatically and set to filename_out (if it is not NULL). */
idf_error idf_write_buf(const idf_parser *idf, FILE *f);
char *idf_jbs_name(const idf_parser *idf);
char *idf_exp_name(const idf_parser *idf, size_t i_spectrum);
char *idf_spectrum_out_name(const idf_parser *idf, size_t i_spectrum);
const char *idf_boolean_to_str(int boolean); /* "true", "false", "unset" trinary */
const char *idf_error_code_to_str(idf_error idferr);
#endif // IDFPARSER_H
