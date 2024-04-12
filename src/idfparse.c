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
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "generic.h"
#ifdef WIN32
#include "win_compat.h"
#endif
#include "jabs_debug.h"
#include "idfparse.h"

xmlNode *idf_findnode_deeper(xmlNode *root, const char *path, const char **path_next) {
    DEBUGMSG("findnode_deeper called with path = \"%s\".", path);
    size_t offset = strcspn(path, "/"); /* We use strcspn to tokenize because we don't have to mutilate path by inserting '\0' in it */
    if(offset == 0) {
        fprintf(stderr, "Weirdness. Offset is zero.\n");
        return NULL;
    }
    int last = (path[offset] == '\0');
    if(path_next) {
        *path_next = path + offset + (last ? 0 : 1);
    }
    xmlNode *cur_node = NULL;
    DEBUGMSG("Should find %.*s in %s", (int)offset, path, root->name);
    for (cur_node = root->children; cur_node; cur_node = cur_node->next) {
        if(cur_node->type == XML_ELEMENT_NODE) {
            size_t len = xmlStrlen(cur_node->name);
            if(len != offset) /* Prevents partial matches in strncmp */
                continue;
            if(strncmp((char *)cur_node->name, path, offset) == 0) {
                return cur_node;
            }
        }
    }
    DEBUGSTR("Returning NULL");
    return NULL;
}

xmlNode *idf_findnode(xmlNode *root, const char *path) {
    const char *path_remains = NULL;
    xmlNode *node = root;
    if(!root || !path) {
        return NULL;
    }
    DEBUGMSG("findnode(node name = %s, path = %s) called.", root->name, path);
    do {
        node = idf_findnode_deeper(node, path, &path_remains);
        path = path_remains;
    } while(node && path_remains && *path_remains != '\0');
    return node;
}

size_t idf_foreach(idf_parser *idf, xmlNode *node, const char *name, idf_error (*f)(idf_parser *idf, xmlNode *node)) {
    if(!node) {
        return IDF2JBS_FAILURE;
    }
    xmlNode *cur = NULL;
    size_t n = 0;
    for (cur = node->children; cur; cur = cur->next) {
        if(cur->type == XML_ELEMENT_NODE) {
            if(idf_stringeq(cur->name, name)) {
                DEBUGMSG("foreach found %s.", name);
                if(f && f(idf, cur)) {
                    return IDF2JBS_FAILURE;
                }
                n++;
            }
        }
    }
    return n;
}

char *idf_node_content_to_str(const xmlNode *node) {
    if(!node) {
        return strdup("");
    }
    xmlChar *content = xmlNodeGetContent(node);
    if(!content) {
        return strdup("");
    }
    return (char *) content;
}

const xmlChar *idf_xmlstr(const char *s) {
    return (const xmlChar *)s;
}

double idf_node_content_to_double(const xmlNode *node) {
    xmlChar *content = xmlNodeGetContent(node);
    if(!content) {
        return 0.0;
    }
    double out = strtod((char *)content, NULL); /* TODO: error checking */
    xmlChar *unit = xmlGetProp(node, idf_xmlstr("units"));
    if(unit) {
        out *= idf_unit_string_to_SI(unit);
        free(unit);
    }
    xmlChar *mode = xmlGetProp(node, idf_xmlstr("mode"));
    if(mode) {
        out *= idf_unit_mode(mode);
        free(mode);
    }
    free(content);
    return out;
}

int idf_node_content_to_boolean(const xmlNode *node) {
    xmlChar *content = xmlNodeGetContent(node);
    int out = -1; /* Invalid */
    if(idf_stringeq(content, "true")) {
        out = TRUE;
    } else if(idf_stringeq(content, "false")) {
        out = FALSE;
    }
    free(content);
    return out;
}

double idf_unit_string_to_SI(const xmlChar *unit) {
    for(const idf_unit *u = idf_units; u->unit; u++) {
        if(xmlStrEqual(unit, idf_xmlstr(u->unit))) {
            return u->factor;
        }
    }
    DEBUGMSG("No such unit: \"%s\"!", unit);
    return 0.0;
}

double idf_unit_mode(const xmlChar *mode) {
    if(idf_stringeq(mode, IDF_MODE_FWHM)) {
        return 1.0/C_FWHM;
    }
    return 1.0;
}

int idf_stringeq(const void *a, const void *b) {
    if(!a || !b)
        return -1;
    return (strcmp(a, b) == 0);
}

int idf_stringneq(const void *a, const void *b, size_t n) {
    if(!a || !b)
        return -1;
    return (strncmp(a, b, n) == 0);
}

idf_error idf_write_simple_data_to_file(const char *filename, const char *x, const char *y) {
    FILE *f = fopen(filename, "w");
    if(!f) {
        return IDF2JBS_FAILURE;
    }
    size_t n = 0;
    while(1) {
        char *x_end, *y_end;
        double x_dbl = strtod(x, &x_end);
        double y_dbl = strtod(y, &y_end);
        if(y == y_end || x == x_end)  { /* No conversion */
            break;
        }
        fprintf(f, "%g %g\n", floor(x_dbl), y_dbl);
        y = y_end;
        x = x_end;
        n++;
    }
    fclose(f);
    if(n >= 2) { /* At least two successful conversions should be performed before we can call it a spectrum */
        return IDF2JBS_SUCCESS;
    } else {
        return IDF2JBS_FAILURE;
    }
}

idf_error idf_output_printf(idf_parser *idf, const char *format, ...) {
    va_list argp;
    va_start(argp, format);
    if(idf->pos_write + IDF_BUF_MSG_MAX >= idf->buf_size) {
        idf_buffer_realloc(idf); /* We trust that there is now enough space... */
    }
    if(!idf->buf) {
        return IDF2JBS_FAILURE;
    }
    size_t n = vsnprintf(idf->buf + idf->pos_write, IDF_BUF_MSG_MAX, format, argp);
    if(n >= IDF_BUF_MSG_MAX) {
        DEBUGMSG("Output truncated, because %zu > %i", n, IDF_BUF_MSG_MAX);
        n = IDF_BUF_MSG_MAX - 1;
    }
    idf->pos_write += n;
    va_end(argp);
    return IDF2JBS_SUCCESS;
}

int idf_output_puts(idf_parser *idf, const char *s) {
    if(!idf) {
        return EOF;
    }
    size_t len = strlen(s);
    while(idf->buf && idf->pos_write + len >= idf->buf_size) {
        idf_buffer_realloc(idf);
    }
    if(!idf->buf) {
        return EOF;
    }
    strncat(idf->buf, s, len);
    return 0;
}

idf_error idf_buffer_realloc(idf_parser *idf) {
    size_t size_new;
    if(!idf->buf) {
        idf->buf_size = 0;
        idf->pos_write = 0;
        size_new = IDF_BUF_SIZE_INITIAL;
    } else {
        size_new = idf->buf_size * 2; /* TODO: good strategy? */
    }
    DEBUGMSG("Buffer will be reallocated. New size: %zu, old size: %zu", size_new, idf->buf_size);
    idf->buf = realloc(idf->buf, size_new);
    if(!idf->buf) {
        idf->buf_size = 0;
        return IDF2JBS_FAILURE;
    }
    idf->buf_size = size_new;
    idf->buf[idf->pos_write] = '\0'; /* Terminating the string is necessary only after the initial allocation, but this doesn't hurt much. */
    return IDF2JBS_SUCCESS;
}

idf_parser *idf_file_read(const char *filename) {
    idf_parser *idf = calloc(1, sizeof(idf_parser));
    if(!filename) {
        idf->error = IDF2JBS_FAILURE_COULD_NOT_READ;
        return idf;
    }
    char *fn = strdup(filename);
    char *extension = jabs_file_extension(fn);
    if(!idf_stringeq(extension, ".xml") && !idf_stringeq(extension, ".idf") && !idf_stringeq(extension, ".xnra")) {
        DEBUGMSG("Extension %s is not valid.", extension);
        free(fn);
        idf->error = IDF2JBS_FAILURE_WRONG_EXTENSION;
        return idf;
    }
    idf->doc = xmlReadFile(fn, NULL, 0);
    if(!idf->doc) {
        idf->error = IDF2JBS_FAILURE_COULD_NOT_READ;
        free(fn);
        return idf;
    }
    idf->root_element = xmlDocGetRootElement(idf->doc);
    if(!idf_stringeq(idf->root_element->name, "idf")) {
        idf->error = IDF2JBS_FAILURE_NOT_IDF_FILE;
        free(fn);
        return idf;
    }
    idf->filename = strdup(fn);
    *extension = '\0'; /* extension still points to start of extension, rest of fn can be used for basename */
    idf->basename = strdup(fn);
    free(fn);
    if(strlen(idf->basename) == 0) { /* Without this filename of ".xml" would be valid and this could cause problems. */
        idf->error = IDF2JBS_FAILURE;
        return idf;
    }
    idf->buf = NULL;
    idf_buffer_realloc(idf);
    return idf;
}
void idf_file_free(idf_parser *idf) {
    if(!idf) {
        return;
    }
    xmlFreeDoc(idf->doc);
    free(idf->filename);
    free(idf->basename);
    free(idf->sample_basename);
    free(idf->buf);
    free(idf);
}

idf_error idf_write_buf_to_file(const idf_parser *idf, char **filename_out) {
    char *filename = idf_jbs_name(idf);
    FILE *f = fopen(filename, "w");
    if(!f) {
        free(filename);
        return IDF2JBS_FAILURE_COULD_NOT_WRITE_OUTPUT;
    }
    if(idf_write_buf(idf, f)) {
       return IDF2JBS_FAILURE_COULD_NOT_WRITE_OUTPUT;
    }
    if(filename_out) {
        *filename_out = filename;
    } else {
        free(filename);
    }
    fclose(f);
    return IDF2JBS_SUCCESS;
}

idf_error idf_write_buf(const idf_parser *idf, FILE *f) {
    return (fputs(idf->buf, f) == EOF) ? IDF2JBS_FAILURE : IDF2JBS_SUCCESS;
}

char *idf_jbs_name(const idf_parser *idf) {
    char *out;
    int len = asprintf(&out, "%s%s", idf->basename, JABS_FILE_SUFFIX);
    if(len < 0) {
        return NULL;
    }
    return out;
}

char *idf_exp_name(const idf_parser *idf, size_t i_spectrum) {
    char *out;
    int len;
    if(idf->n_spectra == 1) {
        len = asprintf(&out, "%s_exp%s", idf->sample_basename, EXP_SPECTRUM_FILE_SUFFIX);
    } else {
        len = asprintf(&out, "%s_exp%zu%s", idf->sample_basename, i_spectrum, EXP_SPECTRUM_FILE_SUFFIX);
    }
    if(len < 0) {
        return NULL;
    }
    return out;
}

char *idf_sim_name(const idf_parser *idf, size_t i_spectrum) {
    char *out;
    int len;
    if(idf->n_spectra == 1) {
        len = asprintf(&out, "%s_sim%s", idf->sample_basename, EXP_SPECTRUM_FILE_SUFFIX);
    } else {
        len = asprintf(&out, "%s_sim%zu%s", idf->sample_basename, i_spectrum, EXP_SPECTRUM_FILE_SUFFIX);
    }
    if(len < 0) {
        return NULL;
    }
    return out;
}

char *idf_spectrum_out_name(const idf_parser *idf, size_t i_spectrum) {
    char *out;
    int len;
    if(idf->n_spectra == 1) {
        len = asprintf(&out, "%s_out%s", idf->sample_basename, SAVE_SPECTRUM_FILE_SUFFIX);
    } else {
        len = asprintf(&out, "%s_out%zu%s", idf->sample_basename, i_spectrum, SAVE_SPECTRUM_FILE_SUFFIX);
    }
    if(len < 0) {
        return NULL;
    }
    return out;
}

const char *idf_boolean_to_str(int boolean) {
    if(boolean == TRUE) {
        return "true";
    } else if(boolean == FALSE) {
        return "false";
    } else {
        return "unset";
    }
}

const char *idf_error_code_to_str(idf_error idferr) {
    switch(idferr) {
        case IDF2JBS_SUCCESS:
            return "success";
        case IDF2JBS_FAILURE:
            return "generic failure";
        case IDF2JBS_FAILURE_COULD_NOT_READ:
            return "could not read file";
        case IDF2JBS_FAILURE_NOT_IDF_FILE:
            return "not an IDF file";
        case IDF2JBS_FAILURE_WRONG_EXTENSION:
            return "wrong file extension";
        case IDF2JBS_FAILURE_NO_SAMPLES_DEFINED:
            return "no samples defined in file";
        default:
            return "unknown error code";
    }
}

xmlNodePtr idf_new_node_printf(const xmlChar *name, const char * restrict format, ...) { /* Creates a new xmlNode with content using formatted print */
    va_list argp;
    va_start(argp, format);
    char *s;
    int len = vasprintf(&s, format, argp);
    va_end(argp);
    if(len < 0) {
        return NULL;
    }
    xmlNodePtr n = xmlNewNode(NULL, BAD_CAST name);
    if(s) {
        xmlNodeSetContent(n, BAD_CAST s);
    }
    free(s);
    return n;
}

xmlNodePtr idf_new_node_units(const xmlChar *name, const xmlChar *unit, const xmlChar *mode, double value) {
    xmlNodePtr n = xmlNewNode(NULL, BAD_CAST name);
    if(unit) {
        value /= idf_unit_string_to_SI(unit);
        xmlSetProp(n, BAD_CAST "units", unit);
    }
    if(mode) {
        value /= idf_unit_mode(mode);
        xmlSetProp(n, BAD_CAST "mode", mode);
    }
    char *out;
    int len = asprintf(&out, "%g", value);
    if(len < 0) {
        xmlFreeNode(n);
        return NULL;
    }
    if(out) {
        xmlNodeSetContent(n, BAD_CAST out);
    }
    free(out);
    return n;
}
