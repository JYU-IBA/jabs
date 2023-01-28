#include <string.h>
#ifdef WIN32
#include "win_compat.h"
#endif
#include "idfparse.h"

xmlNode *idf_findnode_deeper(xmlNode *root, const char *path, const char **path_next) {
#ifdef DEBUG
    fprintf(stderr, "findnode_deeper called with path = \"%s\".\n", path);
#endif
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
#ifdef DEBUG
    fprintf(stderr, "Should find %.*s in %s\n", (int)offset, path, root->name);
#endif
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
#ifdef DEBUG
    fprintf(stderr, "Returning NULL\n");
#endif
    return NULL;
}

xmlNode *idf_findnode(xmlNode *root, const char *path) {
    const char *path_remains = NULL;
    xmlNode *node = root;
    if(!root || !path) {
        return NULL;
    }
#ifdef DEBUG
    fprintf(stderr, "findnode(node name = %s, path = %s) called.\n", root->name, path);
#endif
    do {
        node = idf_findnode_deeper(node, path, &path_remains);
        path = path_remains;
    } while(node && path_remains && *path_remains != '\0');
    return node;
}

idf_error idf_foreach(idf_parser *idf, xmlNode *node, const char *name, idf_error (*f)(idf_parser *idf, xmlNode *node)) {
    if(!node) {
        return IDF2JBS_FAILURE;
    }
    xmlNode *cur = NULL;
    int n = 0;
    for (cur = node->children; cur; cur = cur->next) {
        if(cur->type == XML_ELEMENT_NODE) {
            if(idf_stringeq(cur->name, name)) {
#ifdef DEBUG
                fprintf(stderr, "foreach found %s.\n", name);
#endif
                if(f(idf, cur)) {
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
    double out = strtod((char *)content, NULL); /* TODO: error checking */
    xmlChar *unit = xmlGetProp(node, idf_xmlstr("units"));
    if(unit) {
        out *= idf_unit_string_to_SI(unit);
        free(unit);
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

double idf_unit_string_to_SI(xmlChar *unit) {
    for(const idf_unit *u = idf_units; u->unit; u++) {
        if(xmlStrEqual(unit, idf_xmlstr(u->unit))) {
            return u->factor;
        }
    }
    return 0.0;
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
    if(n > 1) { /* At least two successfull conversions should be performed before we can call it a spectrum */
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
        fprintf(stderr, "Output truncated.\n");
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
#ifdef DEBUG
    fprintf(stderr, "Buffer will be reallocated. New size: %zu, old size: %zu\n", size_new, idf->buf_size);
#endif
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
    size_t l = strlen(filename);
    if(l < 4) {
        fprintf(stderr, "Filename %s is suspiciously short.\n", filename);
        idf->error = IDF2JBS_FAILURE_COULD_NOT_READ;
        return idf;
    }
    const char *extension = filename + l - 4;
    if(!idf_stringeq(extension, ".xml") && !idf_stringeq(extension, ".idf")) {
#ifdef DEBUG
        fprintf(stderr, "Extension %s is not valid.\n", extension);
#endif
        idf->error = IDF2JBS_FAILURE_COULD_NOT_READ;
        return idf;
    }
    idf->doc = xmlReadFile(filename, NULL, 0);
    if(!idf->doc) {
        idf->error = IDF2JBS_FAILURE_COULD_NOT_READ;
        return idf;
    }
    idf->root_element = xmlDocGetRootElement(idf->doc);
    if(!idf_stringeq(idf->root_element->name, "idf")) {
        idf->error = IDF2JBS_FAILURE_NOT_IDF_FILE;
        return idf;
    }
    idf->filename = strdup(filename);
    idf->basename = strndup(filename, l - 4);
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
    free(idf->buf);
    free(idf);
}

idf_error idf_write_buf_to_file(const idf_parser *idf, char **filename_out) {
    char *filename = idf_file_name_with_suffix(idf, JABS_FILE_SUFFIX);
    FILE *f = fopen(filename, "w");
    if(!f) {
        free(filename);
        return IDF2JBS_FAILURE;
    }
    if(idf_write_buf(idf, f)) {
       return IDF2JBS_FAILURE;
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

char *idf_file_name_with_suffix(const idf_parser *idf, const char *suffix) {
    size_t len = strlen(idf->basename) + strlen(suffix) + 1;
    char *fn = malloc(sizeof(char) * len);
    strcpy(fn, idf->basename);
    strcat(fn, suffix);
    return fn;
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
