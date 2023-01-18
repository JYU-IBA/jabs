#include <string.h>
#include "idfparse.h"


xmlNode *findnode_deeper(xmlNode *root, const char *path, const char **path_next) {
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

xmlNode *findnode(xmlNode *root, const char *path) {
    const char *path_remains = NULL;
    xmlNode *node = root;
    if(!root || !path) {
        return NULL;
    }
#ifdef DEBUG
    fprintf(stderr, "findnode(node name = %s, path = %s) called.\n", root->name, path);
#endif
    do {
        node = findnode_deeper(node, path, &path_remains);
        path = path_remains;
    } while(node && path_remains && *path_remains != '\0');
    return node;
}

int idf_foreach(idfparser *idf, xmlNode *node, const char *name, int (*f)(idfparser *idf, xmlNode *node)) {
    xmlNode *cur = NULL;
    int n = 0;
    for (cur = node->children; cur; cur = cur->next) {
        if(cur->type == XML_ELEMENT_NODE) {
            if(idf_stringeq(cur->name, name)) {
#ifdef DEBUG
                fprintf(stderr, "Found sample.\n");
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

int idf_nodename_equals(const xmlNode *node, const char *s) {
    if(!node || !s)
        return -1;
    return (strcmp((char *)node->name, s) == 0);
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

double idf_unit_string_to_SI(xmlChar *unit) {
    for(const idfunit *u = idfunits; u->unit; u++) {
        if(xmlStrEqual(unit, idf_xmlstr(u->unit))) {
            return u->factor;
        }
    }
    return 0.0;
}

int idf_stringeq(const void *a, const void *b) {
    if(!a || !b)
        return -1;
    return (strcmp(a,b) == 0);
}

int idf_write_simple_data_to_file(const char *filename, const char *x, const char *y) {
#if 0
    FILE *f = fopen(filename, "w");
    if(!f) {
        return IDF2JBS_FAILURE;
    }
#endif
    while(1) {
        char *x_end, *y_end;
        double x_dbl = strtod(x, &x_end);
        double y_dbl = strtod(y, &y_end);
        if(y == y_end || x == x_end)  { /* No conversion */
            break;
        }
        fprintf(stderr, "%g %g\n", floor(x_dbl), y_dbl);
        y = y_end;
        x = x_end;
    }
#if 0
    fclose(f);
#endif
    return IDF2JBS_SUCCESS;
}
