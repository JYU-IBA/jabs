#ifndef IDFPARSER_H
#define IDFPARSER_H
#include <libxml/tree.h>
#include <jibal_units.h>

#define IDF2JBS_SUCCESS (0)
#define IDF2JBS_FAILURE (-1)
#define IDF2JBS_FAILURE_COULD_NOT_READ (-2)
#define IDF2JBS_FAILURE_NOT_IDF_FILE (-3)

typedef struct {
    const char *unit;
    double factor;
} idfunit;

static const idfunit idfunits[] = {
        {"fraction", 1.0},
        {"degree", C_DEG},
        {"1e15at/cm2", C_TFU},
        {"#particles", 1.0},
        {"amu", C_U},
        {"keV", C_KEV},
        {"keV/channel", C_KEV},
        {"keV/channel^2", C_KEV},
        {"MeV", C_MEV},
        {"mm", C_MM},
        {"us", C_US},
        {0, 0}};



typedef struct idfparser {
    char *filename;
    xmlDoc *doc;
    xmlNode *root_element;
} idfparser;


xmlNode *findnode_deeper(xmlNode *root, const char *path, const char **path_next); /* Called by findnode(), goes one step deeper in path */
xmlNode *findnode(xmlNode *root, const char *path); /* Finds the first node with given path. Path should look "like/this/thing" */
int idf_foreach(idfparser *idf, xmlNode *node, const char *name, int (*f)(idfparser *idf, xmlNode *node)); /* Runs f(parser,child_node) for each child element of node element name "name" */
int idf_nodename_equals(const xmlNode *node, const char *s);
char *idf_node_content_to_str(const xmlNode *node);
const xmlChar *idf_xmlstr(const char *s);
double idf_node_content_to_double(const xmlNode *node);
double idf_unit_string_to_SI(xmlChar *unit);
int idf_stringeq(const void *a, const void *b);
int idf_write_simple_data_to_file(const char *filename, const char *x, const char *y);
#endif // IDFPARSER_H
