#ifndef IDFPARSER_H
#define IDFPARSER_H
#include <libxml/tree.h>
#include <jibal_units.h>

#define IDF_BUF_MSG_MAX (1024) /* Maximum length of input string to idf_output_printf. Longer strings will be truncated. */
#define IDF_BUF_SIZE_INITIAL (8*IDF_BUF_MSG_MAX)
#define JABS_FILE_SUFFIX (".jbs")
#define SPECTRUM_FILE_SUFFIX (".dat")

typedef enum idf_error {
    IDF2JBS_SUCCESS = 0,
    IDF2JBS_FAILURE = -1,
    IDF2JBS_FAILURE_COULD_NOT_READ = -2,
    IDF2JBS_FAILURE_NOT_IDF_FILE = -3
} idf_error;

typedef struct idf_unit {
    const char *unit;
    double factor;
} idf_unit;

static const idf_unit idf_units[] = {
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
        {"msr", C_MSR},
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
} idf_parser;

xmlNode *findnode_deeper(xmlNode *root, const char *path, const char **path_next); /* Called by findnode(), goes one step deeper in path */
xmlNode *findnode(xmlNode *root, const char *path); /* Finds the first node with given path. Path should look "like/this/thing" */
idf_error idf_foreach(idf_parser *idf, xmlNode *node, const char *name, idf_error (*f)(idf_parser *idf, xmlNode *node)); /* Runs f(parser,child_node) for each child element of node element name "name" */
int idf_nodename_equals(const xmlNode *node, const char *s);
char *idf_node_content_to_str(const xmlNode *node);
const xmlChar *idf_xmlstr(const char *s);
double idf_node_content_to_double(const xmlNode *node);
double idf_unit_string_to_SI(xmlChar *unit);
int idf_stringeq(const void *a, const void *b);
idf_error idf_write_simple_data_to_file(const char *filename, const char *x, const char *y);
idf_error idf_output_printf(idf_parser *idf, const char *format, ...);
idf_error idf_buffer_realloc(idf_parser *idf);
idf_parser *idf_file_read(const char *filename);
void idf_file_free(idf_parser *idf);
idf_error idf_write_buf_to_file(const idf_parser *idf, char **filename_out); /* Name is generated automatically and set to filename_out (if it is not NULL). */
idf_error idf_write_buf(const idf_parser *idf, FILE *f);
char *idf_file_name_with_suffix(const idf_parser *idf, const char *suffix);
#endif // IDFPARSER_H
