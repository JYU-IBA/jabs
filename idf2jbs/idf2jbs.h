#include <libxml/tree.h>
#include <jibal_units.h>
#define IDF2JBS_SUCCESS (0)
#define IDF2JBS_FAILURE (1)

#define CALIB_PARAMS_MAX (10)

#define IDF_UNIT_TFU "1e15at/cm2"
#define JABS_UNIT_TFU "tfu"

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
    int tmp;
} idfparser;

int parse_xml(const char *filename);
xmlNode *findnode(xmlNode *root, const char *path);
int nodename_equals(const xmlNode *node, const char *s);
double unit_string_to_SI(xmlChar *unit);
xmlNode *findnode_deeper(xmlNode *root, const char *path, const char **path_next);
xmlNode *findnode(xmlNode *root, const char *path);
double node_content_to_double(const xmlNode *node); /* Performs conversion to SI if possible and unit is defined */
char *node_content_to_str(const xmlNode *node);
const xmlChar *xmlstr(const char *s); /* The purpose of this function is to reduce compiler warnings. It's just a recast. */
int write_simple_data_to_file(const char *filename, const char *x, const char *y);
