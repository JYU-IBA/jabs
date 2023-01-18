#include <stdio.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "idfparse.h"
#include "idfelementparsers.h"
#include "idf2jbs.h"

int idffile_parse(const char *filename) {
    xmlDoc *doc = NULL;
    xmlNode *root_element = NULL;
    doc = xmlReadFile(filename, NULL, 0);
    if(!doc) {
        return IDF2JBS_FAILURE_COULD_NOT_READ;
    }
    root_element = xmlDocGetRootElement(doc);
    if(!idf_stringeq(root_element->name, "idf")) {
        return IDF2JBS_FAILURE_NOT_IDF_FILE;
    }
    idfparser *idf = malloc(sizeof(idfparser));
    idf->doc = doc;
    idf->root_element = root_element;
    idf_foreach(idf, root_element, "sample", idf_parse_sample);
    xmlFreeDoc(doc);
    return IDF2JBS_SUCCESS;
}
