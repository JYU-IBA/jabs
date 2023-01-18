#include <libxml/parser.h>
#include <libxml/tree.h>
#include "idfparse.h"
#include "idfelementparsers.h"
#include "idf2jbs.h"

idf_error idf_parse(const char *filename, char **filename_out) {
    idf_parser *idf = idf_file_read(filename);
    if(!idf) {
        return IDF2JBS_FAILURE;
    }
    if(idf->error) {
        return IDF2JBS_FAILURE;
    }
    idf_foreach(idf, idf->root_element, "sample", idf_parse_sample);
    idf_output_printf(idf, "simulate\n");
    char *fn = NULL;
    idf_write_buf_to_file(idf,  &fn);
    if(filename_out) {
        *filename_out = fn;
    } else {
        free(fn);
    }
    idf_file_free(idf);
    return IDF2JBS_SUCCESS;
}
