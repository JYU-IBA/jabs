#include <libxml/parser.h>
#include <libxml/tree.h>
#include "idfparse.h"
#include "idfelementparsers.h"
#include "idf2jbs.h"

idf_error idffile_parse(const char *filename) {
    idf_parser *idf = idf_file_read(filename);
    if(!idf) {
        return IDF2JBS_FAILURE;
    }
    idf_foreach(idf, idf->root_element, "sample", idf_parse_sample);
    idf_output_printf(idf, "simulate\n");
    char *filename_out = NULL;
    idf_write_buf_to_file(idf,  &filename_out);
    free(filename_out);
#ifdef DEBUG
    if(filename_out) {
        fprintf(stderr, "Wrote file %s\n", filename_out);
    }
#endif
    idf_file_free(idf);
    return IDF2JBS_SUCCESS;
}
