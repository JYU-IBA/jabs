#include <libxml/parser.h>
#include <libxml/tree.h>
#include "idfparse.h"
#include "idfelementparsers.h"
#include "idf2jbs.h"
#include "options.h"

idf_error idf_parse(const char *filename, char **filename_out) {
    idf_parser *idf = idf_file_read(filename);
    if(!idf) {
        return IDF2JBS_FAILURE;
    }
    if(idf->error) {
        return idf->error;
    }
    idf_output_printf(idf, "#Converted from IDF file \"%s\" by JaBS version %s\n", filename, jabs_version());
    idf_foreach(idf, idf_findnode(idf->root_element, "notes"), "note", idf_parse_note);
    idf->n_samples = idf_foreach(idf, idf->root_element, "sample", NULL);

    if(idf->n_samples == 0) {
        idf_file_free(idf);
        return IDF2JBS_FAILURE_NO_SAMPLES_DEFINED;
    }

    idf_foreach(idf, idf->root_element, "sample", idf_parse_sample);

    char *fn = NULL;
    idf->error = idf_write_buf_to_file(idf,  &fn);
    if(filename_out) {
        *filename_out = fn;
    } else {
        free(fn);
    }
    idf_file_free(idf);
    return idf->error;
}
