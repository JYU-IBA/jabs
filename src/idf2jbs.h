#include "idfparse.h"
idf_error idf_parse(const char *filename, char **filename_out);
/* Reads an IDF file, parses it and saves output as an JaBS script file (automatic file name stored in newly allocated filename_out if it is not NULL). */
