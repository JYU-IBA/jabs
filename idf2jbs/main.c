#include <stdlib.h>
#include <stdio.h>

#include "idf2jbs.h"

int main(int argc, char **argv) {
    if(argc < 2) {
        fprintf(stderr, "Not enough arguments! Usage: idf2jbs <filename>\n");
        return EXIT_FAILURE;
    }
    const char *filename = argv[1];
    int ret = idffile_parse(filename);
    return ret;
}
