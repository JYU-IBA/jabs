#include <stdlib.h>
#include <stdio.h>
#include <jibal.h>

int main(int argc, char **argv) {
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return 1;
    }
    if(argc < 2) {
        fprintf(stderr, "Not enough arguments! Usage: %s: isotope\n", argv[0]);
        return 1;
    }
    const jibal_isotope *incident = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    if(incident) {
        fprintf(stderr, "Mass %.3lf u\n", incident->mass/C_U);
    }
    jibal_free(jibal);
    return 0;
}
