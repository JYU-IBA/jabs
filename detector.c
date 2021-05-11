#include <stdio.h>
#include <stdlib.h>
#include <jibal_config.h>
#include "defaults.h"
#include "detector.h"
extern inline double detector_calibrated(const detector *det, size_t ch);

int detector_sanity_check(const detector *det) {
    if(!det) {
        fprintf(stderr, "No detector!\n");
        return -1;
    }
    if (det->resolution <= 0.0) {
        fprintf(stderr, "Warning: detector resolution (%g) is negative.\n", det->resolution);
        return -1;
    }
    if (det->slope <= 0.0) {
        fprintf(stderr, "Warning: detector slope (%g) is negative.\n", det->slope);
        return -1;
    }
    return 0;
}

detector *detector_from_file(const jibal *jibal, const char *filename) {
    detector *det = detector_default();
    if(!det)
        return NULL;
    FILE *f = fopen(filename, "r");
    if(!f) {
        fprintf(stderr, "Could not read detector from file \"%s\".\n", filename);
        return det;
    }
    det->resolution = C_FWHM * sqrt(det->resolution); /* Convert resolution to FWHM from variance for the duration of input parsing */
    jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT,   "slope",      &det->slope,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "offset",     &det->offset,     NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "resolution", &det->resolution, NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "theta",      &det->theta,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "phi",        &det->phi,        NULL},
            {JIBAL_CONFIG_VAR_INT,    "number",     &det->number,     NULL},
            {JIBAL_CONFIG_VAR_INT,    "channels",   &det->channels,   NULL},
            {JIBAL_CONFIG_VAR_INT,    "compress",   &det->compress,   NULL},
            {JIBAL_CONFIG_VAR_STRING, "foil",       &det->foil_description,       NULL},
            {JIBAL_CONFIG_VAR_NONE,NULL,NULL,NULL}
    };
    jibal_config_var_read(jibal->units, f, filename, vars);
    det->resolution /= C_FWHM;
    det->resolution *= det->resolution;
    if(det->foil_description) {
        det->foil = sample_from_sample_model(sample_model_from_string(jibal, det->foil_description));
        if(!det->foil) {
            fprintf(stderr, "Error: detector foil description %s was not parsed successfully!\n", det->foil_description);

        }
    }
#ifdef DEBUG
    fprintf(stderr, "Read detector from \"%s\":\n", filename);
    detector_print(stderr, det);
#endif
    return det;
}

detector *detector_default() {
    detector *det = malloc(sizeof(detector));
    det->theta = THETA;
    det->phi = 0.0;
    det->resolution = (DETECTOR_RESOLUTION*DETECTOR_RESOLUTION);
    det->slope = ENERGY_SLOPE;
    det->offset = 0.0*C_KEV;
    det->number = 1; /* This implies default file format has channel numbers. Values are in the second column (number 1). */
    det->channels = 16384;
    det->compress = 1;
    det->foil = NULL;
    det->foil_description = NULL;
    return det;
}

void detector_free(detector *det) {
    if(!det)
        return;
    sample_free(det->foil);
    free(det);
}

void detector_print(FILE *f, const detector *det) {
    fprintf(f, "slope = %g keV\n", det->slope/C_KEV);
    fprintf(f, "offset = %g keV\n", det->offset/C_KEV);
    fprintf(f, "resolution = %g keV\n", C_FWHM*sqrt(det->resolution)/C_KEV);
    fprintf(f, "theta = %g deg\n", det->theta/C_DEG);
    fprintf(f, "phi = %g deg\n", det->phi/C_DEG);
    fprintf(f, "number = %zu\n", det->number);
    fprintf(f, "channels = %zu\n", det->channels);
    if(det->foil_description) {
        fprintf(f, "foil = %s\n", det->foil_description);
    }
}
