#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jibal_config.h>

#include "generic.h"
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
    FILE *f = fopen(filename, "r");
    if(!f) {
        fprintf(stderr, "Could not read detector from file \"%s\".\n", filename);
        return NULL;
    }
    detector *det = detector_default();
    if(!det)
        return NULL;
    det->resolution = C_FWHM * sqrt(det->resolution); /* Convert resolution to FWHM from variance for the duration of input parsing */
    jibal_config_file *cf = jibal_config_file_init(jibal->units);
    jibal_config_var *vars = detector_make_vars(det); /* Will be freed when config is free'd */
    jibal_config_file_set_vars(cf, vars);
    if(jibal_config_file_read(cf, filename)) {
        fprintf(stderr, "Could not read detector from \"%s\"\n", filename);
        detector_free(det);
        jibal_config_file_free(cf);
        return NULL;
    }
    jibal_config_file_free(cf);
    det->resolution /= C_FWHM;
    det->resolution *= det->resolution;
    detector_update_foil(jibal, det);
#ifdef DEBUG
    fprintf(stderr, "Read detector from \"%s\":\n", filename);
    detector_print(NULL, det);
#endif
    return det;
}

detector *detector_default() {
    detector *det = malloc(sizeof(detector));
    det->theta = DETECTOR_THETA;
    det->phi = DETECTOR_PHI;
    det->solid = DETECTOR_SOLID;
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

int detector_print(const char *filename, const detector *det) {
    if(!det)
        return EXIT_FAILURE;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return EXIT_FAILURE;
    fprintf(f, "slope = %g keV\n", det->slope/C_KEV);
    fprintf(f, "offset = %g keV\n", det->offset/C_KEV);
    fprintf(f, "resolution = %g keV\n", C_FWHM*sqrt(det->resolution)/C_KEV);
    fprintf(f, "theta = %g deg\n", det->theta/C_DEG);
    fprintf(f, "phi = %g deg\n", det->phi/C_DEG);
    fprintf(f, "solid = %g msr\n", det->solid/C_MSR);
    fprintf(f, "number = %zu\n", det->number);
    fprintf(f, "channels = %zu\n", det->channels);
    if(det->foil_description) {
        fprintf(f, "foil = %s\n", det->foil_description);
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

int detector_update_foil(const jibal *jibal, detector *det) {
    if(!det)
        return EXIT_FAILURE;
    sample_free(det->foil);
    det->foil = NULL;
    if(!det->foil_description)
        return EXIT_SUCCESS; /* Successfully cleared foil */
    sample_model *sm = sample_model_from_string(jibal, det->foil_description);
    det->foil = sample_from_sample_model(sm);
    sample_model_free(sm);
    if(!det->foil) {
        fprintf(stderr, "Error: detector foil description %s was not parsed successfully!\n", det->foil_description);
        free(det->foil_description);
        det->foil_description = NULL;
    }
    return EXIT_SUCCESS;
}

int detector_set_var(const jibal *jibal, detector *det, const char *var_str, const char *val_str) {
    if(!jibal || !det || !var_str || !val_str)
        return EXIT_FAILURE;
    jibal_config_var *vars = detector_make_vars(det);
    if(!vars)
        return EXIT_FAILURE;
    int found = FALSE;
    int error = FALSE;
    for(jibal_config_var *var = vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
        if(strcmp(var_str, var->name) == 0) {
            jibal_config_var_set(jibal->units, var, val_str, NULL);
            if(var->variable == &det->foil_description) {
                if(detector_update_foil(jibal, det)) {
                    error = TRUE;
                }
            } else if(var->variable == &det->resolution) {
                det->resolution /= C_FWHM;
                det->resolution *= det->resolution;
            }
            found = TRUE;
            break;
        }
    }
    free(vars);
    if(found && !error) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}

jibal_config_var *detector_make_vars(detector *det) {
    if(!det)
        return NULL;
    jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT, "slope",      &det->slope,          NULL},
            {JIBAL_CONFIG_VAR_UNIT, "offset",     &det->offset,         NULL},
            {JIBAL_CONFIG_VAR_UNIT, "resolution", &det->resolution,     NULL},
            {JIBAL_CONFIG_VAR_UNIT, "theta",      &det->theta,          NULL},
            {JIBAL_CONFIG_VAR_UNIT, "phi",        &det->phi,            NULL},
            {JIBAL_CONFIG_VAR_UNIT, "solid",      &det->solid,          NULL},
            {JIBAL_CONFIG_VAR_INT,  "number",     &det->number,         NULL},
            {JIBAL_CONFIG_VAR_INT,  "channels",   &det->channels,       NULL},
            {JIBAL_CONFIG_VAR_INT,  "compress",   &det->compress,       NULL},
            {JIBAL_CONFIG_VAR_STRING,"foil",       &det->foil_description,NULL},
            {JIBAL_CONFIG_VAR_NONE, NULL, NULL, NULL}
    };
    int n_vars;
    for(n_vars = 0; vars[n_vars].type != 0; n_vars++);
    size_t var_size = sizeof(jibal_config_var)*(n_vars + 1); /* +1 because the null termination didn't count */
    jibal_config_var *vars_out = malloc(var_size);
    if(vars_out) {
        memcpy(vars_out, vars, var_size); /* Note that this does not make a deep copy */
    }
    return vars_out;
}

