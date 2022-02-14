#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jibal_config.h>

#include "generic.h"
#include "defaults.h"
#include "detector.h"
#include "message.h"
#include "win_compat.h"

extern inline double detector_calibrated(const detector *det, size_t ch);

char *detector_calibration_to_string(const detector *det) {
    if(!det)
        return NULL;
    const calibration *c = det->calibration;
    char *out = NULL;
    asprintf_append(&out, "%s", calibration_name(c));
    if(!c)
        return out;
    switch(c->type) {
        case CALIBRATION_LINEAR:
            asprintf_append(&out, " slope %g%s offset %g%s",
#ifdef DETECTOR_NATIVE_SPECTRA /* TODO: when simulating ToF spectra with a ToF detector, we want to use the code below, otherwise slope and offset are in energy units (see else-branch) */
                            calibration_get_param(c, CALIBRATION_PARAM_SLOPE)/detector_param_unit_factor(det), detector_param_unit(det),
                            calibration_get_param(c, CALIBRATION_PARAM_OFFSET)/detector_param_unit_factor(det), detector_param_unit(det)
#else
                            calibration_get_param(c, CALIBRATION_PARAM_SLOPE)/C_KEV, "keV",
                            calibration_get_param(c, CALIBRATION_PARAM_OFFSET)/C_KEV, "keV"
#endif
                            );
            break;
        case CALIBRATION_POLY:
            for(size_t i = 0; i < calibration_get_number_of_params(c); i++) {
                asprintf_append(&out, " %g%s", calibration_get_param(c, i)/C_KEV, "keV");
            }
            break;
        default:
            break;
    }
    return out;
}


const char *detector_type_name(const detector *det) {
    return detector_option[det->type].s;
}

int detector_sanity_check(const detector *det) {
    if(!det) {
        jabs_message(MSG_ERROR, stderr, "No detector!\n");
        return -1;
    }
    if(det->resolution <= 0.0) {
        jabs_message(MSG_ERROR, stderr, "Warning: detector resolution (%g) is negative.\n", det->resolution);
        return -1;
    }
    double slope = calibration_get_param(det->calibration, CALIBRATION_PARAM_SLOPE);
    if(slope < 0.0) {
        jabs_message(MSG_ERROR, stderr, "Warning: detector slope (%g) is negative.\n", slope);
        return -1;
    }
    if(det->type == DETECTOR_TOF && det->length < 1 * C_MM) {
        jabs_message(MSG_ERROR, stderr, "Warning: length (%g) is small (%g mm)\n", det->length/C_MM);
        return -1;
    }
    return 0;
}

detector *detector_from_file(const jibal *jibal, const char *filename) {
    FILE *f = fopen(filename, "r");
    if(!f) {
        jabs_message(MSG_ERROR, stderr, "Could not read detector from file \"%s\".\n", filename);
        return NULL;
    }
    detector *det = detector_default(NULL);
    if(!det)
        return NULL;
    jibal_config_file *cf = jibal_config_file_init(jibal->units);
    jibal_config_var *vars = detector_make_vars(det); /* Will be freed when config is free'd */
    jibal_config_file_set_vars(cf, vars);
    if(jibal_config_file_read(cf, filename)) {
        jabs_message(MSG_ERROR, stderr, "Could not read detector from \"%s\"\n", filename);
        detector_free(det);
        jibal_config_file_free(cf);
        return NULL;
    }
    jibal_config_file_free(cf);
    detector_update_foil(det); /* TODO: foil and aperture are not read! */
#ifdef DEBUG
    fprintf(stderr, "Read detector from \"%s\":\n", filename);
    detector_print(NULL, det);
#endif
    return det;
}

detector **detectors_from_file(const jibal *jibal, const char *filename, size_t *n_detectors_out) {
    char *line = NULL;
    size_t line_size = 0;
    size_t lineno = 0;
    *n_detectors_out = 0;

    size_t n_detectors = 0;
    detector **detectors = NULL;
    detector *det = NULL;

    FILE *in = fopen_file_or_stream(filename, "r");
    if(!in)
        return NULL;

    char **header_strings = NULL;
    int n_header_strings = 0;
    if(getline(&line, &line_size, in) > 0) {
        header_strings = string_to_argv(line, &n_header_strings);
    } else {
        return NULL;
    }

    while(getline(&line, &line_size, in) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strips all kinds of newlines! */
        if(strlen(line) >= 1 && *line == '#') /* Comment */
            continue;
        char *line_split = line;
        char *col;
        size_t n = 0; /* Number of columns on this row */
        while((col = strsep(&line_split, " \t")) != NULL) {
            if(*col == '\0') {/* Multiple separators are treated as one */
                continue;
            }
            if(n == 0) { /* First column, (re)allocate space for a new detector */
                n_detectors++;
                detectors = realloc(detectors, sizeof(detector) * (n_detectors));
                det = detector_default(detectors[n_detectors-1]);
            }
            if(n < (size_t) n_header_strings) {
                detector_set_var(jibal, det, header_strings[n], col);
            }
            n++;
        }
    }
    *n_detectors_out = n_detectors;
    argv_free(header_strings, n_header_strings);
    return detectors;
}

detector *detector_default(detector *det) {
    if(!det) {
        det = malloc(sizeof(detector));
    } else {
        calibration_free(det->calibration);
    }
    det->type = DETECTOR_ENERGY;
    det->theta = DETECTOR_THETA;
    det->phi = DETECTOR_PHI;
    det->solid = DETECTOR_SOLID;
    det->aperture = NULL;
    det->distance = DETECTOR_DISTANCE;
    det->length = DETECTOR_LENGTH;
    det->resolution = DETECTOR_RESOLUTION;
    det->resolution_variance = 0.0; /* Calculated before needed. */
    det->column = 1; /* This implies default file format has channel numbers. Values are in the second column (number 1). */
    det->channels = 16384;
    det->compress = 1;
    det->foil = NULL;
    det->foil_sm = NULL;
    det->calibration = calibration_init_linear();
    calibration_set_param(det->calibration, CALIBRATION_PARAM_SLOPE, ENERGY_SLOPE);
    return det;
}

void detector_free(detector *det) {
    if(!det)
        return;
    sample_model_free(det->foil_sm);
    sample_free(det->foil);
    aperture_free(det->aperture);
    calibration_free(det->calibration);
    free(det);
}

int detector_print(const char *filename, const detector *det) {
    if(!det)
        return EXIT_FAILURE;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return EXIT_FAILURE;
    jabs_message(MSG_INFO, f, "type = %s\n", detector_type_name(det));
    char *calib_str = detector_calibration_to_string(det);
    jabs_message(MSG_INFO, f, "calibration = %s\n", calib_str);
    free(calib_str);
    if(det->type == DETECTOR_ENERGY) {
        jabs_message(MSG_INFO, f, "resolution = %g keV\n", det->resolution/C_KEV);
    } else if(det->type == DETECTOR_TOF) {
        jabs_message(MSG_INFO, f, "length = %g mm\n", det->length/C_MM);
        jabs_message(MSG_INFO, f, "resolution = %g ps\n", det->resolution/C_PS);
    } else if(det->type == DETECTOR_ELECTROSTATIC) {
        jabs_message(MSG_INFO, f, "resolution = %g\n", det->resolution);
    }
    jabs_message(MSG_INFO, f, "theta = %g deg\n", det->theta/C_DEG);
    jabs_message(MSG_INFO, f, "phi = %g deg\n", det->phi/C_DEG);
#if 0
    jabs_message(MSG_INFO, f, "angle from horizontal = %.3lf deg\n", detector_angle(det, 'x')/C_DEG);
    jabs_message(MSG_INFO, f, "angle from vertical = %.3lf deg\n", detector_angle(det, 'y')/C_DEG);
#endif
    jabs_message(MSG_INFO, f, "solid = %g msr\n", det->solid/C_MSR);
    if(det->aperture) {
        char *s = aperture_to_string(det->aperture);
        jabs_message(MSG_INFO, f, "aperture = %s\n", s);
        free(s);
    }
    jabs_message(MSG_INFO, f, "distance = %g mm\n", det->distance/C_MM);
#if 0
    if(det->distance > 1.0*C_MM) {
        jabs_message(MSG_INFO, f, "solid angle (calculated) = %.4lf msr\n", detector_solid_angle_calc(det)/C_MSR);
    }
#endif
    jabs_message(MSG_INFO, f, "column = %zu\n", det->column);
    jabs_message(MSG_INFO, f, "channels = %zu\n", det->channels);
    if(det->foil) {
        char *foil_str = sample_model_to_string(det->foil_sm);
        if(foil_str) {
            jabs_message(MSG_INFO, f, "foil = %s\n", foil_str);
            free(foil_str);
        }
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

int detector_aperture_set_from_argv(const jibal *jibal, detector *det, int *argc, char * const **argv) {
    aperture *a = aperture_from_argv(jibal, argc, argv);
    if(a) {
        aperture_free(det->aperture);
        det->aperture = a;
    } else {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int detector_foil_set_from_argv(const jibal *jibal, detector *det, int *argc, char * const **argv) {
    sample_model *sm = sample_model_from_argv(jibal, argc, argv);
    if(sm) {
        free(det->foil_sm);
        det->foil_sm = sm;
        detector_update_foil(det);
    } else {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int detector_update_foil(detector *det) {
    det->foil = sample_from_sample_model(det->foil_sm);
    return 0;
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
            {JIBAL_CONFIG_VAR_OPTION, "type",       &det->type,          detector_option},
            {JIBAL_CONFIG_VAR_UNIT,   "length",     &det->length,            NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "resolution", &det->resolution,        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "theta",      &det->theta,             NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "phi",        &det->phi,               NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "solid",      &det->solid,             NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "distance",   &det->distance,          NULL},
            {JIBAL_CONFIG_VAR_INT,    "column",     &det->column,            NULL},
            {JIBAL_CONFIG_VAR_INT,    "channels",   &det->channels,          NULL},
            {JIBAL_CONFIG_VAR_INT,    "compress",   &det->compress,          NULL},
            {JIBAL_CONFIG_VAR_NONE, NULL, NULL,                              NULL}
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

double detector_angle(const detector *det, const char direction) { /* Gives detector angle (to an axis, see angle_tilt()) */
    double angle = C_PI - angle_tilt(det->theta, det->phi, direction); /* The pi is here because our detector angles are defined oddly */
    angle = fmod(angle, C_2PI);
    if(angle > C_PI)
        angle -= C_2PI;
    return angle;
}

double detector_theta_deriv(const detector *det, const char direction) { /* We don't need to use this. Probably doesn't work when det->theta = 180 deg TODO: remove. */
    static const double delta = 0.001*C_DEG;
    double theta, phi;

    theta = delta;
    if(direction == 'x') {
        phi = 0.0;
    } else if(direction == 'y') {
        phi = C_PI_2;
    } else {
        return 0.0;
    }
    rotate(theta, phi, det->theta, det->phi, &theta, &phi);
    double result = (theta - det->theta)/delta;
#ifdef DEBUG
    fprintf(stderr, "(%.7lf deg - %.7lf deg)/(%g deg) = %g\n", theta/C_DEG, det->theta/C_DEG, delta/C_DEG, result);
#endif
    return result; /* TODO: sign of result? */
}

double detector_solid_angle_calc(const detector *det) {
    if(det->distance < 1.0*C_MM)
        return 0.0;
    if(!det->aperture)
        return 0.0;
    if(det->aperture->type == APERTURE_CIRCLE) {
        return 2.0*C_PI*(1.0 - 1.0/(sqrt(pow2(det->aperture->diameter/2.0/det->distance)+1.0)));
    }
    if(det->aperture->type == APERTURE_RECTANGLE) {
        double alpha = det->aperture->width / det->distance / 2.0;
        double beta =  det->aperture->height / det->distance / 2.0;
        return 4.0 * atan(alpha * beta / sqrt(1 + pow2(alpha) + pow2(beta)));
    }
    return 0.0;
}

double detector_resolution(const detector *det, const jibal_isotope *isotope, double E) { /* Returns variance! */
    switch(det->type) {
        case DETECTOR_NONE:
            return 0.0;
        case DETECTOR_ENERGY:
            return det->resolution_variance;
        case DETECTOR_TOF:
            return det->resolution_variance * pow3(2*E)/(pow2(det->length)*isotope->mass);
        case DETECTOR_ELECTROSTATIC:
            return det->resolution_variance * pow2(E);
    }
    return 0.0; /* Never reached */
}

void detector_update(detector *det) {
    det->resolution_variance = pow2(det->resolution/C_FWHM);
#ifdef DEBUG
    fprintf(stderr, "Updated detector, resolution = %g keV FWHM, variance = %g\n", det->resolution/C_KEV, det->resolution_variance);
#endif
}

const char *detector_param_unit(const detector *det) { /* Match these with detector_param_unit_factor() */
    switch(det->type) {
        case DETECTOR_TOF:
            return "ns";
        case DETECTOR_ENERGY:
            return "keV";
        default:
            return "";
    }
}

double detector_param_unit_factor(const detector *det) {
    switch(det->type) {
        case DETECTOR_TOF:
            return C_NS;
        case DETECTOR_ENERGY:
            return C_KEV;
        default:
            return 1.0;
    }
}
