#define _GNU_SOURCE /* asprintf on Linux */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <jibal_config.h>

#include "generic.h"
#include "defaults.h"
#include "detector.h"
#include "message.h"
#include "win_compat.h"

extern inline double detector_calibrated(const detector *det, int Z, size_t ch);

calibration *detector_get_calibration(const detector *det, int Z) {
    if(det == NULL || Z < JIBAL_ANY_Z)
        return NULL;
    if(Z == JIBAL_ANY_Z || Z > det->cal_Z_max) {
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Z = %i is either any Z (-1) or larger than %i. Returning default calibration (%p)\n", Z, det->cal_Z_max, (void *)det->calibration);
#endif
        return det->calibration;
    }
    assert(det->cal_Z_max > 0 && det->calibration);
    if(det->calibration_Z[Z]) {
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "There is a Z-specific calibration (%p)\n", (void *)det->calibration_Z[Z]);
#endif
        return det->calibration_Z[Z];
    } else {
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "There is a table of Z-specific calibrations, but nothing in there for Z = %i. Returning default calibration (%p)\n", Z, (void *)det->calibration);
#endif
        return det->calibration; /* Fallback */
    }
}

int detector_set_calibration_Z(const jibal_config *jibal_config, detector *det, calibration *cal, int Z) { /* TODO: test! */
    if(!det || !cal)
        return EXIT_FAILURE;
    if(Z == JIBAL_ANY_Z) {
        calibration_free(det->calibration);
        det->calibration = cal;
        return EXIT_SUCCESS;
    }
    if(Z < 1) {
        return EXIT_FAILURE;
    }
    if(Z > jibal_config->Z_max) { /* This is not strictly necessary, but makes things easier later on. */
        return EXIT_FAILURE;
    }
    if(Z > det->cal_Z_max) {
        det->calibration_Z = realloc(det->calibration_Z, sizeof(calibration *) * (Z+1)); /* Allocate more space */
        if(!det->calibration_Z) {
            det->cal_Z_max = 0;
            return EXIT_FAILURE;
        }
        for(int i = (int) det->cal_Z_max + 1; i <= Z; i++) { /* Initialize */
            det->calibration_Z[i] = NULL;
        }
        det->cal_Z_max = Z;
#ifdef DEBUG
        fprintf(stderr, "More space allocated for detector = %p calibrations. Z_max = %i.\n", (void *) det, det->cal_Z_max);
#endif
    }
    calibration_free(det->calibration_Z[Z]);
    det->calibration_Z[Z] = cal;
#ifdef DEBUG
    fprintf(stderr, "Calibration initialized (Z = %i) with this: %p.", Z, (void *)cal);
#endif

    return EXIT_SUCCESS;
}




const char *detector_type_name(const detector *det) {
    return detector_option[det->type].s;
}

int detector_sanity_check(const detector *det) {
    if(!det) {
        jabs_message(MSG_ERROR, stderr, "No detector!\n");
        return -1;
    }
    if(det->calibration->resolution <= 0.0) {
        jabs_message(MSG_ERROR, stderr, "Warning: detector resolution (%g) is negative.\n", det->calibration->resolution);
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
    det->column = 1; /* This implies default file format has channel numbers. Values are in the second column (number 1). */
    det->channels = 16384;
    det->compress = 1;
    det->foil = NULL;
    det->foil_sm = NULL;
    det->calibration = calibration_init_linear();
    det->cal_Z_max = -1;
    det->calibration_Z = NULL;
    calibration_set_param(det->calibration, CALIBRATION_PARAM_SLOPE, ENERGY_SLOPE);
    calibration_set_param(det->calibration, CALIBRATION_PARAM_RESOLUTION, DETECTOR_RESOLUTION);
    return det;
}

void detector_free(detector *det) {
    if(!det)
        return;
    sample_model_free(det->foil_sm);
    sample_free(det->foil);
    aperture_free(det->aperture);
    detector_calibrations_free(det);
    free(det);
}

void detector_calibrations_free(detector *det) {
    calibration_free(det->calibration);
    if(det->calibration_Z) {
        for(int Z = 0;  Z <= det->cal_Z_max; Z++) {
            calibration_free(det->calibration_Z[Z]);
        }
        free(det->calibration_Z);
    }
    det->cal_Z_max = -1;
}

int detector_print(const jibal *jibal, const char *filename, const detector *det) {
    if(!det)
        return EXIT_FAILURE;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return EXIT_FAILURE;
    jabs_message(MSG_INFO, f, "type = %s\n", detector_type_name(det));
    char *calib_str = calibration_to_string(det->calibration);
    jabs_message(MSG_INFO, f, "calibration = %s\n", calib_str);
    free(calib_str);
    char *reso_str = detector_resolution_to_string(det, JIBAL_ANY_Z);
    jabs_message(MSG_INFO, f, "resolution = %s\n", reso_str);

    for(int Z = 0; Z <= (int)det->cal_Z_max; Z++) {
        const calibration *c = detector_get_calibration(det, Z);
        if(det->calibration == c || !c) /* Z calibration is same as default (fallback) or NULL (shouldn't happen) */
            continue;
        const char *elem_name = jibal_element_name(jibal->elements, Z);
        char *Z_calib_str = calibration_to_string(c);
        jabs_message(MSG_INFO, f, "calibration(%s) = %s\n", elem_name, Z_calib_str);
        free(Z_calib_str);
        char *Z_reso_str = detector_resolution_to_string(det, Z);
        jabs_message(MSG_INFO, f, "resolution(%s) = %s\n", elem_name, Z_reso_str);
        free(Z_reso_str);
    }
    if(det->type == DETECTOR_TOF) {
        jabs_message(MSG_INFO, f, "length = %g mm\n", det->length/C_MM);
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
    det->aperture = aperture_set_from_argv(jibal, det->aperture, argc, argv);
    if(!det->aperture) {
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
    calibration *c = detector_get_calibration(det, isotope->Z);
    if(!c) {
        return 0.0;
    }
    switch(det->type) {
        case DETECTOR_NONE:
            return 0.0;
        case DETECTOR_ENERGY:
            return c->resolution_variance;
        case DETECTOR_TOF:
            return c->resolution_variance * pow3(2*E)/(pow2(det->length)*isotope->mass);
        case DETECTOR_ELECTROSTATIC:
            return c->resolution_variance  * pow2(E);
    }
    return 0.0; /* Never reached */
}

char *detector_resolution_to_string(const detector *det, int Z) {
    const calibration *c = detector_get_calibration(det, Z);
    double resolution = calibration_get_param(c, CALIBRATION_PARAM_RESOLUTION);
    char *out = NULL;
    asprintf(&out, "%g %s", resolution/detector_param_unit_factor(det), detector_param_unit(det));
    return out;
}

void detector_update(detector *det) {
    det->calibration->resolution_variance = pow2(det->calibration->resolution/C_FWHM);
#ifdef DEBUG
    fprintf(stderr, "Updated detector, resolution = %g keV FWHM, variance = %g\n", det->calibration->resolution/C_KEV, det->calibration->resolution_variance);
#endif
    for(int Z = 1; Z  <= (int)det->cal_Z_max; Z++) {
        calibration *c = det->calibration_Z[Z];
        if(!c || c == det->calibration) {
            continue;
        }
        c->resolution_variance = pow2(c->resolution/C_FWHM);
#ifdef DEBUG
        fprintf(stderr, "Updated detector, Z = %i, resolution = %g keV FWHM, variance = %g\n", Z, det->calibration->resolution/C_KEV, det->calibration->resolution_variance);
#endif
    }

}

const char *detector_param_unit(const detector *det) { /* Units to match detector type (at the moment for resolution). Match these with detector_param_unit_factor() */
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
