#define _GNU_SOURCE /* asprintf on Linux */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <jibal_config.h>

#include "jabs_debug.h"
#include "generic.h"
#include "defaults.h"
#include "detector.h"
#include "message.h"
#include "win_compat.h"

extern inline double detector_calibrated(const detector *det, int Z, size_t ch);

calibration *detector_get_calibration(const detector *det, int Z) {
    assert(det != NULL && Z >= JIBAL_ANY_Z);
    if(Z == JIBAL_ANY_Z || Z > det->cal_Z_max) {
        DEBUGVERBOSEMSG("Z = %i is either any Z (-1) or larger than %i. Returning default calibration (%p)", Z, det->cal_Z_max, (void *)det->calibration);
        return det->calibration;
    } else if(det->calibration_Z[Z]) {
        DEBUGVERBOSEMSG("There is a Z-specific calibration (%p)", (void *)det->calibration_Z[Z]);
        return det->calibration_Z[Z];
    } else {
        DEBUGVERBOSEMSG("There is a table of Z-specific calibrations, but nothing in there for Z = %i. Returning default calibration (%p)", Z, (void *)det->calibration);
        return det->calibration; /* Fallback */
    }
}

int detector_set_calibration_Z(const jibal_config *jibal_config, detector *det, calibration *cal, int Z) {
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
        DEBUGMSG("More space allocated for detector = %p calibrations. Z_max = %i.", (void *) det, det->cal_Z_max);
    }
    calibration_free(det->calibration_Z[Z]);
    det->calibration_Z[Z] = cal;
    DEBUGMSG("Calibration initialized (Z = %i) with this: %p.", Z, (void *)cal);
    return EXIT_SUCCESS;
}




const char *detector_type_name(const detector *det) {
    return detector_option[det->type].s;
}

int detector_sanity_check(const detector *det, size_t n_channels) {
    if(!det) {
        jabs_message(MSG_ERROR, stderr, "No detector!\n");
        return -1;
    }
    if(!isfinite(det->calibration->resolution) || det->calibration->resolution <= 0.0) {
        jabs_message(MSG_ERROR, stderr, "Warning: detector resolution (%g) is negative or not finite.\n", det->calibration->resolution);
        return -1;
    }
    if(calibration_is_monotonically_increasing(det->calibration, n_channels) == FALSE) {
        jabs_message(MSG_ERROR, stderr, "Warning: detector default calibration is not monotonically increasing! (is slope, %g keV/ch, negative? are %zu channels too much?)\n", calibration_get_param(det->calibration, CALIBRATION_PARAM_SLOPE) / C_KEV, n_channels);
        return -1;
    }
    for(int Z = 0; Z <= (int)det->cal_Z_max; Z++) {
        const calibration *c = detector_get_calibration(det, Z);
        if(det->calibration == c || !c) /* Z calibration is same as default (fallback) or NULL (shouldn't happen) */
            continue;
        if(calibration_is_monotonically_increasing(det->calibration, det->channels) == FALSE) {
            jabs_message(MSG_ERROR, stderr, "Warning: detector calibration for Z = %i is not monotonically increasing!\n", Z);
            return -1;
        }
    }
    if(det->type == DETECTOR_TOF && det->length < 1 * C_MM) {
        jabs_message(MSG_ERROR, stderr, "Warning: ToF length is small (%g mm)\n", det->length/C_MM);
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
    det->name = NULL; /* When added, default name is given */
    det->type = DETECTOR_ENERGY;
    det->theta = DETECTOR_THETA_DEFAULT;
    det->phi = DETECTOR_PHI_DEFAULT;
    det->solid = DETECTOR_SOLID_DEFAULT;
    det->aperture = NULL;
    det->distance = DETECTOR_DISTANCE_DEFAULT;
    det->length = DETECTOR_LENGTH_DEFAULT;
    det->column = 1; /* This implies default file format has channel numbers. Values are in the second column (number 1). */
    det->channels = CHANNELS_MAX_DEFAULT;
    det->compress = 1;
    det->foil = NULL;
    det->foil_sm = NULL;
    det->calibration = calibration_init_linear();
    det->cal_Z_max = -1;
    det->calibration_Z = NULL;
    calibration_set_param(det->calibration, CALIBRATION_PARAM_SLOPE, ENERGY_SLOPE_DEFAULT);
    calibration_set_param(det->calibration, CALIBRATION_PARAM_RESOLUTION, DETECTOR_RESOLUTION_DEFAULT);
    return det;
}

void detector_free(detector *det) {
    if(!det)
        return;
    sample_model_free(det->foil_sm);
    sample_free(det->foil);
    aperture_free(det->aperture);
    detector_calibrations_free(det);
    free(det->name);
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

int detector_print(const jibal *jibal, const detector *det) {
    if(!det) {
        DEBUGSTR("Tried to print detector that does not exist.\n");
        return EXIT_FAILURE;
    }
    jabs_message(MSG_INFO, stderr, "name = %s\n", det->name);
    jabs_message(MSG_INFO, stderr, "type = %s\n", detector_type_name(det));
    char *calib_str = calibration_to_string(det->calibration);
    jabs_message(MSG_INFO, stderr, "calibration = %s\n", calib_str);
    free(calib_str);
    char *reso_str = detector_resolution_to_string(det, JIBAL_ANY_Z);
    jabs_message(MSG_INFO, stderr, "resolution = %s\n", reso_str);
    free(reso_str);

    for(int Z = 0; Z <= (int)det->cal_Z_max; Z++) {
        const calibration *c = detector_get_calibration(det, Z);
        if(det->calibration == c || !c) /* Z calibration is same as default (fallback) or NULL (shouldn't happen) */
            continue;
        const char *elem_name = jibal_element_name(jibal->elements, Z);
        char *Z_calib_str = calibration_to_string(c);
        jabs_message(MSG_INFO, stderr, "calibration(%s) = %s\n", elem_name, Z_calib_str);
        free(Z_calib_str);
        char *Z_reso_str = detector_resolution_to_string(det, Z);
        jabs_message(MSG_INFO, stderr, "resolution(%s) = %s\n", elem_name, Z_reso_str);
        free(Z_reso_str);
    }
    if(det->type == DETECTOR_TOF) {
        jabs_message(MSG_INFO, stderr, "length = %g mm\n", det->length/C_MM);
    }
    jabs_message(MSG_INFO, stderr, "theta = %g deg\n", det->theta/C_DEG);
    jabs_message(MSG_INFO, stderr, "phi = %g deg\n", det->phi/C_DEG);
#if 0
    jabs_message(MSG_INFO, stderr, "angle from horizontal = %.3lf deg\n", detector_angle(det, 'x')/C_DEG);
    jabs_message(MSG_INFO, stderr, "angle from vertical = %.3lf deg\n", detector_angle(det, 'y')/C_DEG);
#endif

    jabs_message(MSG_INFO, stderr, "distance = %g mm\n", det->distance/C_MM);
    jabs_message(MSG_INFO, stderr, "solid = %g msr\n", det->solid/C_MSR);
    if(det->aperture) {
        char *s = aperture_to_string(det->aperture);
        jabs_message(MSG_INFO, stderr, "aperture = %s\n", s);
        free(s);
    }
    if(det->distance > 1.0*C_MM && det->aperture) {
        jabs_message(MSG_INFO, stderr, "solid angle (calculated) = %.4lf msr\n", detector_solid_angle_calc(det)/C_MSR);
    }
    jabs_message(MSG_INFO, stderr, "column = %zu\n", det->column);
    jabs_message(MSG_INFO, stderr, "channels = %zu\n", det->channels);
    if(det->foil) {
        char *foil_str = sample_model_to_string(det->foil_sm);
        if(foil_str) {
            jabs_message(MSG_INFO, stderr, "foil = %s\n", foil_str);
            free(foil_str);
        }
    }
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

double detector_angle(const detector *det, const char direction) { /* Gives detector angle (to an axis, see angle_tilt()) */
    double angle = C_PI - angle_tilt(det->theta, det->phi, direction); /* The pi is here because our detector angles are defined oddly */
    angle = fmod(angle, C_2PI);
    if(angle > C_PI)
        angle -= C_2PI;
    return angle;
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
    if(asprintf(&out, "%.4g%s", resolution/detector_param_unit_factor(det), detector_param_unit(det)) < 0) {
        return NULL;
    }
    return out;
}

void detector_update(detector *det) {
    if(!det) {
        return;
    }
    det->calibration->resolution_variance = pow2(det->calibration->resolution/C_FWHM);
    DEBUGMSG("Updating detector, resolution = %g keV FWHM, variance = %g", det->calibration->resolution/C_KEV, det->calibration->resolution_variance);
    for(int Z = 1; Z  <= (int)det->cal_Z_max; Z++) {
        calibration *c = det->calibration_Z[Z];
        if(!c || c == det->calibration) {
            continue;
        }
        c->resolution_variance = pow2(c->resolution/C_FWHM);
        DEBUGMSG("Updating detector, Z = %i, resolution = %g keV FWHM, variance = %g", Z, det->calibration->resolution/C_KEV, det->calibration->resolution_variance);
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

int detector_set_name(detector *det, const char *name) {
    if(!name) {
        return EXIT_FAILURE;
    }
    if(isdigit(*name)) {
        return EXIT_FAILURE;
    }
    static const char *keywords[] = {"first", "last", "unnamed", "none",
                                     "aperture", "beta", "calibration", "channels",
                                     "column", "compress", "distance", "foil", "length",
                                     "name", "offset", "phi", "resolution", "slope",
                                     "solid", "theta", "type", 0}; /* Not allowed names */


    for(const char **k = keywords; *k != 0; k++) {
        if(strcmp(*k, name) == 0) {
            return EXIT_FAILURE;
        }
    }
    free(det->name);
    det->name = strdup(name);
    if(!det->name) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

const char *detector_name(detector *det) {
    if(det){
        if(det->name) {
            return det->name;
        } else {
            return "unnamed";
        }
    }
}
