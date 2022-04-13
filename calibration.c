#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include "string.h"
#include <jibal_units.h>
#include "generic.h"
#include "calibration.h"

extern inline double calibration_eval(const calibration *c, double x);

calibration *calibration_init() {
    calibration *c = malloc(sizeof(calibration));
    if(!c)
        return NULL;
    c->f = calibration_none;
    c->type = CALIBRATION_NONE;
    c->params = NULL;
    return c;
}

void calibration_free(calibration *c) {
    if(!c)
        return;
    if(c->type == CALIBRATION_POLY) {
        calibration_params_poly *p = (calibration_params_poly *)c->params;
        free(p->a);
    }
    free(c->params);
    free(c);
}

calibration *calibration_init_linear() {
    calibration *c = calibration_init();
    if(!c)
        return NULL;
    c->f = calibration_linear;
    c->type = CALIBRATION_LINEAR;
    calibration_params_linear *p = malloc(sizeof(calibration_params_linear));
    p->offset = 0.0;
    p->slope = 0.0;
    c->params = p;
    c->resolution = 0.0;
    c->resolution_variance = 0.0;
    return c;
}

calibration *calibration_init_poly(size_t n) {
    calibration *c = calibration_init();
    if(!c)
        return NULL;
    c->f = calibration_poly;
    c->type = CALIBRATION_POLY;
    calibration_params_poly *p = malloc(sizeof(calibration_params_poly));
    c->params = p;
    p->n = n;
    p->a = calloc(p->n + 1, sizeof(double));
    return c;
}

double calibration_linear(const void *params, double x) {
    calibration_params_linear *p = (calibration_params_linear *)params;
    return p->offset + p->slope * x;
}

double calibration_poly(const void *params, double x) {
    calibration_params_poly *p = (calibration_params_poly *)params;
    double mul = x;
    double sum = p->a[0];
    for(size_t i = 1; i <= p->n; i++) {
        sum += mul * p->a[i];
        mul *= x;
    }
    return sum;
}

double calibration_none(const void *params, double x) {
    (void) params;
    return x;
}

int calibration_set_param(calibration *c, int i, double value) {
    if(!c || !c->params)
        return EXIT_FAILURE;
    if(i < 0) {
        if(i == CALIBRATION_PARAM_RESOLUTION) {
            c->resolution = value;
            return EXIT_SUCCESS;
        }
        return EXIT_FAILURE;
    }
    if(c->type == CALIBRATION_LINEAR) {
        calibration_params_linear *p = (calibration_params_linear *) c->params;
        switch(i) {
            case CALIBRATION_PARAM_OFFSET:
                p->offset = value;
                break;
            case CALIBRATION_PARAM_SLOPE:
                p->slope = value;
                break;
            default:
                return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    }
    if(c->type == CALIBRATION_POLY) {
        calibration_params_poly *p = (calibration_params_poly *) c->params;
        if(i <= (int)p->n) {
            p->a[i] = value;
            return EXIT_SUCCESS;
        } else {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

double calibration_get_param(const calibration *c, int i) {
    if(!c || !c->params)
        return 0.0;
    const double *val = calibration_get_param_ref((calibration *) c, i);
    if(val) {
        return *val;
    } else {
        return 0.0;
    }
}

size_t calibration_get_number_of_params(const calibration *c) { /* Number does not include "resolution" parameter! */
    if(c->type == CALIBRATION_LINEAR) {
        return 2;
    }
    if(c->type == CALIBRATION_POLY) {
        const calibration_params_poly *pp = (const calibration_params_poly *) c->params;
        return pp->n + 1;
    }
    return 0;
}

double *calibration_get_param_ref(calibration *c, int i) {
    if(!c || !c->params)
        return NULL;
    if(i == CALIBRATION_PARAM_RESOLUTION) {
        return &(c->resolution);
    }
    if(c->type == CALIBRATION_LINEAR) {
        calibration_params_linear *p = (calibration_params_linear *) c->params;
        switch(i) {
            case CALIBRATION_PARAM_OFFSET:
                return &(p->offset);
            case CALIBRATION_PARAM_SLOPE:
                return &(p->slope);
            default:
                break;
        }
        return NULL;
    }
    if(c->type == CALIBRATION_POLY) {
        calibration_params_poly *p = (calibration_params_poly*) c->params;
        if(i <= (int) p->n) {
            return &(p->a[i]);
        } else {
            return NULL;
        }
    }
    return NULL;
}

int calibration_copy_params(calibration *dst, calibration *src) {
    if(!src || !dst)
        return EXIT_FAILURE;
    size_t n_dst = calibration_get_number_of_params(dst);
    size_t n_src = calibration_get_number_of_params(src);
    size_t n = n_dst < n_src ? n_dst : n_src;
    for(int i = CALIBRATION_PARAM_RESOLUTION; i < (int)n; i++) { /* Loop includes resolution. Note that the parameters are copied over as-is. This works with linear and poly, but maybe not with others. */
        calibration_set_param(dst, i , calibration_get_param(src, i));
    }
    return EXIT_SUCCESS;
}

const char *calibration_name(const calibration *c) {
    if(!c)
        return calibration_option[CALIBRATION_NONE].s;
    return calibration_option[c->type].s;
}

char *calibration_to_string(const calibration *c) { /* Note that this does not include "resolution" although resolution parameter is considered part of a calibration. The reason is that calibrations don't know the type of detector. */
    char *out = NULL;
    asprintf_append(&out, "%s", calibration_name(c));
    if(!c)
        return out;
    switch(c->type) {
        case CALIBRATION_LINEAR:
            asprintf_append(&out, " slope %g%s offset %g%s",
                            calibration_get_param(c, CALIBRATION_PARAM_SLOPE)/C_KEV, "keV",
                            calibration_get_param(c, CALIBRATION_PARAM_OFFSET)/C_KEV, "keV"
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

char *calibration_param_name(calibration_type type, calibration_param_type i) {
    char *s = NULL;
    switch(i) {
        case CALIBRATION_PARAM_RESOLUTION:
            s = "resolution";
            break;
        case CALIBRATION_PARAM_OFFSET:
            s = "offset";
            break;
        case CALIBRATION_PARAM_SLOPE:
            s = "slope";
            break;
        case CALIBRATION_PARAM_QUAD:
            s = "quad";
            break;
        default:
            break;
    }
    if(s) {
        return strdup(s);
    }
    if(asprintf(&s, "calib_p%i", i) < 0) {
        return NULL;
    }
    return s;
}
