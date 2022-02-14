#include <assert.h>
#include "generic.h"
#include "calibration.h"

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
    return x;
}

int calibration_set_param(calibration *c, calibration_param_type type, double value) {
    if(!c || !c->params)
        return EXIT_FAILURE;
    if(c->type != CALIBRATION_LINEAR)
        return EXIT_FAILURE;
    calibration_params_linear *p = (calibration_params_linear *) c->params;
    switch(type) {
        case CALIBRATION_PARAM_NONE:
            break;
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

double calibration_get_param(const calibration *c, calibration_param_type type) {
    if(!c || !c->params)
        return 0.0;
    if(c->type != CALIBRATION_LINEAR)
        return 0.0;
    calibration_params_linear *p = (calibration_params_linear *) c->params;
    switch(type) {
        case CALIBRATION_PARAM_OFFSET:
            return p->offset;
        case CALIBRATION_PARAM_SLOPE:
            return p->slope;
        default:
            break;
    }
    return 0.0;
}

size_t calibration_get_number_of_params(const calibration *c) {
    if(c->type == CALIBRATION_LINEAR) {
        return 2;
    }
    if(c->type == CALIBRATION_POLY) {
        const calibration_params_poly *pp = (const calibration_params_poly *) c->params;
        return pp->n + 1;
    }
    return 0;
}

double calibration_get_param_number(const calibration *c, size_t i) {
    if(c->type == CALIBRATION_LINEAR) {
        const calibration_params_linear *pl = (const calibration_params_linear *) c->params;
        switch(i) {
            case 0:
                return pl->offset;
            case 1:
                return pl->slope;
            default:
                return 0.0;
        }
    }
    if(c->type == CALIBRATION_POLY) {
        const calibration_params_poly *pp = (const calibration_params_poly *) c->params;
        if(i <= pp->n) {
            return pp->a[i];
        } else {
            return 0.0;
        }
    }
    return 0.0;
}

double *calibration_get_param_ref(calibration *c, calibration_param_type type) {
    if(!c || !c->params)
        return NULL;
    if(c->type != CALIBRATION_LINEAR)
        return NULL;
    calibration_params_linear *p = (calibration_params_linear *) c->params;
    switch(type) {
        case CALIBRATION_PARAM_OFFSET:
            return &(p->offset);
        case CALIBRATION_PARAM_SLOPE:
            return &(p->slope);
        default:
            break;
    }
    return NULL;
}

const char *calibration_name(const calibration *c) {
    if(!c)
        return calibration_option[CALIBRATION_NONE].s;
    return calibration_option[c->type].s;
}
