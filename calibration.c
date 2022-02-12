#include <assert.h>
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
    free(c->params);
    free(c);
}

calibration *calibration_init_linear(double offset, double slope) {
    calibration *c = calibration_init();
    if(!c)
        return NULL;
    c->f = calibration_linear;
    c->type = CALIBRATION_LINEAR;
    calibration_params_linear *p = malloc(sizeof(calibration_params_linear));
    p->offset = offset;
    p->slope = slope;
    c->params = p;
    return c;
}



double calibration_linear(const void *params, double x) {
    calibration_params_linear *p = (calibration_params_linear *)params;
    return p->offset + p->slope * x;
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
