#ifndef CALIB_CALIBRATION_H
#define CALIB_CALIBRATION_H
#include <inttypes.h>
#include <stdlib.h>
#include <jibal_option.h>

typedef enum calibration_type {
    CALIBRATION_NONE = 0,
    CALIBRATION_ARB = 1,
    CALIBRATION_LINEAR = 2
} calibration_type;

static const jibal_option calibration_option[] = {
        {JIBAL_OPTION_STR_NONE, CALIBRATION_NONE},
        {"arbitrary", CALIBRATION_ARB},
        {"linear", CALIBRATION_LINEAR},
        {NULL, 0}
};

typedef enum calibration_param_type {
    CALIBRATION_PARAM_NONE = 0,
    CALIBRATION_PARAM_OFFSET = 1,
    CALIBRATION_PARAM_SLOPE = 2
} calibration_param_type;

typedef struct calibration {
    calibration_type type;
    double (*f)(const void *, double);
    void *params;
} calibration;

typedef struct calibration_params_linear {
    double offset;
    double slope;
} calibration_params_linear;

calibration *calibration_init();
void calibration_free(calibration *c);
calibration *calibration_init_linear(double offset, double slope);
double calibration_linear(const void *params, double x);
double calibration_eval(const calibration *c, double x);
int calibration_set_param(calibration *c, calibration_param_type type, double value);
double calibration_get_param(const calibration *c, calibration_param_type type);
const char *calibration_name(const calibration *c);
#endif //CALIB_CALIBRATION_H
