#ifndef CALIB_CALIBRATION_H
#define CALIB_CALIBRATION_H
#include <inttypes.h>
#include <stdlib.h>

typedef enum calibration_type {
    CALIBRATION_NONE = 0,
    CALIBRATION_ARB = 1,
    CALIBRATION_LINEAR = 2
} calibration_type;

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
#endif //CALIB_CALIBRATION_H
