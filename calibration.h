#ifndef CALIB_CALIBRATION_H
#define CALIB_CALIBRATION_H
#include <inttypes.h>
#include <stdlib.h>
#include <jibal_option.h>

typedef enum calibration_type {
    CALIBRATION_NONE = 0,
    CALIBRATION_ARB = 1,
    CALIBRATION_LINEAR = 2,
    CALIBRATION_POLY = 3
} calibration_type;

static const jibal_option calibration_option[] = {
        {JIBAL_OPTION_STR_NONE, CALIBRATION_NONE},
        {"arbitrary", CALIBRATION_ARB}, /* Not implemented */
        {"linear", CALIBRATION_LINEAR},
        {"polynomial", CALIBRATION_POLY},
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

typedef struct calibration_params_poly {
    size_t n;
    double *a; /* Array, size n+1, polynomial is f(x) = a[0] + a[1] * x + a[2] * x^2 + ... + a[n] * x^n */
} calibration_params_poly;

calibration *calibration_init();
void calibration_free(calibration *c);
calibration *calibration_init_linear();
calibration *calibration_init_poly(size_t n); /* n is the degree of the polynomial */
double calibration_linear(const void *params, double x);
double calibration_poly(const void *params, double x);
double calibration_none(const void *params, double x);
inline double calibration_eval(const calibration *c, double x) {return c->f(c->params, x);}
int calibration_set_param(calibration *c, calibration_param_type type, double value);
double calibration_get_param(const calibration *c, calibration_param_type type);
size_t calibration_get_number_of_params(const calibration *c);
double calibration_get_param_number(const calibration *c, size_t i); /* get i'th param in range [0..n-1], get n by  calibration_get_number_of_params()*/
double *calibration_get_param_ref(calibration *c, calibration_param_type type);
const char *calibration_name(const calibration *c);
#endif //CALIB_CALIBRATION_H
