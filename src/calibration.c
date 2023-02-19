#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <jibal_units.h>
#include <jibal_phys.h>
#include "generic.h"
#include "calibration.h"
#include "jabs_debug.h"

extern inline double calibration_eval(const calibration *c, size_t ch);

calibration *calibration_init() {
    calibration *c = calloc(1, sizeof(calibration));
    if(!c)
        return NULL;
    c->f = calibration_none;
    c->type = CALIBRATION_NONE;
    return c;
}

calibration *calibration_clone(const calibration *c_orig) {
    if(!c_orig) {
        return NULL;
    }
    calibration *c;
    switch(c_orig->type) {
        case CALIBRATION_POLY:
            c = calibration_init_poly(calibration_get_number_of_params(c_orig));
            calibration_copy_params(c, c_orig);
            break;
        case CALIBRATION_LINEAR:
            c = calibration_init_linear();
            calibration_copy_params(c, c_orig);
            break;
        default:
            c = NULL;
            break;
    }
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
    calibration_params_linear *p = calloc(1, sizeof(calibration_params_linear));
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

double calibration_linear(const void *params, size_t ch) {
    calibration_params_linear *p = (calibration_params_linear *)params;
    return p->offset + p->slope * ch;
}

double calibration_poly(const void *params, size_t x) {
    calibration_params_poly *p = (calibration_params_poly *)params;
    double mul = x;
    double sum = p->a[0];
    for(size_t i = 1; i <= p->n; i++) {
        sum += mul * p->a[i];
        mul *= x;
    }
    return sum;
}

double calibration_none(const void *params, size_t ch) {
    (void) params;
    (void) ch;
    return 0.0;
}

size_t calibration_inverse(const calibration *cal, double E, size_t ch_max) {
    if(cal->type == CALIBRATION_LINEAR) {
        calibration_params_linear *p = (calibration_params_linear *) cal->params;
        double ch = (E - p->offset) / p->slope;
        if(ch >= 0.0 && ch < ch_max) {
            return (size_t) ch;
        } else {
            return 0;
        }
    } else if(cal->type == CALIBRATION_POLY) {
        if(calibration_get_number_of_params(cal) == 3) {
            const double a = calibration_get_param(cal, CALIBRATION_PARAM_OFFSET);
            const double b = calibration_get_param(cal, CALIBRATION_PARAM_SLOPE);
            const double c = calibration_get_param(cal, CALIBRATION_PARAM_QUAD);
            double discriminant = pow2(b) - 4 * (a - E) * c;
            if(c == 0.0) { /* This is annoying. */
                return 0;
            }
            if(discriminant < 0) {
                return 0;
            } else {
                return -b * sqrt(discriminant) / (2.0 * c);
            }
        }
        size_t lo = 0, mi, hi = ch_max;
        while(hi - lo > 1) {
            mi = (hi + lo) / 2;
            if(E >= calibration_eval(cal, mi)) {
                lo = mi;
            } else {
                hi = mi;
            }
        }
        return lo;
    } else {
        return 0;
    }
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
    assert(c && c->params);
    const double *val = calibration_get_param_ref((calibration *) c, i);
    assert(val);
    return *val;
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

int calibration_copy_params(calibration *dst, const calibration *src) {
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
    (void) type; /* We could use the calibration type to give different names for different parameters */
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
    if(asprintf(&s, "p%i", i) < 0) {
        return NULL;
    }
    return s;
}

int calibration_is_monotonically_increasing(const calibration *cal, size_t n_channels) {
    if(!cal) {
        return FALSE;
    }
    size_t n_params = calibration_get_number_of_params(cal);
    if(n_params < 2) { /* Calibration either constant, type is CALIBRATION_NONE or something else. */
        return FALSE;
    } else if(n_params == 2 && calibration_get_param(cal, CALIBRATION_PARAM_SLOPE) > 0.0) { /* Linear calibration (also CALIBRATION_POLY with n == 2) */
        return TRUE;
    } else if(n_params == 3) { /* Quadratic polynomial, a + b * x + c * x * x, easy enough to solve */
        double b = calibration_get_param(cal, CALIBRATION_PARAM_SLOPE);
        if(b <= 0.0) { /* Derivative (b + 2 * c *x) is going to be negative at least at x = 0 if b is negative. Also in practice b is always > 0! */
            return FALSE;
        }
        double c = calibration_get_param(cal, CALIBRATION_PARAM_QUAD);
        if(c >= 0.0) { /* For all positive x, the derivative (b + 2 * c * x) of this parabola is positive, since b > 0 && c >= 0 */
            return TRUE;
        }
        double peak_x = -b/(2.0*c); /* Parabola "top" or "bottom" at this energy, since (b + 2 * c * peak_x) = 0.0, we will have non-monotonous behaviour if this is somewhere inside [0, n_channels]. Now b is > 0 && c < 0, so peak_x > 0 */
        if(peak_x > (n_channels + 1)) { /* +1 for GSL histogram range safety */
            return TRUE;
        } else {
            return FALSE;
        }
    } else { /* Higher order, use brute force! */
        double y_old = calibration_eval(cal, 0);
        for(size_t i = 1; i <= n_channels; i++) { /* n+1, because histograms have n+1 ranges for n channels */
            double y = calibration_eval(cal, i);
            if(y < y_old) {
                DEBUGMSG("Calibration fails monotonicity check at %zu, where %g keV < %g keV.", i, y / C_KEV, y_old / C_KEV);
                return FALSE;
            }
            y_old = y;
        }
        return TRUE;
    }
}
