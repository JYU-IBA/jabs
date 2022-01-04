#include <stdlib.h>
#include <jibal_units.h>
#include "aperture.h"

const char *aperture_name(const aperture *a) {
    return aperture_option[a->type].s;
}

aperture aperture_default() {
    aperture a;
    a.type = APERTURE_NONE;
    a.diameter = 0.0;
    a.width = 0.0;
    a.height = 0.0;
    return a;
}

double aperture_shape_parameter(const aperture *a) {
    static const double shape_circle = 0.5 * C_FWHM/2.0;
    static const double shape_rect = 0.5 * C_FWHM/1.7320508075688772;
    if(a->type == APERTURE_CIRCLE) {
        return shape_circle;
    }
    if(a->type == APERTURE_RECTANGLE) {
        return shape_rect;
    }
    return 0.0;
}

double aperture_width_shape_product(const aperture *a, const char direction) {
    double shape = aperture_shape_parameter(a);
    if(a->type == APERTURE_CIRCLE) {
        return shape * a->diameter;
    }
    if(a->type == APERTURE_RECTANGLE) {
        if(direction == 'x') {
            return shape * a->width;
        }
        if(direction == 'y') {
            return shape * a->height;
        }
        return 0.0;
    }
    return 0.0;
}
