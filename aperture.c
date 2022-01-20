#include <stdlib.h>
#include <jibal_units.h>
#include <string.h>
#include "generic.h"
#include "aperture.h"

const char *aperture_name(const aperture *a) {
    if(!a)
        return aperture_option[APERTURE_NONE].s;
    return aperture_option[a->type].s;
}

void aperture_free(aperture *a) {
    free(a);
}

aperture *aperture_default() {
    aperture *a = malloc(sizeof(aperture));
    a->type = APERTURE_NONE;
    a->diameter = 0.0;
    a->width = 0.0;
    a->height = 0.0;
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
    if(!a)
        return 0.0;
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

aperture *aperture_from_argv(const jibal *jibal, int * const argc, char * const ** const argv) {
    aperture *a = aperture_default();
    int found = 0;
    for(const jibal_option *o = aperture_option; o->s; o++) { /* look for a matching aperture_option keyword (circle, rectangle) */
        if(strcmp(o->s, (*argv)[0]) == 0) {
            a->type = o->val;
            (*argv)++;
            (*argc)--;
            found = 1;
            break;
        }
    }

    if(!found) {
        aperture_free(a);
        return NULL;
    }

    while((*argc) >= 2) {
        if(strcmp((*argv)[0], "width") == 0) {
            a->width = jibal_get_val(jibal->units, UNIT_TYPE_DISTANCE, (*argv)[1]);
        } else if(strcmp((*argv)[0], "height") == 0) {
            a->height = jibal_get_val(jibal->units, UNIT_TYPE_DISTANCE, (*argv)[1]);
        } else if(strcmp((*argv)[0], "diameter") == 0) {
            a->diameter = jibal_get_val(jibal->units, UNIT_TYPE_DISTANCE, (*argv)[1]);
        } else {
            break;
        }
        (*argc) -= 2;
        (*argv) += 2;
    }
    return a; /* Does not guarantee the aperture makes sense, unparsed arguments may remain! */
}

aperture *aperture_from_string(const jibal *jibal, const char *str) {
    int argc_orig = 0;
    char **argv_orig = string_to_argv(str, &argc_orig);
    char **argv = argv_orig;
    int argc = argc_orig;
    aperture *a = aperture_from_argv(jibal, &argc, (char *const **) &argv);
    argv_free(argv_orig, argc_orig);
    return a;
}
