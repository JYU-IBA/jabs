#include "ion.h"



void ion_set_isotope(ion *ion, const jibal_isotope *isotope) {
    if(!isotope)
        return;
    ion->isotope = isotope;
    ion->mass = isotope->mass;
    ion->Z = isotope->Z;
}

void ion_set_angle(ion *ion, double angle) {
    if(!ion)
        return;
    ion->angle = angle;
    ion->cosine = cos(angle);
    ion->inverse_cosine = 1.0/ion->cosine;
}

