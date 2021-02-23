/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <assert.h>
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

double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, const jibal_isotope *isotopes) {
    size_t i = (isotope-isotopes);
    assert(i < ion->nucl_stop_isotopes);
    const double epsilon = ion->nucl_stop[i].eps0 * ion->E;
#ifdef ACCURATE_NUCLEAR_STOPPING
    if(epsilon <= 30.0) {
        return ion->nucl_stop[i].k*log(1+1.1383*epsilon)/(2*(epsilon+0.01321*pow(epsilon, 0.21226)+0.19593*pow(epsilon, 0.5)));
#else
    if(epsilon <= M_E) {
        return ion->nucl_stop[i].k/(2.0*M_E); /* Below the maximum (of the else-branch) return maximum. The approximation is bad, but better than nothing. */
#endif
    } else {
        return ion->nucl_stop[i].k*log(epsilon)/(2.0*epsilon);
    }
}
void ion_nuclear_stop_fill_params(ion *ion, const jibal_isotope *isotopes, int n_isotopes) {
    const jibal_isotope *isotope;
    int i = 0;
    ion->nucl_stop = malloc(n_isotopes * sizeof (nucl_stop_pair ));
    for(isotope=isotopes; isotope->A != 0 && i < n_isotopes; isotope++) {
        int Z1 = ion->isotope->Z;
        int Z2 = isotope->Z;
        double m1 = ion->isotope->mass;
        double m2 = isotope->mass;
        double a_u=0.8854*C_BOHR_RADIUS/(pow(Z1, 0.23)+pow(Z2, 0.23));
        double gamma = 4.0*m1*m2/pow(m1+m2, 2.0);
        double eps0 = (1.0/C_KEV)*32.53*m2/(Z1*Z2*(m1+m2)*(pow(Z1, 0.23)+pow(Z2, 0.23)));
        double k = C_PI * (a_u * a_u) * gamma/eps0;
        ion->nucl_stop[i].eps0 = eps0;
        ion->nucl_stop[i].k = k;
        i++;
    }
    ion->nucl_stop_isotopes = n_isotopes; /* This is somewhat redundant, but it is there to remind us. */
}
