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
#include <math.h>
#include "ion.h"
#include "rotate.h"

void ion_reset(ion *ion) {
    ion->isotope = NULL;
    ion->E = 0.0;
    ion->S = 0.0;
    ion->Z = 0;
    ion->mass = 0.0;
    ion->theta = 0.0;
    ion->phi = 0.0;
    ion->cosine_theta = 1.0;
    ion->cosine_phi = 1.0;
    ion->inverse_cosine_theta = 1.0; /* These need to be set matching to the angles */
    ion->inverse_cosine_phi = 1.0;
    ion->nucl_stop = NULL; /* Be careful, memory could leak */
}

void ion_set_isotope(ion *ion, const jibal_isotope *isotope) {
    if(!isotope)
        return;
    ion->isotope = isotope;
    ion->mass = isotope->mass;
    ion->mass_inverse = 1.0/isotope->mass;
    ion->Z = isotope->Z;
}

void ion_set_angle(ion *ion, double theta, double phi) {
    if(!ion)
        return;
    if(theta == ion->theta && phi == ion->phi) {
        return;
    }
    ion->theta = fmod(theta, C_2PI);
    if(ion->theta > C_PI) { /* Limit theta to [0, pi] */
        ion->theta = C_2PI - ion->theta;
        phi += C_PI;
    }
    ion->phi = fmod(phi, C_2PI);
    ion->cosine_theta = cos(ion->theta);
    ion->cosine_phi = cos(ion->phi);
    ion->inverse_cosine_theta = 1.0/ion->cosine_theta;
    ion->inverse_cosine_phi = 1.0/ion->cosine_phi;
}

double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, int accurate) {
#ifdef DEBUG
    assert(isotope->i < ion->nucl_stop->n_isotopes);
    assert(ion->nucl_stop->t[isotope->i].target == isotope);
#endif
    nucl_stop_pair x = ion->nucl_stop->t[isotope->i];
    const double epsilon = x.eps0 * ion->E;
    if(accurate) {
        if(epsilon <= 30.0) {
            return x.k * log(1 + 1.1383 * epsilon) / (2 * (epsilon + 0.01321 * pow(epsilon, 0.21226) + 0.19593 * sqrt(epsilon)));
        }
    } else {
        if (epsilon <= M_E) {
            return x.k / (2.0 * M_E); /* Below the maximum (of the else-branch) return maximum. The approximation is bad, but better than nothing. */
        }
    }
    return (x.k * log(epsilon) / (2.0 * epsilon));
}

void ion_rotate(ion *ion, double theta2, double phi2) { /* Wrapper for rotate() */
    double theta, phi;
    rotate(theta2, phi2, ion->theta, ion->phi, &theta, &phi);
    ion_set_angle(ion, theta, phi);
}

void ion_print(FILE *f, const ion *ion) {
    fprintf(f, "ion %s (Z=%i, mass=%.3lf u), E = %.3lf keV, theta = %.3lf deg (cos_theta = %.3lf, 1/cos_theta = %.3lf), phi = %.3lf deg. Nuclear stopping isotopes %zu calculated for %s.\n",
            ion->isotope->name, ion->Z, ion->mass/C_U, ion->E/C_KEV,  ion->theta/C_DEG, ion->cosine_theta, ion->inverse_cosine_theta, ion->phi/C_DEG, ion->nucl_stop->n_isotopes, ion->nucl_stop->incident->name);
}
