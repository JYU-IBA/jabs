/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <assert.h>
#include <math.h>
#include "jabs_debug.h"
#include "nuclear_stopping.h"

extern inline double nuclear_stopping_universal(const nuclear_stopping *ns, double epsilon);

#define NUCLEAR_STOPPING_EPSILON_HIGH (30.0) /* Above this, a high energy formula is used */

static void nuclear_stopping_universal_low(double epsilon, double *f, double *df) { /* Universal nuclear stopping and its derivative, low energy formula */
    const double L = log(1.0 + 1.1383 * epsilon);
    const double p = pow(epsilon, 0.21226);
    const double s = sqrt(epsilon);
    const double D = epsilon + 0.01321 * p + 0.19593 * s;
    const double dD = 1.0 + 0.01321 * 0.21226 * p / epsilon + 0.19593 * 0.5 / s;
    *f = L / (2.0 * D);
    *df = (1.1383 / (1.0 + 1.1383 * epsilon) * D - L * dD) / (2.0 * D * D);
}

static void nuclear_stopping_universal_high(double epsilon, double *f, double *df) { /* Universal nuclear stopping and its derivative, high energy formula */
    *f = log(epsilon) / (2.0 * epsilon);
    *df = (1.0 - log(epsilon)) / (2.0 * epsilon * epsilon);
}

double nuclear_stopping_universal_exact(double epsilon) {
    double f, df;
    if(epsilon <= NUCLEAR_STOPPING_EPSILON_HIGH) {
        nuclear_stopping_universal_low(epsilon, &f, &df);
    } else {
        nuclear_stopping_universal_high(epsilon, &f, &df);
    }
    return f;
}

static double *nuclear_stopping_universal_table_new(void) { /* Cubic Hermite polynomials for each segment. Segment ends are exactly representable, NUCLEAR_STOPPING_EPSILON_HIGH is a segment end. */
    double *table = malloc(4 * NUCLEAR_STOPPING_TABLE_SEGMENTS * sizeof(double));
    if(!table) {
        return NULL;
    }
    const int n_per_octave = 1 << NUCLEAR_STOPPING_TABLE_BITS;
    for(int e = NUCLEAR_STOPPING_TABLE_EXP_MIN; e < NUCLEAR_STOPPING_TABLE_EXP_MAX; e++) {
        const double h = ldexp(1.0, e - NUCLEAR_STOPPING_TABLE_BITS); /* Segment width */
        for(int j = 0; j < n_per_octave; j++) {
            const double x0 = ldexp(1.0, e) + j * h;
            const double x1 = x0 + h;
            double f0, df0, f1, df1;
            if(x0 + h / 2.0 > NUCLEAR_STOPPING_EPSILON_HIGH) { /* Both ends of segment are evaluated with the same formula */
                nuclear_stopping_universal_high(x0, &f0, &df0);
                nuclear_stopping_universal_high(x1, &f1, &df1);
            } else {
                nuclear_stopping_universal_low(x0, &f0, &df0);
                nuclear_stopping_universal_low(x1, &f1, &df1);
            }
            double *c = table + 4 * (((e - NUCLEAR_STOPPING_TABLE_EXP_MIN) << NUCLEAR_STOPPING_TABLE_BITS) + j);
            c[0] = f0;
            c[1] = h * df0;
            c[2] = 3.0 * (f1 - f0) - h * (2.0 * df0 + df1);
            c[3] = 2.0 * (f0 - f1) + h * (df0 + df1);
        }
    }
    return table;
}

nuclear_stopping *nuclear_stopping_new(const jibal_isotope *incident, const jibal_isotope *isotopes) {
    nuclear_stopping *ns = malloc(sizeof(nuclear_stopping));
    if(!ns) {
        return NULL;
    }
    size_t n_isotopes = jibal_isotopes_n(isotopes);
    if(n_isotopes == 0) {
        free(ns);
        return NULL;
    }
    ns->t = calloc(n_isotopes, sizeof(nucl_stop_pair));
    ns->universal = nuclear_stopping_universal_table_new();
    if(!ns->t || !ns->universal) {
        free(ns->t);
        free(ns->universal);
        free(ns);
        return NULL;
    }
    size_t i = 0;
    ns->incident = incident;
    int Z1 = incident->Z;
    double m1 = incident->mass;
    for(const jibal_isotope *target = isotopes; target->A != 0 && i < n_isotopes; target++) {

        int Z2 = target->Z;
        double m2 = target->mass;
        double zpower = 1.0/(pow(Z1, 0.23) + pow(Z2, 0.23));
        double a_u = 0.8854*C_BOHR_RADIUS * zpower;
        double gamma = 4.0 * m1 * m2/pow(m1 + m2, 2.0);
        double eps0 = ((1.0/C_KEV) * 32.53 * m2/(Z1 * Z2 * (m1 + m2))) * zpower;
        double k = C_PI * (a_u * a_u) * gamma/eps0;
        ns->t[i].eps0 = eps0;
        ns->t[i].k = k;
#ifdef DEBUG
        ns->t[i].target = target;
#endif
        i++;
    }
    ns->n_isotopes = n_isotopes;
    ns->refcount = 1;
    return ns;
}

nuclear_stopping *nuclear_stopping_shared_copy(nuclear_stopping *ns) {
    assert(ns);
    if(!ns) {
        return NULL;
    }
    DEBUGMSG("Shallow copy of ns = %p made (incident: %s), refcount is %i", (void *)ns, ns->incident->name, ns->refcount);
    assert(ns->refcount > 0);
    ns->refcount += 1;
    return ns;
}

void nuclear_stopping_free(nuclear_stopping *ns) {
    if(!ns) {
        return;
    }
    ns->refcount--;
    DEBUGVERBOSEMSG("nuclear_stopping_free(ns = %p), incident %s, called, refcount is now %i", (void *)ns, ns->incident->name, ns->refcount);
    assert(ns->refcount >= 0);
    if(ns->refcount == 0) {
        free(ns->t);
        free(ns->universal);
        free(ns);
    }
}
