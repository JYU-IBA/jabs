/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <assert.h>
#include "jabs_debug.h"
#include "nuclear_stopping.h"

nuclear_stopping *nuclear_stopping_new(const jibal_isotope *incident, const jibal_isotope *isotopes) {
    nuclear_stopping *ns = malloc(sizeof(nuclear_stopping));
    if(!ns) {
        return NULL;
    }
    size_t n_isotopes = jibal_isotopes_n(isotopes);
    assert(n_isotopes > 0);
    ns->t = calloc(n_isotopes, sizeof(nucl_stop_pair));
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
    if(!ns) {
        return NULL;
    }
    DEBUGMSG("Shallow copy of ns = %p made (incident: %s), refcount is %i", (void *)ns, ns->incident->name, ns->refcount);
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
        free(ns);
    }
}
