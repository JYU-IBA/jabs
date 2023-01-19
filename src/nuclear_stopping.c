/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include "nuclear_stopping.h"

nuclear_stopping *nuclear_stopping_new(const jibal_isotope *incident, const jibal_isotope *isotopes) {
    nuclear_stopping *ns = malloc(sizeof(nuclear_stopping));
    if(!ns) {
        return NULL;
    }
    size_t n_isotopes = jibal_isotopes_n(isotopes);
    size_t i = 0;
    for(const jibal_isotope *target = isotopes; target->A != 0 && i < n_isotopes; target++) {
        int Z1 = incident->Z;
        int Z2 = target->Z;
        double m1 = incident->mass;
        double m2 = target->mass;
        double a_u=0.8854*C_BOHR_RADIUS/(pow(Z1, 0.23)+pow(Z2, 0.23));
        double gamma = 4.0*m1*m2/pow(m1+m2, 2.0);
        double eps0 = (1.0/C_KEV)*32.53*m2/(Z1*Z2*(m1+m2)*(pow(Z1, 0.23)+pow(Z2, 0.23)));
        double k = C_PI * (a_u * a_u) * gamma/eps0;
        ns->nucl_stop[i].eps0 = eps0;
        ns->nucl_stop[i].k = k;
        i++;
    }
    ns->nucl_stop_isotopes = n_isotopes;
    ns->refcount = 1;
    return ns;
}

nuclear_stopping *nuclear_stopping_shallow_copy(nuclear_stopping *ns) {
    if(!ns) {
        return NULL;
    }
#ifdef DEBUG
    fprintf(stderr, "Shallow copy of ns = %p made, refcount is %zu\n", (void *)ns, ns->refcount);
#endif
    ns->refcount += 1;
    return ns;
}
