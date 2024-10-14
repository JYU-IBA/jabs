/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_NUCLEAR_STOPPING_H
#define JABS_NUCLEAR_STOPPING_H
#include <stdlib.h>
#include <jibal_masses.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double k;
    double eps0;
#ifdef DEBUG
    const jibal_isotope *target; /* target is not needed, but it is here for debugging purposes */
#endif
} nucl_stop_pair;

typedef struct {
    nucl_stop_pair *t; /* table, indexing refers to ith element in jibal_isotopes (as given to nuclear_stopping_new()) */
    size_t n_isotopes;
    int refcount; /* nuclear_stopping_shared_copy() increases refcount, free() decreases. Actual free on zero. */
    const jibal_isotope *incident;
} nuclear_stopping;

nuclear_stopping *nuclear_stopping_new(const jibal_isotope *incident, const jibal_isotope *isotopes);
nuclear_stopping *nuclear_stopping_shared_copy(nuclear_stopping *ns);
void nuclear_stopping_free(nuclear_stopping *ns);
#ifdef __cplusplus
}
#endif
#endif //JABS_NUCLEAR_STOPPING_H
