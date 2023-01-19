/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

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

typedef struct {
    double k;
    double eps0;
} nucl_stop_pair;

typedef struct {
    nucl_stop_pair *nucl_stop;
    size_t nucl_stop_isotopes;
    int refcount;
} nuclear_stopping;

nuclear_stopping *nuclear_stopping_new(const jibal_isotope *incident, const jibal_isotope *isotopes);
#endif //JABS_NUCLEAR_STOPPING_H
