/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_PROB_DIST_H
#define JABS_PROB_DIST_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef struct {
    double x;
    double p;
} prob_point;

typedef struct {
    size_t n;
    prob_point *points;
} prob_dist; /* Discrete probability distribution */


prob_dist *prob_dist_alloc(size_t n);
void prob_dist_free(prob_dist *pd);
prob_dist *prob_dist_gaussian(size_t n);
#endif //JABS_PROB_DIST_H
