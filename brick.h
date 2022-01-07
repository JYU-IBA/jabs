#ifndef JABS_BRICK_H
#define JABS_BRICK_H

#include <gsl/gsl_histogram.h>
#include "sample.h"

typedef struct {
    depth d;
    double E_0;
    double E; /* Energy */
    double S; /* Energy loss straggling (variance) */
    double S_geo_x; /* Geometric straggling */
    double S_geo_y; /* Geometric straggling */
    double Q; /* Counts */
} brick;

void brick_int(double sigma_low, double  sigma_high, double E_low, double E_high, gsl_histogram *h, double Q);
void brick_int2(gsl_histogram *h, const brick *bricks, size_t n_bricks, double S, double scale);

#endif // JABS_BRICK_H
