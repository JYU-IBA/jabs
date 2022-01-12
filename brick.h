#ifndef JABS_BRICK_H
#define JABS_BRICK_H

#include <gsl/gsl_histogram.h>
#include "sample.h"

typedef struct {
    depth d; /* Depth from sample surface */
    double thick; /* Thickness of brick (difference of depth compared to previous brick) */
    double E_0; /* Incident ion energy at depth */
    double S_0;  /* Incident ion energy loss straggling at depth */
    double E; /* Energy of reaction product (as detected, after detector foil) */
    double S; /* Energy loss straggling (variance) of reaction product */
    double S_geo_x; /* Geometric straggling in "width" direction */
    double S_geo_y; /* Geometric straggling in "height" direction */
    double Q; /* Counts in brick (concentration, cross section, fluence, solid angle etc) */
    double sc; /* Concentration cross-section product. Possibly weighted by straggling and affected by concentration gradients. */
} brick;

void brick_int(double sigma_low, double  sigma_high, double E_low, double E_high, gsl_histogram *h, double Q);
void brick_int2(gsl_histogram *h, const brick *bricks, size_t n_bricks, double S, double scale);

#endif // JABS_BRICK_H
