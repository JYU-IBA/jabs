#ifndef JABS_BRICK_H
#define JABS_BRICK_H

#include <jibal_masses.h>
#include <gsl/gsl_histogram.h>
#include "sample.h"
#include "detector.h"

typedef struct {
    depth d; /* Depth from sample surface */
    double thick; /* Thickness of brick (difference of depth compared to previous brick) */
    double E_0; /* Incident ion energy at depth */
    double S_0;  /* Incident ion energy loss straggling at depth */
    double E_r; /* Energy after reaction */
    double S_r; /* Energy loss straggling after reaction */
    double E; /* Energy of reaction product (as detected, after detector foil) */
    double S; /* Energy loss straggling (variance) of reaction product */
    double S_geo_x; /* Geometric straggling in "width" direction */
    double S_geo_y; /* Geometric straggling in "height" direction */
    double Q; /* Counts in brick (concentration, cross section, fluence, solid angle etc) */
    double sc; /* Concentration cross-section product. Possibly weighted by straggling and affected by concentration gradients. */
    double sigma; /* Sum (quadrature) of all broadening. */
} brick;

void bricks_calculate_sigma(const detector *det, const jibal_isotope *isotope, brick *bricks, size_t last_brick); /* Sums up all the contributions to sigma (including detector resolution) */
void bricks_convolute(gsl_histogram *h, const brick *bricks, size_t last_brick, double scale, double sigmas_cutoff);
#endif // JABS_BRICK_H
