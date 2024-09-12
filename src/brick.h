/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_BRICK_H
#define JABS_BRICK_H

#include <jibal_masses.h>
#include "sample.h"
#include "detector.h"
#include "histogram.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    int valid;
    depth d; /* Depth from sample surface */
    double thick; /* Thickness of brick (difference of depth compared to previous brick) */
    double E_0; /* Incident ion energy at depth */
    double S_0;  /* Incident ion energy loss straggling at depth */
    double E_r; /* Energy after reaction */
    double S_r; /* Energy loss straggling after reaction */
    double E_s; /* Energy at surface (after reaction) */
    double S_s; /* Energy loss straggling (after reaction) */
    double E; /* Energy of reaction product (as detected, after detector foil) */
    double S; /* Energy loss straggling (variance) of reaction product */
    double S_geo_x; /* Geometric straggling in "width" direction */
    double S_geo_y; /* Geometric straggling in "height" direction */
    double Q; /* Counts in brick (concentration, cross section, fluence, solid angle etc) */
    double sc; /* Concentration cross-section product. Possibly weighted by straggling and affected by concentration gradients. */
    double S_sum; /* Sum (quadrature) of all broadening. */
    double deriv; /* dE_0/dE_r, estimate. Multiply this by wanted E_r change to get E_0 change  */
    double dE; /* Energy change in brick */
    double effective_stopping;
} brick;

void bricks_calculate_sigma(const detector *det, const jibal_isotope *isotope, brick *bricks, size_t last_brick); /* Sums up all the contributions to sigma (including detector resolution) */
void bricks_convolute(jabs_histogram *h, const calibration *c, const brick *bricks, size_t last_brick, double scale, double sigmas_cutoff, double emin, int accurate);
#ifdef __cplusplus
}
#endif
#endif // JABS_BRICK_H
