/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_PILEUP_H
#define JABS_PILEUP_H

#include <stdio.h>
#include <jibal_option.h>
#include "histogram.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PILEUP_GAUSSIAN_ORDER (7) /* Semi-Gaussian pulse shape is CR-RC^n with this n, 0.1 % width is 2.97 times the peaking time (Ortec: 2.9) */
#define PILEUP_GAUSSIAN_PEAKING_FACTOR (2.2) /* Peaking time of Gaussian shaper is this times shaping time */
#define PILEUP_GAUSSIAN_TAIL (1.0e-3) /* Pulses closer than the time it takes for a pulse to fall to this fraction of its maximum are merged */

typedef enum {
    PILEUP_SHAPER_NONE = 0,
    PILEUP_SHAPER_TRAPEZOID = 1,
    PILEUP_SHAPER_GAUSSIAN = 2
} pileup_shaper;

static const jibal_option pileup_shaper_option[] = {
        {JIBAL_OPTION_STR_NONE, PILEUP_SHAPER_NONE},
        {"trapezoid",           PILEUP_SHAPER_TRAPEZOID},
        {"gaussian",            PILEUP_SHAPER_GAUSSIAN},
        {NULL,                  0}
};

typedef struct pileup_params {
    pileup_shaper shaper;
    double rate; /* Count rate of pulses (1/s) */
    double rise; /* Trapezoid: rise time (s) */
    double flat; /* Trapezoid: flat top (s), amplitude is sampled at the end of the flat top, i.e. rise + flat after the pulse */
    double shaping; /* Gaussian: shaping time (s), peaking time is PILEUP_GAUSSIAN_PEAKING_FACTOR times this */
    double pair_resolution; /* Pair resolution time of pile-up rejector (s), zero if there is no rejector */
} pileup_params;

typedef struct pileup_result {
    double scale; /* Pulses that are not piled up are scaled by this */
    double fraction; /* Recorded pile-up events per recorded event */
} pileup_result;

int pileup_calculate(const jabs_histogram *h, jabs_histogram *h_pileup, const pileup_params *p, pileup_result *result);
/* Pile-up of pulse pairs (first order). Pulse amplitudes are proportional to channel (bin i is amplitude [i, i+1)).
 * Pile-up events (pairs) are added to h_pileup (which may have more channels than h), so that the spectrum with pile-up is scale * h + h_pileup.
 * Pulses arrive randomly (Poisson) at the given rate (counts in live time divided by live time). A trailing pulse arriving dt after the leading pulse
 * within the merging time T is merged with it (one recorded event). The trailing pulse arrives while the system is busy, so it is a dead time loss,
 * which is corrected by live time. Every pulse in live time is therefore a leading pulse and without a pile-up rejector the number of counts is
 * conserved. With a pile-up rejector pairs with pair_resolution <= dt < T are rejected, these losses are not corrected.
 *
 * Trapezoid: amplitude is sampled at the end of the flat top, T = rise + flat. Amplitude A1 + g * A2, where g = 1 for dt <= flat and
 * g = (rise + flat - dt) / rise during the rise of the trailing pulse.
 *
 * Gaussian (analog semi-Gaussian shaping, peak detecting ADC, as in Wielopolski and Gardner, NIM 133 (1976) 303): pulse shape is unipolar CR-RC^n,
 * amplitude is the maximum of the sum of the two pulses. T is the time from the start of a pulse until it falls to PILEUP_GAUSSIAN_TAIL of maximum. */

double pileup_merge_time(const pileup_params *p); /* Merging time T, see above */

#ifdef __cplusplus
}
#endif
#endif // JABS_PILEUP_H
