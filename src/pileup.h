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

#define PILEUP_GAUSSIAN_ORDER_DEFAULT (7) /* Semi-Gaussian pulse shape is CR-RC^n with this n by default, 0.1 % width is 2.97 times the peaking time (Ortec: 2.9) */
#define PILEUP_GAUSSIAN_ORDER_MAX (16)
#define PILEUP_GAUSSIAN_PEAKING_FACTOR (2.2) /* Peaking time of Gaussian shaper is this times shaping time (any order) */
#define PILEUP_GAUSSIAN_TAIL (1.0e-3) /* Pulses closer than the time it takes for a pulse to fall to this fraction of its maximum are a pair */

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

typedef enum {
    PILEUP_TAIL_SUM = 0, /* One event per pair (trailing pulse is lost), Gaussian: amplitude is the (highest) maximum of the sum (Wielopolski and Gardner) */
    PILEUP_TAIL_FIRST = 1, /* One event per pair, Gaussian: amplitude is held from the first maximum until conversion starts */
    PILEUP_TAIL_BOTH = 2 /* As PILEUP_TAIL_FIRST, but a trailing pulse resolved from the leading one is also recorded (Molodtsov and Gurbich) */
} pileup_tail;

static const jibal_option pileup_tail_option[] = {
        {"sum",   PILEUP_TAIL_SUM},
        {"first", PILEUP_TAIL_FIRST},
        {"both",  PILEUP_TAIL_BOTH},
        {NULL,    0}
};

typedef struct pileup_params {
    pileup_shaper shaper;
    pileup_tail tail;
    double rate; /* Count rate of pulses (1/s), i.e. in real time */
    double rise; /* Trapezoid: rise time (s) */
    double flat; /* Trapezoid: flat top (s), amplitude is sampled at the end of the flat top, i.e. rise + flat after the pulse */
    double shaping; /* Gaussian: shaping time (s), peaking time is PILEUP_GAUSSIAN_PEAKING_FACTOR times this */
    size_t order; /* Gaussian: pulse shape is CR-RC^order */
    double conversion_delay; /* Gaussian: conversion starts this long after the first maximum (peak detection) of the pulse (s) */
    double pair_resolution; /* Pulses closer than this are not resolved (pile-up rejector or trigger filter), zero if not known */
    int reject; /* Pairs resolved by pile-up inspection (pair_resolution > 0) are rejected */
} pileup_params;

typedef struct pileup_result {
    double scale; /* Pulses that are not piled up are scaled by this */
    double fraction; /* Recorded pile-up events per recorded event */
} pileup_result;

int pileup_calculate(const jabs_histogram *h, jabs_histogram *h_pileup, const pileup_params *p, pileup_result *result);
/* Pile-up of pulse pairs (first order). Pulse amplitudes are proportional to channel (bin i is amplitude [i, i+1)).
 * Events from pairs are added to h_pileup (which may have more channels than h), so that the spectrum with pile-up is scale * h + h_pileup.
 * Pulses arrive randomly (Poisson) at the given rate (counts in real time divided by real time). The system is dead (live time clock is stopped)
 * from the arrival of a pulse until its signal has decayed and it has been converted. A pulse arriving in live time is a leading pulse. A trailing
 * pulse arriving dt after the leading pulse within time T forms a pair with it. Fraction exp(-rate * T) of leading pulses are not piled up. Pulses
 * arriving in dead time are lost (dead time losses, corrected by live time) unless they pile up with the leading pulse. A pair gives one event
 * (tail sum or first, or unresolved pulses) or two (tail both). A recorded trailing pulse may be followed by another pulse within T and so on, so
 * without a rejector second events are weighted by 1 / (1 - p), where p = 1 - exp(-rate * T) (the chain is approximated by pairs). With a pile-up
 * rejector (reject and pair_resolution > 0) pairs with pair_resolution <= dt < T are rejected and recorded pairs are weighted by exp(-rate * T), the
 * probability that the pair is not rejected due to a third pulse. The result corresponds to as many leading pulses
 * as there are counts in h, if h is the spectrum in real time the caller should scale the result by live time / real time.
 *
 * Trapezoid: amplitude is sampled at the end of the flat top, T = rise + flat. Leading pulse is A1 + g * A2, where g = 1 for dt <= flat and
 * g = (rise + flat - dt) / rise during the rise of the trailing pulse. Pairs closer than pair_resolution (trigger filter) are one event. With tail both
 * and no rejection, a trailing pulse at least pair_resolution after the leading one triggers its own event, sampled on the falling edge of the
 * leading pulse, A2 + max(1 - dt / rise, 0) * A1.
 *
 * Gaussian (analog semi-Gaussian shaping, peak detecting ADC, as in Wielopolski and Gardner, NIM 133 (1976) 303 and Molodtsov and Gurbich,
 * NIMB 267 (2009) 3484): pulse shape is unipolar CR-RC^n, T is the time from the start of a pulse until it falls to PILEUP_GAUSSIAN_TAIL of maximum.
 * The sum of two pulses has one or two maxima. Tail sum: amplitude is the highest maximum. Tail first: the peak detector holds the first maximum,
 * but follows the signal until conversion starts (conversion_delay after the first maximum), amplitude is the maximum before that. Tail both: the
 * second maximum, if it comes after the start of the conversion, is recorded as another event (tail pile-up, conversion is fast enough not to lose it). */

double pileup_merge_time(const pileup_params *p); /* Time T, see above */

#ifdef __cplusplus
}
#endif
#endif // JABS_PILEUP_H
