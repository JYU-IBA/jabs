/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <math.h>
#include <gsl/gsl_math.h>
#include <jibal_phys.h>
#include "pileup.h"

#define PILEUP_GAUSSIAN_TABLE_NQ (64) /* Steps of fraction q = A2 / (A1 + A2) of amplitude in trailing pulse */
#define PILEUP_GAUSSIAN_TABLE_DU (0.0125) /* Step of time between pulses (units of peaking time) */
#define PILEUP_GAUSSIAN_SEARCH_DT (0.01) /* Step when searching for maxima of sum of pulses (units of peaking time) */
#define PILEUP_GAUSSIAN_TAIL_START (2.0) /* Time between pulses (units of peaking time) after which nodes are PILEUP_GAUSSIAN_TAIL_STEP grid steps apart, sum changes slowly */
#define PILEUP_GAUSSIAN_TAIL_STEP (4)
#define PILEUP_DEPOSIT_DK_MIN (0.25) /* Amplitude (channels) of a pair changes at least this much between nodes, shorter intervals are merged */
#define PILEUP_DEPOSIT_WIDTH_MIN (1.0e-6) /* Narrower uniform distributions of amplitude are deposited as points */

typedef struct pileup_gaussian_table { /* Maxima of (1 - q) * pulse(t) + q * pulse(t - u), index i_q * n_u + i_u. Times are in units of peaking time. */
    size_t order;
    double t_tail; /* Time from start of pulse until it falls to PILEUP_GAUSSIAN_TAIL of maximum */
    size_t n_u;
    double *t1; /* Time of first maximum */
    double *m1; /* First maximum */
    double *t2; /* Time of second maximum, negative if there is only one */
    double *m2; /* Second maximum */
} pileup_gaussian_table;

static pileup_gaussian_table *gaussian_tables[PILEUP_GAUSSIAN_ORDER_MAX + 1];

static double pileup_gaussian_pulse(size_t order, double t) { /* CR-RC^n, t in units of peaking time, maximum is 1 at t = 1 */
    if(t <= 0.0) {
        return 0.0;
    }
    return exp((double) order * (log(t) + 1.0 - t));
}

static double pileup_gaussian_sum(size_t order, double q, double u, double t) {
    return (1.0 - q) * pileup_gaussian_pulse(order, t) + q * pileup_gaussian_pulse(order, t - u);
}

static void pileup_gaussian_sum_maxima(size_t order, double q, double u, double *t1, double *m1, double *t2, double *m2) {
    const size_t n = (size_t) ceil((u + 1.0) / PILEUP_GAUSSIAN_SEARCH_DT) + 2; /* Both pulses are falling after u + 1 */
    size_t n_max = 0;
    double s_prev = 0.0;
    double s = pileup_gaussian_sum(order, q, u, PILEUP_GAUSSIAN_SEARCH_DT);
    *t1 = 0.0;
    *m1 = 0.0;
    *t2 = -1.0;
    *m2 = 0.0;
    for(size_t i = 1; i < n && n_max < 2; i++) {
        const double t = i * PILEUP_GAUSSIAN_SEARCH_DT;
        const double s_next = pileup_gaussian_sum(order, q, u, t + PILEUP_GAUSSIAN_SEARCH_DT);
        if(s > s_prev && s >= s_next) { /* Local maximum, parabola through the sample and its neighbours */
            const double curv = s_prev - 2.0 * s + s_next;
            double t_max = t, s_max = s;
            if(curv < 0.0) {
                t_max = t - PILEUP_GAUSSIAN_SEARCH_DT * (s_next - s_prev) / (2.0 * curv);
                s_max = s - pow2(s_next - s_prev) / (8.0 * curv);
            }
            if(n_max == 0) {
                *t1 = t_max;
                *m1 = s_max;
            } else {
                *t2 = t_max;
                *m2 = s_max;
            }
            n_max++;
        }
        s_prev = s;
        s = s_next;
    }
}

static pileup_gaussian_table *pileup_gaussian_table_init(size_t order) {
    pileup_gaussian_table *tab = calloc(1, sizeof(pileup_gaussian_table));
    if(!tab) {
        return NULL;
    }
    tab->order = order;
    double low = 1.0, high = 100.0;
    for(int i = 0; i < 100; i++) { /* Bisection, pulse is decreasing after peak */
        const double mid = (low + high) / 2.0;
        if(pileup_gaussian_pulse(order, mid) > PILEUP_GAUSSIAN_TAIL) {
            low = mid;
        } else {
            high = mid;
        }
    }
    tab->t_tail = (low + high) / 2.0;
    tab->n_u = (size_t) ceil(tab->t_tail / PILEUP_GAUSSIAN_TABLE_DU) + 2;
    const size_t n = (PILEUP_GAUSSIAN_TABLE_NQ + 1) * tab->n_u;
    tab->t1 = malloc(n * sizeof(double));
    tab->m1 = malloc(n * sizeof(double));
    tab->t2 = malloc(n * sizeof(double));
    tab->m2 = malloc(n * sizeof(double));
    if(!tab->t1 || !tab->m1 || !tab->t2 || !tab->m2) {
        free(tab->t1);
        free(tab->m1);
        free(tab->t2);
        free(tab->m2);
        free(tab);
        return NULL;
    }
    for(size_t i_q = 0; i_q <= PILEUP_GAUSSIAN_TABLE_NQ; i_q++) {
        const double q = (1.0 * i_q) / PILEUP_GAUSSIAN_TABLE_NQ;
        for(size_t i_u = 0; i_u < tab->n_u; i_u++) {
            const size_t k = i_q * tab->n_u + i_u;
            pileup_gaussian_sum_maxima(order, q, i_u * PILEUP_GAUSSIAN_TABLE_DU, &tab->t1[k], &tab->m1[k], &tab->t2[k], &tab->m2[k]);
        }
    }
    return tab;
}

static const pileup_gaussian_table *pileup_gaussian_table_get(size_t order) {
    if(order < 1 || order > PILEUP_GAUSSIAN_ORDER_MAX) {
        return NULL;
    }
    const pileup_gaussian_table *tab;
#pragma omp critical(pileup_gaussian_table)
    {
        if(!gaussian_tables[order]) {
            gaussian_tables[order] = pileup_gaussian_table_init(order);
        }
        tab = gaussian_tables[order];
    }
    return tab;
}

double pileup_merge_time(const pileup_params *p) {
    switch(p->shaper) {
        case PILEUP_SHAPER_TRAPEZOID:
            return p->rise + p->flat;
        case PILEUP_SHAPER_GAUSSIAN: {
            const pileup_gaussian_table *tab = pileup_gaussian_table_get(p->order);
            return tab ? tab->t_tail * PILEUP_GAUSSIAN_PEAKING_FACTOR * p->shaping : 0.0;
        }
        default:
            return 0.0;
    }
}

typedef struct pileup_deposit { /* Uniform distributions of amplitude are added as ramps to the cumulative distribution F, F(x) is the probability of amplitude below x */
    size_t n; /* Number of bins */
    double *jump; /* F(x) - F(x - 1) gets a contribution jump[x]... */
    double *slope; /* ...and a contribution of slope[x'] for all x' <= x */
} pileup_deposit;

static void pileup_deposit_ramp(pileup_deposit *d, double a, double c) { /* Adds c * max(x - a, 0) to F */
    const size_t x0 = a <= 0.0 ? 0 : (size_t) ceil(a);
    if(x0 > d->n) {
        return;
    }
    d->jump[x0] += c * (x0 - a);
    d->slope[x0 + 1] += c;
}

static void pileup_deposit_uniform(pileup_deposit *d, double low, double high, double w) { /* Adds probability w distributed uniformly between low and high */
    if(high - low < PILEUP_DEPOSIT_WIDTH_MIN) { /* Point at low, to bin floor(low) */
        const size_t x0 = low < 0.0 ? 0 : (size_t) floor(low) + 1;
        if(x0 <= d->n) {
            d->jump[x0] += w;
        }
        return;
    }
    const double c = w / (high - low);
    pileup_deposit_ramp(d, low, c);
    pileup_deposit_ramp(d, high, -c);
}

static void pileup_deposit_curve(pileup_deposit *d, double S, const double *v, const double *w, size_t n_nodes, double c) {
    /* Amplitude S * v[l] at nodes, linear between nodes, probability c * w[l] between nodes l and l + 1. Negative v: no event. Intervals are merged
     * until amplitude changes by PILEUP_DEPOSIT_DK_MIN. Amplitudes are spread within bins, this is approximated by two copies shifted by +-0.25. */
    size_t l_start = 0;
    double w_sum = 0.0;
    int open = FALSE;
    for(size_t l = 0; l + 1 < n_nodes; l++) {
        if(v[l] < 0.0 || v[l + 1] < 0.0) {
            continue;
        }
        if(!open) {
            l_start = l;
            w_sum = 0.0;
            open = TRUE;
        }
        w_sum += w[l];
        const int last = (l + 2 == n_nodes) || v[l + 2] < 0.0;
        if(!last && S * fabs(v[l + 1] - v[l_start]) < PILEUP_DEPOSIT_DK_MIN) {
            continue;
        }
        const double low = S * GSL_MIN_DBL(v[l_start], v[l + 1]);
        const double high = S * GSL_MAX_DBL(v[l_start], v[l + 1]);
        pileup_deposit_uniform(d, low - 0.25, high - 0.25, 0.5 * c * w_sum);
        pileup_deposit_uniform(d, low + 0.25, high + 0.25, 0.5 * c * w_sum);
        open = FALSE;
    }
}

static void pileup_trapezoid_fraction(double *Q, size_t n, const jabs_histogram *h, double N, double r, double dt0, double dt1, double f0, double f1, double weight) {
    /* Adds distribution of f * A, where A is the amplitude of a pulse, to Q (multiplied by weight). Pulse arrives between dt0 and dt1, f changes linearly from f0 to f1. */
    if(dt1 <= dt0) {
        return;
    }
    const double f_low = GSL_MIN_DBL(f0, f1), f_high = GSL_MAX_DBL(f0, f1);
    for(size_t j = 0; j < n; j++) {
        const double m = weight * h->bin[j] / N;
        if(m == 0.0) {
            continue;
        }
        const double x = j + 0.5;
        if((f_high - f_low) * x < PILEUP_DEPOSIT_WIDTH_MIN) {
            Q[GSL_MIN((size_t) floor(f0 * x), n - 1)] += m * exp(-r * dt0) * -expm1(-r * (dt1 - dt0));
            continue;
        }
        for(size_t k = (size_t) floor(f_low * x); k < n && k <= f_high * x; k++) {
            const double fa = GSL_MAX_DBL(k / x, f_low), fb = GSL_MIN_DBL((k + 1) / x, f_high);
            if(fb <= fa) {
                continue;
            }
            const double ta = dt0 + (fa - f0) / (f1 - f0) * (dt1 - dt0), tb = dt0 + (fb - f0) / (f1 - f0) * (dt1 - dt0);
            const double t_low = GSL_MIN_DBL(ta, tb), t_high = GSL_MAX_DBL(ta, tb);
            Q[k] += m * exp(-r * t_low) * -expm1(-r * (t_high - t_low)); /* Probability that pulse arrives between t_low and t_high */
        }
    }
}

static int pileup_trapezoid(const jabs_histogram *h, size_t n, double N, jabs_histogram *h_pileup, const pileup_params *p, double dt_max, double weight, double weight2) {
    const double r = p->rate;
    const double t_merge = p->rise + p->flat;
    const double t_res = p->pair_resolution > 0.0 ? GSL_MIN_DBL(p->pair_resolution, t_merge) : t_merge; /* Pairs closer than this are not resolved */
    const int reject = dt_max < t_merge; /* Pairs further apart than dt_max are rejected, leading pulse is recorded if trailing pulse arrives before dt_max */
    double *Q = calloc(n, sizeof(double)); /* Distribution of amplitude added to a pulse by the other pulse of the pair, normalized to probability */
    if(!Q) {
        return EXIT_FAILURE;
    }
    /* Leading pulse, A1 + g * A2. Trailing pulse within flat top (g = 1) or during its rise (g decreases linearly to zero at t_merge). */
    pileup_trapezoid_fraction(Q, n, h, N, r, 0.0, GSL_MIN_DBL(p->flat, dt_max), 1.0, 1.0, weight);
    if(p->rise > 0.0 && dt_max > p->flat) {
        pileup_trapezoid_fraction(Q, n, h, N, r, p->flat, dt_max, 1.0, (t_merge - dt_max) / p->rise, weight);
    }
    if(p->tail == PILEUP_TAIL_BOTH && !reject && t_res < t_merge) { /* Resolved trailing pulse, A2 + max(1 - dt / rise, 0) * A1 */
        if(p->rise > 0.0 && p->rise > t_res) {
            const double dt_end = GSL_MIN_DBL(p->rise, t_merge);
            pileup_trapezoid_fraction(Q, n, h, N, r, t_res, dt_end, 1.0 - t_res / p->rise, 1.0 - dt_end / p->rise, weight * weight2);
        }
        const double dt_clean = GSL_MAX_DBL(t_res, p->rise); /* Leading pulse has returned to zero, trailing pulse is recorded as is */
        if(dt_clean < t_merge) {
            const double w = weight * weight2 * exp(-r * dt_clean) * -expm1(-r * (t_merge - dt_clean));
            for(size_t i = 0; i < n && i < h_pileup->n; i++) {
                h_pileup->bin[i] += w * h->bin[i];
            }
        }
    }
    for(size_t i = 0; i < n; i++) {
        const double h_i = 0.5 * h->bin[i];
        if(h_i == 0.0) {
            continue;
        }
        for(size_t k = 0; k < n && i + k < h_pileup->n; k++) { /* Sum of uniform amplitudes in bins i and k is a triangle centered at i + k + 1, half in each bin */
            const double c = h_i * Q[k];
            h_pileup->bin[i + k] += c;
            if(i + k + 1 < h_pileup->n) {
                h_pileup->bin[i + k + 1] += c;
            }
        }
    }
    free(Q);
    return EXIT_SUCCESS;
}

static double pileup_gaussian_v2(const double *v_low, const double *v_high, size_t i_u, double f_u, double f_q) {
    /* Interpolates second event amplitude (negative: no event) in u and q. If only one row of q has the event, the nearest row is used. */
    double a = v_low[i_u], b = v_high[i_u];
    if(f_u > 0.0) {
        a = (a < 0.0 || v_low[i_u + 1] < 0.0) ? -1.0 : a + f_u * (v_low[i_u + 1] - a);
        b = (b < 0.0 || v_high[i_u + 1] < 0.0) ? -1.0 : b + f_u * (v_high[i_u + 1] - b);
    }
    if(a >= 0.0 && b >= 0.0) {
        return a + f_q * (b - a);
    }
    if(a >= 0.0 && f_q < 0.5) {
        return a;
    }
    if(b >= 0.0 && f_q >= 0.5) {
        return b;
    }
    return -1.0;
}

static int pileup_gaussian(const jabs_histogram *h, size_t n, double N, jabs_histogram *h_pileup, const pileup_params *p, double dt_max, double weight, double weight2) {
    const pileup_gaussian_table *tab = pileup_gaussian_table_get(p->order);
    if(!tab) {
        return EXIT_FAILURE;
    }
    const double t_peak = PILEUP_GAUSSIAN_PEAKING_FACTOR * p->shaping;
    const double u_max = dt_max / t_peak; /* Pairs closer than this (units of peaking time) are recorded */
    const double delay = p->conversion_delay / t_peak;
    const size_t i_u_max = (size_t) floor(u_max / PILEUP_GAUSSIAN_TABLE_DU);
    if(i_u_max + 2 > tab->n_u) {
        return EXIT_FAILURE;
    }
    const size_t n_tab = (PILEUP_GAUSSIAN_TABLE_NQ + 1) * tab->n_u;
    const size_t n_nodes_max = i_u_max + 2;
    double *v1_tab = malloc(n_tab * sizeof(double)); /* Amplitude of first event (units of A1 + A2) */
    double *v2_tab = malloc(n_tab * sizeof(double)); /* Amplitude of second event, negative if there is no second event */
    size_t *i_u = malloc(n_nodes_max * sizeof(size_t)); /* Nodes are on the grid of the table (except the last), amplitudes are linear in time between nodes */
    double *w = malloc(n_nodes_max * sizeof(double)); /* Probability of trailing pulse between nodes l and l + 1 */
    double *v1 = malloc(n_nodes_max * sizeof(double));
    double *v2 = malloc(n_nodes_max * sizeof(double));
    size_t *nonzero = malloc(n * sizeof(size_t));
    pileup_deposit d = {.n = h_pileup->n, .jump = calloc(h_pileup->n + 2, sizeof(double)), .slope = calloc(h_pileup->n + 2, sizeof(double))};
    int status = EXIT_FAILURE;
    if(!v1_tab || !v2_tab || !i_u || !w || !v1 || !v2 || !nonzero || !d.jump || !d.slope) {
        goto cleanup;
    }
    for(size_t i_q = 0; i_q <= PILEUP_GAUSSIAN_TABLE_NQ; i_q++) {
        const double q = (1.0 * i_q) / PILEUP_GAUSSIAN_TABLE_NQ;
        for(size_t k_u = 0; k_u < tab->n_u; k_u++) {
            const size_t k = i_q * tab->n_u + k_u;
            v2_tab[k] = -1.0;
            if(tab->t2[k] < 0.0 || p->tail == PILEUP_TAIL_SUM || tab->t1[k] + delay >= tab->t2[k]) { /* Highest maximum */
                v1_tab[k] = GSL_MAX_DBL(tab->m1[k], tab->m2[k]);
                continue;
            }
            /* Two maxima, conversion starts before the second one. Peak detector holds the first maximum or follows the rising second pulse. */
            v1_tab[k] = GSL_MAX_DBL(tab->m1[k], pileup_gaussian_sum(tab->order, q, k_u * PILEUP_GAUSSIAN_TABLE_DU, tab->t1[k] + delay));
            if(p->tail == PILEUP_TAIL_BOTH) {
                v2_tab[k] = tab->m2[k];
            }
        }
    }
    size_t n_nodes = 0;
    for(size_t l = 0; l <= i_u_max; l += (l * PILEUP_GAUSSIAN_TABLE_DU < PILEUP_GAUSSIAN_TAIL_START) ? 1 : PILEUP_GAUSSIAN_TAIL_STEP) {
        i_u[n_nodes++] = l;
    }
    const double u_last = u_max / PILEUP_GAUSSIAN_TABLE_DU; /* Last node, generally not on the grid */
    const double f_u_last = u_last - i_u_max;
    if(f_u_last > 0.0 || i_u[n_nodes - 1] != i_u_max) {
        i_u[n_nodes++] = i_u_max; /* Last node, value is interpolated from grid points i_u_max and i_u_max + 1 */
    }
    for(size_t l = 0; l + 1 < n_nodes; l++) {
        const double u_low = i_u[l] * PILEUP_GAUSSIAN_TABLE_DU;
        const double u_high = (l + 2 == n_nodes) ? u_max : i_u[l + 1] * PILEUP_GAUSSIAN_TABLE_DU;
        w[l] = exp(-p->rate * t_peak * u_low) * -expm1(-p->rate * t_peak * (u_high - u_low));
    }
    size_t n_nonzero = 0;
    for(size_t i = 0; i < n; i++) {
        if(h->bin[i] != 0.0) {
            nonzero[n_nonzero++] = i;
        }
    }
    for(size_t a = 0; a < n_nonzero; a++) {
        const size_t i = nonzero[a]; /* Leading pulse */
        const double c_i = weight * h->bin[i] / N;
        for(size_t b = 0; b < n_nonzero; b++) {
            const size_t j = nonzero[b]; /* Trailing pulse */
            const double c = c_i * h->bin[j];
            const double S = i + j + 1.0; /* Sum of amplitudes at centers of bins */
            const double q = (j + 0.5) / S * PILEUP_GAUSSIAN_TABLE_NQ;
            const size_t i_q = GSL_MIN((size_t) q, PILEUP_GAUSSIAN_TABLE_NQ - 1);
            const double f_q = q - i_q;
            const double *v1_low = v1_tab + i_q * tab->n_u, *v1_high = v1_low + tab->n_u;
            const double *v2_low = v2_tab + i_q * tab->n_u, *v2_high = v2_low + tab->n_u;
            int second = FALSE;
            for(size_t l = 0; l + 1 < n_nodes; l++) {
                v1[l] = v1_low[i_u[l]] + f_q * (v1_high[i_u[l]] - v1_low[i_u[l]]);
                v2[l] = pileup_gaussian_v2(v2_low, v2_high, i_u[l], 0.0, f_q);
                second |= (v2[l] >= 0.0);
            }
            const double v1_a = v1_low[i_u_max] + f_u_last * (v1_low[i_u_max + 1] - v1_low[i_u_max]);
            const double v1_b = v1_high[i_u_max] + f_u_last * (v1_high[i_u_max + 1] - v1_high[i_u_max]);
            v1[n_nodes - 1] = v1_a + f_q * (v1_b - v1_a);
            v2[n_nodes - 1] = pileup_gaussian_v2(v2_low, v2_high, i_u_max, f_u_last, f_q);
            pileup_deposit_curve(&d, S, v1, w, n_nodes, c);
            if(second) {
                pileup_deposit_curve(&d, S, v2, w, n_nodes, c * weight2);
            }
        }
    }
    double slope = 0.0;
    for(size_t k = 0; k < h_pileup->n; k++) { /* Bin k is F(k + 1) - F(k) */
        slope += d.slope[k + 1];
        h_pileup->bin[k] += d.jump[k + 1] + slope;
    }
    status = EXIT_SUCCESS;
cleanup:
    free(v1_tab);
    free(v2_tab);
    free(i_u);
    free(w);
    free(v1);
    free(v2);
    free(nonzero);
    free(d.jump);
    free(d.slope);
    return status;
}

int pileup_calculate(const jabs_histogram *h, jabs_histogram *h_pileup, const pileup_params *p, pileup_result *result) {
    result->scale = 1.0;
    result->fraction = 0.0;
    if(!h || !h_pileup || p->rate < 0.0 || p->rise < 0.0 || p->flat < 0.0 || p->shaping < 0.0 || p->pair_resolution < 0.0 || p->conversion_delay < 0.0) {
        return EXIT_FAILURE;
    }
    const double t_merge = pileup_merge_time(p);
    if(p->rate == 0.0 || t_merge == 0.0) {
        return EXIT_SUCCESS;
    }
    size_t n = 0; /* Channels up to the last non-zero channel */
    double N = 0.0;
    for(size_t i = 0; i < h->n; i++) {
        if(h->bin[i] != 0.0) {
            n = i + 1;
        }
        N += h->bin[i];
    }
    if(n == 0 || N <= 0.0) {
        return EXIT_SUCCESS;
    }
    const double dt_max = (p->reject && p->pair_resolution > 0.0) ? GSL_MIN_DBL(p->pair_resolution, t_merge) : t_merge; /* Pairs closer than this are recorded, otherwise rejected */
    const double clean = exp(-p->rate * t_merge); /* Probability that no pulse follows within T */
    const int reject = dt_max < t_merge;
    const double weight = reject ? clean : 1.0; /* Recorded pairs are rejected if another pulse follows within T, see pileup.h */
    const double weight2 = reject ? 1.0 : 1.0 / clean; /* Recorded trailing pulses may be followed by further pulses, see pileup.h */
    result->scale = clean; /* Every pulse arriving in live time is a leading pulse */
    double N_pileup = 0.0;
    for(size_t i = 0; i < h_pileup->n; i++) {
        N_pileup -= h_pileup->bin[i];
    }
    int status;
    switch(p->shaper) {
        case PILEUP_SHAPER_TRAPEZOID:
            status = pileup_trapezoid(h, n, N, h_pileup, p, dt_max, weight, weight2);
            break;
        case PILEUP_SHAPER_GAUSSIAN:
            status = pileup_gaussian(h, n, N, h_pileup, p, dt_max, weight, weight2);
            break;
        default:
            return EXIT_FAILURE;
    }
    if(status) {
        return status;
    }
    for(size_t i = 0; i < h_pileup->n; i++) {
        N_pileup += h_pileup->bin[i];
    }
    result->fraction = N_pileup / (result->scale * N + N_pileup);
    return EXIT_SUCCESS;
}
