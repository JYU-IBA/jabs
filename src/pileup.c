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
#define PILEUP_GAUSSIAN_TABLE_NU_MAX (512)
#define PILEUP_GAUSSIAN_SEARCH_DT (0.01) /* Step when searching for maximum of sum of pulses (units of peaking time) */
#define PILEUP_GAUSSIAN_TAIL_START (2.0) /* Time between pulses (units of peaking time) after which nodes are PILEUP_GAUSSIAN_TAIL_STEP grid steps apart, sum changes slowly */
#define PILEUP_GAUSSIAN_TAIL_STEP (4)
#define PILEUP_DEPOSIT_DK_MIN (0.25) /* Amplitude (channels) of a pair changes at least this much between nodes, shorter intervals are merged */
#define PILEUP_DEPOSIT_WIDTH_MIN (1.0e-6) /* Narrower uniform distributions of amplitude are deposited as points */

typedef struct pileup_gaussian_table {
    double t_tail; /* Time from start of pulse until it falls to PILEUP_GAUSSIAN_TAIL of maximum (units of peaking time) */
    size_t n_u;
    double m[PILEUP_GAUSSIAN_TABLE_NQ + 1][PILEUP_GAUSSIAN_TABLE_NU_MAX]; /* Maximum of (1 - q) * pulse(t) + q * pulse(t - u), m[i_q][i_u] */
} pileup_gaussian_table;

static pileup_gaussian_table gaussian_table;
static int gaussian_table_ready = 0;

static double pileup_gaussian_pulse(double t) { /* CR-RC^n, t in units of peaking time, maximum is 1 at t = 1 */
    if(t <= 0.0) {
        return 0.0;
    }
    return exp(PILEUP_GAUSSIAN_ORDER * (log(t) + 1.0 - t));
}

static double pileup_gaussian_sum_max(double q, double u) {
    const double t_max = u + 2.0;
    const size_t n = (size_t) ceil(t_max / PILEUP_GAUSSIAN_SEARCH_DT);
    size_t i_best = 0;
    double best = 0.0;
    for(size_t i = 1; i < n; i++) {
        const double t = i * PILEUP_GAUSSIAN_SEARCH_DT;
        const double s = (1.0 - q) * pileup_gaussian_pulse(t) + q * pileup_gaussian_pulse(t - u);
        if(s > best) {
            best = s;
            i_best = i;
        }
    }
    if(i_best == 0) {
        return best;
    }
    const double t = i_best * PILEUP_GAUSSIAN_SEARCH_DT; /* Parabola through the best sample and its neighbours */
    const double s_minus = (1.0 - q) * pileup_gaussian_pulse(t - PILEUP_GAUSSIAN_SEARCH_DT) + q * pileup_gaussian_pulse(t - PILEUP_GAUSSIAN_SEARCH_DT - u);
    const double s_plus = (1.0 - q) * pileup_gaussian_pulse(t + PILEUP_GAUSSIAN_SEARCH_DT) + q * pileup_gaussian_pulse(t + PILEUP_GAUSSIAN_SEARCH_DT - u);
    const double curv = s_minus - 2.0 * best + s_plus;
    if(curv >= 0.0) {
        return best;
    }
    return best - pow2(s_plus - s_minus) / (8.0 * curv);
}

static void pileup_gaussian_table_init(pileup_gaussian_table *tab) {
    double low = 1.0, high = 10.0;
    for(int i = 0; i < 100; i++) { /* Bisection, pulse is decreasing after peak */
        const double mid = (low + high) / 2.0;
        if(pileup_gaussian_pulse(mid) > PILEUP_GAUSSIAN_TAIL) {
            low = mid;
        } else {
            high = mid;
        }
    }
    tab->t_tail = (low + high) / 2.0;
    tab->n_u = GSL_MIN((size_t) ceil(tab->t_tail / PILEUP_GAUSSIAN_TABLE_DU) + 2, PILEUP_GAUSSIAN_TABLE_NU_MAX);
    for(size_t i_q = 0; i_q <= PILEUP_GAUSSIAN_TABLE_NQ; i_q++) {
        const double q = (1.0 * i_q) / PILEUP_GAUSSIAN_TABLE_NQ;
        for(size_t i_u = 0; i_u < tab->n_u; i_u++) {
            tab->m[i_q][i_u] = pileup_gaussian_sum_max(q, i_u * PILEUP_GAUSSIAN_TABLE_DU);
        }
    }
}

static const pileup_gaussian_table *pileup_gaussian_table_get(void) {
#pragma omp critical(pileup_gaussian_table)
    {
        if(!gaussian_table_ready) {
            pileup_gaussian_table_init(&gaussian_table);
            gaussian_table_ready = 1;
        }
    }
    return &gaussian_table;
}

double pileup_merge_time(const pileup_params *p) {
    switch(p->shaper) {
        case PILEUP_SHAPER_TRAPEZOID:
            return p->rise + p->flat;
        case PILEUP_SHAPER_GAUSSIAN:
            return pileup_gaussian_table_get()->t_tail * PILEUP_GAUSSIAN_PEAKING_FACTOR * p->shaping;
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

static int pileup_trapezoid(const jabs_histogram *h, size_t n, double N, jabs_histogram *h_pileup, const pileup_params *p, double dt_max) {
    const double r = p->rate;
    const double t_sample = p->rise + p->flat;
    const double w_full = -expm1(-r * GSL_MIN_DBL(p->flat, dt_max)); /* Trailing pulse within flat top, full amplitude is added */
    const int partial = p->rise > 0.0 && dt_max > p->flat; /* Trailing pulse during its rise, fraction g of amplitude is added */
    const double g_low = partial ? (t_sample - dt_max) / p->rise : 1.0;
    const double c_partial = exp(-r * dt_max);
    const double r_rise = r * p->rise;
    double *Q = calloc(n, sizeof(double)); /* Distribution of amplitude added by the trailing pulse, normalized to probability */
    if(!Q) {
        return EXIT_FAILURE;
    }
    for(size_t j = 0; j < n; j++) {
        const double m = h->bin[j] / N;
        if(m == 0.0) {
            continue;
        }
        Q[j] += w_full * m;
        if(!partial) {
            continue;
        }
        const double x = j + 0.5; /* Amplitude of trailing pulse, added amplitude g * x for g in [g_low, 1], density of g is r * rise * exp(-r * (t_sample - g * rise)) */
        double C_prev = 0.0;
        for(size_t k = (size_t) floor(g_low * x); k <= j; k++) {
            const double b = GSL_MIN_DBL(k + 1.0, x);
            const double C = c_partial * expm1(r_rise * (b / x - g_low)); /* Cumulative probability from g_low to b / x */
            Q[k] += m * (C - C_prev);
            C_prev = C;
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

static int pileup_gaussian(const jabs_histogram *h, size_t n, double N, jabs_histogram *h_pileup, const pileup_params *p, double dt_max) {
    const pileup_gaussian_table *tab = pileup_gaussian_table_get();
    const double t_peak = PILEUP_GAUSSIAN_PEAKING_FACTOR * p->shaping;
    const double u_max = dt_max / t_peak; /* Pairs closer than this (units of peaking time) are recorded */
    const size_t i_u_max = (size_t) floor(u_max / PILEUP_GAUSSIAN_TABLE_DU);
    if(i_u_max + 2 > tab->n_u) {
        return EXIT_FAILURE;
    }
    const size_t n_nodes_max = i_u_max + 2;
    size_t *i_u = malloc(n_nodes_max * sizeof(size_t)); /* Nodes are on the grid of the table (except the last), sum amplitude is linear in time between nodes */
    double *w = malloc(n_nodes_max * sizeof(double)); /* Probability of trailing pulse between nodes l and l + 1 */
    double *m = malloc(n_nodes_max * sizeof(double));
    size_t *nonzero = malloc(n * sizeof(size_t));
    pileup_deposit d = {.n = h_pileup->n, .jump = calloc(h_pileup->n + 2, sizeof(double)), .slope = calloc(h_pileup->n + 2, sizeof(double))};
    if(!i_u || !w || !m || !nonzero || !d.jump || !d.slope) {
        free(i_u);
        free(w);
        free(m);
        free(nonzero);
        free(d.jump);
        free(d.slope);
        return EXIT_FAILURE;
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
        const double c_i = h->bin[i] / N;
        for(size_t b = 0; b < n_nonzero; b++) {
            const size_t j = nonzero[b]; /* Trailing pulse */
            const double c = c_i * h->bin[j];
            const double S = i + j + 1.0; /* Sum of amplitudes at centers of bins */
            const double q = (j + 0.5) / S * PILEUP_GAUSSIAN_TABLE_NQ;
            const size_t i_q = GSL_MIN((size_t) q, PILEUP_GAUSSIAN_TABLE_NQ - 1);
            const double f_q = q - i_q;
            const double *m_low = tab->m[i_q], *m_high = tab->m[i_q + 1];
            for(size_t l = 0; l + 1 < n_nodes; l++) {
                m[l] = m_low[i_u[l]] + f_q * (m_high[i_u[l]] - m_low[i_u[l]]);
            }
            const double m_a = m_low[i_u_max] + f_u_last * (m_low[i_u_max + 1] - m_low[i_u_max]);
            const double m_b = m_high[i_u_max] + f_u_last * (m_high[i_u_max + 1] - m_high[i_u_max]);
            m[n_nodes - 1] = m_a + f_q * (m_b - m_a);
            size_t l_start = 0;
            double w_sum = 0.0;
            for(size_t l = 0; l + 1 < n_nodes; l++) { /* Amplitude S * m. Intervals are merged until amplitude changes by PILEUP_DEPOSIT_DK_MIN. Amplitudes are spread within bins, this is approximated by two copies shifted by +-0.25. */
                w_sum += w[l];
                if(l + 2 < n_nodes && S * fabs(m[l + 1] - m[l_start]) < PILEUP_DEPOSIT_DK_MIN) {
                    continue;
                }
                const double low = S * GSL_MIN_DBL(m[l_start], m[l + 1]);
                const double high = S * GSL_MAX_DBL(m[l_start], m[l + 1]);
                pileup_deposit_uniform(&d, low - 0.25, high - 0.25, 0.5 * c * w_sum);
                pileup_deposit_uniform(&d, low + 0.25, high + 0.25, 0.5 * c * w_sum);
                l_start = l + 1;
                w_sum = 0.0;
            }
        }
    }
    double slope = 0.0;
    for(size_t k = 0; k < h_pileup->n; k++) { /* Bin k is F(k + 1) - F(k) */
        slope += d.slope[k + 1];
        h_pileup->bin[k] += d.jump[k + 1] + slope;
    }
    free(i_u);
    free(w);
    free(m);
    free(nonzero);
    free(d.jump);
    free(d.slope);
    return EXIT_SUCCESS;
}

int pileup_calculate(const jabs_histogram *h, jabs_histogram *h_pileup, const pileup_params *p, pileup_result *result) {
    result->scale = 1.0;
    result->fraction = 0.0;
    if(!h || !h_pileup || p->rate < 0.0 || p->rise < 0.0 || p->flat < 0.0 || p->shaping < 0.0 || p->pair_resolution < 0.0) {
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
    const double dt_max = p->pair_resolution > 0.0 ? GSL_MIN_DBL(p->pair_resolution, t_merge) : t_merge; /* Pairs closer than this are recorded as pile-up events, otherwise rejected */
    const double clean = exp(-p->rate * t_merge); /* Probability that no pulse follows within merging time */
    int status;
    switch(p->shaper) {
        case PILEUP_SHAPER_TRAPEZOID:
            status = pileup_trapezoid(h, n, N, h_pileup, p, dt_max);
            break;
        case PILEUP_SHAPER_GAUSSIAN:
            status = pileup_gaussian(h, n, N, h_pileup, p, dt_max);
            break;
        default:
            return EXIT_FAILURE;
    }
    if(status) {
        return status;
    }
    const double piled = -expm1(-p->rate * dt_max);
    result->scale = clean; /* Every pulse arriving in live time is a leading pulse, see pileup.h */
    result->fraction = piled / (clean + piled);
    return EXIT_SUCCESS;
}
