/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_MS_H
#define JABS_MS_H

/* Small angle multiple scattering (MS).
 *
 * Angular spread is calculated using the tabulated universal projected angular distributions of
 * G. Amsel, G. Battistig, A. L'Hoir, Nucl. Instr. and Meth. B 201 (2003) 325 (Tables 1 and 3, TF and LJ potentials).
 * Thick multielemental targets are handled with the "transformation method" of
 * E. Szilagyi, F. Paszti, G. Amsel, Nucl. Instr. and Meth. B 100 (1995) 103, Section 4.4.1.
 * Energy spread is obtained from correlated angular and lateral (chord angle) spread, Amsel et al. Section 4.3 and 8.1.
 *
 * The shape of MS distributions is approximated by Pearson VII (Student's t) distributions, shape parameter is
 * determined from the peak height and HWHM product (M*w) given by Amsel et al. The Pearson VII distributions are
 * represented as scale mixtures of Gaussians (ms_mix), which allows convolution with Gaussian straggling. */

#include <jibal_option.h>
#include "ion.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MS_MIX_N 10 /* Number of Gaussian components in scale mixture */

typedef enum ms_model_type {
    MS_MODEL_NONE = 0,
    MS_MODEL_AMSEL_TF = 1,
    MS_MODEL_AMSEL_LJ = 2
} ms_model_type;

static const jibal_option ms_model_option[] = {
        {JIBAL_OPTION_STR_NONE, MS_MODEL_NONE},
        {"TF", MS_MODEL_AMSEL_TF},
        {"LJ", MS_MODEL_AMSEL_LJ},
        {NULL, 0}
};

typedef struct ms_mix {
    double w[MS_MIX_N]; /* Weights, sum is one */
    double v[MS_MIX_N]; /* Variances, increasing */
} ms_mix;

typedef struct ms_table_point {
    double log_tau; /* ln of reduced thickness */
    double log_w; /* ln of reduced HWHM */
    double m_w; /* Characteristic exponent of HWHM, w ~ tau^(1/m_w) */
    double C_b; /* Ratio of HWHMs of angular and chord angle distributions */
    double m_b; /* Exponent for linear combination of angular and chord angle distributions */
    double nu; /* Shape (Student's t degrees of freedom, also "m" of Szilagyi et al. Eq. (41)) */
} ms_table_point;

typedef struct ms_model {
    ms_model_type type;
    size_t n;
    ms_table_point *t;
    size_t n_nu;
    double log_nu_min;
    double log_nu_step;
    ms_mix *unit; /* Table of n_nu scale mixtures with unit HWHM */
} ms_model;

ms_model *ms_model_new(ms_model_type type);
void ms_model_free(ms_model *ms);
const char *ms_model_name(ms_model_type type);
double ms_screening_length(int Z2); /* Thomas-Fermi screening length, assuming bare projectile */
double ms_scaling_factor(int Z1, int Z2, double a, double E); /* mu, Amsel et al. Eq. (22), reduced angle is mu times lab angle */
double ms_w_reduced(const ms_model *ms, double tau); /* Reduced HWHM of projected angular distribution at reduced thickness tau */
double ms_tau(const ms_model *ms, double w_reduced); /* Inverse of ms_w_reduced() */
void ms_table_params(const ms_model *ms, double tau, double *C_b, double *m_b, double *nu);
void ms_step_element(const ms_model *ms, ion_ms *s, int Z1, int Z2, double c, double E, double h); /* Updates angular spread after path length h (at/m2) in element Z2 with concentration c */
double ms_tau_star(const ms_model *ms, const ion_ms *s); /* Equivalent reduced thickness for determining shape (Szilagyi et al. Eq. (53)) */
double ms_sz(double k, double m); /* HWHM scaling factor of alpha + k * beta, compared to alpha (Amsel et al. Eq. (72)) */
double ms_leg_hwhm(const ms_model *ms, const ion_ms *s, const double a[2], const double b[2], double *nu); /* HWHM of energy spread a * alpha + b * beta in two orthogonal directions, shape returned in nu */
void ms_mix_t(const ms_model *ms, double hwhm, double nu, ms_mix *out); /* Student's t (Pearson VII) distribution with given HWHM as scale mixture */
void ms_mix_convolve(const ms_mix *a, const ms_mix *b, ms_mix *out);
void ms_mix_zero(ms_mix *m);
double ms_mix_hwhm(const ms_mix *m);
#ifdef __cplusplus
}
#endif
#endif // JABS_MS_H
