/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <jibal_units.h>
#include <jibal_phys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "jabs_debug.h"
#include "defaults.h"
#include "ms.h"

#define MS_TABLE_N 23 /* n = log2(tau) from -10 to 12 */
#define MS_TABLE_N_MIN (-10)
#define MS_NU_AVERAGING 1.65 /* Szilagyi et al. Eq. (52), central value of characteristic exponent */
#define MS_MIX_NU_MIN 0.1
#define MS_MIX_NU_MAX 1000.0
#define MS_MIX_N_NU 81
#define MS_MIX_P_LOW 0.02 /* Probability of variance multiplier in the lowest bin */
#define MS_MIX_P_HIGH 1.0e-5 /* Probability of variance multiplier in the highest bin, if not limited by MS_MIX_SIGMA_CAP */
#define MS_MIX_SIGMA_CAP 30.0 /* Largest standard deviation of a component in units of HWHM. Truncates the extreme tails. */

/* Amsel et al. NIMB 201 (2003) 325, Table 1 (TF) and Table 3 (LJ). Columns: w (reduced HWHM), M*w, m_w, C_b, m_b. Row i corresponds to tau = 2^(i - 10) */
static const double ms_table_tf[MS_TABLE_N][5] = {
        {2.894e-05, 0.200, 0.626, 2.176, 0.622},
        {8.797e-05, 0.203, 0.626, 2.177, 0.620},
        {2.665e-04, 0.209, 0.627, 2.184, 0.608},
        {8.005e-04, 0.218, 0.634, 2.175, 0.623},
        {2.362e-03, 0.234, 0.649, 2.177, 0.620},
        {6.723e-03, 0.255, 0.680, 2.172, 0.629},
        {1.799e-02, 0.282, 0.734, 2.157, 0.656},
        {4.414e-02, 0.312, 0.817, 2.126, 0.716},
        {9.811e-02, 0.340, 0.923, 2.081, 0.809},
        {1.988e-01, 0.364, 1.042, 2.030, 0.925},
        {3.730e-01, 0.383, 1.160, 1.982, 1.049},
        {6.599e-01, 0.398, 1.269, 1.940, 1.168},
        {1.117e+00, 0.409, 1.363, 1.905, 1.275},
        {1.830e+00, 0.418, 1.445, 1.877, 1.369},
        {2.923e+00, 0.426, 1.513, 1.855, 1.449},
        {4.581e+00, 0.431, 1.570, 1.828, 1.517},
        {7.075e+00, 0.436, 1.618, 1.824, 1.573},
        {1.080e+01, 0.440, 1.658, 1.812, 1.620},
        {1.633e+01, 0.443, 1.692, 1.803, 1.660},
        {2.451e+01, 0.445, 1.721, 1.795, 1.693},
        {3.657e+01, 0.448, 1.745, 1.789, 1.722},
        {5.427e+01, 0.449, 1.766, 1.784, 1.746},
        {8.020e+01, 0.451, 1.784, 1.779, 1.766}
};

static const double ms_table_lj[MS_TABLE_N][5] = {
        {1.710e-07, 0.086, 0.382, 2.332, 0.382},
        {1.049e-06, 0.089, 0.382, 2.332, 0.382},
        {6.433e-06, 0.093, 0.382, 2.332, 0.382},
        {3.930e-05, 0.101, 0.384, 2.332, 0.382},
        {2.365e-04, 0.118, 0.390, 2.332, 0.382},
        {1.349e-03, 0.149, 0.410, 2.332, 0.382},
        {6.663e-03, 0.199, 0.467, 2.327, 0.390},
        {2.553e-02, 0.258, 0.579, 2.297, 0.431},
        {7.377e-02, 0.310, 0.739, 2.205, 0.572},
        {1.712e-01, 0.347, 0.913, 2.108, 0.752},
        {3.439e-01, 0.374, 1.074, 2.028, 0.930},
        {6.303e-01, 0.392, 1.212, 1.967, 1.089},
        {1.087e+00, 0.406, 1.327, 1.921, 1.224},
        {1.800e+00, 0.416, 1.421, 1.887, 1.336},
        {2.894e+00, 0.424, 1.497, 1.861, 1.427},
        {4.553e+00, 0.431, 1.560, 1.841, 1.502},
        {7.048e+00, 0.436, 1.611, 1.826, 1.563},
        {1.077e+01, 0.439, 1.654, 1.814, 1.614},
        {1.631e+01, 0.443, 1.689, 1.804, 1.655},
        {2.449e+01, 0.445, 1.719, 1.796, 1.690},
        {3.655e+01, 0.447, 1.744, 1.789, 1.719},
        {5.426e+01, 0.449, 1.765, 1.784, 1.743},
        {8.019e+01, 0.451, 1.783, 1.779, 1.763}
};

static double ms_t_peak_hwhm_product(double nu) { /* Peak height times HWHM of Student's t distribution. Gaussian limit is sqrt(ln 2/pi) = 0.4697 */
    return exp(gsl_sf_lngamma((nu + 1.0) / 2.0) - gsl_sf_lngamma(nu / 2.0)) / sqrt(C_PI) * sqrt(pow(2.0, 2.0 / (nu + 1.0)) - 1.0);
}

static double ms_t_nu_from_peak_hwhm_product(double mw) { /* Inverse of the above, by bisection in log(nu) */
    double lo = log(1e-3), hi = log(1e5);
    if(mw >= ms_t_peak_hwhm_product(exp(hi))) {
        return exp(hi);
    }
    if(mw <= ms_t_peak_hwhm_product(exp(lo))) {
        return exp(lo);
    }
    for(int i = 0; i < 100; i++) {
        double mid = (lo + hi) / 2.0;
        if(ms_t_peak_hwhm_product(exp(mid)) < mw) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return exp((lo + hi) / 2.0);
}

static double ms_mix_pdf(const ms_mix *m, double x) {
    double sum = 0.0;
    for(int k = 0; k < MS_MIX_N; k++) {
        if(m->w[k] > 0.0) {
            sum += m->w[k] * exp(-0.5 * x * x / m->v[k]) / sqrt(2.0 * C_PI * m->v[k]);
        }
    }
    return sum;
}

double ms_mix_hwhm(const ms_mix *m) {
    double half = ms_mix_pdf(m, 0.0) / 2.0;
    if(half <= 0.0) {
        return 0.0;
    }
    double lo = 0.0, hi = 1.0;
    while(ms_mix_pdf(m, hi) > half) {
        hi *= 2.0;
    }
    for(int i = 0; i < 100; i++) {
        double mid = (lo + hi) / 2.0;
        if(ms_mix_pdf(m, mid) > half) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / 2.0;
}

void ms_mix_zero(ms_mix *m) {
    for(int k = 0; k < MS_MIX_N; k++) {
        m->w[k] = (k == 0) ? 1.0 : 0.0;
        m->v[k] = 0.0;
    }
}

static void ms_mix_unit(double nu, ms_mix *m) { /* Student's t with nu degrees of freedom and unit HWHM. t is a Gaussian with variance s^2 * lambda, where lambda ~ InvGamma(nu/2, nu/2). The distribution of lambda is binned logarithmically. */
    const double a = nu / 2.0, b = 2.0 / nu; /* Shape and scale of gamma distributed 1/lambda */
    const double s2 = 1.0 / (nu * (pow(2.0, 2.0 / (nu + 1.0)) - 1.0)); /* Squared scale of t giving unit HWHM */
    const double v_cap = pow2(MS_MIX_SIGMA_CAP / sqrt(2.0 * log(2.0)));
    double lo = 1.0 / gsl_cdf_gamma_Qinv(MS_MIX_P_LOW, a, b);
    double hi = 1.0 / gsl_cdf_gamma_Qinv(1.0 - MS_MIX_P_HIGH, a, b);
    hi = GSL_MIN_DBL(hi, v_cap / s2);
    if(!(hi > lo * 1.0001)) {
        hi = lo * 1.0001;
    }
    double r = pow(hi / lo, 1.0 / MS_MIX_N);
    double cdf_prev = 0.0;
    for(int k = 0; k < MS_MIX_N; k++) {
        double e_low = lo * pow(r, k);
        double e_high = e_low * r;
        double cdf = (k == MS_MIX_N - 1) ? 1.0 : gsl_cdf_gamma_Q(1.0 / e_high, a, b); /* P(lambda < e_high) */
        m->w[k] = cdf - cdf_prev;
        m->v[k] = GSL_MIN_DBL(sqrt(e_low * e_high) * s2, v_cap);
        cdf_prev = cdf;
    }
    double h = ms_mix_hwhm(m);
    for(int k = 0; k < MS_MIX_N; k++) {
        m->v[k] /= h * h;
    }
}

ms_model *ms_model_new(ms_model_type type) {
    const double (*table)[5];
    switch(type) {
        case MS_MODEL_AMSEL_TF:
            table = ms_table_tf;
            break;
        case MS_MODEL_AMSEL_LJ:
            table = ms_table_lj;
            break;
        default:
            return NULL;
    }
    ms_model *ms = malloc(sizeof(ms_model));
    if(!ms) {
        return NULL;
    }
    ms->type = type;
    ms->n = MS_TABLE_N;
    ms->t = malloc(ms->n * sizeof(ms_table_point));
    ms->n_nu = MS_MIX_N_NU;
    ms->unit = malloc(ms->n_nu * sizeof(ms_mix));
    if(!ms->t || !ms->unit) {
        ms_model_free(ms);
        return NULL;
    }
    for(size_t i = 0; i < ms->n; i++) {
        ms_table_point *p = &ms->t[i];
        p->log_tau = (MS_TABLE_N_MIN + (int)i) * M_LN2;
        p->log_w = log(table[i][0]);
        p->nu = ms_t_nu_from_peak_hwhm_product(table[i][1]);
        p->m_w = table[i][2];
        p->C_b = table[i][3];
        p->m_b = table[i][4];
    }
    ms->log_nu_min = log(MS_MIX_NU_MIN);
    ms->log_nu_step = (log(MS_MIX_NU_MAX) - ms->log_nu_min) / (ms->n_nu - 1);
    for(size_t i = 0; i < ms->n_nu; i++) {
        ms_mix_unit(exp(ms->log_nu_min + i * ms->log_nu_step), &ms->unit[i]);
    }
    return ms;
}

void ms_model_free(ms_model *ms) {
    if(!ms) {
        return;
    }
    free(ms->t);
    free(ms->unit);
    free(ms);
}

const char *ms_model_name(ms_model_type type) {
    return jibal_option_get_string(ms_model_option, type);
}

double ms_screening_length(int Z2) {
    return 0.8853 * C_BOHR_RADIUS * pow(Z2, -1.0 / 3.0);
}

double ms_scaling_factor(int Z1, int Z2, double a, double E) {
    return E * a / (2.0 * Z1 * Z2 * C_E * C_E * C_COULOMB);
}

double ms_w_reduced(const ms_model *ms, double tau) { /* Piecewise linear in log-log scale, extrapolation with local exponents */
    if(tau <= 0.0) {
        return 0.0;
    }
    double lt = log(tau);
    const ms_table_point *first = &ms->t[0], *last = &ms->t[ms->n - 1];
    if(lt <= first->log_tau) {
        return exp(first->log_w + (lt - first->log_tau) / first->m_w);
    }
    if(lt >= last->log_tau) {
        return exp(last->log_w + (lt - last->log_tau) / last->m_w);
    }
    size_t i = (size_t)((lt - first->log_tau) / M_LN2);
    if(i >= ms->n - 1) {
        i = ms->n - 2;
    }
    const ms_table_point *p = &ms->t[i];
    const ms_table_point *p2 = &ms->t[i + 1];
    double frac = (lt - p->log_tau) / (p2->log_tau - p->log_tau);
    return exp(p->log_w + frac * (p2->log_w - p->log_w));
}

double ms_tau(const ms_model *ms, double w_reduced) {
    if(w_reduced <= 0.0) {
        return 0.0;
    }
    double lw = log(w_reduced);
    const ms_table_point *first = &ms->t[0], *last = &ms->t[ms->n - 1];
    if(lw <= first->log_w) {
        return exp(first->log_tau + (lw - first->log_w) * first->m_w);
    }
    if(lw >= last->log_w) {
        return exp(last->log_tau + (lw - last->log_w) * last->m_w);
    }
    size_t lo = 0, hi = ms->n - 1;
    while(hi - lo > 1) {
        size_t mi = (lo + hi) / 2;
        if(lw >= ms->t[mi].log_w) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    const ms_table_point *p = &ms->t[lo];
    const ms_table_point *p2 = &ms->t[hi];
    double frac = (lw - p->log_w) / (p2->log_w - p->log_w);
    return exp(p->log_tau + frac * (p2->log_tau - p->log_tau));
}

void ms_table_params(const ms_model *ms, double tau, double *C_b, double *m_b, double *nu) {
    const ms_table_point *first = &ms->t[0], *last = &ms->t[ms->n - 1];
    double lt = tau > 0.0 ? log(tau) : first->log_tau;
    const ms_table_point *p, *p2;
    double frac;
    if(lt <= first->log_tau) {
        p = p2 = first;
        frac = 0.0;
    } else if(lt >= last->log_tau) {
        p = p2 = last;
        frac = 0.0;
    } else {
        size_t i = (size_t)((lt - first->log_tau) / M_LN2);
        if(i >= ms->n - 1) {
            i = ms->n - 2;
        }
        p = &ms->t[i];
        p2 = &ms->t[i + 1];
        frac = (lt - p->log_tau) / (p2->log_tau - p->log_tau);
    }
    if(C_b) {
        *C_b = p->C_b + frac * (p2->C_b - p->C_b);
    }
    if(m_b) {
        *m_b = p->m_b + frac * (p2->m_b - p->m_b);
    }
    if(nu) {
        *nu = exp(log(p->nu) + frac * (log(p2->nu) - log(p->nu)));
    }
}

void ms_step_element(const ms_model *ms, ion_ms *s, int Z1, int Z2, double c, double E, double h) {
    double a = ms_screening_length(Z2);
    double tau = C_PI * a * a * c * h;
    if(tau <= 0.0) {
        return;
    }
    double mu = ms_scaling_factor(Z1, Z2, a, E);
    double tau_eq = ms_tau(ms, mu * s->w); /* Szilagyi et al. Eq. (48), reduced thickness of this element (at this energy) that would produce the current angular spread */
    s->w = ms_w_reduced(ms, tau_eq + tau) / mu; /* Eq. (49) */
    double mu_pow = pow(mu, -MS_NU_AVERAGING); /* Eq. (50) */
    s->sum_a += tau * mu_pow * mu;
    s->sum_b += tau * mu_pow;
}

double ms_tau_star(const ms_model *ms, const ion_ms *s) {
    if(s->w <= 0.0 || s->sum_b <= 0.0) {
        return 0.0;
    }
    double mu_mean = s->sum_a / s->sum_b;
    return ms_tau(ms, mu_mean * s->w);
}

double ms_sz(double k, double m) {
    if(fabs(k) < 1e-9) {
        return 1.0;
    }
    double q = (k < -1.0) ? 1.0 : -1.0;
    double val = fabs((pow(fabs(1.0 + k), m + 1.0) + q) / (k * (m + 1.0)));
    return pow(val, 1.0 / m);
}

double ms_leg_hwhm(const ms_model *ms, const ion_ms *s, const double a[2], const double b[2], double *nu) {
    if(nu) {
        *nu = MS_MIX_NU_MAX;
    }
    if(!ms || s->w <= 0.0) {
        return 0.0;
    }
    double C_b, m_b;
    ms_table_params(ms, ms_tau_star(ms, s), &C_b, &m_b, nu);
    double sum = 0.0;
    for(int i = 0; i < 2; i++) {
        double h;
        if(fabs(a[i]) <= 1e-12 * fabs(b[i])) { /* Lateral spread (chord angle) only, also handles a = b = 0 */
            h = fabs(b[i]) * s->w / C_b;
        } else {
            h = fabs(a[i]) * s->w * ms_sz(b[i] / a[i], m_b);
        }
        sum += h * h; /* Projections to orthogonal directions add quadratically, Amsel et al. Eq. (A.10) */
    }
    return sqrt(sum);
}

void ms_mix_t(const ms_model *ms, double hwhm, double nu, ms_mix *out) {
    if(hwhm <= 0.0) {
        ms_mix_zero(out);
        return;
    }
    double x = (log(nu) - ms->log_nu_min) / ms->log_nu_step;
    size_t i;
    double frac;
    if(x <= 0.0) {
        i = 0;
        frac = 0.0;
    } else if(x >= ms->n_nu - 1) {
        i = ms->n_nu - 2;
        frac = 1.0;
    } else {
        i = (size_t)x;
        frac = x - i;
    }
    const ms_mix *lo = &ms->unit[i], *hi = &ms->unit[i + 1];
    double w_sum = 0.0;
    for(int k = 0; k < MS_MIX_N; k++) {
        out->w[k] = lo->w[k] + frac * (hi->w[k] - lo->w[k]);
        out->v[k] = (lo->v[k] + frac * (hi->v[k] - lo->v[k])) * hwhm * hwhm;
        w_sum += out->w[k];
    }
    for(int k = 0; k < MS_MIX_N; k++) {
        out->w[k] /= w_sum;
    }
}

void ms_mix_convolve(const ms_mix *a, const ms_mix *b, ms_mix *out) { /* Variances of all component pairs are summed, then components are combined into MS_MIX_N logarithmic bins */
    double v_min = DBL_MAX, v_max = 0.0;
    for(int i = 0; i < MS_MIX_N; i++) {
        for(int j = 0; j < MS_MIX_N; j++) {
            if(a->w[i] <= 0.0 || b->w[j] <= 0.0) {
                continue;
            }
            double v = a->v[i] + b->v[j];
            v_min = GSL_MIN_DBL(v_min, v);
            v_max = GSL_MAX_DBL(v_max, v);
        }
    }
    if(v_max <= 0.0) {
        ms_mix_zero(out);
        return;
    }
    double l_min = log(v_min);
    double l_step = (log(v_max) - l_min) / MS_MIX_N;
    double wv[MS_MIX_N];
    for(int k = 0; k < MS_MIX_N; k++) {
        out->w[k] = 0.0;
        wv[k] = 0.0;
    }
    for(int i = 0; i < MS_MIX_N; i++) {
        for(int j = 0; j < MS_MIX_N; j++) {
            double w = a->w[i] * b->w[j];
            if(w <= 0.0) {
                continue;
            }
            double v = a->v[i] + b->v[j];
            int k = l_step > 0.0 ? (int)((log(v) - l_min) / l_step) : 0;
            k = GSL_MAX_INT(0, GSL_MIN_INT(k, MS_MIX_N - 1));
            out->w[k] += w;
            wv[k] += w * v;
        }
    }
    for(int k = 0; k < MS_MIX_N; k++) {
        if(out->w[k] > 0.0) {
            out->v[k] = wv[k] / out->w[k];
        } else { /* Empty bin. Variance is set to bin center so that it is sensible to pair it with another mixture. */
            out->v[k] = exp(l_min + (k + 0.5) * l_step);
        }
    }
}
