/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

    Accuracy test of stopping (energy loss and straggling) integration. Trajectories computed with stop_step() using the
    default step parameters and various fixed steps are compared to a reference computed with very small fixed steps. The reference
    is independent of stop_step(): it integrates energy and energy loss variance directly (see trajectory_reference()).
    Errors are evaluated at every step of the tested trajectory by interpolating the reference trajectory in depth.

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <jibal.h>
#include <jibal_units.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include "defaults.h"
#include "ion.h"
#include "sample.h"
#include "stop.h"
#include "sim_calc_params.h"

#define REF_STEP (0.1 * C_KEV)
#define REF_CHECK_STEP (0.4 * C_KEV) /* RK4 error scales as h^4, error of reference is roughly (difference to this) / 256 */
#define REF_DERIVATIVE_EPS (1e-5) /* Relative energy difference for numerical derivative of stopping in reference calculation */
#define PASSES_RK4 (5) /* Passes over sample isotopes per step in stop_step(): four stages (stopping and straggling together) and stopping at the end */
#define PASSES_EULER (3) /* Stopping and straggling at start and midpoint, stopping at the end */
#define PASSES_REF (16) /* Four stages, each with stopping at three energies and straggling */
#define SIGMA_REL_THRESHOLD (1.0 * C_KEV) /* Relative straggling errors are evaluated only where reference sigma exceeds this */

typedef struct scenario {
    const char *name;
    const char *ion;
    double E;
    double sigma0; /* Initial straggling (sigma) */
    const char *sample; /* Sample model string, or "file:" followed by sample model file name */
    double theta; /* Angle with respect to sample normal, > 90 deg means travelling towards the surface */
    double x0; /* Starting depth */
    int exiting; /* Use exiting ion step parameters as defaults */
    double emin;
} scenario;

static const scenario scenarios[] = {
    {"1H 2 MeV into Si",                 "1H",    2.0 * C_MEV,  0.0,          "Si 500000tfu",                     0.0,             0.0,            FALSE, 50.0 * C_KEV},
    {"4He 2 MeV into Au",                "4He",   2.0 * C_MEV,  0.0,          "Au 50000tfu",                      7.0 * C_DEG,     0.0,            FALSE, 50.0 * C_KEV},
    {"4He 2 MeV into TiN/SiO2 gradient", "4He",   2.0 * C_MEV,  0.0,          "file:../example/tests/point_by_point_profile.txt", 7.0 * C_DEG, 0.0, FALSE, 50.0 * C_KEV},
    {"35Cl 8.5 MeV into Al2O3/Si 70deg", "35Cl",  8.5 * C_MEV,  0.0,          "Au 1tfu Al2O3 1000tfu Si 10000tfu", 70.0 * C_DEG,   0.0,            FALSE, 10.0 * C_KEV},
    {"127I 50 MeV into SiO2 75deg",      "127I",  50.0 * C_MEV, 0.0,          "SiO2 30000tfu",                    75.0 * C_DEG,    0.0,            FALSE, 500.0 * C_KEV},
    {"4He 1.6 MeV out of Au 165deg",     "4He",   1.6 * C_MEV,  5.0 * C_KEV,  "Au 20000tfu",                      165.0 * C_DEG,   5000.0 * C_TFU, TRUE,  10.0 * C_KEV},
    {"4He 1.2 MeV out of Si under Au",   "4He",   1.2 * C_MEV,  5.0 * C_KEV,  "Au 1000tfu Si 30000tfu",           130.0 * C_DEG,   15000.0 * C_TFU, TRUE, 10.0 * C_KEV},
    {"4He 600 keV out of Si 170deg",     "4He",   0.6 * C_MEV,  5.0 * C_KEV,  "Si 10000tfu",                      170.0 * C_DEG,   3000.0 * C_TFU, TRUE,  10.0 * C_KEV},
    {"16O 4 MeV recoil out of Al2O3",    "16O",   4.0 * C_MEV,  20.0 * C_KEV, "Au 1tfu Al2O3 1000tfu Si 10000tfu", 110.0 * C_DEG,  900.0 * C_TFU,  TRUE,  10.0 * C_KEV},
    {"1H 1 MeV recoil out of Al2O3",     "1H",    1.0 * C_MEV,  5.0 * C_KEV,  "Au 1tfu Al2O3 1000tfu Si 10000tfu", 110.0 * C_DEG,  900.0 * C_TFU,  TRUE,  10.0 * C_KEV},
    {"35Cl 300 keV out of Si 140deg",   "35Cl",  0.3 * C_MEV,  5.0 * C_KEV,  "Au 1tfu Al2O3 1000tfu Si 10000tfu", 140.0 * C_DEG,  1300.0 * C_TFU, TRUE,  10.0 * C_KEV},
    {"28Si 300 keV recoil out of Si",    "28Si",  0.3 * C_MEV,  5.0 * C_KEV,  "Au 1tfu Al2O3 1000tfu Si 10000tfu", 110.0 * C_DEG,  1300.0 * C_TFU, TRUE,  10.0 * C_KEV},
    {"4He 2 MeV into Au/Si grazing 80deg", "4He",  2.0 * C_MEV,  0.0,          "Au 200tfu Si 20000tfu",            80.0 * C_DEG,    0.0,            FALSE, 50.0 * C_KEV},
    {"4He 1.8 MeV out of Si grazing 85deg", "4He",  1.8 * C_MEV,  5.0 * C_KEV,  "Au 100tfu Si 20000tfu",            95.0 * C_DEG,    1000.0 * C_TFU, TRUE,  10.0 * C_KEV},
    {"4He 1.5 MeV out of 20x(Ni/Ti)",      "4He",  1.5 * C_MEV,  5.0 * C_KEV,  "Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Ni 50tfu Ti 50tfu Si 10000tfu", 170.0 * C_DEG, 1950.0 * C_TFU, TRUE, 10.0 * C_KEV},
    {"4He 1.4 MeV out of TiN/SiO2 gradient", "4He", 1.4 * C_MEV, 5.0 * C_KEV,  "file:../example/tests/point_by_point_profile.txt", 165.0 * C_DEG, 1900.0 * C_TFU, TRUE, 10.0 * C_KEV},
    {"4He 900 keV out of Au, stops",       "4He",  0.9 * C_MEV,  5.0 * C_KEV,  "Au 20000tfu",                      170.0 * C_DEG,   15000.0 * C_TFU, TRUE, 10.0 * C_KEV},
    {"4He 3.9 MeV NRA alpha out of TiN",   "4He",  3.9 * C_MEV,  10.0 * C_KEV, "TiN 20000tfu Si 20000tfu",         150.0 * C_DEG,   15000.0 * C_TFU, TRUE, 10.0 * C_KEV},
    {"1H 3 MeV NRA proton out of thick Si", "1H",  3.0 * C_MEV,  10.0 * C_KEV, "Si 300000tfu",                     135.0 * C_DEG,   100000.0 * C_TFU, TRUE, 10.0 * C_KEV},
    {"48Ti 6 MeV recoil out of TiN",       "48Ti", 6.0 * C_MEV,  20.0 * C_KEV, "TiN 2000tfu Si 10000tfu",          115.0 * C_DEG,   1500.0 * C_TFU, TRUE,  10.0 * C_KEV},
    {"1H 1.5 MeV recoil out of Au/C grazing", "1H", 1.5 * C_MEV, 5.0 * C_KEV,  "Au 100tfu C 50000tfu",             100.0 * C_DEG,   3000.0 * C_TFU, TRUE,  10.0 * C_KEV},
    {"12C 3 MeV recoil out of C grazing",  "12C",  3.0 * C_MEV,  10.0 * C_KEV, "Au 100tfu C 50000tfu",             100.0 * C_DEG,   1500.0 * C_TFU, TRUE,  10.0 * C_KEV},
};

typedef struct trajectory {
    double *x;
    double *E;
    double *S;
    size_t n;
    size_t n_alloc;
    int stopped; /* Ion stopped (E below emin) inside sample */
    size_t n_evals; /* Passes over sample isotopes */
    size_t n_rejected; /* Rejected adaptive steps */
    double t; /* CPU time (s) */
} trajectory;

static void trajectory_add(trajectory *tr, const depth *d, const ion *ion) {
    if(tr->n == tr->n_alloc) {
        tr->n_alloc = tr->n_alloc ? tr->n_alloc * 2 : 1024;
        tr->x = realloc(tr->x, tr->n_alloc * sizeof(double));
        tr->E = realloc(tr->E, tr->n_alloc * sizeof(double));
        tr->S = realloc(tr->S, tr->n_alloc * sizeof(double));
    }
    tr->x[tr->n] = d->x;
    tr->E[tr->n] = ion->E;
    tr->S[tr->n] = ion->S;
    tr->n++;
}

static void trajectory_free(trajectory *tr) {
    free(tr->x);
    free(tr->E);
    free(tr->S);
}

/* Same logic as stop_sample_exit(), but works in both directions and records the trajectory. */
static trajectory trajectory_compute(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params, ion ion, depth d, const sample *sample, double emin) {
    trajectory tr = {0};
    int last = FALSE;
    const int adaptive = stop->rk4 && stop->accuracy > 0.0;
    stop_adaptive state = {0};
    emin = GSL_MAX_DBL(emin, ion.ion_gsto->emin);
    clock_t start = clock();
    trajectory_add(&tr, &d, &ion);
    while(1) {
        if(ion.inverse_cosine_theta > 0.0 && d.x >= (sample->thickness - DEPTH_TOLERANCE)) {
            break;
        }
        if(ion.inverse_cosine_theta < 0.0 && d.x <= DEPTH_TOLERANCE) {
            break;
        }
        if(last) {
            tr.stopped = TRUE;
            break;
        }
        double E_step = stop_step_calc(params, &ion);
        if(adaptive && state.h > 0.0) {
            E_step = stop_step_max(params, &ion);
        }
        if(ion.E - E_step <= emin) {
            if(ion.E <= emin) {
                tr.stopped = TRUE;
                break;
            }
            E_step = ion.E - emin;
            last = TRUE;
        }
        d = adaptive ? stop_step_adaptive(stop, stragg, &ion, sample, d, E_step, &state) : stop_step(stop, stragg, &ion, sample, d, E_step);
        if(last && ion.E - emin > STOP_LAST_STEP_REMAINDER * E_step) { /* Last step was truncated at range boundary or energy loss was smaller than expected */
            last = FALSE;
        }
        trajectory_add(&tr, &d, &ion);
    }
    tr.t = (double)(clock() - start) / CLOCKS_PER_SEC;
    tr.n_evals = adaptive ? state.n_evals : (tr.n - 1) * (stop->rk4 ? PASSES_RK4 : PASSES_EULER);
    tr.n_rejected = state.n_rejected;
    return tr;
}

/* Derivatives of energy and energy loss variance along the path (per unit path length). Two ions at the same depth with energy difference
 * delta obey d(delta)/dl = -(dS/dE) delta, therefore d(variance)/dl = W - 2 (dS/dE) variance. dS/dE is calculated numerically. */
static void reference_derivatives(const jabs_stop *stop, const jabs_stop *stragg, const ion *ion, const sample *sample, depth d, double E, double V, double *dE, double *dV) {
    const double eps = REF_DERIVATIVE_EPS * E;
    double S = stop_sample(stop, ion, sample, d, E);
    double S_E = (stop_sample(stop, ion, sample, d, E + eps) - stop_sample(stop, ion, sample, d, E - eps)) / (2.0 * eps);
    *dE = -S;
    *dV = stop_sample(stragg, ion, sample, d, E) - 2.0 * S_E * V;
}

/* Reference trajectory, independent of stop_step(). Fourth order Runge-Kutta for the system (E, variance), step taken so that
 * energy changes by approximately E_step. Steps end at range boundaries. */
static trajectory trajectory_reference(const jabs_stop *stop, const jabs_stop *stragg, double E_step, ion ion, depth d, const sample *sample, double emin) {
    trajectory tr = {0};
    emin = GSL_MAX_DBL(emin, ion.ion_gsto->emin);
    clock_t start = clock();
    trajectory_add(&tr, &d, &ion);
    while(1) {
        if((ion.inverse_cosine_theta > 0.0 && d.x >= (sample->thickness - DEPTH_TOLERANCE)) || (ion.inverse_cosine_theta < 0.0 && d.x <= DEPTH_TOLERANCE)) {
            break;
        }
        if(ion.E <= emin) {
            tr.stopped = TRUE;
            break;
        }
        depth d_next = stop_next_crossing(&ion, sample, &d);
        d.i = d_next.i;
        double h_max = (d_next.x - d.x) * ion.inverse_cosine_theta; /* Path length to next range boundary, always positive */
        double h = E_step / stop_sample(stop, &ion, sample, d, ion.E);
        depth d_half = d, d_full = d;
        if(h >= h_max) {
            h = h_max;
            d_full = d_next;
        } else {
            d_full.x = d.x + h * ion.cosine_theta;
        }
        d_half.x = (d.x + d_full.x) / 2.0;
        double E = ion.E, V = ion.S;
        double kE1, kV1, kE2, kV2, kE3, kV3, kE4, kV4;
        reference_derivatives(stop, stragg, &ion, sample, d, E, V, &kE1, &kV1);
        reference_derivatives(stop, stragg, &ion, sample, d_half, E + h / 2.0 * kE1, V + h / 2.0 * kV1, &kE2, &kV2);
        reference_derivatives(stop, stragg, &ion, sample, d_half, E + h / 2.0 * kE2, V + h / 2.0 * kV2, &kE3, &kV3);
        reference_derivatives(stop, stragg, &ion, sample, d_full, E + h * kE3, V + h * kV3, &kE4, &kV4);
        ion.E += h / 6.0 * (kE1 + 2.0 * kE2 + 2.0 * kE3 + kE4);
        ion.S += h / 6.0 * (kV1 + 2.0 * kV2 + 2.0 * kV3 + kV4);
        d = d_full;
        trajectory_add(&tr, &d, &ion);
    }
    tr.t = (double)(clock() - start) / CLOCKS_PER_SEC;
    tr.n_evals = (tr.n - 1) * PASSES_REF;
    return tr;
}

/* Linear interpolation of reference E and S at depth x. Returns FALSE if x is outside the reference trajectory. */
static int trajectory_interpolate(const trajectory *ref, double x, double *E, double *S) {
    int increasing = ref->x[ref->n - 1] > ref->x[0];
    size_t lo = 0, hi = ref->n - 1;
    if(increasing ? (x < ref->x[lo] || x > ref->x[hi]) : (x > ref->x[lo] || x < ref->x[hi])) {
        return FALSE;
    }
    while(hi - lo > 1) {
        size_t mi = (lo + hi) / 2;
        if((ref->x[mi] <= x) == increasing) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    double dx = ref->x[hi] - ref->x[lo];
    double frac = dx == 0.0 ? 0.0 : (x - ref->x[lo]) / dx;
    *E = ref->E[lo] + frac * (ref->E[hi] - ref->E[lo]);
    *S = ref->S[lo] + frac * (ref->S[hi] - ref->S[lo]);
    return TRUE;
}

typedef struct trajectory_error {
    double E_max; /* Maximum absolute energy error along path */
    double E_end; /* Energy error at the end of path */
    double sigma_max; /* Maximum absolute sigma error */
    double sigma_rel_max; /* Maximum relative sigma error (where sigma_ref > SIGMA_REL_THRESHOLD) */
    double sigma_rel_end;
} trajectory_error;

static trajectory_error trajectory_compare(const trajectory *ref, const trajectory *tr) {
    trajectory_error err = {0};
    for(size_t i = 0; i < tr->n; i++) {
        double E_ref, S_ref;
        if(!trajectory_interpolate(ref, tr->x[i], &E_ref, &S_ref)) {
            continue;
        }
        double dE = tr->E[i] - E_ref;
        double sigma_ref = sqrt(S_ref);
        double dsigma = sqrt(tr->S[i]) - sigma_ref;
        if(fabs(dE) > fabs(err.E_max)) {
            err.E_max = dE;
        }
        if(fabs(dsigma) > fabs(err.sigma_max)) {
            err.sigma_max = dsigma;
        }
        if(sigma_ref > SIGMA_REL_THRESHOLD) {
            double rel = dsigma / sigma_ref;
            if(fabs(rel) > fabs(err.sigma_rel_max)) {
                err.sigma_rel_max = rel;
            }
            err.sigma_rel_end = rel;
        }
        err.E_end = dE;
    }
    return err;
}

static int assign_stopping_sample(jibal_gsto *gsto, int Z1, const sample *sample) {
    for(size_t i = 0; i < sample->n_isotopes; i++) {
        int Z2 = sample->isotopes[i]->Z;
        if(!jibal_gsto_auto_assign(gsto, Z1, Z2)) {
            fprintf(stderr, "Could not assign stopping for Z1 = %i, Z2 = %i.\n", Z1, Z2);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

static void print_result(const char *label, double step, const trajectory *tr, const trajectory_error *err) {
    char step_str[32];
    if(step > 0.0) {
        snprintf(step_str, sizeof(step_str), "%.1f", step / C_KEV);
    } else {
        snprintf(step_str, sizeof(step_str), "auto");
    }
    fprintf(stdout, "  %-10s %8s %8zu %9zu %11.2f %11.2f %10.2f %9.4f %9.4f %9.2f  %s at %.1f tfu, E = %.3f keV", label, step_str, tr->n - 1, tr->n_evals,
            err->E_max / C_EV, err->E_end / C_EV, err->sigma_max / C_EV, 100.0 * err->sigma_rel_max, 100.0 * err->sigma_rel_end, 1e3 * tr->t,
            tr->stopped ? "stopped" : "exited", tr->x[tr->n - 1] / C_TFU, tr->E[tr->n - 1] / C_KEV);
    if(tr->n_rejected) {
        fprintf(stdout, ", %zu rejected", tr->n_rejected);
    }
    fprintf(stdout, "\n");
}

static int scenario_run(jibal *jibal, const scenario *sc) {
    sample_model *sm;
    if(strncmp(sc->sample, "file:", 5) == 0) {
        sm = sample_model_from_file(jibal, sc->sample + 5);
    } else {
        sm = sample_model_from_string(jibal, sc->sample);
    }
    if(!sm) {
        fprintf(stderr, "Could not make sample model \"%s\".\n", sc->sample);
        return EXIT_FAILURE;
    }
    sample *sample = sample_from_sample_model(sm);
    sample_model_free(sm);
    const jibal_isotope *isotope = jibal_isotope_find(jibal->isotopes, sc->ion, 0, 0);
    if(!sample || !isotope) {
        fprintf(stderr, "Could not set up scenario \"%s\".\n", sc->name);
        return EXIT_FAILURE;
    }
    if(assign_stopping_sample(jibal->gsto, isotope->Z, sample)) {
        return EXIT_FAILURE;
    }
    jibal_gsto_load_all(jibal->gsto);

    ion ion;
    ion_reset(&ion);
    ion_set_isotope(&ion, isotope);
    ion_set_angle(&ion, sc->theta, 0.0);
    ion.E = sc->E;
    ion.S = pow2(sc->sigma0);
    ion.nucl_stop = nuclear_stopping_new(isotope, jibal->isotopes);
    ion.ion_gsto = ion_gsto_new(isotope, jibal->gsto);
    depth d = depth_seek(sample, sc->x0);

    sim_calc_params *p = sim_calc_params_defaults(NULL);
    sim_calc_params *p_fast = sim_calc_params_defaults_fast(sim_calc_params_defaults(NULL));
    sim_calc_params *p_acc = sim_calc_params_defaults_accurate(sim_calc_params_defaults(NULL));
    jabs_stop stop = {.gsto = jibal->gsto, .type = GSTO_STO_TOT, .rk4 = TRUE, .emin = sc->emin, .ms = NULL};
    jabs_stop stragg = stop;
    stragg.type = GSTO_STO_STRAGG;
    jabs_stop_step_params fixed = {.step = REF_STEP, .min = 0.0, .max = 0.0, .sigmas = 0.0, .max_rel = 0.0};

    fprintf(stdout, "\n%s (sample %s, theta = %g deg, x0 = %g tfu, E0 = %g keV, sigma0 = %g keV)\n", sc->name, sc->sample, sc->theta / C_DEG, sc->x0 / C_TFU, sc->E / C_KEV, sc->sigma0 / C_KEV);
    trajectory ref = trajectory_reference(&stop, &stragg, REF_STEP, ion, d, sample, sc->emin);
    double E_lost = ref.E[0] - ref.E[ref.n - 1];
    fprintf(stdout, "  reference: %zu steps of %g keV, path ends at %g tfu with E = %.3f keV (%s), energy loss %.3f keV, sigma %.3f keV\n",
            ref.n - 1, REF_STEP / C_KEV, ref.x[ref.n - 1] / C_TFU, ref.E[ref.n - 1] / C_KEV, ref.stopped ? "stopped" : "exited", E_lost / C_KEV, sqrt(ref.S[ref.n - 1]) / C_KEV);
    fprintf(stdout, "  %-10s %8s %8s %9s %11s %11s %10s %9s %9s %9s\n", "method", "step/keV", "steps", "evals", "dE_max/eV", "dE_end/eV", "dsig_max/eV", "dsig_max%", "dsig_end%", "time/ms");

    trajectory tr = trajectory_reference(&stop, &stragg, REF_CHECK_STEP, ion, d, sample, sc->emin);
    trajectory_error err = trajectory_compare(&ref, &tr);
    print_result("ref check", REF_CHECK_STEP, &tr, &err);
    trajectory_free(&tr);

    tr = trajectory_compute(&stop, &stragg, &fixed, ion, d, sample, sc->emin); /* What stop_step() converges to */
    err = trajectory_compare(&ref, &tr);
    print_result("rk4 fine", fixed.step, &tr, &err);
    trajectory_free(&tr);

    const sim_calc_params *presets[] = {p, p_fast, p_acc};
    static const char *preset_names[] = {"default", "fast", "accurate"};
    for(size_t i = 0; i < sizeof(presets) / sizeof(presets[0]); i++) { /* Presets, adaptive steps (exiting_stop_accuracy) are used for exiting ions only, like in stop_sample_exit() */
        jabs_stop s = stop, sg = stragg;
        s.rk4 = sg.rk4 = presets[i]->rk4;
        if(sc->exiting) {
            s.accuracy = sg.accuracy = presets[i]->exiting_stop_accuracy;
            s.stragg_accuracy = sg.stragg_accuracy = EXITING_STOP_STRAGG_ACCURACY;
        }
        tr = trajectory_compute(&s, &sg, sc->exiting ? &presets[i]->exiting_stop_params : &presets[i]->incident_stop_params, ion, d, sample, sc->emin);
        err = trajectory_compare(&ref, &tr);
        print_result(preset_names[i], 0.0, &tr, &err);
        trajectory_free(&tr);
    }

    static const double steps_keV[] = {1.0, 5.0, 20.0, 50.0, 100.0, 250.0, 500.0};
    for(size_t i = 0; i < sizeof(steps_keV) / sizeof(steps_keV[0]); i++) {
        jabs_stop_step_params f = fixed;
        f.step = steps_keV[i] * C_KEV;
        if(f.step > E_lost / 2.0) {
            break;
        }
        for(int rk4 = TRUE; rk4 >= FALSE; rk4--) {
            jabs_stop s = stop, sg = stragg;
            s.rk4 = sg.rk4 = rk4;
            tr = trajectory_compute(&s, &sg, &f, ion, d, sample, sc->emin);
            err = trajectory_compare(&ref, &tr);
            print_result(rk4 ? "rk4" : "euler", f.step, &tr, &err);
            trajectory_free(&tr);
        }
    }
    trajectory_free(&ref);
    sim_calc_params_free(p);
    sim_calc_params_free(p_fast);
    sim_calc_params_free(p_acc);
    nuclear_stopping_free(ion.nucl_stop);
    ion_gsto_free(ion.ion_gsto);
    sample_free(sample);
    return EXIT_SUCCESS;
}

/* Checks spline interpolation (see ion_gsto_table_init()) of all stopping and straggling tables against linear interpolation in em: the largest relative difference,
 * overshoot (spline leaves the range of the two nodes of a segment where the data is monotonic over four nodes) and negative values. Interpolation error is estimated by
 * a leave-out test: table of every other node is used to estimate the left-out nodes. */
static void tables_check(jibal_gsto *gsto) {
    fprintf(stdout, "%-12s %5s %6s %10s %16s %9s %10s %9s %6s | leave-out rel. error (median, max): %-19s %-19s\n", "file", "type", "points", "combos", "max |spl-lin|/lin", "at Z1-Z2", "overshoots", "max ovs", "neg", "linear", "spline");
    for(size_t i_file = 0; i_file < gsto->n_files; i_file++) {
        gsto_file_t *file = &gsto->files[i_file];
        if(file->xpoints < 3) {
            continue;
        }
        jibal_gsto_assign_range(gsto, file->Z1_min, file->Z1_max, file->Z2_min, file->Z2_max, file); /* Data is loaded for assigned combinations only */
        if(!jibal_gsto_load(gsto, FALSE, file)) {
            fprintf(stderr, "Could not load %s.\n", file->name);
            continue;
        }
        double max_rel = 0.0;
        int max_Z1 = 0, max_Z2 = 0;
        size_t n_combos = 0, n_overshoot = 0, n_negative = 0;
        double max_overshoot = 0.0; /* Relative to segment value */
        size_t n_half = (file->xpoints + 1) / 2;
        double *err_lin = NULL, *err_spl = NULL;
        size_t n_err = 0;
        gsto_file_t half = *file;
        half.xpoints = n_half;
        half.xscale = GSTO_XSCALE_ARBITRARY; /* jibal_gsto_em_to_index() must not use the speedups of the original file */
        half.em = malloc(n_half * sizeof(double));
        double *half_data = malloc(n_half * sizeof(double));
        for(size_t i = 0; i < n_half; i++) {
            half.em[i] = file->em[2 * i];
        }
        for(int Z1 = file->Z1_min; Z1 <= file->Z1_max; Z1++) {
            for(int Z2 = file->Z2_min; Z2 <= file->Z2_max; Z2++) {
                const double *data = jibal_gsto_file_get_data(file, Z1, Z2);
                if(!data || data[file->xpoints / 2] <= 0.0) {
                    continue;
                }
                n_combos++;
                jabs_gsto_table t;
                ion_gsto_table_init(&t, file, data, 1.0);
                if(!t.deriv) {
                    continue;
                }
                for(size_t i = 0; i + 1 < file->xpoints; i++) {
                    int monotonic = i > 0 && i + 2 < file->xpoints && ((data[i - 1] <= data[i] && data[i] <= data[i + 1] && data[i + 1] <= data[i + 2]) || (data[i - 1] >= data[i] && data[i] >= data[i + 1] && data[i + 1] >= data[i + 2]));
                    double lo = GSL_MIN_DBL(data[i], data[i + 1]), hi = GSL_MAX_DBL(data[i], data[i + 1]);
                    int overshoot = FALSE;
                    for(int k = 1; k < 16; k++) {
                        double em = exp(t.lnem[i] + (t.lnem[i + 1] - t.lnem[i]) * k / 16.0);
                        double spl = ion_gsto_table_eval(&t, em);
                        double lin = jibal_linear_interpolation(file->em[i], file->em[i + 1], data[i], data[i + 1], em);
                        if(spl < 0.0) {
                            n_negative++;
                        }
                        if(monotonic && (spl < lo - 1e-9 * hi || spl > hi + 1e-9 * hi)) {
                            overshoot = TRUE;
                            max_overshoot = GSL_MAX_DBL(max_overshoot, GSL_MAX_DBL(lo - spl, spl - hi) / hi);
                        }
                        double rel = fabs(spl - lin) / GSL_MAX_DBL(fabs(lin), 1e-6 * fabs(data[file->xpoints / 2]));
                        if(rel > max_rel) {
                            max_rel = rel;
                            max_Z1 = Z1;
                            max_Z2 = Z2;
                        }
                    }
                    n_overshoot += overshoot;
                }
                ion_gsto_table_free(&t);
                for(size_t i = 0; i < n_half; i++) {
                    half_data[i] = data[2 * i];
                }
                jabs_gsto_table th;
                ion_gsto_table_init(&th, &half, half_data, 1.0);
                err_lin = realloc(err_lin, (n_err + file->xpoints) * sizeof(double));
                err_spl = realloc(err_spl, (n_err + file->xpoints) * sizeof(double));
                for(size_t i = 1; i + 1 < n_half * 2 - 1 && i < file->xpoints; i += 2) {
                    double em = file->em[i];
                    if(data[i] <= 0.0) {
                        continue;
                    }
                    size_t j = i / 2;
                    double lin = jibal_linear_interpolation(half.em[j], half.em[j + 1], half_data[j], half_data[j + 1], em);
                    err_lin[n_err] = fabs(lin / data[i] - 1.0);
                    err_spl[n_err] = th.deriv ? fabs(ion_gsto_table_eval(&th, em) / data[i] - 1.0) : err_lin[n_err];
                    n_err++;
                }
                ion_gsto_table_free(&th);
            }
        }
        free(half.em);
        free(half_data);
        double med_lin = 0.0, max_lin = 0.0, med_spl = 0.0, max_spl = 0.0;
        if(n_err) {
            gsl_sort(err_lin, 1, n_err);
            gsl_sort(err_spl, 1, n_err);
            med_lin = err_lin[n_err / 2];
            max_lin = err_lin[n_err - 1];
            med_spl = err_spl[n_err / 2];
            max_spl = err_spl[n_err - 1];
        }
        free(err_lin);
        free(err_spl);
        fprintf(stdout, "%-12s %5s %6zu %10zu %16.2e %5i-%-3i %10zu %9.2e %6zu | %9.2e %9.2e %9.2e %9.2e\n", file->name, file->type == GSTO_STO_STRAGG ? "strg" : "stop", file->xpoints, n_combos,
                max_rel, max_Z1, max_Z2, n_overshoot, max_overshoot, n_negative, med_lin, max_lin, med_spl, max_spl);
    }
}

int main(int argc, char **argv) {
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return EXIT_FAILURE;
    }
    if(argc > 1 && strcmp(argv[1], "--tables") == 0) {
        tables_check(jibal->gsto);
        jibal_free(jibal);
        return EXIT_SUCCESS;
    }
    fprintf(stdout, "Errors are (test - reference), evaluated at each step of the tested trajectory. Relative sigma errors where reference sigma > %g keV.\n"
                    "evals = number of passes over sample isotopes (stop_sample() calls or equivalent).\n", SIGMA_REL_THRESHOLD / C_KEV);
    size_t n = sizeof(scenarios) / sizeof(scenarios[0]);
    for(size_t i = 0; i < n; i++) {
        if(argc > 1 && !strstr(scenarios[i].name, argv[1])) { /* Optional filter by scenario name */
            continue;
        }
        if(scenario_run(jibal, &scenarios[i])) {
            fprintf(stderr, "Scenario \"%s\" failed.\n", scenarios[i].name);
        }
    }
    jibal_free(jibal);
    return EXIT_SUCCESS;
}
