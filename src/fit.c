/*
 *  Jaakko's Backscattering Simulator (JaBS)
 *  Copyright (C) 2021 - 2022 Jaakko Julin
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *   See LICENSE.txt for the full license.
 */

/*
 * This file is based a mostly around the GSL examples for nonlinear least-squares fitting. See e.g.
 * https://www.gnu.org/software/gsl/doc/html/nls.html . Please note that GSL (https://www.gnu.org/software/gsl/) is
 * is distributed under the terms of the GNU General Public License (GPL).
 */

#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fit.h"
#include "jabs.h"
#include "defaults.h"
#include "spectrum.h"
#include "message.h"

int fit_function(const gsl_vector *x, void *params, gsl_vector *f) {
    clock_t start, end;
    struct fit_data *fit_data = (struct fit_data *) params;
#ifdef DEBUG
    fprintf(stderr, "Fit iteration %zu, fit function evaluation %zu\n", fit_data->stats.iter, fit_data->stats.n_evals_iter);
#endif
    size_t i_v_active = fit_data->stats.iter_call - 1; /* Which variable is being varied by the fit algorithm in this particular call inside an iteration (compared to first call (== 0)) */
    fit_variable *var_active = fit_data->stats.iter_call && i_v_active < fit_data->fit_params->n_active ? &fit_data->fit_params->vars[i_v_active] : NULL;
    fit_data->stats.iter_call++;
#ifdef DEBUG
    fprintf(stderr, "Iter %zu (call %zu)\n", fit_data->stats.iter, fit_data->stats.iter_call);
#endif

    for(size_t i = 0; i < fit_data->fit_params->n; i++) {
        fit_variable *var = &fit_data->fit_params->vars[i];
        if(var->active) {
            *(var->value) = gsl_vector_get(x, var->i_v);
            if(fit_data->stats.iter_call == 1) { /* Store the value of fit parameters at the first function evaluation */
                var->value_iter = *(var->value);
            }
#ifdef DEBUG
            fprintf(stderr, "  %zu %12.10lf %12.10lf %s\n", var->i_v, *(var->value)/var->value_orig, *(var->value)/var->value_iter, var->i_v == i_v_active ? "THIS ONE":"");
#endif
        }
    }
#ifdef DEBUG
    fprintf(stderr, "\n");
#endif


    if(var_active && var_active->value == &fit_data->sim->fluence) {
        double scale = *(var_active->value) / var_active->value_iter;
#ifdef DEBUG
        fprintf(stderr, "Varying fluence (iter %zu, call %zu) by %12.10lf.\n", fit_data->stats.iter, fit_data->stats.iter_call, scale);
#endif
        for(size_t i = 0; i < fit_data->n_ws; i++) {
            sim_workspace *ws = fit_data_ws(fit_data, i);
            sim_workspace_histograms_scale(ws, scale);
            sim_workspace_calculate_sum_spectra(ws);
        }
        fit_set_residuals(fit_data, f);
        fit_data->stats.n_evals_iter++;
        return GSL_SUCCESS;
    }

    sample_free(fit_data->sim->sample);
    fit_data->sim->sample = NULL;
    fit_data_workspaces_free(fit_data);

    if(sample_model_sanity_check(fit_data->sm)) {
        fit_data->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }
    if(sim_sanity_check(fit_data->sim)) {
        fit_data->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }
    fit_data->sim->sample = sample_from_sample_model(fit_data->sm);
    sample_renormalize(fit_data->sim->sample);

    if(fit_data_workspaces_init(fit_data)) {
        fit_data_workspaces_free(fit_data); /* Some workspaces may have already been allocated */
        fit_data->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }

#if 0 /* If and when this is enabled, make sure to calculate one complete spectrum with original emin at the end of fit */
    for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) { /* Sets the lowest energy in each simulation according to fit ranges */
        double emin = fit_emin(fit_data, i_det);
        sim_workspace *ws = fit_data_ws(fit_data, i_det);
        if(emin > ws->emin) { /* Only increase the emin, never reduce it. */
            ws->emin = emin;
        }
    }
#endif

    start = clock();
    for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) {
        simulate_with_ds(fit_data->ws[i_det]);
    }
    end = clock();
    fit_data->stats.cputime_iter += (((double) (end - start)) / CLOCKS_PER_SEC);
    fit_data->stats.n_evals_iter++;

    if(fit_data->stats.iter_call == 1) { /* First call of iter, store sum histograms to fit_data */
        fit_data_histo_sum_free(fit_data);
        fit_data->histo_sum_iter = calloc(fit_data->n_ws, sizeof(gsl_histogram *));
        for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) {
            sim_workspace *ws = fit_data_ws(fit_data, i_det);
            fit_data->histo_sum_iter[i_det] = gsl_histogram_clone(ws->histo_sum);
        }
        fit_data->n_histo_sum = fit_data->n_ws;
    }
    fit_data->stats.error = fit_set_residuals(fit_data, f);
    if(fit_data->stats.error)
        return GSL_FAILURE;
    return GSL_SUCCESS;
}

int fit_set_residuals(const struct fit_data *fit_data, gsl_vector *f) {
    size_t i_vec = 0;
    for(size_t i_range = 0; i_range < fit_data->n_fit_ranges; i_range++) {
        /* TODO: actually calculate which spectra we need to simulate when fitting. Calculate emin based on lowest energy in lowest range. */
        if(i_vec >= f->size) {
            jabs_message(MSG_ERROR, stderr, "Too many channels in fits for the residuals vector. This shouldn't happen.\n");
            return FIT_ERROR_IMPOSSIBLE;
        }

        roi *range = &fit_data->fit_ranges[i_range];
        if(range->i_det >= fit_data->sim->n_det) {
            jabs_message(MSG_ERROR, stderr, "Fit range %zu has detector %zu, but we're only supposed to have %zu detectors!\n", i_range + 1, range->i_det, fit_data->sim->n_det);
            return FIT_ERROR_IMPOSSIBLE;
        }
        sim_workspace *ws = fit_data->ws[range->i_det];
        assert(ws);
        gsl_histogram *exp = fit_data->exp[range->i_det];
        assert(exp);
        for(size_t i = range->low; i <= range->high; i++) {
            if(i >= ws->n_channels) { /* Outside range of simulated spectrum */
                gsl_vector_set(f, i_vec, exp->bin[i]);
            } else {
                gsl_vector_set(f, i_vec, exp->bin[i] - ws->histo_sum->bin[i]);
            }
            i_vec++;
        }
    }
    if(i_vec != f->size) {
        jabs_message(MSG_ERROR, stderr, "Not enough channels in fits for the residuals vector. This shouldn't happen.\n");
        return FIT_ERROR_IMPOSSIBLE;
    }
    return FIT_ERROR_NONE;
}

void fit_iter_stats_update(struct fit_data *fit_data, const gsl_multifit_nlinear_workspace *w) {
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&fit_data->stats.rcond, w);
    gsl_blas_ddot(f, f, &fit_data->stats.chisq);
    fit_data->stats.norm = gsl_blas_dnrm2(f);
    fit_data->stats.chisq_dof = fit_data->stats.chisq / fit_data->dof;
    fit_data->stats.n_evals += fit_data->stats.n_evals_iter;
    fit_data->stats.cputime_cumul += fit_data->stats.cputime_iter;
}

void fit_iter_stats_print(const struct fit_stats *stats) {
    jabs_message(MSG_INFO, stderr, "%4zu | %12.6e | %14.8e | %12.7lf | %4zu | %13.3lf | %12.1lf |\n",
                 stats->iter, 1.0 / stats->rcond, stats->norm,
                 stats->chisq_dof, stats->n_evals, stats->cputime_cumul,
                 1000.0 * stats->cputime_iter / stats->n_evals_iter);
}

fit_params *fit_params_new() {
    fit_params *p = malloc(sizeof(fit_params));
    p->n = 0;
    p->n_active = 0;
    p->vars = NULL;
    return p;
}

int fit_params_add_parameter(fit_params *p, double *value, const char *name, const char *unit, double unit_factor) {
    if(!value || !name) {
#ifdef DEBUG
        fprintf(stderr, "Didn't add a fit parameter since a NULL pointer was passed. Value = %p, name = %p (%s).\n", (void *)value, (void *)name, name);
#endif
        return EXIT_FAILURE;
    }
    for(size_t i = 0; i < p->n; i++) {
        if(p->vars[i].value == value) {
#ifdef DEBUG
            fprintf(stderr, "Didn't add fit parameter %s that points to value %p since it already exists.\n", name, (void *)value);
#endif
            return EXIT_SUCCESS; /* Parameter already exists, don't add. */ /* TODO: maybe this requirement could be relaxed? */
        }
    }
    p->n++;
    p->vars = realloc(p->vars, sizeof(fit_variable) * p->n);
    fit_variable *var = &(p->vars[p->n - 1]);
    var->value = value;
    var->name = strdup(name);
    var->unit = unit;
    var->unit_factor = unit_factor;
    var->active = FALSE;
#ifdef DEBUG
    fprintf(stderr, "Fit parameter %s added successfully (total %zu).\n", var->name, p->n);
#endif
    return EXIT_SUCCESS;
}

void fit_params_free(fit_params *p) {
    if(!p)
        return;
    for(size_t i = 0; i < p->n; i++) {
        free(p->vars[i].name); /* This should be fit_variable_free(), but then p->vars should probably be an array of pointers, too). */
    }
    free(p->vars);
    free(p);
}

void fit_params_update(fit_params *p) {
    if(!p)
        return;

    /* TODO: check if some ACTIVE parameters are duplicates. We don't care about inactive ones. */

    p->n_active = 0;
    for(size_t i = 0; i < p->n; i++) {
        if(p->vars[i].active) {
            p->vars[i].i_v = p->n_active;
            p->n_active++;
        } else {
            p->vars[i].i_v = p->n; /* Intentionally invalid index */
        }
    }
}

void fit_stats_print(FILE *f, const struct fit_stats *stats) {
    if(stats->chisq_dof > 0.0) {
        jabs_message(MSG_INFO, f, "Final chisq/dof = %.7lf\n", stats->chisq_dof);
    }
}

int fit_data_fit_range_add(struct fit_data *fit_data, const struct roi *range) { /* Makes a deep copy */
    if(range->low == 0 && range->high == 0) {
        return EXIT_FAILURE;
    }
    if(range->high < range->low) {
        jabs_message(MSG_ERROR, stderr, "Range from %zu to %zu is not valid!\n", range->low, range->high);
        return EXIT_FAILURE;
    }
    fit_data->n_fit_ranges++;
    fit_data->fit_ranges = realloc(fit_data->fit_ranges, fit_data->n_fit_ranges * sizeof(roi));
    if(!fit_data->fit_ranges) {
        fit_data->n_fit_ranges = 0;
        return EXIT_FAILURE;
    }
    fit_data->fit_ranges[fit_data->n_fit_ranges - 1] = *range;
    return EXIT_SUCCESS;
}

void fit_data_fit_ranges_free(struct fit_data *fit_data) {
    if(!fit_data)
        return;
    free(fit_data->fit_ranges);
    fit_data->fit_ranges = NULL;
    fit_data->n_fit_ranges = 0;
}

fit_data *fit_data_new(const jibal *jibal, simulation *sim) {
    struct fit_data *f = malloc(sizeof(struct fit_data));
    f->jibal = jibal;
    f->sim = sim;
    f->exp = calloc(sim->n_det, sizeof(gsl_histogram *)); /* Allocating based on initial number of detectors. */
    f->sm = NULL; /* Can be set later. */
    f->ws = NULL; /* Initialized later */
    f->n_ws = 0; /* Number of allocated workspaces, initially same as number of detectors. */
    f->fit_params = NULL; /* Holds fit parameters AFTER a fit. */
    f->histo_sum_iter = NULL; /* Initialized later */
    f->fit_iter_callback = NULL; /* Optional */
    f->n_fit_ranges = 0;
    f->fit_ranges = NULL;
    fit_data_defaults(f);
    return f;
}

void fit_data_defaults(fit_data *f) {
    f->n_iters_max = FIT_ITERS_MAX;
    f->xtol = FIT_XTOL;
    f->gtol = FIT_GTOL;
    f->ftol = FIT_FTOL;
    f->chisq_tol = FIT_CHISQ_TOL;
    f->chisq_fast_tol = FIT_FAST_CHISQ_TOL;
    f->phase_start = FIT_PHASE_FAST;
    f->phase_stop = FIT_PHASE_SLOW;
}

void fit_data_free(fit_data *fit) {
    if(!fit)
        return;
    fit_params_free(fit->fit_params);
    fit_data_fit_ranges_free(fit);
    fit_data_histo_sum_free(fit);
    free(fit);
}

void fit_data_roi_print(FILE *f, const struct fit_data *fit_data, const struct roi *roi) {
    if(!fit_data) {
        return;
    }
    sim_workspace *ws = fit_data_ws(fit_data, roi->i_det);
    gsl_histogram *exp = fit_data_exp(fit_data, roi->i_det);
    if(!ws)
        return;
    size_t n_exp = spectrum_channels_in_range(exp, roi->low, roi->high);
    size_t n_sim = spectrum_channels_in_range(ws->histo_sum, roi->low, roi->high);
    double exp_cts = spectrum_roi(exp, roi->low, roi->high);
    double sim_cts = spectrum_roi(ws->histo_sum, roi->low, roi->high);

    jabs_message(MSG_INFO, f, "          low = %12zu\n", roi->low);
    jabs_message(MSG_INFO, f, "         high = %12zu\n", roi->high);
    jabs_message(MSG_INFO, f, "        E_low = %12.3lf keV (low energy edge of bin)\n", detector_calibrated(ws->det, JIBAL_ANY_Z, roi->low) / C_KEV);
    jabs_message(MSG_INFO, f, "       E_high = %12.3lf keV (high energy edge of bin)\n", detector_calibrated(ws->det, JIBAL_ANY_Z, roi->high + 1) / C_KEV);
    jabs_message(MSG_INFO, f, "        n_sim = %12zu\n", n_sim);
    jabs_message(MSG_INFO, f, "          sim  = %12g\n", sim_cts);
    if(exp) {
        jabs_message(MSG_INFO, f, "        n_exp = %12zu\n", n_exp);
        jabs_message(MSG_INFO, f, "          exp  = %12g\n", exp_cts);
        jabs_message(MSG_INFO, f, "      exp-sim  = %12g\n", exp_cts - sim_cts);
        jabs_message(MSG_INFO, f, "    sqrt(exp)  = %12.5lf\n", sqrt(exp_cts));
        jabs_message(MSG_INFO, f, "      sim/exp  = %12.5lf\n", sim_cts / exp_cts);
        jabs_message(MSG_INFO, f, "      exp/sim  = %12.5lf\n", exp_cts / sim_cts);
        jabs_message(MSG_INFO, f, "  1/sqrt(exp)  = %12.5lf%%\n", 100.0 / sqrt(exp_cts));
        jabs_message(MSG_INFO, f, "(exp-sim)/exp  = %12.5lf%%\n", 100.0 * (exp_cts - sim_cts) / exp_cts);
    }
}

gsl_histogram *fit_data_exp(const struct fit_data *fit_data, size_t i_det) {
    if(!fit_data || !fit_data->exp)
        return NULL;
    if(i_det >= fit_data->sim->n_det)
        return NULL;
    return fit_data->exp[i_det];
}

gsl_histogram *fit_data_sim(const struct fit_data *fit_data, size_t i_det) {
    sim_workspace *ws = fit_data_ws(fit_data, i_det);
    if(ws)
        return ws->histo_sum;
    else
        return NULL;
}

void fit_data_exp_free(struct fit_data *fit_data) {
    if(!fit_data->exp)
        return;
    for(size_t i_det = 0; i_det < fit_data->sim->n_det; i_det++) {
        if(fit_data->exp[i_det]) {
            gsl_histogram_free(fit_data->exp[i_det]);
            fit_data->exp[i_det] = NULL;
        }
    }
    free(fit_data->exp);
    fit_data->exp = NULL;
}

int fit_data_load_exp(struct fit_data *fit, size_t i_det, const char *filename) {
    gsl_histogram *h = spectrum_read_detector(filename, sim_det(fit->sim, i_det));
    if(!h) {
        jabs_message(MSG_ERROR, stderr, "Reading spectrum from file \"%s\" was not successful.\n", filename);
        return EXIT_FAILURE;
    }
    if(fit->exp[i_det]) {
        gsl_histogram_free(fit->exp[i_det]);
    }
    fit->exp[i_det] = h;
    return EXIT_SUCCESS;
}

void fit_data_histo_sum_free(struct fit_data *fit_data) {
    if(!fit_data->histo_sum_iter)
        return;
    for(size_t i = 0; i < fit_data->n_histo_sum; i++) {
        if(fit_data->histo_sum_iter[i]) {
            gsl_histogram_free(fit_data->histo_sum_iter[i]);
        }
    }
    free(fit_data->histo_sum_iter);
    fit_data->histo_sum_iter = NULL;
    fit_data->n_histo_sum = 0;
}

int fit_data_add_det(struct fit_data *fit_data, detector *det) {
    if(!fit_data || !det)
        return EXIT_FAILURE;
    if(sim_det_add(fit_data->sim, det)) {
        return EXIT_FAILURE;
    }
    fit_data->exp = realloc(fit_data->exp, sizeof(gsl_histogram *) * fit_data->sim->n_det);
    fit_data->exp[fit_data->sim->n_det - 1] = NULL;
    return EXIT_SUCCESS;
}

sim_workspace *fit_data_ws(const struct fit_data *fit_data, size_t i_det) {
    if(!fit_data || !fit_data->ws)
        return NULL;
    if(i_det >= fit_data->n_ws)
        return NULL;
    return fit_data->ws[i_det];
}


size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data) {
    size_t sum = 0;
    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *r = &fit_data->fit_ranges[i];
        detector *det = sim_det(fit_data->sim, r->i_det);
        if(!det) {
            continue;
        }
        if(r->high >= det->channels) { /* Limited by detector */
            sum += (det->channels - r->low);
        } else {
            sum += (r->high - r->low) + 1;
        }
    }
    return sum;
}

int fit_data_workspaces_init(struct fit_data *fit_data) {
    int status = EXIT_SUCCESS;
    fit_data_workspaces_free(fit_data);
    fit_data->n_ws = fit_data->sim->n_det;
    fit_data->ws = calloc(fit_data->n_ws, sizeof(sim_workspace *));
    if(!fit_data->ws)
        return EXIT_FAILURE;
    for(size_t i = 0; i < fit_data->n_ws; i++) {
        detector *det = sim_det(fit_data->sim, i);
        if(detector_sanity_check(det)) {
            jabs_message(MSG_ERROR, stderr, "Detector %zu failed sanity check!\n", i);
            status = EXIT_FAILURE;
            break;
        }
        detector_update(det);
        sim_workspace *ws = sim_workspace_init(fit_data->jibal, fit_data->sim, det);
        if(!ws) {
            jabs_message(MSG_ERROR, stderr, "Workspace %zu failed to initialize!\n", i);
            status = EXIT_FAILURE;
        }
        fit_data->ws[i] = ws;
    }
    if(status == EXIT_FAILURE) {
        fit_data_workspaces_free(fit_data);
    }
    return status;
}

void fit_data_workspaces_free(struct fit_data *fit_data) {
    if(!fit_data->ws) {
        return;
    }
    for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) {
        sim_workspace_free(fit_data->ws[i_det]);
        fit_data->ws[i_det] = NULL;
    }
    free(fit_data->ws);
    fit_data->ws = NULL;
    fit_data->n_ws = 0;
}

struct fit_stats fit_stats_init() {
    struct fit_stats s;
    s.n_evals = 0;
    s.n_evals_iter = 0;
    s.cputime_cumul = 0.0;
    s.cputime_iter = 0.0;
    s.chisq0 = 0.0;
    s.chisq = 0.0;
    s.chisq_dof = 0.0;
    s.rcond = 0.0;
    s.iter = 0;
    s.iter_call = 0;
    s.error = FIT_ERROR_NONE;
    return s;
}

void fit_data_print(FILE *f, const struct fit_data *fit_data) {
    if(!fit_data) {
        return;
    }
    if(fit_data->n_fit_ranges == 0) {
        jabs_message(MSG_ERROR, f, "No fit ranges.\n");
        return;
    }
    jabs_message(MSG_INFO, f, "%zu fit ranges:\n", fit_data->n_fit_ranges);
    jabs_message(MSG_INFO, f, "  i |    low |   high |    exp cts |    sim cts | sim/exp |  sigmas |\n");
    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *range = &fit_data->fit_ranges[i];
        double exp_cts = spectrum_roi(fit_data_exp(fit_data, range->i_det), range->low, range->high);
        double sim_cts = spectrum_roi(fit_data_sim(fit_data, range->i_det), range->low, range->high);
        if(exp_cts == 0.0) {
            jabs_message(MSG_INFO, f, "%3zu | %6lu | %6lu | %10.0lf | %10.1lf |         |         |\n", i + 1, range->low, range->high, exp_cts, sim_cts);
        } else {
            double ratio = sim_cts / exp_cts;
            double sigmas = (sim_cts - exp_cts) / sqrt(exp_cts);
            jabs_message(MSG_INFO, f, "%3zu | %6lu | %6lu | %10.0lf | %10.1lf | %7.5lf | %7.2lf |\n", i + 1, range->low, range->high, exp_cts, sim_cts, ratio, sigmas);
        }
    }
    jabs_message(MSG_INFO, f, "\nFit has %zu channels total.\n", fit_data_ranges_calculate_number_of_channels(fit_data));
}

int jabs_test_delta(const gsl_vector *dx, const gsl_vector *x, double epsabs, double epsrel) { /* test_delta() copied from GSL convergence.c and modified */
    int ok = TRUE;
    for(size_t i = 0; i < x->size; i++) {
        double xi = gsl_vector_get(x, i);
        double dxi = gsl_vector_get(dx, i);
        double tolerance = epsabs + epsrel * fabs(xi);
        double rel = fabs(dxi) / tolerance; /* "How many times over the acceptable tolerance are we */
#ifdef DEBUG
        fprintf(stderr, "Test delta: i %zu, xi %g, dxi %g, tolerance %g, rel %g\n", i, xi, dxi, tolerance, rel);
#endif
        if(rel >= 1.0) {
            ok = FALSE;
            break;
        }
    }
    if(ok)
        return GSL_SUCCESS;
    return GSL_CONTINUE;
}

int jabs_gsl_multifit_nlinear_driver(const size_t maxiter, const double xtol, const double chisq_tol, struct fit_data *fit_data, gsl_multifit_nlinear_workspace *w) {
    int status = 0;
    size_t iter;
    double chisq_dof_old;
    jabs_message(MSG_INFO, stderr, "iter |    cond(J)   |     |f(x)|     |   chisq/dof  | eval | cpu cumul (s) | cpu/eval (ms)|\n");
    for(iter = 0; iter <= maxiter; iter++) {
        fit_data->stats.iter_call = 0;
        fit_data->stats.iter = iter;
        if(iter) {
            chisq_dof_old = fit_data->stats.chisq_dof;
            fit_data->stats.cputime_iter = 0.0;
            fit_data->stats.n_evals_iter = 0;
            status = gsl_multifit_nlinear_iterate(w);
#ifdef DEBUG
            fprintf(stderr, "Iteration status %i (%s)\n", status, gsl_strerror(status));
#endif
        }
        if(fit_data->stats.error) {
            return fit_data->stats.error;
        }
        if(status == GSL_ENOPROG && iter == 1) {
            return FIT_ERROR_NO_PROGRESS;
        }
        fit_iter_stats_update(fit_data, w);
        fit_iter_stats_print(&fit_data->stats);
        if(fit_data->fit_iter_callback) {
            if(fit_data->fit_iter_callback(fit_data->stats)) {
                return FIT_ERROR_ABORTED;
            }
        }

        if(iter == 0)
            continue;
        /* test for convergence */
        status = jabs_test_delta(w->dx, w->x, xtol * xtol, xtol);
        if(status == GSL_SUCCESS) {
            return FIT_SUCCESS_DELTA;
        }
        double chisq_change = 1.0 - fit_data->stats.chisq_dof / chisq_dof_old;
        if(fit_data->stats.chisq_dof > chisq_dof_old) {
            jabs_message(MSG_WARNING, stderr, "Chisq increased, this probably shouldn't happen.\n");
        }
        if(chisq_change < chisq_tol) {
            return FIT_SUCCESS_CHISQ;
        }
    }
    return FIT_ERROR_MAXITER;
}

void fit_report_results(const fit_data *fit, const gsl_multifit_nlinear_workspace *w, const gsl_multifit_nlinear_fdf *fdf) {
    jabs_message(MSG_INFO, stderr, "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    jabs_message(MSG_INFO, stderr, "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    jabs_message(MSG_INFO, stderr, "function evaluations: %zu\n", fit->stats.n_evals);
#ifdef DEBUG
    jabs_message(MSG_INFO, stderr, "function evaluations (GSL): %zu\n", fdf->nevalf);
#endif
    jabs_message(MSG_INFO, stderr, "Jacobian evaluations: %zu\n", fdf->nevaldf);
    jabs_message(MSG_INFO, stderr, "reason for stopping: %s\n", fit_error_str(fit->stats.error));
    jabs_message(MSG_INFO, stderr, "initial |f(x)| = %f\n", sqrt(fit->stats.chisq0));
    jabs_message(MSG_INFO, stderr, "final   |f(x)| = %f\n", sqrt(fit->stats.chisq));
}


void fit_parameters_update(const fit_data *fit, const gsl_multifit_nlinear_workspace *w, const gsl_matrix *covar) {
    const fit_params *fit_params = fit->fit_params;
    double c = GSL_MAX_DBL(1, sqrt(fit->stats.chisq_dof));
    for(size_t i = 0; i < fit_params->n; i++) { /* Update final fitted values to the table (same as used for initial guess) */
        fit_variable *var = &(fit_params->vars[i]);
        if(!var->active)
            continue;
        assert(var->i_v < fit_params->n_active);
        var->value_final = gsl_vector_get(w->x, var->i_v);
        *(var->value) = var->value_final;
        var->err = c * sqrt(gsl_matrix_get(covar, var->i_v, var->i_v));
        var->err_rel = fabs(var->err / var->value_final);
        var->sigmas = fabs(var->value_final - var->value_orig) / var->value_orig / var->err_rel;
    }
}

void fit_parameters_update_changed(const fit_data *fit) {
    const fit_params *fit_params = fit->fit_params;
    for(size_t i = 0; i < fit_params->n; i++) {
        fit_variable *var = &(fit_params->vars[i]);
        if(!var->active)
            continue;
        if(*(var->value) != var->value_final) { /* Values changed by something (renormalization) */
            double scale = *(var->value) / var->value_final;
            var->value_final = (*var->value);
            var->err *= scale;
            /* var->err_rel stays the same */
        }
    }
}

void fit_covar_print(const gsl_matrix *covar) {
    jabs_message(MSG_INFO, stderr, "\nCorrelation coefficients matrix:\n       | ");
    for(size_t i = 0; i < covar->size1; i++) {
        jabs_message(MSG_INFO, stderr, " %4zu  ", i + 1);
    }
    jabs_message(MSG_INFO, stderr, "\n");
    for(size_t i = 0; i < covar->size1; i++) {
        jabs_message(MSG_INFO, stderr, "%6zu | ", i + 1);
        for(size_t j = 0; j <= i && j < covar->size2; j++) {
            jabs_message(MSG_INFO, stderr, " %6.3f", gsl_matrix_get(covar, i, j) / sqrt(gsl_matrix_get(covar, i, i) * gsl_matrix_get(covar, j, j)));
        }
        jabs_message(MSG_INFO, stderr, "\n");
    }
}

int fit(struct fit_data *fit_data) {
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    fdf_params.trs = gsl_multifit_nlinear_trs_lm;
    fdf_params.solver = gsl_multifit_nlinear_solver_svd; /* Robust? */
    struct fit_params *fit_params = fit_data->fit_params;
    if(!fit_params || fit_params->n_active == 0) {
        jabs_message(MSG_ERROR, stderr, "No parameters to fit.\n");
        return EXIT_FAILURE;
    }
    if(!fit_data->exp) {
        jabs_message(MSG_ERROR, stderr, "No experimental spectrum to fit.\n");
        return EXIT_FAILURE;
    }
    gsl_multifit_nlinear_fdf fdf;
    fdf.params = fit_data;
    if(!fit_data->exp) {
        jabs_message(MSG_ERROR, stderr, "No experimental data, can not fit.\n");
        return EXIT_FAILURE;
    }
    if(!fit_data->n_fit_ranges) {
        jabs_message(MSG_ERROR, stderr, "No fit range(s) given, can not fit.\n");
        return EXIT_FAILURE;
    }

    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *range = &fit_data->fit_ranges[i];
        jabs_message(MSG_INFO, stderr, "Fit range %zu [%lu:%lu]\n", i + 1, range->low, range->high);
    }

    fdf.f = &fit_function;
    fdf.df = NULL; /* Jacobian, with NULL using finite difference. TODO: this could be implemented */
    fdf.fvv = NULL; /* No geodesic acceleration */
    fdf.n = fit_data_ranges_calculate_number_of_channels(fit_data);
    fdf.p = fit_params->n_active;
    if(fdf.n < fdf.p) {
        jabs_message(MSG_ERROR, stderr, "Not enough data (%zu points) for given number of free parameters (%zu)\n", fdf.n, fdf.p);
        return -1;
    } else {
        jabs_message(MSG_INFO, stderr, "%zu channels and %zu parameters in fit, %zu degrees of freedom.\n", fdf.n, fdf.p, fdf.n - fdf.p);
    }
    gsl_vector *f;
    gsl_matrix *J;
    fit_data->dof = fdf.n - fdf.p;
    int status;

    double *weights = malloc(sizeof(double) * fdf.n);
    if(!weights)
        return 1;
    size_t i_w = 0;
    for(size_t i_range = 0; i_range < fit_data->n_fit_ranges; i_range++) {
        roi *range = &fit_data->fit_ranges[i_range];
        assert(range);
        detector *det = sim_det(fit_data->sim, range->i_det);
        gsl_histogram *exp = fit_data_exp(fit_data, range->i_det);
        if(!det) {
            jabs_message(MSG_ERROR, stderr, "Detector %zu (fit range %zu) does not exist.\n", range->i_det + 1, i_range + 1);
            free(weights);
            return 1;
        }
        if(!exp) {
            jabs_message(MSG_ERROR, stderr, "Experimental spectrum for detector %zu (fit range %zu) does not exist.\n", range->i_det + 1, i_range + 1);
            free(weights);
            return 1;
        }
        for(size_t i = range->low; i <= range->high && i < det->channels; i++) {
            if(exp->bin[i] > 1.0) {
                weights[i_w] = 1.0 / (exp->bin[i]);
            } else {
                weights[i_w] = 1.0; /* TODO: ?*/
            }
            i_w++;
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Set %zu weights.\n", i_w);
#endif
    assert(i_w == fdf.n);

    gsl_matrix *covar = gsl_matrix_alloc(fit_params->n_active, fit_params->n_active);
    gsl_vector *x = gsl_vector_alloc(fit_params->n_active);
    for(size_t i = 0; i < fit_params->n; i++) { /* Update all (including inactives) */
        fit_variable *var = &(fit_params->vars[i]);
        var->err = 0.0;
        var->value_orig = *(var->value);
    }

    gsl_vector_view wts = gsl_vector_view_array(weights, i_w);

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc(T, &fdf_params, fdf.n, fdf.p);

    sim_calc_params p_orig = *fit_data->sim->params; /* Store original values (will be used in final stage of fitting) */
    for(int phase = fit_data->phase_start; phase <= fit_data->phase_stop; phase++) { /* Phase 1 is "fast", phase 2 normal. */
        assert(phase >= FIT_PHASE_FAST && phase <= FIT_PHASE_SLOW);
        double xtol = fit_data->xtol;
        double chisq_tol = fit_data->chisq_tol;
        /* initialize solver with starting point and weights */
        fit_data->stats = fit_stats_init();
        fit_data->stats.phase = phase;
        if(fit_data->fit_iter_callback) { /* First call to callback quickly (before most initialization) */
            if(fit_data->fit_iter_callback(fit_data->stats)) {
                fit_data->stats.error = FIT_ERROR_ABORTED;
                break;
            }
        }
        if(phase == FIT_PHASE_FAST) {
            sim_calc_params_fast(fit_data->sim->params, TRUE); /* Set current parameters to be faster in phase 0. */
            xtol *= FIT_FAST_XTOL_MULTIPLIER;
            chisq_tol = fit_data->chisq_fast_tol;
        } else if(phase == FIT_PHASE_SLOW) {
            *fit_data->sim->params = p_orig; /* Restore original parameters in phase 1 */
        }
        sim_calc_params_update(fit_data->sim->params);
        for(size_t i = 0; i < fit_params->n; i++) { /* Set active variables to vector */
            fit_variable *var = &(fit_params->vars[i]);
            if(var->active) {
                gsl_vector_set(x, var->i_v, *(var->value)); /* Start phase from initial or fitted (previous phase) results */
            }
        }
        jabs_message(MSG_INFO, stderr, "\nInitializing fit phase %i. Xtol = %e, chisq_tol %e\n", phase, xtol, chisq_tol);
        gsl_multifit_nlinear_winit(x, &wts.vector, &fdf, w);

        /* compute initial cost function */
        f = gsl_multifit_nlinear_residual(w);
        gsl_blas_ddot(f, f, &fit_data->stats.chisq0);

        status = jabs_gsl_multifit_nlinear_driver(fit_data->n_iters_max, xtol, chisq_tol, fit_data, w); /* Fit */
        fit_data->stats.error = status;
        if(status < 0) {
            jabs_message(MSG_ERROR, stderr, "Fit aborted in phase %i, reason: %s.\n", phase, fit_error_str(fit_data->stats.error));
            break;
        }
        jabs_message(MSG_INFO, stderr, "Phase %i finished. CPU time used for actual simulation so far: %.3lf s.\n", phase, fit_data->stats.cputime_cumul);
        fit_report_results(fit_data, w, &fdf);
    }

    if(fit_data->stats.error < 0) { /* Revert changes on error */
        for(size_t i = 0; i < fit_params->n; i++) {
            fit_variable *var = &(fit_params->vars[i]);
            *(var->value) = var->value_orig;
        }
    } else { /* Do final calculations when fit was successful */
        /* compute covariance of best fit parameters */
        J = gsl_multifit_nlinear_jac(w);
        gsl_multifit_nlinear_covar(J, 0.0, covar);

        /* compute final cost */
        gsl_blas_ddot(f, f, &fit_data->stats.chisq);
        fit_data->stats.chisq_dof = fit_data->stats.chisq / fit_data->dof;

        fit_parameters_update(fit_data, w, covar);
        sample_model_renormalize(fit_data->sm);
        fit_parameters_update_changed(fit_data); /* sample_model_renormalize() can and will change concentration values, this will recompute error (assuming relative error stays the same) */
        fit_params_print_final(fit_params);
        fit_covar_print(covar);

        for(size_t i = 0; i < fit_data->n_ws; i++) {
            spectrum_set_calibration(fit_data_exp(fit_data, i), sim_det(fit_data->sim, i), JIBAL_ANY_Z); /* Update the experimental spectra to final calibration (using default calibration) */
        }
    }
    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);
    gsl_vector_free(x);
    free(weights);
    return fit_data->stats.error;
}

int fit_set_roi_from_string(roi *r, const char *str) {
    const char *str_orig = str;
    if(!str)
        return EXIT_FAILURE;
    if(*str != '[') { /* Silent failure, intentionally, since this can signal end of valid roi range arguments. */
#ifdef DEBUG
        fprintf(stderr, "fit_set_roi_from_string() fails silently.\n");
#endif
        return EXIT_FAILURE;
    }
    str++; /* Skipping '[' */
    char *end;
    r->low = strtoull(str, &end, 10);
    str = end;
    if(*str != ':') {
        jabs_message(MSG_ERROR, stderr, "Can not parse range from \"%s\". Is ':' missing?\n", str_orig);
        return EXIT_FAILURE;
    }
    str++; /* Skipping ':' */
    r->high = strtoull(str, &end, 10);
    str = end;
    if(*str != ']') {
        jabs_message(MSG_ERROR, stderr, "Can not parse range from \"%s\". Is ']' missing near \"%s\"?\n", str_orig, str);
        return EXIT_FAILURE;
    }
    str++;
    if(*str != '\0') {
        jabs_message(MSG_ERROR, stderr, "Unexpected input when parsing a range, \"%s\" at end of \"%s\"\n", str, str_orig);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

double fit_emin(struct fit_data *fit, size_t i_det) {
    double emin = fit->sim->beam_E;
    for(size_t i_range = 0; i_range < fit->n_fit_ranges; i_range++) {
        const roi *r = &(fit->fit_ranges[i_range]);
        if(r->i_det != i_det)
            continue;
        const detector *det = sim_det(fit->sim, i_det);
        for(int Z = JIBAL_ANY_Z; Z <= det->cal_Z_max; Z++) {
            double E = detector_calibrated(det, Z, r->low); /* TODO: assumes calibration is increasing monotonously. Is it guaranteed in all cases? */
            E -= 3.0 * det->calibration->resolution; /* TODO: bad approximation. */
            E *= 0.95;
            E -= 10.0 * C_KEV;
            if(E < emin) {
                emin = E;
            }
        }
    }
    return emin;
}

const char *fit_error_str(int error) {
    switch(error) {
        case FIT_SUCCESS_CHISQ:
            return "chi squared change below tolerance";
        case FIT_SUCCESS_DELTA:
            return "step size below tolerance";
        case FIT_ERROR_NONE:
            return "success";
        case FIT_ERROR_GENERIC:
            return "generic error";
        case FIT_ERROR_MAXITER:
            return "maximum number of iterations reached";
        case FIT_ERROR_NO_PROGRESS:
            return "iteration is not making progress";
        case FIT_ERROR_SANITY:
            return "simulation failed sanity check";
        case FIT_ERROR_IMPOSSIBLE:
            return "an impossible thing has happened";
        case FIT_ERROR_ABORTED:
            return "user requested abort";
        default:
            return "unknown";
    }
}
