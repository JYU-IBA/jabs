/*
 *  Jaakko's Backscattering Simulator (JaBS)
 *  Copyright (C) 2021 Jaakko Julin
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

int fit_function(const gsl_vector *x, void *params, gsl_vector * f)
{
    clock_t start, end;
    struct fit_data *fit_data = (struct fit_data *) params;

    for(size_t i = 0; i < fit_data->fit_params->n; i++) {
        *(fit_data->fit_params->func_params[i]) = gsl_vector_get(x, i);
    }

    sample_free(fit_data->sim->sample);
    sample_model_renormalize(fit_data->sm); /* TODO: only necessary if sample concentrations are fitted */
    fit_data->sim->sample = sample_from_sample_model(fit_data->sm);
    fit_data_workspaces_free(fit_data);
    if(sim_sanity_check(fit_data->sim) || fit_data_workspaces_init(fit_data)) { /* Either fails: clean up and return */
        fit_data_workspaces_free(fit_data); /* Some workspaces may have already been allocated */
        gsl_vector_set_all(f, 0.0);
        return GSL_FAILURE;
    }
    for(size_t i_det = 0; i_det < fit_data->sim->n_det; i_det++) {
        start = clock();
        simulate_with_ds(fit_data->ws[i_det]);
        end = clock();
        fit_data->stats.cputime_actual += (((double) (end - start)) / CLOCKS_PER_SEC);
    }
    size_t i_vec = 0;
    for(size_t i_range = 0; i_range < fit_data->n_fit_ranges; i_range++) {
        if(i_vec >= f->size) {
            jabs_message(MSG_ERROR, stderr, "Too many channels in fits for the residuals vector. This shouldn't happen.\n");
            return GSL_FAILURE;
        }

        roi *range = &fit_data->fit_ranges[i_range];
        if(range->i_det >= fit_data->sim->n_det) {
            jabs_message(MSG_ERROR, stderr, "Fit range %zu has detector %zu, but we're only supposed to have %zu detectors!\n", i_range+1, range->i_det, fit_data->sim->n_det);
            return GSL_FAILURE;
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
        return GSL_FAILURE;
    }
    return GSL_SUCCESS;
}

void fit_callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w) {
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    struct fit_data *fit_data = (struct fit_data *)params;
    double rcond;

    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&rcond, w);

    jabs_message(MSG_INFO, stderr, "iter %2zu: cond(J) = %12.6e, |f(x)| = %14.8e", iter, 1.0 / rcond, gsl_blas_dnrm2(f));
#ifndef NO_CHISQ
    double chisq;
    gsl_blas_ddot(f, f, &chisq);
    jabs_message(MSG_INFO, stderr, ", chisq/dof = %10.7lf", chisq/fit_data->dof);
#endif
#ifdef FIT_PRINT_PARAMS
    size_t i;
    for(i = 0; i < fit_data->fit_params->n; i++) {
        fprintf(stderr, ", prob[%zu] = %12.6e", i, gsl_vector_get(w->x, i));
    }
#endif
    jabs_message(MSG_INFO, stderr, "\n");
}

fit_params *fit_params_new() {
    fit_params *p = malloc(sizeof(fit_params));
    p->n = 0;
    p->func_params = NULL;
    p->func_params_err = NULL;
    return p;
}
int fit_params_add_parameter(fit_params *p, double *value) {
    if(!value) {
#ifdef DEBUG
        fprintf(stderr, "Didn't add a fit parameter since a NULL pointer was passed.\n");
#endif
        return EXIT_FAILURE;
    }
    for(size_t i = 0; i < p->n; i++) {
        if(p->func_params[i] == value) {
#ifdef DEBUG
            fprintf(stderr, "Didn't add fit parameter that points to value %p\n", (void *)value);
#endif
            return EXIT_SUCCESS; /* Parameter already exists, don't add. */
        }
    }
    p->n++;
    p->func_params = realloc(p->func_params, sizeof(double *)*p->n);
    p->func_params_err = realloc(p->func_params_err, sizeof(double)*p->n);
    p->func_params[p->n-1] = value;
    return EXIT_SUCCESS;
}
void fit_params_free(fit_params *p) {
    if(!p)
        return;
    free(p->func_params);
    free(p->func_params_err);
    free(p);
}

void fit_stats_print(FILE *f, const struct fit_stats *stats) {
    jabs_message(MSG_INFO, f,"CPU time used for actual simulation: %.3lf s.\n", stats->cputime_actual);
    if(stats->n_evals > 0) {
        jabs_message(MSG_INFO, f, "One simulation: %.3lf ms.\n", 1000.0 * stats->cputime_actual / stats->n_evals);
    }
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
    f->n_iters_max = FIT_ITERS_MAX;
    f->n_fit_ranges = 0;
    f->xtol = FIT_XTOL;
    f->gtol = FIT_GTOL;
    f->ftol = FIT_FTOL;
    f->fit_ranges = NULL;
    f->jibal = jibal;
    f->sim = sim;
    f->exp = calloc(sim->n_det, sizeof(gsl_histogram *)); /* Allocating based on initial number of detectors. */
    f->sm = NULL; /* Can be set later. */
    f->ws = NULL; /* Initialized later */
    f->fit_params = fit_params_new();
    return f;
}

void fit_data_free(fit_data *fit) {
    if(!fit)
        return;
    fit_params_free(fit->fit_params);
    fit_data_fit_ranges_free(fit);
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

    jabs_message(MSG_INFO, f,  "          low = %12zu\n", roi->low);
    jabs_message(MSG_INFO, f,  "         high = %12zu\n", roi->high);
    jabs_message(MSG_INFO, f,  "        E_low = %12.3lf keV (low energy edge of bin)\n", detector_calibrated(ws->det, JIBAL_ANY_Z, roi->low)/C_KEV);
    jabs_message(MSG_INFO, f,  "       E_high = %12.3lf keV (high energy edge of bin)\n", detector_calibrated(ws->det, JIBAL_ANY_Z, roi->high+1)/C_KEV);
    jabs_message(MSG_INFO, f,  "        n_sim = %12zu\n", n_sim);
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
    if(i_det >= fit_data->sim->n_det)
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
    size_t n = fit_data->sim->n_det;
    fit_data->ws = calloc(n, sizeof(sim_workspace *));
    if(!fit_data->ws)
        return EXIT_FAILURE;
    for(size_t i = 0; i < n; i++) {
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
    for(size_t i_det = 0; i_det < fit_data->sim->n_det; i_det++) {
        sim_workspace_free(fit_data->ws[i_det]);
        fit_data->ws[i_det] = NULL;
    }
    free(fit_data->ws);
    fit_data->ws = NULL;
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
    jabs_message(MSG_INFO, f, "  i |    low |   high |    exp cts |    sim cts | sim/exp |\n");
    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *range = &fit_data->fit_ranges[i];
        double exp_cts = spectrum_roi(fit_data->exp[range->i_det], range->low, range->high);
        double sim_cts = fit_data->ws?spectrum_roi(fit_data->ws[range->i_det]->histo_sum, range->low, range->high):0.0;
        if(exp_cts == 0.0) {
            jabs_message(MSG_INFO, f, "%3zu | %6lu | %6lu | %10.0lf | %10.1lf |         |\n", i + 1, range->low, range->high, exp_cts, sim_cts);
        } else {
            double ratio = sim_cts/exp_cts;
            jabs_message(MSG_INFO, f, "%3zu | %6lu | %6lu | %10.0lf | %10.1lf | %7.5lf |\n", i + 1, range->low, range->high, exp_cts, sim_cts, ratio);
        }
    }
    jabs_message(MSG_INFO, f, "\nFit has %zu channels total.\n", fit_data_ranges_calculate_number_of_channels(fit_data));
}

int fit(struct fit_data *fit_data) {
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    fdf_params.trs = gsl_multifit_nlinear_trs_lm;
    fdf_params.solver = gsl_multifit_nlinear_solver_svd; /* Robust? */
    struct fit_params *fit_params = fit_data->fit_params;
    if(!fit_params || fit_params->n == 0) {
        jabs_message(MSG_ERROR, stderr, "No parameters to fit.\n");
        return -1;
    }
    if(!fit_data->exp) {
        jabs_message(MSG_ERROR, stderr,  "No experimental spectrum to fit.\n");
        return -1;
    }
    gsl_multifit_nlinear_fdf fdf;
    fdf.params = fit_data;
    fit_data->stats.cputime_actual = 0.0;
    fit_data->stats.n_evals = 0;
    fit_data->stats.n_iters = 0;
    fit_data->stats.chisq_dof = 0.0;
    if(!fit_data->exp) {
        jabs_message(MSG_ERROR, stderr, "No experimental data, can not fit.\n");
        return -1;
    }
    if(!fit_data->n_fit_ranges) {
        jabs_message(MSG_ERROR, stderr, "No fit range(s) given, can not fit.\n");
        return -1;
    }

    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *range = &fit_data->fit_ranges[i];
        jabs_message(MSG_INFO,  stderr, "Fit range %zu [%lu, %lu]\n", i+1, range->low, range->high);
    }

    fdf.f = &fit_function;
    fdf.df = NULL; /* Jacobian, with NULL using finite difference. TODO: this could be implemented */
    fdf.fvv = NULL; /* No geodesic acceleration */
    fdf.n = fit_data_ranges_calculate_number_of_channels(fit_data);
    fdf.p = fit_params->n;
    if(fdf.n < fdf.p) {
        jabs_message(MSG_ERROR, stderr,"Not enough data (%zu points) for given number of free parameters (%zu)\n", fdf.n, fdf.p);
        return -1;
    } else {
        jabs_message(MSG_INFO,  stderr, "%zu channels and %zu parameters in fit.\n", fdf.n, fdf.p);
    }
    gsl_vector *f;
    gsl_matrix *J;
    double chisq, chisq0;
    fit_data->dof = fdf.n - fdf.p;
    int status, info;
    size_t i, j;

    double *weights = malloc(sizeof(double) * fdf.n);
    if(!weights)
        return 1;
    size_t i_w = 0;
    for(size_t  i_range = 0; i_range < fit_data->n_fit_ranges; i_range++) {
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
            jabs_message(MSG_ERROR, stderr,  "Experimental spectrum for detector %zu (fit range %zu) does not exist.\n", range->i_det + 1, i_range + 1);
            free(weights);
            return 1;
        }
        for(i = range->low; i <= range->high && i < det->channels; i++) {
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

    gsl_matrix *covar = gsl_matrix_alloc (fit_params->n, fit_params->n);
    gsl_vector *x = gsl_vector_alloc(fit_params->n);
    for(i = 0; i < fit_params->n; i++) {
        gsl_vector_set(x, i, *fit_params->func_params[i]); /* Initial values of fitted parameters from function parameter array */
    }

    gsl_vector_view wts = gsl_vector_view_array(weights, i_w);

/* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, fdf.n, fdf.p);

/* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (x, &wts.vector, &fdf, w);

/* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

/* solve the system with a maximum of n_iters_max iterations */
    status = gsl_multifit_nlinear_driver(fit_data->n_iters_max, fit_data->xtol, fit_data->gtol, fit_data->ftol, fit_callback, fit_data, &info, w);

/* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

/* compute final cost */
    gsl_blas_ddot(f, f, &chisq);
    jabs_message(MSG_INFO, stderr,  "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    jabs_message(MSG_INFO, stderr,  "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    fit_data->stats.n_iters = gsl_multifit_nlinear_niter(w);
    jabs_message(MSG_INFO, stderr,"function evaluations: %zu\n", fdf.nevalf);
    fit_data->stats.n_evals = fdf.nevalf;
    jabs_message(MSG_INFO, stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    jabs_message(MSG_INFO, stderr, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    jabs_message(MSG_INFO, stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
    jabs_message(MSG_INFO, stderr, "final   |f(x)| = %f\n", sqrt(chisq));
    jabs_message(MSG_INFO, stderr, "status = %s\n", gsl_strerror (status));

    double c = GSL_MAX_DBL(1, sqrt(chisq / fit_data->dof));
    jabs_message(MSG_INFO, stderr, "Final fitted parameters\n");
    for(i = 0; i < fit_params->n; i++) {
        jabs_message(MSG_INFO, stderr, "    p[%zu] = %g +- %g (%.2lf%%)\n", i,
                gsl_vector_get(w->x, i),
                c * sqrt(gsl_matrix_get(covar, i, i)),
                100.0*(c * sqrt(gsl_matrix_get(covar, i, i)))/gsl_vector_get(w->x, i)
        );
    }
    jabs_message(MSG_INFO, stderr, "Correlation coefficients table:\n     ");
    for(j = 0; j < fit_params->n; j++) {
        jabs_message(MSG_INFO, stderr, "  prob[%zu]   ", j);
    }
    jabs_message(MSG_INFO, stderr, "\n");
    for(i = 0; i < fit_params->n; i++) {
        jabs_message(MSG_INFO, stderr, "prob[%zu] ", i);
        for (j = 0; j <= i; j++) {
            jabs_message(MSG_INFO, stderr, " %8.5f", gsl_matrix_get(covar, i, j)/sqrt(gsl_matrix_get(covar, i, i)*gsl_matrix_get(covar, j, j)));
        }
        jabs_message(MSG_INFO, stderr, "\n");
    }
    fit_data->stats.chisq_dof = chisq / fit_data->dof;
    for(i = 0; i < fit_params->n; i++) { /* Clear all err values */
        fit_params->func_params_err[i] = 0.0;
    }
    for(i = 0; i < fit_params->n; i++) { /* Update final fitted values to the table (same as used for initial guess) */
        *(fit_params->func_params[i]) = gsl_vector_get(w->x, i);
        fit_params->func_params_err[i] = c * sqrt(gsl_matrix_get(covar, i, i));
    }
    for(i = 0; i < fit_data->sim->n_det; i++) {
        spectrum_set_calibration(fit_data_exp(fit_data, i), sim_det(fit_data->sim, i), JIBAL_ANY_Z); /* Update the experimental spectra to final calibration (using default calibration) */
    }
    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);
    gsl_vector_free(x);
    free(weights);
    return 0;
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
        jabs_message(MSG_ERROR, stderr,"Can not parse range from \"%s\". Is ':' missing?\n", str_orig);
        return EXIT_FAILURE;
    }
    str++; /* Skipping ':' */
    r->high = strtoull(str, &end, 10);
    str = end;
    if(*str != ']') {
        jabs_message(MSG_ERROR, stderr,"Can not parse range from \"%s\". Is ']' missing near \"%s\"?\n", str_orig, str);
        return EXIT_FAILURE;
    }
    str++;
    if(*str != '\0') {
        jabs_message(MSG_ERROR, stderr, "Unexpected input when parsing a range, \"%s\" at end of \"%s\"\n", str, str_orig);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
