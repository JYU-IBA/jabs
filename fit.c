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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fit.h"
#include "jabs.h"
#include "defaults.h"

int fit_function(const gsl_vector *x, void *params, gsl_vector * f)
{
    clock_t start, end;
    struct fit_data *p = (struct fit_data *) params;
    gsl_histogram *exp = p->exp;
    for(size_t i = 0; i < p->fit_params->n; i++) {
        *(p->fit_params->func_params[i]) = gsl_vector_get(x, i); /* TODO: set some parameters to fit! */
    }
    sim_workspace_free(p->ws);
    sample_free(p->sample);
    p->sample = sample_from_sample_model(p->sm);
    if(sim_sanity_check(p->sim)) {
        p->ws = NULL; /* Workspace has not been allocated yet. */
    } else {
        p->ws = sim_workspace_init(p->sim, p->reactions, p->sample, p->jibal); /* We intentionally "leak" this */
    }
    if(!p->ws) {
        gsl_vector_set_all(f, 0.0);
        return GSL_FAILURE;
    }
    start = clock();
    simulate_with_ds(p->ws);
    end = clock();
    p->stats->cputime_actual += (((double) (end - start)) / CLOCKS_PER_SEC);
    for(size_t i = p->low_ch; i <= p->high_ch; i++) {
        double sum = 0.0;
        if(i >= p->ws->n_channels) { /* Outside range of simulated spectrum */
            gsl_vector_set(f, i-p->low_ch, exp->bin[i]);
        } else {
            for (size_t j = 0; j < p->ws->n_reactions; j++) { /* Sum comes always first, which means we have to compute it first. */
                if(p->ws->reactions[j].histo && i < p->ws->reactions[j].histo->n)
                    sum += p->ws->reactions[j].histo->bin[i];
            }
            gsl_vector_set(f, i-p->low_ch, exp->bin[i] - sum);
        }
    }
    return GSL_SUCCESS;
}

void fit_callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w) {
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    struct fit_data *fit_data = (struct fit_data *)params;
    double rcond;

    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stderr, "iter %2zu: cond(J) = %12.6e, |f(x)| = %14.8e", iter, 1.0 / rcond, gsl_blas_dnrm2(f));
#ifndef NO_CHISQ
    double chisq;
    gsl_blas_ddot(f, f, &chisq);
    fprintf(stderr, ", chisq/dof = %8.4lf", chisq/fit_data->dof);
#endif
#ifdef FIT_PRINT_PARAMS
    size_t i;
    for(i = 0; i < fit_data->fit_params->n; i++) {
        fprintf(stderr, ", prob[%zu] = %12.6e", i, gsl_vector_get(w->x, i));
    }
#endif
    if(fit_data->print_iters) {
        print_spectra(NULL, fit_data->ws, fit_data->exp);
    }
    fprintf(stderr, "\n");
}

fit_params *fit_params_new() {
    fit_params *p = calloc(1, sizeof(fit_params));
    return p;
}
void fit_params_add_parameter(fit_params *p, double *value) {
    p->n++;
    p->func_params = realloc(p->func_params, sizeof(double *)*p->n);
    p->func_params_err = realloc(p->func_params_err, sizeof(double)*p->n);
    p->func_params[p->n-1] = value;
}
void fit_params_free(fit_params *p) {
    free(p->func_params);
    free(p->func_params_err);
    free(p);
}

fit_data *fit_data_new(const jibal *jibal, simulation *sim, gsl_histogram *exp, sample_model *sm,  reaction * const *reactions, const char *fit_vars, int fit_low, int fit_high, int print_iters) {
    struct fit_data *f = malloc(sizeof(struct fit_data));
    f->n_iters_max = FIT_ITERS_MAX;
    f->low_ch = fit_low;
    if(f->low_ch <= 0)
        f->low_ch = (int)(exp->n*0.1);
    f->high_ch = fit_high;
    if(f->high_ch <= 0 || f->high_ch >= exp->n)
        f->high_ch = exp->n - 1;
    f->jibal = jibal;
    f->sim = sim;
    f->exp = exp;
    f->sample = NULL;
    f->sm = sm;
    f->reactions = reactions;
    f->ws = NULL;
    f->fit_params = fit_params_new();
    f->print_iters = print_iters;
    add_fit_params(sim, sm, f->fit_params, fit_vars);
    if(f->fit_params->n == 0) {
        fit_params_free(f->fit_params);
        free(f);
        return NULL;
    }
    return f;
}

void fit_data_free(fit_data *f) {
    if(!f)
        return;
    fit_params_free(f->fit_params);
    free(f);
}

struct fit_stats fit(gsl_histogram *exp, struct fit_data *fit_data) {
    struct fit_stats stats = {.n_iters = 0, .n_evals = 0, .cputime_actual = 0.0};
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
    struct fit_params *fit_params = fit_data->fit_params;
    gsl_multifit_nlinear_fdf fdf;
    fdf.params = fit_data;
    fit_data->stats = &stats;
    if(!fit_data->exp) {
        fprintf(stderr, "No experimental data, can not fit.\n");
        return stats;
    }
    fdf.f = &fit_function;
    fdf.df = NULL; /* Jacobian, with NULL using finite difference. TODO: this could be implemented */
    fdf.fvv = NULL; /* No geodesic acceleration */
    fdf.n = (fit_data->high_ch-fit_data->low_ch+1);
    fdf.p = fit_params->n;
    if(fdf.n < fdf.p) /* insufficient data points */
        return stats;
    gsl_vector *f;
    gsl_matrix *J;
    double chisq, chisq0;
    fit_data->dof = fdf.n-fdf.p;
    int status, info;
    size_t i, j;
    gsl_matrix *covar = gsl_matrix_alloc (fit_params->n, fit_params->n);
    gsl_vector *x = gsl_vector_alloc(fit_params->n);
    for(i = 0; i < fit_params->n; i++) {
        gsl_vector_set(x, i, *fit_params->func_params[i]); /* Initial values of fitted parameters from function parameter array */
    }
    double *weights = malloc(sizeof(double)*(fit_data->high_ch-fit_data->low_ch+1));
    for(i = fit_data->low_ch; i <= fit_data->high_ch; i++) {
        if(exp->bin[i] > 1.0) {
            weights[i-fit_data->low_ch] = 1.0 / (exp->bin[i]);
        } else {
            weights[i-fit_data->low_ch] = 1.0; /* TODO: ?*/
        }
    }
    gsl_vector_view wts = gsl_vector_view_array(weights, fit_data->high_ch-fit_data->low_ch+1);

    const double xtol = FIT_XTOL;
    const double gtol = FIT_GTOL;
    const double ftol = FIT_FTOL;

/* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, fdf.n, fdf.p);

/* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (x, &wts.vector, &fdf, w);

/* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

/* solve the system with a maximum of n_iters_max iterations */
    status = gsl_multifit_nlinear_driver(fit_data->n_iters_max, xtol, gtol, ftol,
                                         fit_callback, fit_data, &info, w);

/* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

/* compute final cost */
    gsl_blas_ddot(f, f, &chisq);
    fprintf(stderr, "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    stats.n_iters = gsl_multifit_nlinear_niter(w);
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    stats.n_evals = fdf.nevalf;
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));
    fprintf (stderr, "status = %s\n", gsl_strerror (status));

    double c = GSL_MAX_DBL(1, sqrt(chisq / fit_data->dof));
    fprintf(stderr, "Final fitted parameters\n");
    for(i = 0; i < fit_params->n; i++) {
        fprintf(stderr, "    prob[%zu] = %g +- %g (%.2lf%%)\n", i,
                gsl_vector_get(w->x, i),
                c * sqrt(gsl_matrix_get(covar, i, i)),
                100.0*(c * sqrt(gsl_matrix_get(covar, i, i)))/gsl_vector_get(w->x, i)
        );
    }
    fprintf(stderr, "Correlation coefficients table:\n     ");
    for(j = 0; j < fit_params->n; j++) {
        fprintf(stderr, "  prob[%zu]   ", j);
    }
    fprintf(stderr, "\n");
    for(i = 0; i < fit_params->n; i++) {
        fprintf(stderr, "prob[%zu] ", i);
        for (j = 0; j <= i; j++) {
            fprintf(stderr, " %8.5f", gsl_matrix_get(covar, i, j)/sqrt(gsl_matrix_get(covar, i, i)*gsl_matrix_get(covar, j, j)));
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "chisq/dof = %g\n", chisq / fit_data->dof);
    for(i = 0; i < fit_params->n; i++) { /* Clear all err values */
        fit_params->func_params_err[i] = 0.0;
    }
    for(i = 0; i < fit_params->n; i++) { /* Update final fitted values to the table (same as used for initial guess) */
        *(fit_params->func_params[i]) = gsl_vector_get(w->x, i);
        fit_params->func_params_err[i] = c * sqrt(gsl_matrix_get(covar, i, i));
    }
    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);
    gsl_vector_free(x);
    free(weights);
    return stats;
}
