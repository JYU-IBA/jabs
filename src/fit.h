/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_FIT_H
#define JABS_FIT_H

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>
#include <jibal.h>

#include "simulation.h"
#include "simulation_workspace.h"
#include "detector.h"
#include "sample.h"
#include "fit_params.h"
#include "simulation.h"
#include "reaction.h"
#include "sample.h"

#define FIT_SUCCESS_CHISQ (2)
#define FIT_SUCCESS_DELTA (1)
#define FIT_SUCCESS (0)
#define FIT_ERROR_NONE (0)
#define FIT_ERROR_GENERIC (-1)
#define FIT_ERROR_MAXITER (-2)
#define FIT_ERROR_NO_PROGRESS (-3)
#define FIT_ERROR_SANITY (-4)
#define FIT_ERROR_WORKSPACE_INITIALIZATION (-5)
#define FIT_ERROR_IMPOSSIBLE (-6)
#define FIT_ERROR_ABORTED (-7)

#define FIT_PHASE_FAST 1
#define FIT_PHASE_SLOW 2

struct fit_stats {
    int phase;
    size_t n_evals;
    size_t n_evals_iter; /* Number of function evaluations per iteration */
    size_t n_speedup_evals;
    size_t n_speedup_evals_iter;
    double cputime_cumul;
    double cputime_iter;
    double chisq0;
    double chisq;
    double chisq_dof;
    double rcond;
    double norm;
    size_t iter;
    size_t iter_call;
    int error;
};

typedef struct roi {
    size_t i_det; /* detector sim->det[i_det] */
    size_t low;
    size_t high;
} roi;

typedef struct fit_data {
    gsl_histogram **exp; /* experimental data to be fitted, array of n_exp elements */
    size_t n_exp; /* same as sim->n_det, but this keeps track on how many spectra we have allocated in exp and ref */
    gsl_histogram *ref; /* reference spectra, exactly one */
    simulation *sim;
    const jibal *jibal; /* This shouldn't be here, but it is the only place I can think of */
    sample_model *sm;
    fit_params *fit_params; /* Allocated with fit_data_new() and freed by fit_data_free() */
    sim_workspace **ws; /* Allocated and leaked by fitting function! An array of n_ws. */
    size_t n_ws;
    struct roi *fit_ranges; /* Array of fit_range, size n_fit_ranges, freed by fit_data_free() */
    size_t n_fit_ranges;
    size_t n_iters_max; /* Maximum number of iterations, in each fit phase */
    double xtol; /* Tolerance of step size */
    double chisq_tol; /* Chi squared relative change tolerance */
    double chisq_fast_tol; /* Chi squared relative change tolerance (fast phase) */
    double dof; /* Degrees of freedom (calculated) */
    struct fit_stats stats; /* Fit statistics, updated as we iterate */
    int phase_start; /* Fit phase to start from (see FIT_PHASE -defines) */
    int phase_stop; /* Inclusive */
    gsl_histogram **histo_sum_iter; /* Array of histograms, updates every iter. */
    size_t n_histo_sum;
    int (*fit_iter_callback)(struct fit_stats stats);
    int magic_bricks;
} fit_data;

fit_data *fit_data_new(const jibal *jibal, simulation *sim);
void fit_data_defaults(fit_data *f);
void fit_data_free(struct fit_data *fit); /* Doesn't free everything in fit_data. Does free fit_params and fit_ranges */
void fit_data_print(FILE *f, const struct fit_data *fit_data);
void fit_data_roi_print(FILE *f, const struct fit_data *fit_data, const struct roi *roi);
gsl_histogram *fit_data_exp(const struct fit_data *fit_data, size_t i_det);
gsl_histogram *fit_data_sim(const struct fit_data *fit_data, size_t i_det); /* You should probably use fit_data_histo_sum() to get the simulated sum spectra */
gsl_histogram *fit_data_ref(const struct fit_data *fit_data);
fit_params *fit_params_all(fit_data *fit);
void fit_data_exp_alloc(fit_data *fit);
void fit_data_exp_free(fit_data *fit);
int fit_data_load_exp(struct fit_data *fit, size_t i_det, const char *filename);
void fit_data_histo_sum_free(struct fit_data *fit_data);
void fit_data_histo_sum_store(struct fit_data *fit_data);
gsl_histogram *fit_data_histo_sum(const struct fit_data *fit_data, size_t i_det);
int fit_data_add_det(struct fit_data *fit, detector *det);
sim_workspace *fit_data_ws(const struct fit_data *fit_data, size_t i_det);
size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data);
sim_workspace *fit_data_workspace_init(fit_data *fit, size_t i_ws);
int fit_data_workspaces_init(fit_data *fit);
void fit_data_workspaces_free(struct fit_data *fit_data); /* Also sets workspace pointers to NULL */
struct fit_stats fit_stats_init();
int fit(struct fit_data *fit_data);
void fit_covar_print(const gsl_matrix *covar);

int fit_parameters_set_from_vector(struct fit_data *fit, const gsl_vector *x); /* Updates values in fit params as they are varied by the fit algorithm. */
int fit_function(const gsl_vector *x, void *params, gsl_vector *f);
int fit_sanity_check(const fit_data *fit);
int fit_speedup(fit_data *fit);
int fit_speedup_fluence(struct fit_data *fit, const fit_variable *var);
int fit_set_residuals(const struct fit_data *fit_data, gsl_vector *f);
void fit_iter_stats_update(struct fit_data *params, const gsl_multifit_nlinear_workspace *w);
void fit_iter_stats_print(const struct fit_stats *stats);
void fit_stats_print(FILE *f, const struct fit_stats *stats);
int fit_data_fit_range_add(struct fit_data *fit_data, const struct roi *range); /* Makes a deep copy */
void fit_data_fit_ranges_free(struct fit_data *fit_data);
int fit_set_roi_from_string(roi *r, const char *str); /* Parses only low and high from "[low:high]". */
double fit_emin(struct fit_data *fit, size_t i_det); /* Returns lowest energy of fit ranges for detector i_det. Detectors, calibrations and fit ranges must be set before calling. */
const char *fit_error_str(int error);
#endif // JABS_FIT_H
