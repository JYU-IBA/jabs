/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

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
#include "detector.h"
#include "sample.h"

typedef struct fit_variable {
    double *value; /* Pointer to a value. This is not allocated or free'd by fitting related methods. */
    double value_orig;
    double value_final;
    double value_iter;
    double err; /* Error estimate of fit will be stored here. */
    double err_rel; /* Relative error */
    double sigmas; /* Change, relative to error */
    char *name;
    const char *unit;
    double unit_factor;
    int active; /* Set to FALSE by default, if this variable is to be used it should be set to TRUE */
    size_t i_v; /* Index in fit */
} fit_variable;

typedef struct fit_params {
    size_t n; /* Number of function parameters */
    size_t n_active; /* Recalculated by fit_params_update() */
    fit_variable *vars; /* Note that this is NOT an array of pointers. It has n elements. */
} fit_params;

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
#define FIT_ERROR_IMPOSSIBLE (-5)

#define FIT_PHASE_FAST 1
#define FIT_PHASE_SLOW 2

struct fit_stats {
    size_t n_evals;
    size_t n_iters;
    size_t n_evals_iter; /* Number of function evaluations per iteration */
    double cputime_cumul;
    double cputime_iter;
    double chisq0;
    double chisq;
    double chisq_dof;
    double rcond;
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
    gsl_histogram **exp; /* experimental data to be fitted, array of sim->n_det */
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
    double gtol; /* Not used */
    double ftol; /* Not used */
    double chisq_tol; /* Chi squared relative change tolerance */
    double chisq_fast_tol; /* Chi squared relative change tolerance (fast phase) */
    double dof; /* Degrees of freedom (calculated) */
    struct fit_stats stats; /* Fit statistics, updated as we iterate */
    int phase_start; /* Fit phase to start from (see FIT_PHASE -defines) */
    int phase_stop; /* Inclusive */
    gsl_histogram **histo_sum_iter; /* Array of histograms, updates every iter. */
    size_t n_histo_sum;
} fit_data;


fit_data *fit_data_new(const jibal *jibal, simulation *sim);
void fit_data_free(struct fit_data *fit); /* Doesn't free everything in fit_data. Does free fit_params and fit_ranges */
void fit_data_print(FILE *f, const struct fit_data *fit_data);
void fit_data_roi_print(FILE *f, const struct fit_data *fit_data, const struct roi *roi);
gsl_histogram *fit_data_exp(const struct fit_data *fit_data, size_t i_det);
gsl_histogram *fit_data_sim(const struct fit_data *fit_data, size_t i_det);
void fit_data_exp_free(struct fit_data *fit_data);
int fit_data_load_exp(struct fit_data *fit, size_t i_det, const char *filename);
void fit_data_histo_sum_free(struct fit_data *fit_data);
int fit_data_add_det(struct fit_data *fit_data, detector *det);
sim_workspace *fit_data_ws(const struct fit_data *fit_data, size_t i_det);
size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data);
int fit_data_workspaces_init(struct fit_data *fit_data); /* If any workspace initialization fails, this frees allocated memory (function below) and returns non-zero */
void fit_data_workspaces_free(struct fit_data *fit_data); /* Also sets workspace pointers to NULL */
struct fit_stats fit_stats_init();
int fit(struct fit_data *fit_data);
void fit_covar_print(const gsl_matrix *covar);
void fit_parameters_update(const fit_data *fit, const gsl_multifit_nlinear_workspace *w, const gsl_matrix *covar); /* Updates values in fit_params, computes errors */
void fit_parameters_update_changed(const fit_data *fit); /* Checks if values have changed since fit_parameters_update(), computes new error */
int fit_function(const gsl_vector *x, void *params, gsl_vector *f);
void fit_callback(size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
fit_params *fit_params_new();
int fit_params_add_parameter(fit_params *p, double *value, const char *name, const char *unit, double unit_factor); /* Pointer to parameter to be fitted (value) is accessed during fitting (read, write). No guarantees that it stays accessible after the fit is over and user decides to change something! */
void fit_params_free(fit_params *p);
void fit_params_update(fit_params *p);
void fit_stats_print(FILE *f, const struct fit_stats *stats);
int fit_data_fit_range_add(struct fit_data *fit_data, const struct roi *range); /* Makes a deep copy */
void fit_data_fit_ranges_free(struct fit_data *fit_data);
int fit_set_roi_from_string(roi *r, const char *str); /* Parses only low and high from "[low:high]". */
double fit_emin(struct fit_data *fit, size_t i_det); /* Returns lowest energy of fit ranges for detector i_det. Detectors, calibrations and fit ranges must be set before calling. */
const char *fit_error_str(int error);
#endif // JABS_FIT_H
