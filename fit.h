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
    double err; /* Error estimate will be stored here. */
    double err_rel;
    char *name;
} fit_variable;

typedef struct fit_params {
    size_t n; /* Number of function parameters */
    fit_variable *vars; /* Note that this is NOT an array of pointers. */
} fit_params;

#include "simulation.h"
#include "reaction.h"
#include "sample.h"

#define FIT_ERROR_NONE 0
#define FIT_ERROR_MAXITER 1
#define FIT_ERROR_NO_PROGRESS 2
#define FIT_ERROR_SANITY 3
#define FIT_ERROR_IMPOSSIBLE 4

struct fit_stats {
    size_t n_evals;
    size_t n_iters;
    size_t n_evals_iter; /* Number of function evaluations per iteration */
    double cputime_cumul;
    double cputime_iter;
    double chisq0;
    double chisq;
    double chisq_dof;
    size_t iter;
    double rel; /* This is updated as we iterate */
    int error;
    int info; /* From GSL */
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
    sim_workspace **ws; /* Allocated and leaked by fitting function! An array of sim->n_det */
    struct roi *fit_ranges; /* Array of fit_range, size n_fit_ranges, freed by fit_data_free() */
    size_t n_fit_ranges;
    size_t n_iters_max;
    double xtol;
    double gtol;
    double ftol;
    double dof;
    //int lm_accel;
    struct fit_stats stats;
} fit_data;


fit_data *fit_data_new(const jibal *jibal, simulation *sim);
void fit_data_free(struct fit_data *fit); /* Doesn't free everything in fit_data. Does free fit_params and fit_ranges */
void fit_data_print(FILE *f, const struct fit_data *fit_data);
void fit_data_roi_print(FILE *f, const struct fit_data *fit_data, const struct roi *roi);
gsl_histogram *fit_data_exp(const struct fit_data *fit_data, size_t i_det);
gsl_histogram *fit_data_sim(const struct fit_data *fit_data, size_t i_det);
void fit_data_exp_free(struct fit_data *fit_data);
int fit_data_load_exp(struct fit_data *fit, size_t i_det, const char *filename);
int fit_data_add_det(struct fit_data *fit_data, detector *det);
sim_workspace *fit_data_ws(const struct fit_data *fit_data, size_t i_det);
size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data);
int fit_data_workspaces_init(struct fit_data *fit_data); /* If any workspace initialization fails, this frees allocated memory (function below) and returns non-zero */
void fit_data_workspaces_free(struct fit_data *fit_data); /* Also sets workspace pointers to NULL */
struct fit_stats fit_stats_init();
int fit(struct fit_data *fit_data);
int fit_function(const gsl_vector *x, void *params, gsl_vector *f);
void fit_callback(size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
fit_params *fit_params_new();
int fit_params_add_parameter(fit_params *p, double *var, const char *name);
void fit_params_free(fit_params *p);
void fit_stats_print(FILE *f, const struct fit_stats *stats);
int fit_data_fit_range_add(struct fit_data *fit_data, const struct roi *range); /* Makes a deep copy */
void fit_data_fit_ranges_free(struct fit_data *fit_data);
int fit_set_roi_from_string(roi *r, const char *str); /* Parses only low and high from "[low:high]". */
const char *fit_error_str(int error);
#endif // JABS_FIT_H
