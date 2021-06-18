/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

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

typedef struct fit_params {
    size_t n; /* Number of function parameters */
    double **func_params; /* Function parameters array. Size is n. The contents can be modified. */
    double *func_params_err; /* Error estimates of parameters after fit. */
} fit_params;

#include "simulation.h"
#include "reaction.h"
#include "sample.h"

struct fit_stats {
    size_t n_evals;
    size_t n_iters;
    double cputime_actual;
    double chisq_dof;
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
    double dof;
    int print_iters;
    struct fit_stats stats;
} fit_data;


fit_data *fit_data_new(const jibal *jibal, simulation *sim);
void fit_data_free(struct fit_data *f); /* Doesn't free everything in fit_data. Does free fit_params and fit_ranges */
void fit_data_print(FILE *f, const struct fit_data *fit_data);
void fit_data_roi_print(FILE *f, const struct fit_data *fit_data, const struct roi *roi);
gsl_histogram *fit_data_exp(const struct fit_data *fit_data, size_t i_det);
sim_workspace *fit_data_ws(const struct fit_data *fit_data, size_t i_det);
size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data);
int fit_data_workspaces_init(struct fit_data *fit_data); /* If any workspace initialization fails, this frees allocated memory (function below) and returns non-zero */
void fit_data_workspaces_free(struct fit_data *fit_data); /* Also sets workspace pointers to NULL */
void fit_data_exp_free(struct fit_data *fit_data);
int fit(struct fit_data *fit_data);
int fit_function(const gsl_vector *x, void *params, gsl_vector *f);
void fit_callback(size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
fit_params *fit_params_new();
void fit_params_add_parameter(fit_params *p, double *value);
void fit_params_free(fit_params *p);
void fit_stats_print(FILE *f, const struct fit_stats *stats);
void fit_range_add(struct fit_data *fit, const struct roi *range); /* Makes a deep copy */
#endif // JABS_FIT_H
