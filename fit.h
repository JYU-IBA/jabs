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
};

typedef struct fit_data {
    gsl_histogram *exp; /* experimental data to be fitted */
    const simulation *sim;
    reaction * const *reactions;
    const jibal *jibal;
    sample *sample;
    sample_model *sm;
    fit_params *fit_params; /* Allocated with fit_data_new() and freed by fit_data_free() */
    sim_workspace *ws; /* Allocated and leaked by fitting function! */
    size_t low_ch;
    size_t high_ch;
    size_t n_iters_max;
    double dof;
    int print_iters;
    struct fit_stats *stats;
} fit_data;

fit_data *fit_data_new(const jibal *jibal, simulation *sim, gsl_histogram *exp, sample_model *sm,  reaction * const *reactions, const char *fit_vars, int fit_low, int fit_high, int print_iters);
void fit_data_free(struct fit_data *f);

struct fit_stats fit(gsl_histogram *exp, struct fit_data *fit_data);
int fit_function(const gsl_vector *x, void *params, gsl_vector *f);
void fit_callback(size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
fit_params *fit_params_new();
void fit_params_add_parameter(fit_params *p, double *value);
void fit_params_free(fit_params *p);

#endif // JABS_FIT_H
