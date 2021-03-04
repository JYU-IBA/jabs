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



struct fit_data {
    gsl_histogram *exp; /* experimental data to be fitted */
    const simulation *sim;
    reaction *reactions;
    const jibal *jibal;
    sample *sample;
    jibal_layer * const *layers;
    size_t n_layers;
    fit_params *fit_params;
    sim_workspace *ws; /* Handled by fitting function! */
    size_t low_ch;
    size_t high_ch;
    size_t n_iters_max;
    double cputime_actual; /* Statistics! */
    double dof;
};

struct fit_stats {
    size_t n_evals;
    size_t n_iters;
};

struct fit_stats fit(gsl_histogram *exp, struct fit_data *fit_data);
int func_f(const gsl_vector *x, void *params, gsl_vector *f);
void fit_callback(size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
fit_params *fit_params_new();
void fit_params_add_parameter(fit_params *p, double *value);
void fit_params_free(fit_params *p);

#endif // JABS_FIT_H
