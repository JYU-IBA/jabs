/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_PARAMS_H
#define JABS_PARAMS_H
#include "fit_variable.h"

typedef struct fit_params {
    size_t n; /* Number of function parameters */
    size_t n_active; /* Recalculated by fit_params_update(), same as fdf->p in multifit_nlinear */
    fit_variable *vars; /* Note that this is NOT an array of pointers. It has n elements. */
} fit_params;

fit_params *fit_params_new();
void fit_params_free(fit_params *p);
int fit_params_add_parameter(fit_params *p, fit_variable_type type, double *value, const char *name, const char *unit, double unit_factor, size_t i_det);  /* Pointer to parameter to be fitted (value) is accessed during fitting (read, write). No guarantees that it stays accessible after the fit is over and user decides to change something! */
int fit_params_update(fit_params *p);
void fit_params_print(const fit_params *params, int active, const char *pattern, jabs_msg_level msg_level); /* if active is TRUE print only active variables. pattern can be NULL to bypass matching. */
void fit_params_print_final(const fit_params *params);
size_t fit_params_enable(fit_params *params, const char *s, int enable); /* Enable/disable one or more variables matching pattern s. Returns number of matches. */
void fit_parameters_update(const fit_params *fit_params, const gsl_multifit_nlinear_workspace *w, const gsl_matrix *covar, double chisq_dof); /* Updates values in fit_params, computes errors */
void fit_parameters_update_changed(const fit_params *fit_params); /* Checks if values have changed since fit_parameters_update(), computes new error */
int fit_params_enable_using_string(fit_params *params, const char *fit_vars);
fit_variable *fit_params_find_active(const fit_params *params, size_t i_v); /* Find active fit parameter with matching index number i_v as set by fit_params_update() */
const char *fit_variable_type_str(const fit_variable *var);
#endif // JABS_PARAMS_H
