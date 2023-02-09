/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_FIT_VARIABLE_H
#define JABS_FIT_VARIABLE_H
#include <stdlib.h>

typedef struct fit_variable {
    double *value; /* Pointer to a value. This is not allocated or free'd by fitting related methods. */
    double value_orig;
    double value_final;
    double value_iter; /* On first call of fit_function(), store *value to this */
    double err; /* Error estimate of fit will be stored here. */
    double err_rel; /* Relative error */
    double sigmas; /* Change, relative to error */
    char *name;
    const char *unit;
    double unit_factor;
    int active; /* Set to FALSE by default, if this variable is to be used it should be set to TRUE */
    int active_iter_call; /* The *value is different from value_iter, i.e. the value is being changed by the fit algorithm */
    size_t i_v; /* Index in fit */
} fit_variable;

#endif // JABS_FIT_VARIABLE_H
