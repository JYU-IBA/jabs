/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_STOP_H
#define JABS_STOP_H
#include "jibal_gsto.h"
#include "ion.h"
#include "sample.h"
#include "ms.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct stop_range_table { /* Ion in material of constant composition (a range of a sample): path length and straggling integrals as functions of energy, see stop_range_table_make() */
    const jibal_isotope *isotope; /* Ion, stopping data, stopping and straggling types, emin and material are the key of the table */
    const jabs_ion_gsto *ion_gsto;
    gsto_stopping_type stop_type;
    gsto_stopping_type stragg_type;
    double emin;
    size_t n_isotopes;
    const jibal_isotope **isotopes;
    double *conc;
    double bragg;
    double stragg;
    int valid; /* FALSE if table could not be made (e.g. no stopping), then ion is stepped */
    size_t n; /* Number of nodes, equally spaced in ln(E) starting from emin */
    double emax;
    double ln_emin;
    double inv_dlne;
    double *E;
    double *R; /* Path length from E to emin, i.e. integral of 1/S from emin to E */
    double *S; /* Stopping at E */
    double *Q; /* Integral of W/S^3 from emin to E, where W is straggling (energy loss variance per unit path length) */
    double *dQ; /* W/S^3 at E */
} stop_range_table;

typedef struct stop_range_tables { /* Cache of range tables, tables are made when needed */
    stop_range_table **t;
    size_t n;
} stop_range_tables;

typedef struct jabs_stop {
    jibal_gsto *gsto;
    gsto_stopping_type type;
    int rk4;
    double accuracy; /* Exiting ions: absolute energy error tolerance per step of adaptive (Dormand-Prince 5(4)) steps, zero for fixed (RK4) steps. Requires rk4. */
    double stragg_accuracy; /* Exiting ions: relative straggling (variance) error tolerance per step of adaptive steps */
    double emin;
    const ms_model *ms; /* Small angle multiple scattering model, NULL if multiple scattering is not calculated */
    stop_range_tables *range_tables; /* Exiting ions cross ranges of constant composition using these tables (when multiple scattering is not calculated) instead of steps, see stop_sample_exit_remaining(). NULL if not used. */
} jabs_stop;

typedef struct jabs_stop_step_params {
    double step; /* This is energy. Zero is automatic. */
    double min;
    double max;
    double sigmas;
    double max_rel; /* Maximum step relative to energy (zero: no limit). Not applied to fixed steps. */
} jabs_stop_step_params;

typedef struct stop_sw {
    double S; /* Stopping */
    double W; /* Straggling (energy loss variance per unit path length) */
    double S_x; /* Derivative of stopping with respect to (perpendicular) depth at constant energy, i.e. due to concentration gradients */
} stop_sw;

typedef struct stop_adaptive { /* State of adaptive stepping, carried from one step to the next */
    double h; /* Suggested next path length step, zero if not known */
    stop_sw sw_end; /* Stopping and straggling at the end of the previous step, reused as first stage of the next step if valid */
    depth d_end;
    double E_end;
    int valid;
    size_t n_rejected; /* Statistics: steps retried because energy would have gone below zero */
    size_t n_evals; /* Statistics: evaluations of stopping and straggling */
} stop_adaptive;

depth stop_next_crossing(const ion *incident, const sample *sample, const depth *d_from);
void stop_ms_step(const ms_model *ms, ion *ion, const sample *sample, depth d, double E, double h); /* Updates multiple scattering (angular spread) of ion after path length h (at/m2) */
depth stop_step(const jabs_stop *stop, const jabs_stop *stragg, ion *incident, const sample *sample, depth depth_before, double step);
depth stop_step_adaptive(const jabs_stop *stop, const jabs_stop *stragg, ion *incident, const sample *sample, depth depth_before, double step, stop_adaptive *state);
double stop_sample(const jabs_stop *stop, const ion *incident, const sample *sample, depth depth, double E);
double stop_step_calc(const jabs_stop_step_params *params, const ion *ion);
double stop_step_max(const jabs_stop_step_params *params, const ion *ion);
stop_range_tables *stop_range_tables_new(void);
void stop_range_tables_free(stop_range_tables *rt);
int stop_sample_exit(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params_exiting, ion *p, depth depth_start, const sample *sample);
int stop_sample_exit_remaining(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params_exiting, ion *p, depth depth_start, const sample *sample, double *remaining); /* Like stop_sample_exit(), if ion stops, path length it still had to go to the surface is stored in remaining (if not NULL) */
#ifdef __cplusplus
}
#endif
#endif // JABS_STOP_H
