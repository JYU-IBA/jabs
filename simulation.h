/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SIMULATION_H
#define JABS_SIMULATION_H

#include <gsl/gsl_histogram.h>
#include <jibal.h>
#include <time.h>

#include "ion.h"
#include "reaction.h"
#include "brick.h"

typedef struct {
    //int n_channels;
    int n_reactions;
    double energy_slope;
    double energy_offset;
    double energy_resolution; /* Stored as variance, i.e. energy squared (in SI-units J^2) */
    double stop_step_incident;
    double stop_step_exiting;
    double p_sr;
    double alpha;
    double beta;
    double theta;
    //ion ion;
    const jibal_isotope *beam_isotope;
    double beam_E;
    double beam_E_broad; /* Variance */
    int fast;
    double emin;
    size_t depthsteps_max;
} simulation;

typedef struct {
    const reaction *r;
    ion p; /* Reaction product */
    gsl_histogram *histo;
    brick *bricks;
    size_t n_bricks;
    int stop;
} sim_reaction; /* Workspace for a single reaction. Yes, the naming is confusing. */


typedef struct {
    simulation sim;
    int n_reactions; /* Same as sim->n_reactions, but we want to keep it here too, as it is used for allocations */
    double *c; /* Concentrations for n_isotopes at some arbitrary x */
    double c_x; /* at this x */
    int i_range_accel;
    jibal_gsto *gsto;
    const jibal_config *jibal_config;
    int rk4;
    gsto_stopping_type stopping_type;
    int n_channels; /* histogram */
    ion ion;
    sim_reaction *reactions;
    const jibal_isotope *isotopes;
} sim_workspace;


#include "sample.h"
simulation *sim_init();
void sim_free(simulation *sim);
int sim_sanity_check(const simulation *sim);
sim_workspace *sim_workspace_init(const simulation *sim, const reaction *reactions, const sample *sample, const jibal *jibal);
void sim_workspace_free(sim_workspace *ws);
void sim_workspace_recalculate_calibration(sim_workspace *ws, const simulation *sim);
void simulation_print(FILE *f, const simulation *sim);
void convolute_bricks(sim_workspace *ws);
#endif // JABS_SIMULATION_H
