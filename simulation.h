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
#include "detector.h"

typedef struct {
    detector det;
    double stop_step_incident;
    double stop_step_exiting;
    double p_sr;
    double sample_theta; /* Polar angle. Kind of. Zero is sample perpendicular to beam. */
    double sample_phi; /* Typically one uses a zero here, unless doing channeling stuff. Note that this is an azimuthal angle. */
    double alpha, beta, theta; /*  These are for convenience and follow the SimNRA conventions! Don't use them for any physics if you can avoid them. Well, theta *is* very convenient... */
    //ion ion;
    const jibal_isotope *beam_isotope;
    double beam_E;
    double beam_E_broad; /* Variance */
    int fast;
    double emin;
    size_t depthsteps_max;
    int ds;
    int ds_steps_azi;
    int ds_steps_polar;
    int n_ds;
    double channeling; /* a very ad-hoc channeling yield correction */
} simulation;

typedef struct sim_reaction {
    const reaction *r;
    ion p; /* Reaction product */
    gsl_histogram *histo;
    brick *bricks;
    size_t n_bricks;
    int stop;
    double max_depth;
    size_t i_isotope; /* Number of isotope (r->target) in sample->isotopes */
    double theta;
    double K;
    double (*cross_section)(const struct sim_reaction *, double);
    double theta_cm;
    double mass_ratio; /* m1/m2, these are variables to speed up cross section calculation */
    double E_cm_ratio;  /* m1/(m1+m2) */
    double cs_constant; /* Rutherford cross section term for RBS (non-energy dependent) */
    double sigma_to_lab_factor; /* Factor to convert cross section to lab */
    double r_VE_factor; /* Andersen correction factor r_VE = this / E_cm */
} sim_reaction; /* Workspace for a single reaction. Yes, the naming is confusing. */


typedef struct {
    simulation sim;
    size_t n_reactions; /* Number of reactions is counted on init. */
    double *c; /* Concentrations for n_isotopes at some arbitrary x */
    double c_x; /* at this x */
    jibal_gsto *gsto;
    const jibal_config *jibal_config;
    int rk4;
    int nucl_stop_accurate;
    int mean_conc_and_energy;
    gsto_stopping_type stopping_type;
    size_t n_channels; /* histogram */
    ion ion;
    sim_reaction *reactions;
    const jibal_isotope *isotopes;
} sim_workspace;


#include "sample.h"
simulation *sim_init();
void sim_free(simulation *sim);
void sim_calculate_geometry(simulation *sim);
int sim_sanity_check(const simulation *sim);
sim_workspace *sim_workspace_init(const simulation *sim, const reaction *reactions, const sample *sample, const jibal *jibal);
void sim_workspace_free(sim_workspace *ws);
void sim_workspace_recalculate_n_channels(sim_workspace *ws, const simulation *sim);
void simulation_print(FILE *f, const simulation *sim);
void convolute_bricks(sim_workspace *ws);
void sim_reaction_recalculate_internal_variables(sim_reaction *r);
double sim_reaction_cross_section_rbs(const sim_reaction *r, double E);
double sim_reaction_cross_section_erd(const sim_reaction *r, double E);
double sim_reaction_cross_section_none(const sim_reaction *r, double E);
#endif // JABS_SIMULATION_H
