/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_GEOSTRAGG_H
#define JABS_GEOSTRAGG_H
#include "detector.h"
#include "ion.h"
#include "stop.h"
#include "sim_reaction.h"

typedef struct {
    char direction;
    double theta;
    double phi; /* phi angle of "direction" */
    double theta_plus;
    double theta_minus;
    double delta_beta;
    double theta_deriv;
    double beta_deriv;
    double theta_product_plus;
    double phi_product_plus;
    double theta_product_minus;
    double phi_product_minus;
} geostragg_vars_dir;

typedef struct {
    double scatter_theta;
    double theta_product;
    double phi_product;
    geostragg_vars_dir x;
    geostragg_vars_dir y;
} geostragg_vars;


double scattering_angle(const ion *incident, double sample_theta, double sample_phi, const detector *det); /* Calculate scattering angle necessary for ion (in sample coordinate system) to hit detector *//* Calculate scattering angle necessary for ion (in sample coordinate system) to hit detector */
double scattering_angle_exit_deriv(const ion *incident, double sample_theta, double sample_phi, const detector *det); /* Calculates the dtheta/dbeta derivative for geometrical straggling */
double exit_angle_delta(double sample_theta, double sample_phi, const detector *det, const aperture *beam_aperture, char direction);
geostragg_vars geostragg_vars_calculate(const ion *incident, double sample_theta, double sample_phi, const detector *det, const aperture *beam_aperture, int geostragg_enabled, int beta_manual_enabled);
double geostragg(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params_exiting, const sample *sample, const sim_reaction *r, const geostragg_vars_dir *gd, depth d, double E_0);
double beta_deriv(double sample_theta, double sample_phi, const detector *det, char direction);
double exit_angle(double sample_theta, double sample_phi, double det_theta, double det_phi);
#endif // JABS_GEOSTRAGG_H
