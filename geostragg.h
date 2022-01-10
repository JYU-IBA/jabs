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
#include "simulation.h"

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


double scattering_angle(const ion *incident, const sim_workspace *ws); /* Calculate scattering angle necessary for ion (in sample coordinate system) to hit detector */
double scattering_angle_exit_deriv(const ion *incident, const sim_workspace *ws);
double exit_angle_delta(const sim_workspace *ws, char direction);
geostragg_vars geostragg_vars_calculate(const sim_workspace *ws, const ion *incident);
double geostragg(const sim_workspace *ws, const sample *sample, const sim_reaction *r, const geostragg_vars_dir *gd, depth d, double E_0);
double theta_deriv_beta(const detector *det, char direction);
double beta_deriv(const detector *det, const simulation *sim, char direction);
int geostragg_vars_print(FILE *f, const geostragg_vars *g);
#endif // JABS_GEOSTRAGG_H
