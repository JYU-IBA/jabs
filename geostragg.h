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

double scattering_angle_exit_deriv(const ion *incident, const sim_workspace *ws);
double exit_angle_delta(const sim_workspace *ws, char direction);
double geostragg(const sim_workspace *ws, const sample *sample, const sim_reaction *r, depth d, double E_0, double delta_beta, double theta_deriv);
double theta_deriv_beta(const detector *det, char direction);
double beta_deriv(const detector *det, const simulation *sim, char direction);
#endif // JABS_GEOSTRAGG_H
