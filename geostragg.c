/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <assert.h>
#include "jabs.h"
#include "message.h"
#include "geostragg.h"


double scattering_angle(const ion *incident, const sim_workspace *ws) { /* Calculate scattering angle necessary for ion (in sample coordinate system) to hit detector */
    double theta, phi;
    double scatter_theta, scatter_phi;
    rotate(-1.0*ws->sim->sample_theta, ws->sim->sample_phi, incident->theta, incident->phi, &theta, &phi); /* Move from sample coordinates to lab */
    rotate(ws->det->theta, ws->det->phi, theta, phi, &scatter_theta, &scatter_phi); /* Move from lab to detector */
#ifdef DEBUG
    fprintf(stderr, "theta = %g deg, phi = %g deg.\n", theta/C_DEG, phi/C_DEG);
    fprintf(stderr, "scatter_theta = %g deg, scatter_phi = %g deg.\n", scatter_theta/C_DEG, scatter_phi/C_DEG);
    if(!ws->params.ds && !ws->params.geostragg) {
        assert(fabs(ws->det->theta - scatter_theta) < 0.01 * C_DEG); /* with DS this assert will fail */
    }
#endif
return scatter_theta;
}

double scattering_angle_exit_deriv(const ion *incident, const sim_workspace *ws) { /* Calculates the dtheta/dbeta derivative for geometrical straggling */
    double theta, phi;
    double scatter_theta, scatter_phi;
    double scatter_theta2, scatter_phi2;
    double theta_product, phi_product;
    double theta_product2, phi_product2;
    double epsilon = 0.001*C_DEG;
    rotate( -1.0*ws->sim->sample_theta, ws->sim->sample_phi, incident->theta, incident->phi, &theta, &phi);
    rotate( ws->det->theta, ws->det->phi, theta, phi, &scatter_theta, &scatter_phi);
    rotate(ws->det->theta+epsilon, ws->det->phi, theta, phi, &scatter_theta2, &scatter_phi2);
    rotate(ws->det->theta, ws->det->phi, ws->sim->sample_theta, ws->sim->sample_phi,  &theta_product, &phi_product);
    rotate(ws->det->theta+epsilon, ws->det->phi, ws->sim->sample_theta, ws->sim->sample_phi,  &theta_product2, &phi_product2);
    double deriv = (scatter_theta2-scatter_theta)/(theta_product2-theta_product) * -1.0; /* Since exit angle is pi - theta_product, the derivative dtheta/dbeta is -1.0 times this derivative calculated here */
#ifdef DEBUG
fprintf(stderr, "Scat: %g deg, eps: %g deg, product %g deg, eps: %g deg. Deriv %.8lf (-1.0 for IBM)\n", scatter_theta/C_DEG,scatter_theta2/C_DEG, theta_product/C_DEG, theta_product2/C_DEG, deriv);
#endif
return deriv;
}

double exit_angle_delta(const sim_workspace *ws, const char direction) {
    if(!ws->params.geostragg)
        return 0.0;
    if(ws->det->distance < 0.001 * C_MM)
        return 0.0;
    double sample_tilt = angle_tilt(ws->sim->sample_theta, ws->sim->sample_phi, direction);
    double det_tilt = C_PI - detector_angle(ws->det, direction);
    double exit = C_PI - det_tilt - sample_tilt;
    double geo = cos(exit)/cos(sample_tilt);
    double w = aperture_width_shape_product(ws->sim->beam_aperture, direction);
    double delta_beam = w *  geo / ws->det->distance;
    double delta_detector = aperture_width_shape_product(ws->det->aperture, direction) / ws->det->distance;
    double result = sqrt(pow2(delta_beam) + pow2(delta_detector));
#ifdef DEBUG
    fprintf(stderr, "Direction %c: detector angles %g deg, sample angle %g deg, exit angle %g deg\n", direction, det_tilt/C_DEG, sample_tilt/C_DEG, exit/C_DEG);
    fprintf(stderr, "cos(alpha) = %g, cos(beta) = %g, cos(beta)/cos(alpha) = %g\n", cos(sample_tilt), cos(exit), geo);
    fprintf(stderr, "Spread of exit angle in direction '%c' due to beam %g deg, due to detector %g deg. Combined %g deg FWHM.\n", direction, delta_beam/C_DEG, delta_detector/C_DEG, result/C_DEG);
#endif
    return result / C_FWHM;
}

geostragg_vars geostragg_vars_calculate(const sim_workspace *ws, const ion *incident) {
    geostragg_vars g;
    g.scatter_theta = scattering_angle(incident, ws);
    rotate(ws->det->theta, ws->det->phi, ws->sim->sample_theta, ws->sim->sample_phi, &g.theta_product, &g.phi_product); /* Detector in sample coordinate system */
#ifdef DEBUG
    double beta;
    beta = (C_PI - g.theta_product);
    ion_print(stderr, incident);
    fprintf(stderr, "Reaction product angles (in sample) %g deg and %g deg. Exit angle (beta) %g deg.\n", g.theta_product/C_DEG, g.phi_product/C_DEG, beta/C_DEG);
#endif
    g.x.direction = 'x';
    g.y.direction = 'y';
    geostragg_vars_dir *gds[] = {&g.x, &g.y, NULL};
    for(geostragg_vars_dir **gdp = gds; *gdp; gdp++) {
        geostragg_vars_dir *gd = *gdp;
        gd->delta_beta = exit_angle_delta(ws, gd->direction);
#ifdef UNNECESSARY_NUMERICAL_THINGS
        gd->theta_deriv = theta_deriv_beta(ws->det, gd->direction);
        gd->phi = ws->det->phi;
#else
        if(gd->direction == 'x') {
            gd->theta_deriv = -cos(ws->det->phi); /* TODO: check IBM with phi 180 deg and */
            gd->phi = 0.0;
        } else if (gd->direction == 'y') {
            gd->theta_deriv = -sin(ws->det->phi); /* TODO: check Cornell with phi 270 deg */
            gd->phi = C_PI_2;
        }
#endif
        gd->delta_beta = exit_angle_delta(ws, gd->direction) / 2.0;
        gd->beta_deriv = fabs(beta_deriv(ws->det, ws->sim, gd->direction));

        gd->theta_plus = g.scatter_theta + gd->delta_beta * gd->theta_deriv;
        gd->theta_minus = g.scatter_theta - gd->delta_beta * gd->theta_deriv;

        rotate(-1.0 * gd->delta_beta * gd->beta_deriv, g.phi_product, g.theta_product, g.phi_product, &gd->theta_product_plus, &gd->phi_product_plus); /* -1.0 again because of difference between theta and pi - theta */
        rotate(1.0 * gd->delta_beta * gd->beta_deriv, g.phi_product, g.theta_product, g.phi_product, &gd->theta_product_minus, &gd->phi_product_minus);

#ifdef DEBUG
        fprintf(stderr, "Spread in exit angle ('%c') %g deg\n", gd->direction, gd->delta_beta / C_DEG);
        fprintf(stderr, "dBeta/dBeta_%c: %g\n", gd->direction, gd->beta_deriv); /* TODO: verify, check sign */
        fprintf(stderr, "dTheta/dBeta_%c = %g\n", gd->direction, gd->theta_deriv); /* TODO: this should also be valid when sample_phi is not zero? */
        fprintf(stderr, "Delta Beta_%c %g deg, Delta theta %g deg\n", gd->direction, -1.0 * gd->delta_beta * gd->beta_deriv, gd->delta_beta * gd->theta_deriv);
        fprintf(stderr, "Reaction product (+) angles (in sample) %g deg and %g deg. Exit angle (beta) %g deg.\n", gd->theta_product_plus/C_DEG, gd->phi_product_plus/C_DEG, (C_PI - gd->theta_product_plus)/C_DEG);
        fprintf(stderr, "Reaction product (+) scattering angle %g deg.\n", gd->theta_plus/C_DEG);
        fprintf(stderr, "Reaction product (-) angles (in sample) %g deg and %g deg. Exit angle (beta) %g deg.\n", gd->theta_product_minus/C_DEG, gd->phi_product_minus/C_DEG, (C_PI - gd->theta_product_minus)/C_DEG);
        fprintf(stderr, "Reaction product (-) scattering angle %g deg.\n", gd->theta_minus/C_DEG);
        fprintf(stderr, "\n");
#endif
    }
    return g;
}

double geostragg(const sim_workspace *ws, const sample *sample, const sim_reaction *r, const geostragg_vars_dir *gd, const depth d, const double E_0) {
    if(!ws->params.geostragg) {
        return 0.0;
    }
    ion ion;
    ion = r->p;
    ion_set_angle(&ion, gd->theta_product_plus, gd->phi_product_plus);
    ion.E = reaction_product_energy(r->r, gd->theta_plus, E_0);
    ion.S = 0.0; /* We don't need straggling for anything, might as well reset it */
    post_scatter_exit(&ion, d, ws, sample);
    double Eplus = ion.E;
    ion = r->p;
    ion_set_angle(&ion, gd->theta_product_minus, gd->phi_product_minus);
    ion.E = reaction_product_energy(r->r, gd->theta_minus, E_0);
    ion.S = 0.0; /* We don't need straggling for anything, might as well reset it */
    post_scatter_exit(&ion, d, ws, sample);
    double Eminus = ion.E;
    double result = pow2(Eplus - Eminus);
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "Direction (%c), Eplus %g keV, Eminus %g keV, Straggling %g keV FWHM\n", gd->direction, Eplus/C_KEV, Eminus/C_KEV, C_FWHM*sqrt(result)/C_KEV);
#endif
    return result;
}

double theta_deriv_beta(const detector *det, const char direction) { /* TODO: possibly wrong sign */
    double x = detector_angle(det, 'x');
    double y = detector_angle(det, 'y');
    static const double delta = 0.001*C_DEG;
    double theta_from_rot = C_PI - theta_tilt(x, y); /* We shouldn't need to calculate this, but we do. Numerical issues? */
    if(direction == 'x') {
        x += delta;
    } else if(direction == 'y') {
        y += delta;
    } else {
        return 0.0;
    }
    double theta_eps = C_PI - theta_tilt(x, y);
    double deriv = (theta_eps - theta_from_rot) / delta;
#ifdef DEBUG
    fprintf(stderr, "Theta %.5lf deg, from rotations %.5lf deg. Delta = %g deg, Theta_eps = %.5lf deg\n", det->theta/C_DEG, theta_from_rot/C_DEG, delta/C_DEG, theta_eps/C_DEG);
    //assert(fabs(theta_from_rot - det->theta) < 1e-3); /* TODO: there is some numerical issue here */ */
    fprintf(stderr, "Derivative %.5lf (dir '%c').\n", deriv, direction);
#endif
    return deriv;
}

double beta_deriv(const detector *det, const simulation *sim, const char direction) { /* TODO: replace by analytical formula */
    double theta, phi;
    double theta_product, phi_product;
    static const double delta = 0.001*C_DEG;
    rotate(det->theta, det->phi, sim->sample_theta, sim->sample_phi, &theta_product, &phi_product); /* Detector in sample coordinate system */
    theta = delta;
    if(direction == 'x') {
        phi = 0.0;
    } else if(direction == 'y') {
        phi = C_PI_2;
    } else {
        return 0.0;
    }
    rotate(theta, phi, det->theta, det->phi, &theta, &phi);
    rotate(theta, phi, sim->sample_theta, sim->sample_phi, &theta, &phi); /* "new" Detector in sample coordinate system */
    double result = -1.0*(theta - theta_product)/delta; /* -1.0, since beta = pi - theta */
#ifdef DEBUG
    fprintf(stderr, "Beta deriv ('%c') got theta = %g deg (orig %g deg), diff in beta %g deg, result %g\n", direction, theta/C_DEG, theta_product/C_DEG,  (theta_product - theta)/C_DEG, result);
#endif
    return result;
}

int geostragg_vars_print(FILE *f, const geostragg_vars *g) {
    if(!g)
        return EXIT_FAILURE;
    jabs_message(MSG_INFO, f, "theta_product = %g deg\n", g->theta_product);
    jabs_message(MSG_INFO, f, "phi_product = %g deg\n", g->theta_product);
    return EXIT_SUCCESS;
}
