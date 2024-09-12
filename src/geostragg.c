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
#include "jabs_debug.h"
#include "jabs.h"
#include "message.h"
#include "geostragg.h"
#include "rotate.h"


double scattering_angle(const ion *incident, double sample_theta, double sample_phi, const detector *det) { /* Calculate scattering angle necessary for ion (in sample coordinate system) to hit detector */
    double theta, phi;
    double scatter_theta, scatter_phi;
    rotate(-1.0 * sample_theta, sample_phi, incident->theta, incident->phi, &theta, &phi); /* Move from sample coordinates to lab */
    rotate(det->theta, det->phi, theta, phi, &scatter_theta, &scatter_phi); /* Move from lab to detector */
    DEBUGMSG("Scattering angle for incident_theta = %g deg, incident_phi = %g deg and sample_theta = %g deg, sample_phi = %g deg.",
             incident->theta / C_DEG, incident->phi / C_DEG, sample_theta / C_DEG, sample_phi / C_DEG);
    DEBUGMSG("Sample to lab, theta = %g deg, phi = %g deg.", theta/C_DEG, phi/C_DEG);
    DEBUGMSG("Lab to detector (final): scatter_theta = %g deg, scatter_phi = %g deg.", scatter_theta/C_DEG, scatter_phi/C_DEG);
#ifdef JABS_DEBUG_SCATTERING_ASSERT
    if(!ws->params->ds && !ws->params->geostragg) {
        assert(fabs(ws->det->theta - scatter_theta) < 0.01 * C_DEG); /* with DS this assert will fail */
    }
#endif
return scatter_theta;
}

double scattering_angle_exit_deriv(const ion *incident, double sample_theta, double sample_phi, const detector *det) { /* Calculates the dtheta/dbeta derivative for geometrical straggling */
    double theta, phi;
    double scatter_theta, scatter_phi;
    double scatter_theta2, scatter_phi2;
    double theta_product, phi_product;
    double theta_product2, phi_product2;
    double epsilon = 0.001*C_DEG;
    rotate( -1.0*sample_theta, sample_phi, incident->theta, incident->phi, &theta, &phi);
    rotate(det->theta, det->phi, theta, phi, &scatter_theta, &scatter_phi);
    rotate(det->theta+epsilon, det->phi, theta, phi, &scatter_theta2, &scatter_phi2);
    rotate(det->theta, det->phi, sample_theta, sample_phi,  &theta_product, &phi_product);
    rotate(det->theta+epsilon, det->phi, sample_theta, sample_phi,  &theta_product2, &phi_product2);
    double deriv = (scatter_theta2-scatter_theta)/(theta_product2-theta_product) * -1.0; /* Since exit angle is pi - theta_product, the derivative dtheta/dbeta is -1.0 times this derivative calculated here */
    DEBUGMSG("Scat: %g deg, eps: %g deg, product %g deg, eps: %g deg. Deriv %.8lf (-1.0 for IBM)", scatter_theta/C_DEG,scatter_theta2/C_DEG, theta_product/C_DEG, theta_product2/C_DEG, deriv);
return deriv;
}

double exit_angle_delta(double sample_theta, double sample_phi, const detector *det, const aperture *beam_aperture, const char direction) {
    if(det->distance < 0.001 * C_MM) {
        return 0.0;
    }
    double sample_tilt = angle_tilt(sample_theta, sample_phi, direction);
    double det_tilt = C_PI - detector_angle(det, direction);
    double exit = C_PI - det_tilt - sample_tilt;
    double geo = cos(exit)/cos(sample_tilt);
    double w = aperture_width_shape_product(beam_aperture, direction);
    double delta_beam = w *  geo / det->distance;
    double delta_detector = aperture_width_shape_product(det->aperture, direction) / det->distance;
    double result = sqrt(pow2(delta_beam) + pow2(delta_detector));
    DEBUGMSG("Direction %c: detector angles %g deg, sample angle %g deg, exit angle %g deg", direction, det_tilt/C_DEG, sample_tilt/C_DEG, exit/C_DEG);
    DEBUGMSG("cos(alpha) = %g, cos(beta) = %g, cos(beta)/cos(alpha) = %g", cos(sample_tilt), cos(exit), geo);
    DEBUGMSG("Spread of exit angle in direction '%c' due to beam %g deg, due to detector %g deg. Combined %g deg FWHM.", direction, delta_beam/C_DEG, delta_detector/C_DEG, result/C_DEG);
    return result / C_FWHM;
}
geostragg_vars geostragg_vars_calculate(const ion *incident, double sample_theta, double sample_phi, const detector *det, const aperture *beam_aperture, int geostragg_enabled, int beta_manual_enabled) {
    geostragg_vars g;
    g.scatter_theta = scattering_angle(incident, sample_theta, sample_phi, det);
    if(beta_manual_enabled) {
        g.theta_product = C_PI - det->beta;
        g.phi_product = 0.0;
    } else {
        rotate(det->theta, det->phi, sample_theta, sample_phi, &g.theta_product, &g.phi_product); /* Detector in sample coordinate system */
    }
#ifdef DEBUG
    double beta;
    beta = (C_PI - g.theta_product);
#ifdef DEBUG_VERBOSE
    ion_print(stderr, incident);
#endif
    DEBUGMSG("Reaction product angles (in sample) %g deg and %g deg. Exit angle (beta) %g deg.", g.theta_product/C_DEG, g.phi_product/C_DEG, beta/C_DEG);
#endif
    if(!geostragg_enabled) { /* If geometric straggling is disabled, we have already calculated everything necessary */
        return g;
    }
    g.x.direction = 'x';
    g.y.direction = 'y';
    geostragg_vars_dir *gds[] = {&g.x, &g.y, NULL};
    for(geostragg_vars_dir **gdp = gds; *gdp; gdp++) {
        geostragg_vars_dir *gd = *gdp;
        if(gd->direction == 'x') {
            gd->theta_deriv = -cos(det->phi); /* TODO: check IBM with phi 180 deg and */
            gd->phi = 0.0;
        } else if (gd->direction == 'y') {
            gd->theta_deriv = -sin(det->phi); /* TODO: check Cornell with phi 270 deg */
            gd->phi = C_PI_2;
        }
        gd->delta_beta = exit_angle_delta(sample_theta, sample_phi, det, beam_aperture, gd->direction) / 2.0;
        gd->beta_deriv = fabs(beta_deriv(sample_theta, sample_phi, det, gd->direction));

        gd->theta_plus = g.scatter_theta + gd->delta_beta * gd->theta_deriv;
        gd->theta_minus = g.scatter_theta - gd->delta_beta * gd->theta_deriv;

        rotate(-1.0 * gd->delta_beta * gd->beta_deriv, g.phi_product, g.theta_product, g.phi_product, &gd->theta_product_plus, &gd->phi_product_plus); /* -1.0 again because of difference between theta and pi - theta */
        rotate(1.0 * gd->delta_beta * gd->beta_deriv, g.phi_product, g.theta_product, g.phi_product, &gd->theta_product_minus, &gd->phi_product_minus);

        DEBUGMSG("Spread in exit angle ('%c') %g deg", gd->direction, gd->delta_beta / C_DEG);
        DEBUGMSG("dBeta/dBeta_%c: %g", gd->direction, gd->beta_deriv); /* TODO: verify, check sign */
        DEBUGMSG("dTheta/dBeta_%c = %g", gd->direction, gd->theta_deriv); /* TODO: this should also be valid when sample_phi is not zero? */
        DEBUGMSG("Delta Beta_%c %g deg, Delta theta %g deg", gd->direction, -1.0 * gd->delta_beta * gd->beta_deriv, gd->delta_beta * gd->theta_deriv);
        DEBUGMSG("Reaction product (+) angles (in sample) %g deg and %g deg. Exit angle (beta) %g deg.", gd->theta_product_plus/C_DEG, gd->phi_product_plus/C_DEG, (C_PI - gd->theta_product_plus)/C_DEG);
        DEBUGMSG("Reaction product (+) scattering angle %g deg.", gd->theta_plus/C_DEG);
        DEBUGMSG("Reaction product (-) angles (in sample) %g deg and %g deg. Exit angle (beta) %g deg.", gd->theta_product_minus/C_DEG, gd->phi_product_minus/C_DEG, (C_PI - gd->theta_product_minus)/C_DEG);
        DEBUGMSG("Reaction product (-) scattering angle %g deg.", gd->theta_minus/C_DEG);
    }
    return g;
}

double geostragg(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params_exiting, const sample *sample, const sim_reaction *r, const geostragg_vars_dir *gd, const depth d, const double E_0) {
    ion ion;
    ion = r->p;
    ion_set_angle(&ion, gd->theta_product_plus, gd->phi_product_plus);
    ion.E = reaction_product_energy(r->r, gd->theta_plus, E_0);
    ion.S = 0.0; /* We don't need straggling for anything, might as well reset it */
    stop_sample_exit(stop, stragg, params_exiting, &ion, d, sample);
    double Eplus = ion.E;
    ion = r->p;
    ion_set_angle(&ion, gd->theta_product_minus, gd->phi_product_minus);
    ion.E = reaction_product_energy(r->r, gd->theta_minus, E_0);
    ion.S = 0.0; /* We don't need straggling for anything, might as well reset it */
    stop_sample_exit(stop, stragg, params_exiting, &ion, d, sample);
    double Eminus = ion.E;
    double result = pow2(Eplus - Eminus);
    DEBUGVERBOSEMSG("Direction (%c), Eplus %g keV, Eminus %g keV, Straggling %g keV FWHM", gd->direction, Eplus/C_KEV, Eminus/C_KEV, C_FWHM*sqrt(result)/C_KEV);
    return result;
}

double beta_deriv(double sample_theta, double sample_phi, const detector *det, const char direction) { /* TODO: replace by analytical formula */
    double theta, phi;
    double theta_product, phi_product;
    static const double delta = 0.001*C_DEG;
    rotate(det->theta, det->phi, sample_theta, sample_phi, &theta_product, &phi_product); /* Detector in sample coordinate system */
    theta = delta;
    if(direction == 'x') {
        phi = 0.0;
    } else if(direction == 'y') {
        phi = C_PI_2;
    } else {
        return 0.0;
    }
    rotate(theta, phi, det->theta, det->phi, &theta, &phi);
    rotate(theta, phi, sample_theta, sample_phi, &theta, &phi); /* "new" Detector in sample coordinate system */
    double result = -1.0*(theta - theta_product)/delta; /* -1.0, since beta = pi - theta */
    DEBUGMSG("Beta deriv ('%c') got theta = %g deg (orig %g deg), diff in beta %g deg, result %g", direction, theta/C_DEG, theta_product/C_DEG,  (theta_product - theta)/C_DEG, result);
    return result;
}

double exit_angle(double sample_theta, double sample_phi, double det_theta, double det_phi) {
    double theta, phi;
    rotate(sample_theta, sample_phi, det_theta, det_phi, &theta, &phi);
    return C_PI - theta;
}
