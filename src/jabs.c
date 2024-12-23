/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <jibal_units.h>
#include <jibal_kin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "jabs_debug.h"
#include "geostragg.h"
#include "roughness.h"
#include "jabs.h"
#include "defaults.h"
#include "message.h"
#include "stop.h"
#include "win_compat.h"

double cross_section_straggling_fixed(const sim_reaction *sim_r, const prob_dist *pd, double E, double S) {
    const double std_dev = sqrt(S);
    double cs_sum = 0.0;
    for(size_t i = 0; i < pd->n; i++) {
        prob_point *pp = &(pd->points[i]);
        double E_stragg = E + pp->x * std_dev;
        double cs = pp->p * sim_r->cross_section(sim_r, E_stragg);
        cs_sum += cs;
        DEBUGVERBOSEMSG("%zu %10g %10g %10g %10g", i, pp->x, E_stragg/C_KEV, pp->p, cs/C_MB_SR);
    }
#ifdef DEBUG_VERBOSE
    double unweighted = sim_r->cross_section(sim_r, E);
    double diff = (cs_sum - unweighted)/unweighted;
    DEBUGVERBOSEMSG("Got cs %.7lf mb/sr (unweighted by straggling %.7lf mb/sr) diff %.7lf%%.", cs_sum/C_MB_SR, unweighted/C_MB_SR, 100.0*diff);
#endif
    return cs_sum;
}


double cross_section_straggling(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, const prob_dist *pd, double E, double S) {
    if(w) {
        return cross_section_straggling_adaptive(sim_r, w, accuracy, E, S);
    }
    if(pd) {
        return cross_section_straggling_fixed(sim_r, pd, E, S);
    }
    return sim_r->cross_section(sim_r, E); /* Fallback, no weights applied */
}


struct cs_stragg_int_params {
    const sim_reaction *sim_r;
    double sigma;
    double E_mean;
};

double cs_stragg_function(double x, void *params) {
    struct cs_stragg_int_params *p = (struct cs_stragg_int_params *) params;
    double a = (x - p->E_mean) / p->sigma;
    double result = exp(-0.5 * a * a) * p->sim_r->cross_section(p->sim_r, x); /* Gaussian (not normalized!) times cross section */
    DEBUGVERBOSEMSG("cs_stragg_function(), E = %g keV, S = %g keV FWHM, E_mean = %g keV, a = %g. Result %g mb/sr.", x/C_KEV, p->sigma*C_FWHM/C_KEV, p->E_mean/C_KEV, a, (0.398942280401432703/p->sigma) * result/C_MB_SR );
    return result;
}

struct cs_int_params {
    const depth *d_before;
    const depth *d_after;
    double E_front;
    double S_front;
    const sample *sample;
    const sim_reaction *sim_r;
    double stop_slope;
    double stragg_slope;
    const prob_dist *cs_stragg_pd;
    gsl_integration_workspace *w;
    double stragg_int_accuracy;
    depth d; /* changes between calls */
};

double cross_section_straggling_adaptive( const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, double E, double S) { /* Uses real numerical integration */
    struct cs_stragg_int_params params;
    params.sigma = sqrt(S);
    params.E_mean = E;
    params.sim_r = sim_r;
    gsl_function F;
    F.function = &cs_stragg_function;
    F.params = &params;
    double E_low = E - 4.0*params.sigma;
    if(E_low < 1.0*C_KEV)
        E_low = 1.0*C_KEV;
    double E_high = E + 4.0*params.sigma;
    double result, error;
#ifdef SINGULARITIES
    gsl_integration_qags(&F, E_low, E_high, 0, accuracy, w->limit,w, &result, &error);
#else
    gsl_integration_qag(&F, E_low, E_high, 1e-6 * C_MB_SR, accuracy, w->limit, GSL_INTEG_GAUSS15, w, &result, &error); /* For some reason 0 as epsabs doesn't work anymore so a very small number is used instead */
#endif
    static const double inv_sqrt_2pi = 0.398942280401432703;
    result *= (inv_sqrt_2pi / params.sigma) * 1.000063346496191; /* Normalize gaussian. The 1.00006 accounts for tails outside +- 4 sigmas */
    DEBUGVERBOSEMSG("Integrated from %g keV to %g keV in %zu steps (limit %zu), got (%g +- %g) mb/sr", E_low/C_KEV, E_high/C_KEV, w->size, w->limit, result/C_MB_SR, error/C_MB_SR);
    return result;
}

double cs_function(double x, void * params) {
    struct cs_int_params *p = (struct cs_int_params *) params;
    p->d.x = p->d_before->x + p->stop_slope * (x - p->E_front); /* Depth assuming constant stopping inside brick */
    double c = get_conc(p->sample, p->d, p->sim_r->i_isotope);
    double sigma;
    double S = p->S_front + p->stragg_slope * (x - p->E_front);
    sigma = cross_section_straggling(p->sim_r, p->w, p->stragg_int_accuracy, p->cs_stragg_pd, x, S);
    DEBUGVERBOSEMSG("Depth %g tfu, energy %g keV, stragg %g keV, sigma %g mb/sr, c %g %%", p->d.x/C_TFU, x/C_KEV, C_FWHM * sqrt(S)/C_KEV, sigma/C_MB_SR, c/C_PERCENT);
    return c*sigma;
}

double cross_section_concentration_product_adaptive(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after,
                                                    double S_front, double S_back) {
    double result, error;
    struct cs_int_params params;
    params.sim_r = sim_r;
    params.d_before = d_before;
    params.d_after = d_after;
    params.E_front = E_front;
    params.stop_slope = (d_after->x - d_before->x)/(E_back - E_front);
    params.S_front = S_front;
    params.d.i = d_after->i;
    params.sample = sample;
    params.w = ws->w_int_cs_stragg; /* Can be NULL */
    params.cs_stragg_pd = ws->params->cs_stragg_pd; /* Can be NULL */
    params.stragg_slope = (S_back-S_front)/(E_back-E_front);
    gsl_function F;
    F.function = &cs_function;
    F.params = &params;
    gsl_set_error_handler_off();
#ifdef SINGULARITIES
    gsl_integration_qags(&F, E_back, E_front, 0, ws->params->int_cs_accuracy, ws->w_int_cs->limit,
                          ws->w_int_cs, &result, &error);
#else
    gsl_integration_qag(&F, E_back, E_front, 0, ws->params->int_cs_accuracy, ws->w_int_cs->limit, GSL_INTEG_GAUSS15,
                          ws->w_int_cs, &result, &error);
#endif
    double final = result/(E_front - E_back);
    DEBUGVERBOSEMSG("E from %g keV to %g keV, diff %g keV", E_front/C_KEV, E_back/C_KEV, (E_front - E_back) / C_KEV);
    DEBUGVERBOSEMSG("stop avg %g eV/tfu", params.stop_slope/(C_EV/C_TFU));
    DEBUGVERBOSEMSG("integration result          = % 18g", result);
    DEBUGVERBOSEMSG("integration estimated error = % 18g", error);
    DEBUGVERBOSEMSG("integration intervals       = %zu", ws->w_int_cs->size);
    DEBUGVERBOSEMSG("final result                = %g mb/sr", final/C_MB_SR);
    return final;
}

double cross_section_concentration_product_stepping(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back) {
    const double E_step_max = ws->params->cs_energy_step_max;
    const double depth_step_max = ws->params->cs_depth_step_max;
    const double S_avg_FWHM = ws->params->cs_stragg_step_sigmas * sqrt((S_front + S_back) / 2.0);
    const double E_step_nominal = -1.0 * GSL_MIN_DBL(E_step_max, S_avg_FWHM); /* Actual step is negative */
    const double d_diff = fabs(d_after->x - d_before->x); /* Since depth_step_max is (typically) always absolute number, we should use absolute value here to determine number of steps. */
    assert(d_diff > 0);
    assert(E_step_nominal < 0.0);
    double sigma;
    double c; /* Concentration (but only if constant!) */
    if(sample->no_conc_gradients) {
        c = *sample_conc_bin(sample, d_before->i, sim_r->i_isotope); /* Concentration weighting was not done inside the loop */
        if(c < CONC_TOLERANCE) {
            return 0.0;
        }
    }
    depth d = *d_before;
    double E_diff = E_back - E_front;
    size_t n_steps = ceil(GSL_MAX_DBL(E_diff/E_step_nominal, d_diff/depth_step_max)); /* Choose the bigger of two evils. Should not be possible to take zero steps since d_diff should always be greater than zero. */
    const double frac = 1.0/(1.0*(n_steps+1));
    const double x_step = (d_after->x - d_before->x) * frac;
    const double E_step = (E_back - E_front) * frac;
    const double S_step = (S_back - S_front) * frac;
#ifdef JABS_DEBUG_CS
    fprintf(stderr, "E_step_nominal = %g keV, n_steps = %zu, S_avg %g keV, E = %g ... %g keV, actual E_step = %g keV\n", E_step_nominal / C_KEV, n_steps, S_avg_FWHM / C_KEV, E_front / C_KEV, E_back / C_KEV, E_step / C_KEV);
#endif
    double sum = 0.0;
    for(size_t i = 1; i <= n_steps; i++) { /* Compute cross section and concentration product in several "sub-steps" */
        double E = E_front + E_step * i;
        double S = S_front + S_step * i;
        d.x = d_before->x + x_step * i;
        double sigma_partial = cross_section_straggling(sim_r, ws->w_int_cs_stragg, ws->params->int_cs_stragg_accuracy, ws->params->cs_stragg_pd, E, S);
        if(!sample->no_conc_gradients) { /* We have concentration gradient */
            double c = get_conc(sample, d, sim_r->i_isotope);
            sigma_partial *= c;
        }
        sum += sigma_partial;
#ifdef JABS_DEBUG_CS
        fprintf(stderr, "i = %zu, E = %g keV, S = %g keV, (E_front = %g keV, E_back = %g keV), sigma = %g mb/sr\n", i, E/C_KEV, C_FWHM * sqrt(S)/C_KEV, E_front/C_KEV, E_back/C_KEV, sigma_partial / C_MB_SR);
#endif
    }
    sigma = sum / n_steps;
    if(sample->no_conc_gradients) {
        sigma *= c; /* Concentration weighting was not done inside the loop */
    }
    return sigma;
}


double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back) {
    double sigmaconc;
    if(ws->params->mean_conc_and_energy) {
        double c;
        if(sample->no_conc_gradients) {
            c = *sample_conc_bin(sample, d_before->i, sim_r->i_isotope);
        } else {
            depth d;
            d.i = d_before->i;
            d.x = (d_before->x + d_after->x) / 2.0;
            c = get_conc(sample, d, sim_r->i_isotope);
        }
        if(c < CONC_TOLERANCE) {
            return 0.0;
        }
        sigmaconc = sim_r->cross_section(sim_r, (E_front + E_back)/2.0) * c;
    } else if(ws->params->cs_adaptive) {
        sigmaconc = cross_section_concentration_product_adaptive(ws, sample, sim_r, E_front, E_back, d_before, d_after, S_front, S_back);
    } else {
        sigmaconc = cross_section_concentration_product_stepping(ws, sample, sim_r, E_front, E_back, d_before, d_after, S_front, S_back);
    }
    if(sample->ranges[d_before->i].yield != 1.0 || sample->ranges[d_before->i].yield_slope != 0) {
        double depth = (d_before->x + d_after->x) / 2.0 - sample->ranges[d_before->i].x; /* How deep are we (on average), for yield slope calculation */
        double yield = sample->ranges[d_before->i].yield + sample->ranges[d_before->i].yield_slope * (depth / C_TFU);
        yield = GSL_MAX_DBL(yield, 0.0); /* Don't allow negative yields */
        sigmaconc *= yield;
    }
    return sigmaconc;
}

int simulate_reaction(const ion *incident, const depth depth_start, sim_workspace *ws, const sample *sample, const des_table *dt, const geostragg_vars *g, sim_reaction *sim_r) {
    ion ion1 = *incident; /* Shallow copy */
    if(simulate_init_reaction(sim_r, sample, ws->params, g, ws->emin, ion1.ion_gsto->emin, ion1.E + ws->params->sigmas_cutoff * sqrt(ion1.S))) {
        return EXIT_FAILURE;
    }
    if(sim_r->stop) {
        sim_r->last_brick = 0;
        return EXIT_SUCCESS;
    }
    depth d_before, d_after = depth_start;
    size_t i_des = 0;
    brick *b = NULL, *b_prev = NULL;
    int skipped, last = FALSE;
#ifdef DEBUG_VERBOSE
    int crossed;
#endif
    sim_r->last_brick = 0;
    const des *des_min = des_table_min_energy_bin(dt);
    int product_and_incident_go_in_different_directions = (incident->inverse_cosine_theta * sim_r->p.inverse_cosine_theta < 0.0); /* false when transmission, true usually. */
    /* When the above is true, we can safely assume that once reaction product energy goes below some energy, we can stop calculating. */
    size_t i_brick = 0;
    for(i_brick = 0; i_brick < sim_r->n_bricks; i_brick++) {
        assert(ion1.S >= 0.0);
        skipped = FALSE;
#ifdef DEBUG_VERBOSE
        crossed = FALSE;
#endif
        b_prev = b;
        b = &sim_r->bricks[i_brick];
        b->valid = TRUE; /* Will be invalidated if necessary */
        d_before = d_after;
        DEBUGVERBOSEMSG("E = %g keV (incident), i_brick = %zu", ion1.E / C_KEV, i_brick);
        if(ion1.E < des_min->E) {
            DEBUGMSG("E = %g keV (incident) is below %g keV. Changing energy.", ion1.E / C_KEV, des_min->E / C_KEV);
            des_set_ion(des_min, &ion1);
            d_after = des_min->d; /* des_table_find_depth() might change the depth */
        }
        if(i_brick != 0) {
            d_after = des_table_find_depth(dt, &i_des, d_before, &ion1); /* Does this handle E below min? */
        }
        if(ion1.E <= des_min->E) {
            DEBUGMSG("This will be the last brick due to incident energy hitting DES minimum %g keV. d_after = %g tfu.", des_min->E / C_KEV, d_after.x / C_TFU);
            last = TRUE;
        }
        if(i_brick == 0 || d_after.i != d_before.i) { /* There was a layer (depth range) crossing. If step() took this into account when making DES table the only issue is the .i index. depth (.x) is not changed. */
            if(ws->params->bricks_skip_zero_conc_ranges) {
                double conc_start = *sample_conc_bin(sample, d_after.i, sim_r->i_isotope);
                double conc_stop = *sample_conc_bin(sample, d_after.i + 1, sim_r->i_isotope);
                DEBUGVERBOSEMSG("Brick %zu crosses into range %zu, d_before = %g tfu.", i_brick, d_after.i, d_before.x / C_TFU);
                DEBUGVERBOSEMSG("Concentration varies between %g%% and %g%%", conc_start * 100.0, conc_stop * 100.0);
                if(conc_start < CONC_TOLERANCE && conc_stop < CONC_TOLERANCE) { /* This isotope concentration is zero in this layer, skip to next one */
                    d_after = des_next_range(dt, &ion1, d_after);
                    skipped = TRUE;
                    DEBUGMSG("Skipped to %g.", d_after.x / C_TFU);
                }
            }
            d_before.i = d_after.i;
#ifdef DEBUG_VERBOSE
            crossed = TRUE;
#endif
        }
        assert(ion1.S >= 0.0);
        b->d = d_after;
        b->E_0 = ion1.E;
        b->S_0 = ion1.S;
        if(ws->params->geostragg) {
            b->S_geo_x = geostragg(&ws->stop, &ws->stragg, &ws->params->exiting_stop_params, sample, sim_r, &(g->x), d_after, b->E_0);
            b->S_geo_y = geostragg(&ws->stop, &ws->stragg, &ws->params->exiting_stop_params, sample, sim_r, &(g->y), d_after, b->E_0);
        }
        sim_reaction_product_energy_and_straggling(sim_r, &ion1); /* sets sim_r->p */
        assert(sim_r->p.S >= 0.0);
        b->E_r = sim_r->p.E;
        b->S_r = sim_r->p.S;

        if(stop_sample_exit(&ws->stop, &ws->stragg, &ws->params->exiting_stop_params, &sim_r->p, d_after, sample) != 0) {
            DEBUGMSG("Stop before exit, energy of reaction product after reaction %g keV.", b->E_r / C_KEV);
            if(product_and_incident_go_in_different_directions) { /* Lowering incident energy will lower reaction product energy */
                DEBUGSTR("We can stop calculation, because lowering incident energy will also lower reaction product energy. This will be the last brick.");
                last = TRUE;
            }
        }
        b->E_s = sim_r->p.E;
        b->S_s = sim_r->p.S;

        if(ws->det->foil) { /* Energy loss in detector foil */
            depth d_foil = {.i = 0, .x = 0.0};
            ion ion_foil = *&sim_r->p;
            ion_set_angle(&ion_foil, 0.0, 0.0); /* Foils are not tilted. We use a temporary copy of "p" to do this step. */
            if(stop_sample_exit(&ws->stop, &ws->stragg, &ws->params->exiting_stop_params, &ion_foil, d_foil, ws->det->foil)) {
                DEBUGMSG("Stop in detector foil. Energy after reaction was %g keV.", b->E_r / C_KEV);
                if(product_and_incident_go_in_different_directions) { /* Lowering incident energy will lower reaction product energy */
                    DEBUGSTR("We can stop calculation, because lowering incident energy will also lower reaction product energy. This will be the last brick.");
                    last = TRUE; /* Discards entire brick, may not be the perfect choice! */
                    b->valid = FALSE;
                }
            }
            b->E = ion_foil.E;
            b->S = ion_foil.S;
        } else {
            b->E = sim_r->p.E;
            b->S = sim_r->p.S;
        }

        double E_deriv;
        double sigma_conc;

        if(b_prev) {
            double d_diff = depth_diff(d_before, d_after);
            if(d_diff == 0) { /* Zero thickness brick, so no energy change either */
                b->dE = 0.0;
                b->thick = 0.0;
                E_deriv = 10.0; /* Can't calculate derivative, so we assume it is large */
                b->effective_stopping = 1e12; /* Can't calculate due to zero thickness, so assume it's large */
                DEBUGVERBOSEMSG("Zero thickness brick, setting derivative to %g.", E_deriv);
                sigma_conc = 0.0;
            } else {
                b->thick = d_diff;
                sigma_conc = cross_section_concentration_product(ws, sample, sim_r, b_prev->E_0, b->E_0, &d_before, &d_after, b_prev->S_0, b->S_0); /* Product of concentration and sigma for isotope i_isotope target and this reaction. */
                assert(b_prev->E_0 > b->E_0);
                double E0_diff = b_prev->E_0 - b->E_0;
                b->dE = b_prev->E - b->E;
                E_deriv = fabs(b->dE / E0_diff); /* how many keVs does the reaction product energy change for each keV of incident ion1 energy change */
                assert(E_deriv > 0.0);
            }
            b->effective_stopping = b->dE / b->thick; /* Effective stopping from detector point-of-view, reciprocal is used in depth resolution calculation */
        } else { /* First brick */
            sigma_conc = 0.0;
            b->dE = 0.0;
            b->thick = 0.0;
            double stop_incident = stop_sample(&ws->stop, &ion1, sample, d_after, b->E_0) * ion1.inverse_cosine_theta;
            double stop_exiting = stop_sample(&ws->stop, &sim_r->p, sample, d_after, b->E_r) * sim_r->p.inverse_cosine_theta;
            double K = reaction_product_energy(sim_r->r, sim_r->theta, b->E_0) / b->E_0;
            b->effective_stopping = K * stop_incident - stop_exiting;
            DEBUGVERBOSEMSG("Calculating E_deriv = (%g * (%g eV/tfu) - (%g eV/tfu)) / (%g eV/tfu)", K, stop_incident / C_EV_TFU, stop_exiting / C_EV_TFU, stop_incident / C_EV_TFU);
            if(stop_incident < 0.001  * C_EV_TFU) {
                E_deriv = 10.0; /* Reasonable default if there is no stopping? */
            } else {
                E_deriv = fabs(b->effective_stopping / (stop_incident));
            }
        }
        DEBUGVERBOSEMSG("E_deriv before imposing min/max clamping was %g.", E_deriv);
        E_deriv = GSL_MAX_DBL(E_deriv, ENERGY_DERIVATIVE_MIN);
        E_deriv = GSL_MIN_DBL(E_deriv, ENERGY_DERIVATIVE_MAX);
        b->deriv = E_deriv;
        b->sc = sigma_conc;
        b->Q = ion1.inverse_cosine_theta * sigma_conc * b->thick;
        DEBUGVERBOSEMSG("%s %s %3zu %3zu:%10.3lf %3zu:%10.3lf %3zu %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3e %8.3lf %2i",
                sim_r->r->target->name, reaction_type_to_string(sim_r->r->type), i_brick,
                d_before.i, d_before.x / C_TFU,
                d_after.i, d_after.x / C_TFU,
                i_des,
                b->E_0 / C_KEV, sqrt(b->S_0) / C_KEV,
                b->E_r / C_KEV, sqrt(b->S_r) / C_KEV,
                b->E_s / C_KEV, sqrt(b->S_s) / C_KEV,
                b->E / C_KEV, sqrt(b->S) / C_KEV,
                E_deriv, get_conc(sample, d_after, sim_r->i_isotope) * 100.0, sigma_conc / C_MB_SR, b->Q, b->deriv,
                skipped
                );
        assert(!isnan(ion1.E));
        if(ion1.E < sim_r->emin_incident) {
            DEBUGMSG("E = %g keV sufficiently below reaction emin_incident %g keV.", b->E / C_KEV, sim_r->emin_incident / C_KEV);
            break;
        }
        if(ion1.inverse_cosine_theta > 0.0 && d_after.x >= sim_r->max_depth) {
            DEBUGMSG("Max depth of %g tfu reached (d_after.x is %g).", sim_r->max_depth / C_TFU,  d_after.x / C_TFU);
            break;
        }
        if(ion1.inverse_cosine_theta < 0.0 && d_after.x < DEPTH_TOLERANCE) {
            DEBUGMSG("Surface reached, d_after.x = %g tfu.", d_after.x);
            break;
        }
        if(last) {
            DEBUGMSG("Last brick (earlier decision, because des_min = %g keV or because reaction product stopped in the sample or detector foil).", des_min->E / C_KEV);
            break;
        }
        if(!skipped) {
            double S_sigma = sqrt(detector_resolution(ws->det, sim_r->p.isotope, b->E) + b->S);
            assert(S_sigma > 0.0);
            double E_change = -ws->params->brick_width_sigmas * S_sigma / E_deriv;
            assert(E_change < 0.0);
            ion1.E += E_change;
        }
    }
    if(i_brick == sim_r->n_bricks) {
        DEBUGMSG("Maybe not enough bricks (%zu)!", sim_r->n_bricks);
        sim_r->last_brick = i_brick - 1; /* For loop adds +1 if it runs through completely */
    } else {
        sim_r->last_brick = i_brick; /* There was a break before the for loop completed */
    }
    return EXIT_SUCCESS;
}

int simulate(const ion *incident, const depth depth_start, sim_workspace *ws, const sample *sample) {
    DEBUGMSG("Starting simulate(ion = %s (E = %.3lf keV, S = %.3lf keV, angles = %.3lf deg, %.3lf deg in sample), depth_start = %g tfu (i = %zu), ...)",
            incident->isotope->name, incident->E / C_KEV, C_FWHM * sqrt(incident->S) / C_KEV, incident->theta / C_DEG, incident->phi / C_DEG,
            depth_start.x / C_TFU, depth_start.i);
    geostragg_vars g = geostragg_vars_calculate(incident, ws->sim->sample_theta, ws->sim->sample_phi,
                                                ws->det, ws->sim->beam_aperture,
                                                ws->params->geostragg, ws->params->beta_manual);
    des_table *dt = des_table_compute(&ws->stop, &ws->stragg, ws->params, sample, incident, depth_start, ws->emin); /* Depth, energy and straggling of incident ion */
    if(!dt) {
        jabs_message(MSG_ERROR, "DES table computation failed.\n");
        return -1;
    }
#ifdef DEBUG
    des_table_print(stderr, dt);
#endif
    int error = FALSE;
    for(size_t i_reaction = 0; i_reaction < ws->n_reactions; i_reaction++) {
        sim_reaction *sim_r = ws->reactions[i_reaction];
        DEBUGMSG("Simulating reaction i_reaction = %zu type %s target %s (i_isotope = %zu, i_jibal = %zu) product %s (i_jibal = %zu)",
                i_reaction, reaction_type_to_string(sim_r->r->type),
                sim_r->r->target->name, sim_r->i_isotope, sim_r->r->target->i,
                sim_r->r->product->name, sim_r->r->product->i);
        if(simulate_reaction(incident, depth_start, ws, sample, dt, &g, sim_r)) {
            jabs_message(MSG_ERROR, "Simulating reaction %zu (%s) failed.\n", i_reaction + 1, reaction_name(sim_r->r));
            DEBUGMSG("Simulating reaction i_reaction = %zu failed.", i_reaction);
            error = TRUE;
            break;
        }
        DEBUGMSG("Finished. sim_r->last_brick = %zu (%zu/%zu). depth = %g tfu", sim_r->last_brick, sim_r->last_brick + 1, sim_r->n_bricks, sim_r->bricks[sim_r->last_brick].d.x / C_TFU);
        if(sim_r->last_brick + 1 == sim_r->n_bricks) {
            jabs_message(MSG_WARNING, "Reaction %s may have produced incomplete results. Use set bricks_n to increase number of bricks from current setting (%zu)", reaction_name(sim_r->r), sim_r->n_bricks);
        }
    }
    des_table_free(dt);
    if(error) {
        return -1;
    }
    size_t n_meaningful = sim_workspace_histograms_calculate(ws);
    DEBUGMSG("Finished simulate(), %zu out of %zu reactions managed to actually produce something.", n_meaningful, ws->n_reactions);
    return (int)n_meaningful;
}

int simulate_init_reaction(sim_reaction *sim_r, const sample *sample, const sim_calc_params *params, const geostragg_vars *g, double emin, double emin_incident, double emax_incident) {
    assert(sim_r);
    sim_r->last_brick = 0;
    sim_r->stop = FALSE;
    ion_set_angle(&sim_r->p, g->theta_product, g->phi_product);
    if(sim_reaction_recalculate_internal_variables(sim_r, params, g->scatter_theta, emin, emin_incident, emax_incident)) {
        jabs_message(MSG_ERROR, "Recalculating internal variables of reaction failed.\n");
        return EXIT_FAILURE;
    }
    if(sim_r->stop) {
        sim_r->max_depth = 0.0;
        return EXIT_SUCCESS;
    }
    sim_r->max_depth = sample_isotope_max_depth(sample, sim_r->i_isotope);
    sim_reaction_reset_bricks(sim_r);
    if(sim_r->i_isotope >= sample->n_isotopes) { /* No target isotope for reaction. */
        DEBUGMSG("No target isotope (%zu) for reaction", sim_r->i_isotope);
        sim_r->stop = TRUE;
    }
    DEBUGMSG("Simulation reaction %s initialized. Max depth %g tfu. i_isotope=%zu, stop = %i.",
            reaction_name(sim_r->r), sim_r->max_depth / C_TFU, sim_r->i_isotope, sim_r->stop);
    return EXIT_SUCCESS;
}

int assign_stopping_Z2(jibal_gsto *gsto, const simulation *sim, int Z2) { /* Assigns stopping and straggling (GSTO) for given Z2. Goes through all possible Z1s (beam and reaction products). */
    int fail = FALSE;
    int Z1 = sim->beam_isotope->Z;
    DEBUGMSG("Assigning stopping in Z2 = %i.", Z2);
    if(assign_stopping_Z1_Z2(gsto, Z1, Z2)) {
        jabs_message(MSG_ERROR, "Can not assign stopping or straggling (beam Z1 = %i). Z2 = %i.\n", Z1, Z2);
        fail = TRUE;
    }
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        const reaction *r = sim->reactions[i_reaction];
        if(!r)
            continue;
        if(r->product->Z == Z1) {/* Z1 repeats, skip */
            continue;
        }
        Z1 = r->product->Z;
        if(assign_stopping_Z1_Z2(gsto, Z1, Z2)) {
            jabs_message(MSG_ERROR, "Can not assign stopping or straggling for reaction product. Reaction: %s.\n", reaction_name(r));
            fail = TRUE;
        }
    }
    return fail;
}

int assign_stopping_Z1_Z2(jibal_gsto *gsto, int Z1, int Z2) {
    if(Z1 < 1 || Z2 < 1) { /* Not a proper ion, stopping can't be assigned */
        DEBUGMSG("Assigning stopping Z1 = %i, Z2 = %i is a potential issue.", Z1, Z2);
        if(Z1 != 0) { /* Neutrons are a possible legitimate reaction product, don't warn */
            jabs_message(MSG_WARNING, "Assigning stopping for Z1 = %i in Z2 = %i is not possible.\n", Z1, Z2);
        }
        return EXIT_SUCCESS; /* Skips, doesn't fail */
    }
    if(!jibal_gsto_auto_assign(gsto, Z1, Z2)) {
        jabs_message(MSG_ERROR, "Assigning stopping for Z1 = %i in Z2 = %i fails.\n", Z1, Z2);
        return EXIT_FAILURE;
    }
    if(!jibal_gsto_get_assigned_file(gsto, GSTO_STO_ELE, Z1, Z2)) {
        jabs_message(MSG_ERROR, "No electronic stopping assigned for Z1 = %i in Z2 = %i.\n", Z1, Z2);
        return EXIT_FAILURE;
    }
    if(!jibal_gsto_get_assigned_file(gsto, GSTO_STO_STRAGG, Z1, Z2)) {
        jabs_message(MSG_ERROR, "No energy loss straggling assigned for Z1 = %i in Z2 = %i.\n", Z1, Z2);
        return EXIT_FAILURE;
    }
    DEBUGMSG("Assigning stopping Z1 = %i, Z2 = %i was a great success.", Z1, Z2);
    return EXIT_SUCCESS;
}

int assign_stopping(jibal_gsto *gsto, const simulation *sim) {
    /* TODO: simplify this by finding all possible Z1, Z2 combinations, considering target elements, beam and reactions before attempting to assign stopping/straggling (GSTO) */
    sample *sample = sim->sample;
    if(!sample) {
        jabs_message(MSG_ERROR, "Could not assign stopping, because sample is not set!\n");
        return 1;
    }
    int Z2 = -1;
    int fail = FALSE; /* We don't return immediately on failure since the user might want to know ALL failures at one go. */
    for(size_t i = 0; i < sample->n_isotopes; i++) {
        if(Z2 == sample->isotopes[i]->Z) {/* Z2 repeats, skip */
            continue;
        }
        Z2 = sample->isotopes[i]->Z;
        if(assign_stopping_Z2(gsto, sim, Z2)) {
            fail = TRUE;
        }
    }
    for(size_t i_det = 0; i_det < sim->n_det; i_det++) {
        const detector *det = sim->det[i_det];
        if(!det || !det->foil)
            continue;
        for(size_t i = 0; i < det->foil->n_isotopes; i++) {
            if(Z2 == det->foil->isotopes[i]->Z) {/* Z2 repeats, skip */
                continue;
            }
            Z2 = det->foil->isotopes[i]->Z;
            if(assign_stopping_Z2(gsto, sim, Z2)) {
                fail = TRUE;
            }
        }
    }
    return fail ? EXIT_FAILURE : EXIT_SUCCESS;
}

int simulate_with_roughness(sim_workspace *ws) {
    int status = 0;
    const double fluence_original = ws->fluence;
    size_t n_rl = 0; /* Number of rough layers */
    for(size_t i = 0; i < ws->sample->n_ranges; i++) {
        sample_range *r = &(ws->sample->ranges[i]);
        if(r->rough.model == ROUGHNESS_FILE) {
            n_rl++;
            continue;
        }
        r->rough.n *= ws->params->rough_layer_multiplier;
        if(r->rough.n > ROUGHNESS_SUBSPECTRA_MAXIMUM) { /* Artificial limit to n */
            r->rough.n = ROUGHNESS_SUBSPECTRA_MAXIMUM;
        }
        if(r->rough.n == 0) {
            r->rough.model = ROUGHNESS_NONE;
        }
        if(r->rough.model == ROUGHNESS_GAMMA) {
            n_rl++;
        }
    }
    DEBUGMSG("%zu rough layers", n_rl);
    if(!n_rl) {
        ion_set_angle(&ws->ion, 0.0 * C_DEG, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        return simulate(&ws->ion, depth_seek(ws->sample, 0.0 * C_TFU), ws, ws->sample);
    }
    struct sample *sample_rough = sample_copy(ws->sample);
    size_t i_rl = 0;
    thick_prob_dist **tpds = calloc(n_rl, sizeof(thick_prob_dist *));
    if(!tpds) {
        return -1;
    }
    for(size_t i = 0; i < ws->sample->n_ranges; i++) {
        sample_range *r = &(ws->sample->ranges[i]);
        if(r->rough.model == ROUGHNESS_NONE) {
            continue;
        }
        tpds[i_rl] = NULL;
        if(r->rough.model == ROUGHNESS_GAMMA) {
            DEBUGMSG("Range %zu is rough (gamma), amount %g tfu, n = %zu spectra", i, ws->sample->ranges[i].rough.x/C_TFU, ws->sample->ranges[i].rough.n);
            tpds[i_rl] = thickness_probability_table_gamma(r->x, r->rough.x, r->rough.n);
        } else if(r->rough.model == ROUGHNESS_FILE) {
            DEBUGMSG("Range %zu is rough (FILE), filename = \"%s\"", i, r->rough.file->filename);
            tpds[i_rl] = thickness_probability_table_copy(r->rough.file->tpd);
        }
        if(tpds[i_rl]) {
            tpds[i_rl]->i_range = i;
            DEBUGMSG("TPD (i_range %zu) for depth %zu (%.3lf tfu nominal), roughness %.3lf tfu:", tpds[i_rl]->i_range, i, ws->sample->ranges[i].x/C_TFU, ws->sample->ranges[i].rough.x/C_TFU);
#ifdef DEBUG_VERBOSE
            thickness_probability_table_print(stderr, tpds[i_rl]);
#endif
            i_rl++;
        }
    }
    size_t iter_total = 1;
    for(i_rl = 0; i_rl < n_rl; i_rl++) { /* Calculate cumulative product (number of subspectra) */
        thick_prob_dist *tpd = tpds[i_rl];
        if(!tpd) {
            continue;
        }
        tpd->modulo = iter_total;
        iter_total *= tpd->n;
        DEBUGMSG("TPD %zu/%zu, modulo %zu, total %zu (cumulating)\n", i_rl, n_rl, tpd->modulo, iter_total)
    }
    for(size_t i_iter = 0; i_iter < iter_total; i_iter++) {
        DEBUGMSG("Roughness step %zu/%zu.", i_iter+1, iter_total);
        double p = 1.0;
        for(size_t i_range = 0; i_range < ws->sample->n_ranges; i_range++) { /* Reset ranges for every iter */
            sample_rough->ranges[i_range].x = ws->sample->ranges[i_range].x;
        }
        for(i_rl = 0; i_rl < n_rl; i_rl++) {
            thick_prob_dist *tpd = tpds[i_rl]; /* One particular thickness probability distribution ("i"th one) */
            if(!tpd) {
                continue;
            }
            size_t j = (i_iter / tpd->modulo) % tpd->n; /* "j"th roughness element */
            thick_prob *pj = &tpds[i_rl]->p[j]; /* ..is this one */
            p *= pj->prob; /* Probability is multiplied by the "i"th roughness, element "j" to get the subspectra weight */
            double x_diff = pj->x - ws->sample->ranges[tpd->i_range].x; /* Amount to change thickness of this and all subsequent layers */
            DEBUGMSG("Modifying ranges from %zu to %zu by %g tfu.", tpd->i_range, ws->sample->n_ranges, x_diff/C_TFU);
            for(size_t i_range = tpd->i_range; i_range < ws->sample->n_ranges; i_range++) {
                sample_rough->ranges[i_range].x += x_diff;
            }
        }
#ifdef DEBUG
        jabs_message(MSG_DEBUG, "Roughness subspectrum %zu / %zu", i_iter, iter_total);
        sample_print(sample_rough, FALSE, MSG_DEBUG);
#endif
        ws->fluence = p * fluence_original;
        ion_set_angle(&ws->ion, 0.0, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        sample_thickness_recalculate(sample_rough);
        status = simulate(&ws->ion, depth_seek(ws->sample, 0.0), ws, sample_rough);
        if(status < 0) {
            break;
        }
    }
    for(size_t i = 0; i < n_rl; i++) {
        thickness_probability_table_free(tpds[i]);
    }
    sample_free(sample_rough);
    free(tpds);
    ws->fluence = fluence_original;
    return status;
}

int simulate_with_ds(sim_workspace *ws) {
    if(!ws) {
        jabs_message(MSG_ERROR, "Congratulations, you've found a bug in %s:%i.\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    const double fluence = ws->fluence;
    if(simulate_with_roughness(ws) < 0) {
        return EXIT_FAILURE;
    }
    if(!ws->params->ds) {
        sim_workspace_calculate_sum_spectra(ws);
        return EXIT_SUCCESS;
    }
    ion_set_angle(&ws->ion, 0.0, 0.0);
    ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
    ion ion1 = ws->ion;
    ion ion2 = ion1;
    depth d_before = depth_seek(ws->sample, 0.0);
    int ds_steps_polar = ws->params->ds_steps_polar;
    int ds_steps_azi = ws->params->ds_steps_azi;
    sim_calc_params_defaults_fast(ws->params); /* This makes DS faster. Changes to ws->params are not reverted, but they don't affect original sim settings */
    ws->params->incident_stop_params.min *= 3.0;
    sim_calc_params_update(ws->params);
    const jibal_isotope *incident = ws->sim->beam_isotope;
    int last = FALSE;
    jabs_message(MSG_VERBOSE, "Dual scattering simulation starts.\n\n");
    double emin = ws->emin;

    while(1) {
        if(last) {
            break;
        }
        double E_front = ion1.E;
        if(E_front <= emin)
            break;
        double E_step = ws->params->ds_incident_stop_step_factor * stop_step_calc(&ws->params->incident_stop_params, &ion1);
        if(E_front - E_step <= emin) { /* Fix last step to be "just enough" */
            last = TRUE;
            E_step = E_front - emin;
            if(E_step < emin * 0.2) { /* The step would be really small, let's not take it */
                break;
            }
        }
        depth d_after = stop_step(&ws->stop, &ws->stragg, &ion1, ws->sample, d_before, E_step);
        double thick_step = depth_diff(d_before, d_after);
        const depth d_halfdepth = {.x = (d_before.x + d_after.x) /
                                        2.0, .i = d_after.i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double E_back = ion1.E;
        const double E_mean = (E_front + E_back) / 2.0;

        jabs_message(MSG_VERBOSE, "\rDS depth from %9.3lf to %9.3lf tfu, E from %6.1lf to %6.1lf keV.", d_before.x / C_TFU, d_after.x / C_TFU, E_front / C_KEV,
                     E_back / C_KEV);
        int n_running = 0; /* How many reactions are still producing data, largest number of all simulations for this depth step. */
        for(int i_polar = 0; i_polar < ds_steps_polar; i_polar++) {
            const double ds_polar_min = 20.0 * C_DEG;
            const double ds_polar_max = 180.0 * C_DEG;
            double ds_polar_step = (ds_polar_max - ds_polar_min) / (ds_steps_polar * 1.0);
            double ds_polar = ds_polar_min + i_polar * ds_polar_step;
            double cs_sum = 0.0;
            for(size_t i = 0; i < ws->sample->n_isotopes; i++) {
                double c = get_conc(ws->sample, d_halfdepth, i);
                if(c < ABUNDANCE_THRESHOLD) {
                    continue;
                }
                const jibal_isotope *target = ws->sample->isotopes[i];
                if(incident->mass >= target->mass && ds_polar > asin(target->mass / incident->mass)) { /* Scattering not possible */
                    continue;
                }
                double cs = 0.0;
                for(int polar_substep = 0; polar_substep < DUAL_SCATTER_POLAR_SUBSTEPS; polar_substep++) {
                    double ds_polar_sub = ds_polar_step * (1.0 * (polar_substep - (DUAL_SCATTER_POLAR_SUBSTEPS / 2)) / (DUAL_SCATTER_POLAR_SUBSTEPS * 1.0)) + ds_polar;
                    cs += jibal_cross_section_rbs(incident, target, ds_polar_sub, E_mean, JIBAL_CS_ANDERSEN) * sin(ds_polar_sub) / (DUAL_SCATTER_POLAR_SUBSTEPS * 1.0);
                }
                cs_sum += c * cs;

                double fluence_tot = cs_sum * thick_step * ion1.inverse_cosine_theta * (2.0 * C_PI) * ds_polar_step; /* TODO: check calculation after moving from p_sr to fluence!*/
                double fluence_azi = fluence_tot / (1.0 * (ds_steps_azi));
                for(int i_azi = 0; i_azi < ds_steps_azi; i_azi++) {
                    ion2 = ion1;
                    double K = jibal_kin_rbs(incident->mass, target->mass, ds_polar, '+');
                    ion2.E *= K;
                    ion2.S *= pow2(K);
                    double ds_azi = C_2PI * (1.0 * i_azi) / (ds_steps_azi * 1.0);
                    ion_rotate(&ion2, ds_polar, ds_azi); /* Dual scattering: first scattering to some angle (scattering angle: ds_polar). Note that this does not follow SimNRA conventions. */
                    ws->fluence = fluence_azi * fluence;
                    double scatangle = scattering_angle(&ion2, ws->sim->sample_theta, ws->sim->sample_phi, ws->det);
                    if(d_before.x == 0.0) {
                        DEBUGMSG("DS polar %.3lf, azi %.3lf, scatter %.3lf", ds_polar/C_DEG, ds_azi/C_DEG, scatangle/C_DEG);
                    }
                    if(scatangle > 19.99999 * C_DEG) {
                        int n_ok = simulate(&ion2, d_halfdepth, ws, ws->sample);
                        if(n_ok < 0) {
                            return EXIT_FAILURE;
                        }
                        if(n_ok > n_running) {
                            n_running = n_ok;
                        }
                    }
                }
            }
        }
        if(ws->sample->ranges[ws->sample->n_ranges - 1].x - d_after.x < 0.01 * C_TFU)
            break;
        d_before = d_after;
        jabs_message(MSG_VERBOSE, " %3i reactions ok.", n_running);
        if(n_running == 0) {
            break;
        }
    }
    jabs_message(MSG_VERBOSE, "\nDual scattering simulation complete.\n");
    sim_workspace_calculate_sum_spectra(ws);
    return EXIT_SUCCESS;
}
