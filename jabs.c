/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

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
#include <string.h>
#include <assert.h>
#include <math.h>
#include <jibal_units.h>
#include <jibal_kin.h>
#include <jibal_r33.h>

#include "geostragg.h"
#include "generic.h"
#include "rotate.h"
#include "roughness.h"
#include "jabs.h"
#include "defaults.h"
#include "message.h"
#include "win_compat.h"

extern inline double normal_pdf_std(double x);

double stop_sample(const sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, const depth depth, double E) {
    double em=E/incident->mass;
    double S1 = 0.0;
#ifdef NEUTRONS_EXIST
    if(incident->Z == 0)
        return 0.0; /* Return something neutron-specific (or whatever) */
#endif
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        double c;
        if(sample->no_conc_gradients) {
            c = *sample_conc_bin(sample, depth.i, i_isotope);
        } else {
            c = get_conc(sample, depth, i_isotope);
        }
        if(c < CONC_TOLERANCE)
            continue;
        if (type == GSTO_STO_TOT) {
            S1 += c * (
                    jibal_gsto_get_em(ws->gsto, GSTO_STO_ELE, incident->Z, sample->isotopes[i_isotope]->Z, em)
                    #ifdef NUCLEAR_STOPPING_FROM_JIBAL
                    +jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, sample->isotopes[i_isotope]->Z, sample->isotopes[i_isotope]->mass)
                    #else
                    + ion_nuclear_stop(incident, sample->isotopes[i_isotope], ws->isotopes, ws->params.nucl_stop_accurate)
                    #endif
                    );
        } else {
            S1 += c * (jibal_gsto_get_em(ws->gsto, type, incident->Z, sample->isotopes[i_isotope]->Z, em));
        }
    }
    switch(type) {
        case GSTO_STO_ELE:
        case GSTO_STO_TOT:
            return S1 * sample->ranges[depth.i].bragg;
        case GSTO_STO_STRAGG:
            return S1 * sample->ranges[depth.i].stragg;
        default:
            return S1;
    }
}

depth next_crossing(const ion *incident, const sample *sample, const depth *d_from) {
    depth d = *d_from;
    if(incident->inverse_cosine_theta > 0) { /* Going deeper */
        while(d.i < sample->n_ranges - 1 && d.x >= sample->ranges[d.i + 1].x) {
            d.i++;
        }
        if(d.i >= sample->n_ranges - 1) { /* There is probably a bug elsewhere in the code if you try to go this deep (deeper than last depth bin). */
            d.i = sample->n_ranges - 1;
            d.x = sample->ranges[d.i].x;
            fprintf(stderr, "Warning: probably too deep! This is a bug!\n");
        } else {
            d.x = sample->ranges[d.i + 1].x;
        }
    } else if(incident->inverse_cosine_theta < 0.0) { /* Going towards the surface */
        while(d.i > 0 && d.x <= sample->ranges[d.i].x) {
            d.i--;
        }
        d.x = sample->ranges[d.i].x;
    } else {
        jabs_message(MSG_ERROR, stderr, "WARNING: Inverse cosine is exactly zero. This is an issue!\n");
    }
    return d;
}

depth stop_step(const sim_workspace *ws, ion *incident, const sample *sample, depth depth, double step) {
    double k1, k2, k3, k4, stop, dE, E;
    struct depth depth_next = next_crossing(incident, sample, &depth);
    double h_max_perp = depth_next.x - depth.x;
#ifdef DEBUG_STOP_STEP
    if(depth_next.i != depth.i) {
        fprintf(stderr, "stop_step crossing depth range from %zu to %zu at depth %lf tfu. E = %.3lf keV, Inverse cosine %lf\n", depth.i, depth_next.i, depth.x/C_TFU, incident->E/C_KEV, incident->inverse_cosine_theta);
    } else {
        fprintf(stderr, "stop_step depth %g tfu (i=%zu) distance to next crossing %g tfu.\n", depth.x/C_TFU, depth.i, h_max_perp/C_TFU);
    }
    assert(fabs(h_max_perp) > 0.01 * C_TFU);
#endif
    depth.i = depth_next.i; /* depth may have the right depth (.x), but the index can be old. We need to cross the depth range somewhere, and it happens to be here. All calculations after this take place inside the same depth range (index depth.i). */
    /* k1...k4 are slopes of energy loss (stopping) at various x (depth) and E. Note convention: positive values, i.e. -dE/dx! */
    E = incident->E;
    k1 = stop_sample(ws, incident, sample, ws->stopping_type, depth, E);
    if(k1 < 0.001*C_EV_TFU) { /* Fail on positive values, zeroes (e.g. due to zero concentrations) and too small negative values */
#ifdef DEBUG_STOP_STEP
        fprintf(stderr, "stop_step returns no progress, because k1 = %g eV/tfu (x = %.3lf tfu, E = %.3lg keV)\n", k1/C_EV_TFU, depth.x/C_TFU, E/C_KEV);
#endif
        return depth;
    }
    double h_max = h_max_perp * incident->inverse_cosine_theta; /*  we can take bigger steps since we are going sideways. Note that inverse_cosine_theta can be negative and in this case h_max should also be negative so h_max is always positive! */
    assert(h_max >= 0.0);
    double h = (step / k1); /* (energy) step should always be positive, as well as k1, so depth step h (not perpendicular, but "real" depth) is always positive  */
    assert(h > 0.0);
    struct depth halfdepth;
    struct depth fulldepth;
    halfdepth.i = depth.i;
    if(h >= h_max) { /* Depth step would take us beyond a depth range. We stop exactly on the boundary */
        h = h_max;
        halfdepth.x = depth.x + h_max_perp/2.0;
        fulldepth = depth_next;
        if(h < 0.001*C_TFU) {
            return fulldepth;
        }
    } else {
        double h_perp = h*incident->cosine_theta; /* x + h_perp is the actual perpendicular depth */
        halfdepth.x = depth.x + h_perp/2.0;
        fulldepth.i = depth.i;
        fulldepth.x = depth.x + h_perp;
    }


    if(ws->params.rk4) {
        k2 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k1);
        k3 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k2);
        k4 = stop_sample(ws, incident, sample, ws->stopping_type, fulldepth, E - h * k3);
        stop = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    } else {
        stop = k1;
    }
    assert(stop > 0.0);
    dE =  -1.0* h * stop; /* Energy change in thickness "h". It is always negative! */
#ifndef STATISTICAL_STRAGGLING
    double s_ratio = stop_sample(ws, incident, sample, ws->stopping_type, fulldepth, E + dE) / k1; /* Ratio of stopping for non-statistical broadening. TODO: at x? */
#ifdef DEBUG
    //if((s_ratio)*(s_ratio) < 0.9 || (s_ratio)*(s_ratio) > 1.1) { /* Non-statistical broadening. */
    //   fprintf(stderr, "YIKES, s_ratio = %g, sq= %g\n", s_ratio, (s_ratio)*(s_ratio));
    //}
#endif
#endif
    incident->S *= pow2(s_ratio);
    incident->S += h * stop_sample(ws, incident, sample, GSTO_STO_STRAGG, halfdepth, E + (0.5 * dE)); /* Straggling, calculate at mid-energy */
    incident->E += dE;
    return fulldepth; /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}


double normal_pdf(double x, double mean, double sigma) {
    static const double inv_sqrt_2pi = 0.398942280401432703;
    double a = (x - mean) / sigma;
    return (inv_sqrt_2pi / sigma) * exp(-0.5 * a * a);
}

double cross_section_straggling(const sim_reaction *sim_r, int n_steps, double E, double S) {
    static const double sigmas = 2.0; /* TODO: what is enough or too much? */
    const double std_dev = sqrt(S);
    const int half_n = n_steps/2;
    const double w = sigmas/(half_n);
    double cs_sum = 0.0;
    double prob_sum = 0.0;
    //static const double inv_sqrt_2pi = 0.398942280401432703;
    for(int i = 0; i < n_steps; i++) {
        double x = w*(i-half_n);
        double E_stragg = E + x * std_dev;
        //double prob = normal_pdf(x, 0.0, 1.0) * w; /* TODO: if this is always a normal distribution and n_steps doesn't change, this function call could be replaced by a lookup table */
        double prob = normal_pdf_std(x);
        prob_sum += prob;
        cs_sum += prob * sim_r->cross_section(sim_r, E_stragg);
    }
    cs_sum /= prob_sum;
#ifdef DEBUG_CS_WEIGHT
    double unweighted = sim_r->cross_section(sim_r, E);
    double diff = (cs_sum-unweighted)/unweighted;
    fprintf(stderr, "Got cs %.7lf mb/sr (unweighted by straggling %.7lf mb/sr) diff %.5lf%%. Sum of probs %lf%% (compensated for).\n", cs_sum/C_MB_SR, unweighted/C_MB_SR, 100.0*diff, prob_sum*100.0);
#endif
    return cs_sum;
}


double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back) {
   if(ws->params.mean_conc_and_energy) { /* This if-branch is slightly faster (maybe) and also serves as a testing branch, since it is a lot easier to understand... */
        const depth d_halfdepth = {.x = (d_before->x + d_after->x)/2.0, .i = d_after->i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double c = get_conc(sample, d_halfdepth, sim_r->i_isotope);
        if(c < ABUNDANCE_THRESHOLD)
            return 0.0;
        const double E_mean = (E_front + E_back) / 2.0;
        assert(sim_r->cross_section);
        double sigma = sim_r->cross_section(sim_r, E_mean);
        return sigma * c * sample->ranges[d_halfdepth.i].yield;
    } else {
        depth d;
        d.i = d_after->i;
        const double x_step = (d_after->x - d_before->x) * ws->params.cs_frac;
        const double E_step = (E_back - E_front) * ws->params.cs_frac;
        const double S_step = (S_back - S_front) * ws->params.cs_frac;
        double sum = 0.0;
        for(int i = 1; i <= ws->params.cs_n_steps; i++) { /* Compute cross section and concentration product in several "sub-steps" */
            d.x = d_before->x + x_step * i;
            double E = E_front + E_step * i;
#ifdef DEBUG_VERBOSE
            fprintf(stderr, "i=%i, E = %g keV, (E_front = %g keV, E_back = %g keV)\n", i, E/C_KEV, E_front/C_KEV, E_back/C_KEV);
#endif
            double c = get_conc(sample, d, sim_r->i_isotope);
            double sigma;
            if(ws->params.cs_n_stragg_steps > 1) { /* Further weighted with straggling */
                double S = S_front + S_step * i;
                sigma = cross_section_straggling(sim_r, ws->params.cs_n_stragg_steps, E, S);
            } else {
                sigma = sim_r->cross_section(sim_r, E);
            }
            sum += sigma * c;
        }
       return sample->ranges[d.i].yield * sum/(ws->params.cs_n_steps*1.0);
   }
    return 0.0;
}

void post_scatter_exit(ion *p, const depth depth_start, const sim_workspace *ws, const sample *sample) {
    depth d = depth_start;
    while(1) { /* Exit from sample (hopefully) */
#ifdef DEBUG_REACTION
        fprintf(stderr, "  Exiting... depth = %g tfu (i = %zu)\n", d.x, d.i);
#endif
        if(d.x <= DEPTH_TOLERANCE) {
            break;
        }
        depth d_after = stop_step(ws, p, sample, d, ws->params.stop_step_exiting == 0.0?ws->params.stop_step_fudge_factor*(p->E*0.07+sqrt(p->S)+10.0*C_KEV):ws->params.stop_step_exiting); /* TODO: 7% of energy plus straggling plus 10 keV is a weird rule. Automatic stop size should be based more on required accuracy in stopping. */
        if(p->E < ws->sim->emin) {
#ifdef DEBUG_REACTION
            fprintf(stderr,
                            "  Reaction %lu with %s: Energy below EMIN when surfacing from %.3lf tfu, break break.\n",
                            i, r->r->target->name, d_after.x / C_TFU);
#endif
            return;
        }
        assert(d_after.x <= d.x /*|| (d_exit.x == d.x && d_exit.i != d.i)*/); /* Going towards the surface */
        d = d_after;
    }
}

void foil_traverse(ion *p, const sample *foil, sim_workspace *ws) {
    if(!foil)
        return;
    depth d = {.i = 0, .x = 0.0};
    ion ion_foil = *p;
    ion_set_angle(&ion_foil, 0.0, 0.0); /* Foils are not tilted. We use a temporary copy of "p" to do this step. */
    size_t i = 0;
    while(1) {
        i++;
        if(foil->ranges[foil->n_ranges-1].x - d.x < DEPTH_TOLERANCE) {
            break;
        }
        depth d_after  = stop_step(ws, &ion_foil, foil, d, ws->params.stop_step_exiting == 0.0?ws->params.stop_step_fudge_factor*(p->E*0.05+sqrt(p->S)+1.0*C_KEV):ws->params.stop_step_exiting);
        if(ion_foil.E < ws->sim->emin) {
            break;
        }
        d = d_after;
    }
    p->E = ion_foil.E;
    p->S = ion_foil.S;
}

double stop_step_calculate(const sim_workspace *ws, const ion *ion) { /* Calculate stop step to take */
    if(ws->params.stop_step_incident > 0) {
        return ws->params.stop_step_incident;
    }
    //double E_step = ws->params.stop_step_incident == 0.0?ws->params.stop_step_fudge_factor*sqrt(detector_resolution(ws->det, ion1.isotope, ion1.E)+ion1.S):ws->params.stop_step_incident;
    double broad = sqrt(ion->S) + ws->params.stop_step_add;
    if(broad < ws->params.stop_step_min) {
        return ws->params.stop_step_min;
    }
    return ws->params.stop_step_fudge_factor * broad; /* Fudge factor also affects the minimum stop step */
}

int simulate(const ion *incident, const depth depth_start, sim_workspace *ws, const sample *sample) { /* Ion is expected to be in the sample coordinate system at starting depth. Also note that sample may be slightly different (e.g. due to roughness) to ws->sim->sample */
    assert(sample->n_ranges);
    int warnings = 0;
    double thickness = sample->ranges[sample->n_ranges-1].x;
    size_t i_depth;
    ion ion1 = *incident; /* Shallow copy of the incident ion */
    geostragg_vars g = geostragg_vars_calculate(ws, incident);
    if(g.theta_product < 90.0*C_DEG) {
        jabs_message(MSG_ERROR, stderr, "Transmission geometry not supported, reaction product will not exit sample (angles in sample %g deg, %g deg).\n", g.theta_product/C_DEG, g.phi_product/C_DEG);
        return EXIT_FAILURE;
    }
#ifdef DEBUG
    fprintf(stderr, "Simulate from depth %g tfu (index %zu), detector theta = %g deg, calculated theta = %g deg. %zu reactions.\n", depth_start.x/C_TFU, depth_start.i, ws->det->theta/C_DEG, g.scatter_theta/C_DEG, ws->n_reactions);
    fprintf(stderr, "Ion energy at start %g keV, straggling %g keV FWHM.\n", incident->E/C_KEV, C_FWHM * sqrt(incident->S) / C_KEV);
#endif
    depth d_before = depth_start;
    for(size_t i = 0; i < ws->n_reactions; i++) {
#ifdef DEBUG
        fprintf(stderr, "Initializing reaction %zu\n", i);
#endif
        sim_reaction *r = &ws->reactions[i];
        if(!r->r)
            continue;
        r->last_brick = 0;
        r->stop = FALSE;
        r->theta = g.scatter_theta;
        ion_set_angle(&r->p, g.theta_product, g.phi_product);
        sim_reaction_recalculate_internal_variables(r);
        if(r->stop) {
            r->max_depth = 0.0;
            continue;
        }
        sim_reaction_product_energy_and_straggling(r, &ion1);
        r->max_depth = sample_isotope_max_depth(sample, r->i_isotope);
        sim_reaction_reset_bricks(r);
        brick *b = &r->bricks[0];
        b->E_0 = ion1.E;
        b->S_0 = ion1.S;
        b->E_r = r->p.E;
        b->S_r = r->p.S;
        b->d = d_before;
        b->E_0 = ion1.E;
        b->Q = 0.0;
        if(r->i_isotope >= sample->n_isotopes) { /* No target isotope for reaction. */
            r->stop = TRUE;
        }
        post_scatter_exit(&r->p, b->d, ws, sample); /* Calculates the exit energy if calculation doesn't start from the surface */
        foil_traverse(&r->p, ws->det->foil, ws);
        b->E = r->p.E;
        b->S = r->p.S;
        if(ws->params.geostragg) {
            b->S_geo_x = geostragg(ws, sample, r, &g.x, d_before, ion1.E);
            b->S_geo_y = geostragg(ws, sample, r, &g.y, d_before, ion1.E);
        }
#ifdef DEBUG
        fprintf(stderr, "Simulation reaction %zu: %s %s. Max depth %g tfu. i_isotope=%zu, stop = %i.\nCross section at %g keV is %g mb/sr, exit %g keV\n",
                i, reaction_name(r->r), r->r->target->name, r->max_depth / C_TFU, r->i_isotope, r->stop, ion1.E/C_KEV, r->cross_section(r, ion1.E)/C_MB_SR, r->p.E/C_KEV);
        ion_print(stderr, &r->p);
#endif
    }
    i_depth=1;
    if(fabs(ion1.cosine_theta) < 1e-6) {
#ifdef DEBUG
        fprintf(stderr, "Ion was going sideways in the sample, we nudged it a little.\n");
#endif
        ion1.theta += 0.01 * C_DEG;
        ion1.cosine_theta = cos(ion1.theta);
    }
    while(1) {
        if(warnings > SIMULATE_WARNING_LIMIT) {
            fprintf(stderr, "Warning limit reached. Won't calculate anything.\n");
            break;
        }
        if(d_before.x >= thickness) /* We're in too deep. */
            break;
        if(ion1.inverse_cosine_theta < 0.0 && d_before.x < 0.001*C_TFU) /* We're coming out of the sample and about to exit */
            break;
        if (ion1.E < ws->sim->emin) {
#ifdef DEBUG
            fprintf(stderr, "Break due to low energy (%.3lf keV < %.3lf keV), x = %.3lf, i_range = %lu.\n", ion1.E/C_KEV, ws->sim->emin/C_KEV, d_before.x/C_TFU, d_before.i);
#endif
            break;
        }
        const double E_front = ion1.E;
        const double S_front = ion1.S;
        double E_step = stop_step_calculate(ws, &ion1);
        depth d_after = stop_step(ws, &ion1, sample, d_before, E_step);
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "After:  %g tfu in range %zu\n", d_after.x/C_TFU, d_after.i);
#endif
        const double d_diff = depth_diff(d_before, d_after);
        /* DEPTH BIN [x, x+d_diff) */
        const double E_back = ion1.E;
        const double S_back = ion1.S;
        if(fabs(d_diff) < 0.001*C_TFU && E_front-E_back < 0.001*C_KEV) {
            jabs_message(MSG_WARNING, stderr, "Warning: no or very little progress was made (E step (goal) %g keV, E from %g keV to E = %g keV, depth = %g tfu, d_diff = %g tfu), check stopping or step size.\n", E_step/C_KEV, E_front/C_KEV , E_back/C_KEV, d_before.x/C_TFU, d_diff/C_TFU);
            //sample_print(stderr, sample, FALSE);
            //d_before.x += incident->inverse_cosine_theta*0.0001*C_TFU;
            d_before = d_after;
            d_before.x += 0.002*C_TFU;
            warnings++;
            continue;
        }

#ifdef DEBUG_VERBOSE
        double E_diff = E_front-E_back;
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3lf keV to  %8.3lf keV (diff %6.4lf keV)\n", x/C_TFU, (x+h)/C_TFU, E_front/C_KEV, ws->ion.E/C_KEV, E_diff/C_KEV);
#endif
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(ion1.S) / C_KEV);
#endif
        int all_stop = 1; /* Do we have any reactions left (or are we too deep in the sample) */
        for (size_t i = 0; i < ws->n_reactions; i++) {
            sim_reaction *r = &ws->reactions[i];
            if(r->stop)
                continue;
            if(i_depth >= r->n_bricks) {
                fprintf(stderr, "Too many bricks (%zu max). Data partial.\n", r->n_bricks);
                r->stop = TRUE;
                continue;
            }
            brick *b = &r->bricks[i_depth];
            if(!ws->params.ds && d_before.x >= r->max_depth) { /* Reactions stop when we are too deep in the sample, unless, of course, if DS is enabled. TODO: check optimizations for DS */
#ifdef DEBUG
                fprintf(stderr, "Reaction %lu with %s stops, because maximum depth is reached at x = %.3lf tfu.\n",
                        i, r->r->target->name, d_before.x / C_TFU); /* TODO: give reactions a name */
#endif
                b->Q = 0.0;
                r->stop = TRUE;
                continue;
            }
            all_stop = 0;
            b->d = d_after;
            b->E_0 = ion1.E; /* Sort of energy just before the reaction. */
            b->S_0 = ion1.S;
            assert(r->p.E > 0.0);

#ifdef DEBUG_REACTION
            fprintf(stderr, "Reaction %s (%zu): %s\n", reaction_name(r->r), i, r->r->target->name);
#endif
            if(ws->params.geostragg) {
                b->S_geo_x = geostragg(ws, sample, r, &g.x, d_after, ion1.E);
                b->S_geo_y = geostragg(ws, sample, r, &g.y, d_after, ion1.E);
            }
            sim_reaction_product_energy_and_straggling(r, &ion1);
            b->E_r = r->p.E;
            b->S_r = r->p.S;
            post_scatter_exit(&r->p, d_after, ws, sample);
            foil_traverse(&r->p, ws->det->foil, ws);
            b->E = r->p.E;
            b->S = r->p.S;
            #ifdef DEBUG
            fprintf(stderr, "Reaction %2zu depth from %8.3lf tfu to %8.3lf tfu, E = %8.3lf keV, Straggling: eloss %7.3lf keV, geo %7.3lf keV\n", i, d_before.x/C_TFU, d_after.x/C_TFU, b->E/C_KEV, sqrt(b->S)/C_KEV, sqrt(b->S_geo_x+b->S_geo_y)/C_KEV);
            #endif
            if (r->p.E < ws->sim->emin) {
                r->stop = TRUE;
                b->Q = 0.0;
                continue;
            }
            double sigma_conc = cross_section_concentration_product(ws, sample, r, E_front, E_back, &d_before, &d_after, S_front, S_back); /* Product of concentration and sigma for isotope i_isotope target and this reaction. */
            b->thick = d_diff;
            if(sigma_conc > 0.0) {
                if(d_after.i == sample->n_ranges - 2) {
                    sigma_conc *= ws->sim->channeling_offset + ws->sim->channeling_slope * (E_front + E_back)/2.0;
                }
                b->Q = ion1.inverse_cosine_theta * sigma_conc * d_diff;
                b->sc = sigma_conc;
                assert(b->Q >= 0.0);
#ifdef DEBUG_VERBOSE
                fprintf(stderr, "    %s: type=%i, E_front = %.3lf, E_back = %.3lf, E_out = %.3lf (sigma*conc = %g mb/sr, Q = %g (thickness = %.4lf tfu)\n",
                        r->r->target->name, r->r->type, E_front/C_KEV, E_back/C_KEV, r->p.E/C_KEV, sigma_conc/C_MB_SR, b->Q, d_diff/C_TFU);
#endif
            } else {
                b->Q = 0.0;
                b->sc = 0.0;
            }
            r->last_brick = i_depth;
        }
        d_before = d_after;
        i_depth++;

        if(all_stop) {
#ifdef DEBUG
            fprintf(stderr, "All reactions have ceased by depth %g.\n", d_before.x/C_TFU);
#endif
            break;
        }
    }
    convolute_bricks(ws);
    return EXIT_SUCCESS;
}

int assign_stopping_Z2(jibal_gsto *gsto, const simulation *sim, int Z2) { /* Assigns stopping and straggling (GSTO) for given Z2. Goes through all possible Z1s (beam and reaction products). */
    int fail = FALSE;
    int Z1 = sim->beam_isotope->Z;
#ifdef DEBUG
    fprintf(stderr, "Assigning stopping in Z2 = %i\n", Z2);
#endif
    if(assign_stopping_Z1_Z2(gsto, Z1, Z2)) {
        jabs_message(MSG_ERROR, stderr, "Can not assign stopping or straggling for beam.\n");
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
            jabs_message(MSG_ERROR, stderr, "Can not assign stopping or straggling for reaction product. Reaction: %s.\n", reaction_name(r));
            fail = TRUE;
        }
    }
    return fail;
}

int assign_stopping_Z1_Z2(jibal_gsto *gsto, int Z1, int Z2) {
    if(Z1 < 1 || Z2 < 1) {
        jabs_message(MSG_WARNING, stderr, "Assign of stopping for Z1 = %i in Z2 = %i is not possible.\n");
        return EXIT_SUCCESS; /* Skips, doesn't fail */
    }
    if(!jibal_gsto_auto_assign(gsto, Z1, Z2)) {
        jabs_message(MSG_ERROR, stderr, "Assign of stopping for Z1 = %i in Z2 = %i fails.\n", Z1, Z2);
        return EXIT_FAILURE;
    }
    if(!jibal_gsto_get_assigned_file(gsto, GSTO_STO_ELE, Z1, Z2)) {
        jabs_message(MSG_ERROR, stderr, "No electronic stopping assigned for Z1 = %i in Z2 = %i.\n", Z1, Z2);
        return EXIT_FAILURE;
    }
    if(!jibal_gsto_get_assigned_file(gsto, GSTO_STO_STRAGG, Z1, Z2)) {
        jabs_message(MSG_ERROR, stderr, "No energy loss straggling assigned for Z1 = %i in Z2 = %i.\n", Z1, Z2);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int assign_stopping(jibal_gsto *gsto, const simulation *sim) {
    /* TODO: simplify this by finding all possible Z1, Z2 combinations, considering target elements, beam and reactions before attempting to assign stopping/straggling (GSTO) */
    struct sample *sample = sim->sample;
    if(!sample) {
        jabs_message(MSG_ERROR, stderr, "Could not assign stopping, because sample is not set!\n");
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
    if(fail)
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}

int print_spectra(const char *filename, const sim_workspace *ws, const gsl_histogram *exp) {
    char sep = ' ';
    if(!ws)
        return EXIT_FAILURE;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return EXIT_FAILURE;
    if(filename) {
        size_t l = strlen(filename);
        if(l > 4 && strncmp(filename+l-4, ".csv", 4) == 0) { /* For CSV: print header line */
            sep = ','; /* and set the separator! */
            fprintf(f, "\"Channel\",\"Energy (keV)\",\"Simulated\"");
            if(exp) {
                fprintf(f, ",\"Experimental\"");
            }
            for(size_t j = 0; j < ws->n_reactions; j++) {
                const reaction *r = ws->reactions[j].r;
                fprintf(f, ",\"%s (%s)\"", r->target->name, reaction_name(r));
            }
            fprintf(f, "\n");
        }
    }
    for(size_t i = 0; i < ws->n_channels; i++) {
        fprintf(f,"%lu%c%.3lf%c", i, sep, detector_calibrated(ws->det, JIBAL_ANY_Z, i)/C_KEV, sep); /* Channel, energy. TODO: Z-specific calibration can have different energy (e.g. for a particular reaction). */
        if(ws->histo_sum->bin[i] == 0.0) {
            fprintf(f, "0"); /* Tidier output with a clean zero sum */
        } else {
            fprintf(f, "%e", ws->histo_sum->bin[i]);
        }
        if(exp) {
            if(i < exp->n) {
                fprintf(f, "%c%g", sep, exp->bin[i]);
            } else {
                fprintf(f, "%c0", sep);
            }
        }
        for (size_t j = 0; j < ws->n_reactions; j++) {
            if(i >= ws->reactions[j].histo->n || ws->reactions[j].histo->bin[i] == 0.0) {
                fprintf(f,"%c0", sep);
            } else {
                fprintf(f, "%c%e", sep, ws->reactions[j].histo->bin[i]);
            }
        }
        fprintf(f, "\n");
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

fit_params *fit_params_all(fit_data *fit) {
    simulation *sim = fit->sim;
    sample_model *sm = fit->sm;
    if(!sim)
        return NULL;
    size_t param_name_max_len = 256; /* Laziness. We use a fixed size temporary string. snprintf is used, so no overflows should occur, but very long names may be truncated. */
    char *param_name = malloc(sizeof(char) * param_name_max_len);
    fit_params *params = fit_params_new();
    fit_params_add_parameter(params, &sim->fluence, "fluence", "", 1.0);
    fit_params_add_parameter(params, &sim->channeling_offset, "channeling", "", 1.0);
    fit_params_add_parameter(params, &sim->channeling_slope, "channeling_slope", "1/keV", 1.0/C_KEV);
    fit_params_add_parameter(params, &sim->sample_theta, "alpha", "deg", C_DEG);
    fit_params_add_parameter(params, &sim->beam_E, "energy", "keV", C_KEV);
    for(size_t i_det = 0; i_det < sim->n_det; i_det++) {
        detector *det = sim_det(sim, i_det);
        char *det_name = NULL;
        if(asprintf(&det_name, "det%zu_", i_det + 1) < 0) {
            return NULL;
        }
        snprintf(param_name, param_name_max_len, "%ssolid", det_name);
        fit_params_add_parameter(params, &det->solid, param_name, "msr", C_MSR);

        for(int Z = JIBAL_ANY_Z; Z <= det->cal_Z_max; Z++) {
            calibration *c = detector_get_calibration(det, Z);
            if(Z != JIBAL_ANY_Z && c == det->calibration) /* No Z-specific calibration */
                continue;
            assert(c);
            size_t n = calibration_get_number_of_params(c);
            for(int i = CALIBRATION_PARAM_RESOLUTION; i < (int) n; i++) {
                char *calib_param_name = calibration_param_name(c->type, i);
                snprintf(param_name, param_name_max_len, "%scalib%s%s_%s",
                         det_name,
                         (Z == JIBAL_ANY_Z)?"":"_",
                         (Z == JIBAL_ANY_Z)?"":jibal_element_name(fit->jibal->elements, Z),
                         calib_param_name);
                free(calib_param_name);
                fit_params_add_parameter(params, calibration_get_param_ref(c, i), param_name, "keV", C_KEV);
            }
        }
        free(det_name);
    }
    if(sm) {
        for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) {
            sample_range *r = &(sm->ranges[i_range]);
            snprintf(param_name, param_name_max_len, "thick%zu", i_range + 1);
            fit_params_add_parameter(params, &(r->x), param_name, "tfu", C_TFU);
            if(r->rough.model != ROUGHNESS_NONE) {
                snprintf(param_name, param_name_max_len, "rough%zu", i_range + 1);
                fit_params_add_parameter(params, &(r->rough.x), param_name, "tfu", C_TFU);
            }
            for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
                if(*sample_model_conc_bin(sm, i_range, i_mat) < CONC_TOLERANCE) /* Don't add fit variables for negative, zero or very low concentrations. */
                    continue;
                snprintf(param_name, param_name_max_len, "conc%zu_%s", i_range + 1, sm->materials[i_mat]->name);
                fit_params_add_parameter(params, sample_model_conc_bin(sm, i_range, i_mat), param_name, "%", C_PERCENT);
            }
        }
    }
    free(param_name);
    return params;
}

int print_bricks(const char *filename, const sim_workspace *ws) {
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return EXIT_FAILURE;
    for(size_t i = 0; i < ws->n_reactions; i++) {
        const sim_reaction *r = &ws->reactions[i];
        fprintf(f, "#%s %s\n", reaction_name(r->r), r->r->target->name);
        fprintf(f, "#i  brick  depth    thick      E_0  S_0(el)      E_r  S_r(el)   E(det)    S(el)    S(geo) sigma*conc        Q\n");
        for(size_t j = 0; j <= r->last_brick; j++) {
            brick *b = &r->bricks[j];
            fprintf(f, "%2lu %4lu %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %10.1lf %8e\n",
                    i, j, b->d.x/C_TFU, b->thick/C_TFU,
                    b->E_0/C_KEV, C_FWHM * sqrt(b->S_0)/C_KEV,
                    b->E_r/C_KEV, C_FWHM * sqrt(b->S_r)/C_KEV,
                    b->E/C_KEV, C_FWHM * sqrt(b->S)/C_KEV,
                    C_FWHM * sqrt(b->S_geo_x+b->S_geo_y)/C_KEV, b->sc/C_MB_SR, b->Q * ws->fluence * ws->det->solid);
        }
        fprintf(f, "\n\n");
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}



int simulate_with_roughness(sim_workspace *ws) {
    int status = EXIT_SUCCESS;
    double p_sr = ws->sim->fluence;
    size_t n_rl = 0; /* Number of rough layers */
    for(size_t i = 0; i < ws->sample->n_ranges; i++) {
        sample_range *r = &(ws->sample->ranges[i]);
        r->rough.n *= ws->params.rough_layer_multiplier;
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
#ifdef DEBUG
    fprintf(stderr, "%zu rough layers\n", n_rl);
#endif
    if(!n_rl) {
        ion_set_angle(&ws->ion, 0.0*C_DEG, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        return simulate(&ws->ion, depth_seek(ws->sample, 0.0*C_TFU), ws, ws->sample);
    }
    struct sample *sample_rough = sample_copy(ws->sample);
    size_t *index = malloc(sizeof(size_t) * n_rl);
    size_t *modulos = malloc(sizeof(size_t) * n_rl);
    size_t j = 0;
    thick_prob_dist **tpd = malloc(sizeof(thick_prob_dist *) * n_rl);
    for(size_t i = 0; i < ws->sample->n_ranges; i++) {
        sample_range *r = &(ws->sample->ranges[i]);
        if(r->rough.model == ROUGHNESS_GAMMA) {
#ifdef DEBUG
            fprintf(stderr, "Range %zu is rough (gamma), amount %g tfu, n = %zu spectra\n", i, ws->sample->ranges[i].rough.x/C_TFU, ws->sample->ranges[i].rough.n);
#endif
            tpd[j] = thickness_probability_table_gen(r->x, r->rough.x, r->rough.n);
#ifdef DEBUG
            fprintf(stderr, "TPD for depth %zu (%.3lf tfu nominal), roughness %.3lf tfu:\n", i, ws->sample->ranges[i].x/C_TFU, ws->sample->ranges[i].rough.x/C_TFU);
            for(size_t i_tpd = 0; i_tpd < tpd[j]->n; i_tpd++) {
                fprintf(stderr, "%zu: %.3lf tfu %.3lf%%\n", i_tpd, tpd[j]->p[i_tpd].x/C_TFU, tpd[j]->p[i_tpd].prob/C_PERCENT);
            }
#endif
            index[j] = i;
            if(j)
                modulos[j] = modulos[j-1] *  tpd[j-1]->n;
            else
                modulos[j] = 1;
            j++;
        }
    }
    size_t iter_total = modulos[n_rl-1] * tpd[n_rl-1]->n;
    for(size_t i_iter = 0; i_iter < iter_total; i_iter++) {
#ifdef DEBUG
        fprintf(stderr, "Gamma roughness step %zu/%zu\n", i_iter+1, iter_total);
#endif
        double p = 1.0;
        for(size_t i_range = 0; i_range < ws->sample->n_ranges; i_range++) { /* Reset ranges for every iter */
            sample_rough->ranges[i_range].x = ws->sample->ranges[i_range].x;
        }
        for(size_t i = 0; i < n_rl; i++) {
            j = (i_iter / modulos[i]) % tpd[i]->n; /* "j"th roughness element */
            //fprintf(stderr, " %zu", j);

            size_t i_range = index[i];
            p *= tpd[i]->p[j].prob; /* Probability is multiplied by the "i"th roughness, element "j" */
            double x_diff = tpd[i]->p[j].x - ws->sample->ranges[i_range].x; /* Amount to change thickness of this and and all subsequent layers */
            for(; i_range < ws->sample->n_ranges; i_range++) {
                sample_rough->ranges[i_range].x += x_diff;
            }
#ifdef DEBUG
            fprintf(stderr, "Gamma roughness diff %g tfu (from %g tfu, index i_range=%zu), probability %.3lf%%)\n", x_diff/C_TFU, ws->sample->ranges[i_range].x/C_TFU, i_range, tpd[i]->p[j].prob*100.0);
            fprintf(stderr, "Gamma roughness, ranges (%zu):", ws->sample->n_ranges);
            for(i_range = 0; i_range < ws->sample->n_ranges; i_range++) {
                fprintf(stderr, ", %zu: %g tfu ", i_range, sample_rough->ranges[i_range].x/C_TFU);
            }
            fprintf(stderr, "\n");
#endif
        }
        ws->fluence = p * p_sr;
        ion_set_angle(&ws->ion, 0.0, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        status = simulate(&ws->ion, depth_seek(ws->sample, 0.0), ws, sample_rough);
        if(status != EXIT_SUCCESS)
            break;
    }
    for(size_t i = 0; i < n_rl; i++) {
        thickness_probability_table_free(tpd[i]);
    }
    sample_free(sample_rough);
    free(modulos);
    free(index);
    free(tpd);
    return status;
}

int simulate_with_ds(sim_workspace *ws) {
    if(!ws) {
        jabs_message(MSG_ERROR, stderr, "No workspace, no simulation.\n");
        return EXIT_FAILURE;
    }
    double fluence = ws->fluence;
    if(simulate_with_roughness(ws)) {
        return EXIT_FAILURE;
    }
    if(!ws->params.ds) {
        sim_workspace_calculate_sum_spectra(ws);
        return EXIT_SUCCESS;
    }
    ion_set_angle(&ws->ion, 0.0, 0.0);
    ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
    ion ion1 = ws->ion;
    ion ion2 = ion1;
    depth d_before = depth_seek(ws->sample, 0.0);
    sim_calc_params_fast(&ws->params, TRUE); /* This makes DS faster. Changes to ws->params are not reverted, but they don't affect original sim settings */
    jabs_message(MSG_ERROR, stderr, "\n");
    const jibal_isotope *incident = ws->sim->beam_isotope;
    while(1) {
        double E_front = ion1.E;
        if(E_front < ws->sim->emin)
            break;
        depth d_after = stop_step(ws, &ion1, ws->sample, d_before, ws->params.stop_step_fudge_factor*sqrt(detector_resolution(ws->det, ion1.isotope, ion1.E)+ion1.S)); /* TODO: step? */
        double thick_step = depth_diff(d_before, d_after);
        const depth d_halfdepth = {.x = (d_before.x + d_after.x)/2.0, .i = d_after.i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double E_back = ion1.E;
        const double E_mean = (E_front + E_back) / 2.0;

        jabs_message(MSG_VERBOSE, stderr, "DS depth from %9.3lf tfu to %9.3lf tfu, E from %6.1lf keV to %6.1lf keV. p*sr = %g\n", d_before.x/C_TFU, d_after.x/C_TFU, E_front/C_KEV, E_back/C_KEV, fluence);
        double p_sum = 0.0;
        for(int i_polar = 0; i_polar < ws->params.ds_steps_polar; i_polar++) {
            const double ds_polar_min = 20.0*C_DEG;
            const double ds_polar_max = 180.0*C_DEG;
            double ds_polar_step = (ds_polar_max-ds_polar_min)/(ws->params.ds_steps_polar*1.0);
            double ds_polar = ds_polar_min + i_polar * ds_polar_step;
            double cs_sum = 0.0;
            for(size_t i = 0; i < ws->sample->n_isotopes; i++) {
                double c = get_conc(ws->sample, d_halfdepth, i);
                if(c < ABUNDANCE_THRESHOLD)
                    continue;
                const jibal_isotope *target = ws->sample->isotopes[i];
                if(incident->mass >= target->mass && ds_polar > asin(target->mass / incident->mass)) { /* Scattering not possible */
                    continue;
                }
                double cs = 0.0;
                for(int polar_substep = 0; polar_substep < DUAL_SCATTER_POLAR_SUBSTEPS; polar_substep++) {
                    double ds_polar_sub = ds_polar_step*(1.0*(polar_substep-(DUAL_SCATTER_POLAR_SUBSTEPS/2))/(DUAL_SCATTER_POLAR_SUBSTEPS*1.0)) + ds_polar;
                    cs += jibal_cross_section_rbs(incident, target, ds_polar_sub, E_mean, JIBAL_CS_ANDERSEN) * sin(ds_polar_sub)/(DUAL_SCATTER_POLAR_SUBSTEPS*1.0);
                }
                cs_sum += c * cs;
            }
            double fluence_tot = cs_sum * thick_step * ion1.inverse_cosine_theta * (2.0 * C_PI) * ds_polar_step; /* TODO: check calculation after moving from p_sr to fluence!*/
            p_sum += fluence_tot;
            double fluence_azi = fluence_tot / (1.0 * (ws->params.ds_steps_azi));
            for(int i_azi = 0; i_azi < ws->params.ds_steps_azi; i_azi++) {
                ion2 = ion1;
                double ds_azi = C_2PI * (1.0 * i_azi) / (ws->params.ds_steps_azi * 1.0);
                ion_rotate(&ion2, ds_polar, ds_azi); /* Dual scattering: first scattering to some angle (scattering angle: ds_polar). Note that this does not follow SimNRA conventions. */
                ws->fluence = fluence_azi * fluence;
#ifdef DEBUG
                if(d_before.x == 0.0) {
                    fprintf(stderr, "DS polar %.3lf, azi %.3lf, scatter %.3lf\n", ds_polar/C_DEG, ds_azi/C_DEG, scattering_angle(&ion2, ws)/C_DEG);
                }
#endif
                if(scattering_angle(&ion2, ws) > 19.99999*C_DEG) {
                    if(simulate(&ion2, d_halfdepth, ws, ws->sample)) {
                        return EXIT_FAILURE;
                    }
                }
            }
        }
        fluence -= p_sum * fluence;
        if(ws->sample->ranges[ws->sample->n_ranges-1].x - d_after.x < 0.01*C_TFU)
            break;
        d_before = d_after;
    }
    sim_workspace_calculate_sum_spectra(ws);
    return EXIT_SUCCESS;
}

void fit_params_print(const fit_params *params, int active, const char *pattern) { /* Prints current values of all possible fit variables matching pattern. Pattern can be NULL too. */
    if(!params)
        return;
    if(params->n) {
        if(pattern) {
            jabs_message(MSG_INFO, stderr, "All possible fit variables matching pattern \"%s\":\n", pattern);
        } else {
            jabs_message(MSG_INFO, stderr, "All possible fit variables (use 'show fitvar <pattern>' to see variables matching pattern, wildcards are '*' and '?'):\n");
        }
    } else {
        jabs_message(MSG_INFO, stderr, "No fit variables.\n");
    }
    for(size_t i = 0; i < params->n; i++) {
        const fit_variable *var = &params->vars[i];
        if(active && !var->active)
            continue;
        if(pattern && !is_match(var->name, pattern))
            continue;
        jabs_message(MSG_INFO, stderr, "%s %24s = %g %s\n", var->active ? "X" : " ", var->name, *(var->value) / var->unit_factor, var->unit);
    }
}

void fit_params_print_final(const fit_params *params) { /* Prints final values of active fit variables. */
    if(!params)
        return;
    if(params->n_active) {
        jabs_message(MSG_INFO, stderr, "Final fit variables (%zu/%zu):\n", params->n_active, params->n);
    } else {
        jabs_message(MSG_INFO, stderr, "No fitted variables of total %zu.\n", params->n);
        return;
    }
    jabs_message(MSG_INFO, stderr, "  i |                 variable |  unit |       value |     error | rel % |  orig. value | multipl. | sigmas |\n");
    for(size_t i = 0; i < params->n; i++) {
        const fit_variable *var = &params->vars[i];
        if(!var->active)
            continue;
        //jabs_message(MSG_INFO, stderr, "%24s(%3s) = %12g +- %12g (%.1lf%%)\n", var->name, var->unit, var->value_final/var->unit_factor, var->err/var->unit_factor, var->err/var->value_final*100.0);
        jabs_message(MSG_INFO, stderr, "%3zu | %24s | %5s | %11.6g | %9.4g | %5.1lf | %12.6g | %8.4lf | %6.1lf |\n", var->i_v + 1, var->name, var->unit,
                     var->value_final/var->unit_factor,
                     var->err/var->unit_factor,
                     100.0 * var->err_rel,
                     var->value_orig/var->unit_factor,
                     var->value_final/var->value_orig,
                     var->sigmas
        );
    }
}

size_t fit_params_enable(fit_params *params, const char *s, int enable) {
    size_t n_match = 0;
    for(size_t i = 0; i < params->n; i++) {
        fit_variable *var = &params->vars[i];
        if(is_match(var->name, s)) {
            var->active = enable;
            n_match++;
        }
    }
    return n_match;
}
