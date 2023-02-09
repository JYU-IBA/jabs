/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <assert.h>
#include <gsl/gsl_math.h>
#include <jibal_gsto.h>
#include <jibal_stragg.h>
#include "jabs_debug.h"
#include "stop.h"
#include "defaults.h"

depth stop_next_crossing(const ion *incident, const sample *sample, const depth *d_from) {
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
        abort();
    }
    return d;
}


depth stop_step(const jabs_stop *stop, const jabs_stop *stragg, ion *incident, const sample *sample, depth depth_before, double step) {
    double k1, k2, k3, k4, stopping, dE, E;
    struct depth depth_next = stop_next_crossing(incident, sample, &depth_before);
    double h_max_perp = depth_next.x - depth_before.x;
#ifdef DEBUG_STOP_STEP
    if(depth_next.i != depth_before.i) {
        fprintf(stderr, "step crossing depth range from %zu to %zu at depth %lf tfu. E = %.3lf keV, Inverse cosine %lf\n", depth_before.i, depth_next.i, depth_before.x/C_TFU, incident->E/C_KEV, incident->inverse_cosine_theta);
    } else {
        fprintf(stderr, "step depth_before %g tfu (i=%zu) distance to next crossing %g tfu.\n", depth_before.x/C_TFU, depth_before.i, h_max_perp/C_TFU);
    }
    assert(fabs(h_max_perp) > 0.01 * C_TFU);
#endif
    depth_before.i = depth_next.i; /* depth_before may have the right depth (.x), but the index can be old. We need to cross the depth range somewhere, and it happens to be here. All calculations after this take place inside the same depth_before range (index depth_before.i). */
    /* k1...k4 are slopes of energy loss (stopping) at various x (depth) and E. Note convention: positive values, i.e. -dE/dx! */
    E = incident->E;
    if(incident->Z == 0) {
        k1 = 0.0;
    } else {
        k1 = stop_sample(stop, incident, sample, depth_before, E);
    }
    double h_max = h_max_perp * incident->inverse_cosine_theta; /*  we can take bigger steps since we are going sideways. Note that inverse_cosine_theta can be negative and in this case h_max should also be negative so h_max is always positive! */
    double h = (step / k1); /* (energy) step should always be positive, as well as k1, so depth step h (not perpendicular, but "real" depth) is always positive  */
    if(k1 < STOP_STEP_MINIMUM_STOPPING) { /* We are (almost) dividing by zero if there is no or very little stopping. Assume stopping is zero and make a jump. */
        DEBUGVERBOSEMSG("step returns no progress, because k1 = %g eV/tfu (x = %.3lf tfu, E = %.3lg keV)", k1/C_EV_TFU, depth_before.x/C_TFU, E/C_KEV);
        h = STOP_STEP_DEPTH_FALLBACK;
        if(h > h_max_perp) {
            return depth_next;
        }
        depth_before.x += h;
        return depth_before;
    }
    assert(h_max >= 0.0);
    assert(h > 0.0);
    h = GSL_MAX_DBL(h, STOP_STEP_ABSOLUTE_MINIMUM_STEP);
    depth halfdepth, fulldepth;
    halfdepth.i = depth_before.i;
    if(h >= h_max) { /* Depth step would take us beyond a depth range. We stop exactly on the boundary */
        h = h_max;
        halfdepth.x = depth_before.x + h_max_perp / 2.0;
        fulldepth = depth_next;
        if(h < 0.001 * C_TFU) {
            return fulldepth;
        }
    } else {
        double h_perp = h * incident->cosine_theta; /* x + h_perp is the actual perpendicular depth */
        halfdepth.x = depth_before.x + h_perp / 2.0;
        fulldepth.i = depth_before.i;
        fulldepth.x = depth_before.x + h_perp;
    }


    if(stop->rk4) {
        k2 = stop_sample(stop, incident, sample, halfdepth, E - (h / 2.0) * k1);
        k3 = stop_sample(stop, incident, sample, halfdepth, E - (h / 2.0) * k2);
        k4 = stop_sample(stop, incident, sample, fulldepth, E - h * k3);
        stopping = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    } else {
        stopping = k1;
    }
    assert(stopping > 0.0);
    dE = -1.0 * h * stopping; /* Energy change in thickness "h". It is always negative! */
#ifndef NO_STATISTICAL_STRAGGLING
    double s_ratio = stop_sample(stop, incident, sample, fulldepth, E + dE) / k1; /* Ratio of stopping for non-statistical broadening. TODO: at x? */
    incident->S *= pow2(s_ratio);
#endif
    incident->S += h * stop_sample(stragg, incident, sample, halfdepth, E + (0.5 * dE)); /* Straggling, calculate at mid-energy */
    incident->E += dE;
    return fulldepth; /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}

double stop_sample(const jabs_stop *stop, const ion *incident, const sample *sample, const depth depth, double E) {
    const double em = E * incident->mass_inverse;
    int Z2_old = 0; /* Not valid */
    double out = 0.0;
    int lo, hi = 0;
    const gsto_file_t *file;
    const double *data;
    double unit_factor;
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        double c;
        const jibal_isotope *target = sample->isotopes[i_isotope];
        if(sample->no_conc_gradients) {
            c = *sample_conc_bin(sample, depth.i, i_isotope);
        } else {
            c = get_conc(sample, depth, i_isotope);
        }
        if(c < CONC_TOLERANCE)
            continue;
        int Z2 = target->Z;
        if(Z2 != Z2_old) { /* This saves some computing time, but for electronic stopping we could just sum all isotopic concentrations (and note that those stay constant while we are doing stopping calculations) */
            if(stop->type == GSTO_STO_TOT) {
                file = incident->ion_gsto->gsto_data[Z2].stopfile;
                data = incident->ion_gsto->gsto_data[Z2].stopdata;
                if(file->stounit == GSTO_STO_UNIT_EV15CM2) {
                    unit_factor = C_EV_TFU;
                } else {
                    unit_factor = 1.0;
                }
            } else if(stop->type == GSTO_STO_STRAGG) {
                file = incident->ion_gsto->gsto_data[Z2].straggfile;
                data = incident->ion_gsto->gsto_data[Z2].straggdata;
                if(file->straggunit == GSTO_STRAGG_UNIT_BOHR) {
                    unit_factor = jibal_stragg_bohr(incident->Z, Z2);
                } else {
                    unit_factor = 1.0;
                }
            }
            lo = jibal_gsto_em_to_index(file, em);
            hi = lo + 1;
        }
        double S = unit_factor * jibal_linear_interpolation(file->em[lo], file->em[hi], data[lo], data[hi], em); /* TODO: check boundaries (low energy should be ok, but high maybe not) */

        if(stop->type == GSTO_STO_TOT) {
            S += ion_nuclear_stop(incident, target, stop->nuclear_stopping_accurate);
        }
        out += c * S;
    }
    switch(stop->type) {
        case GSTO_STO_ELE:
        case GSTO_STO_TOT:
            return out * sample->ranges[depth.i].bragg;
        case GSTO_STO_STRAGG:
            return out * sample->ranges[depth.i].stragg;
        default:
            return out;
    }
}

double stop_sample_old(const jabs_stop *stop, const ion *incident, const sample *sample, const depth depth, double E) {
    const double em = E * incident->mass_inverse;
    double S1 = 0.0;
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        double c;
        const jibal_isotope *target = sample->isotopes[i_isotope];
        if(sample->no_conc_gradients) {
            c = *sample_conc_bin(sample, depth.i, i_isotope);
        } else {
            c = get_conc(sample, depth, i_isotope);
        }
        if(c < CONC_TOLERANCE)
            continue;
        if(stop->type == GSTO_STO_TOT) {
            S1 += c * (
                    jibal_gsto_get_em(stop->gsto, GSTO_STO_ELE, incident->Z, target->Z, em)
                    #ifdef NUCLEAR_STOPPING_FROM_JIBAL
                    + jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, target->Z, target->mass)
                    #else
                    + ion_nuclear_stop(incident, target, stop->nuclear_stopping_accurate)
#endif
            );
        } else {
            S1 += c * (jibal_gsto_get_em(stop->gsto, stop->type, incident->Z, target->Z, em));
        }
    }
    switch(stop->type) {
        case GSTO_STO_ELE:
        case GSTO_STO_TOT:
            return S1 * sample->ranges[depth.i].bragg;
        case GSTO_STO_STRAGG:
            return S1 * sample->ranges[depth.i].stragg;
        default:
            return S1;
    }
}

double stop_step_calc(const jabs_stop_step_params *params, const ion *ion) {
    if(params->step > 0.0) {
        return params->step;
    }
    double step = params->sigmas * sqrt(ion->S);
    step = GSL_MIN_DBL(step, params->max);
    step = GSL_MAX_DBL(step, params->min);
    return step;
}

int stop_sample_exit(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params_exiting, ion *p, const depth depth_start, const sample *sample) {
    depth d = depth_start;
    while(1) { /* Exit from sample (hopefully) */
        if(p->inverse_cosine_theta > 0.0 && d.x >= (sample->thickness - DEPTH_TOLERANCE)) { /* Exit through back (transmission) */
            return 0;
        }
        if(p->inverse_cosine_theta < 0.0 && d.x <= DEPTH_TOLERANCE) { /* Exit (surface, front of sample) */
            return 0;
        }
        double E_step = stop_step_calc(params_exiting, p);
        depth d_after = stop_step(stop, stragg, p, sample, d, E_step);
        if(p->E < stop->emin) {
            DEBUGVERBOSEMSG("Energy %g keV below EMIN when surfacing from %.3lf tfu, break break.", p->E / C_KEV, d_after.x / C_TFU);
            return -1;
        }
        d = d_after;
    }
}
