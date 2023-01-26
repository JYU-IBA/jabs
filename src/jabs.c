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
#include <string.h>
#include <assert.h>
#include <math.h>
#include <jibal_units.h>
#include <jibal_kin.h>
#include <jibal_r33.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "geostragg.h"
#include "generic.h"
#include "rotate.h"
#include "roughness.h"
#include "jabs.h"
#include "defaults.h"
#include "message.h"
#include "win_compat.h"

des_table *des_table_init(size_t n) {
    des_table *dt = malloc(sizeof(des_table));
    dt->t = NULL;
    dt->depth_interval_index = NULL;
    dt->n = 0;
    dt->n_ranges = 0;
    if(des_table_realloc(dt, n)) {
        free(dt);
        return NULL;
    }
    return dt;
}
int des_table_realloc(des_table *dt, size_t n) {
#ifdef DEBUG
    fprintf(stderr, "DES table %p realloc to %zu\n", (void *)dt, n);
#endif
    dt->t = realloc(dt->t, n * sizeof(des));
    if(dt->t) {
        dt->n = n;
        return 0;
    } else {
        dt->n = 0;
        return -1;
    }
}

void des_table_free(des_table *dt) {
    if(!dt) {
        return;
    }
    free(dt->t);
    free(dt->depth_interval_index);
    free(dt);
}

size_t des_table_size(const des_table *dt) {
    if(!dt) {
        return 0;
    }
    return dt->n;
}

des *des_table_element(const des_table *dt, size_t i) {
    if(!dt) {
        return NULL;
    }
    if(i > dt->n) {
        return NULL;
    }
    return &(dt->t[i]);
}

void des_table_rebuild_index(des_table *dt) {
    if(!dt->n)
        return;
    size_t n_ranges = dt->t[dt->n-1].d.i + 1; /* Number of ranges is one greater than the last index. Make sure the last element is nice and valid. */
    dt->n_ranges = n_ranges;
    if(dt->depth_interval_index) {
        free(dt->depth_interval_index);
    }
    dt->depth_interval_index = malloc(sizeof(size_t) * (n_ranges + 1));
    size_t i_range_old = 0, i_range_new = 0;
    /* NOTE: depth_start may have not been start of bin! This shouldn't matter if everything is ok...*/
    dt->depth_interval_index[0] = 0;
    for(size_t i = 1; i < dt->n; i++) {
        const des *des = &(dt->t[i]);
        if(des->d.i > i_range_old) {
            i_range_new = des->d.i;
            assert(i_range_new < n_ranges);
            for(size_t i_range = i_range_old + 1; i_range <= i_range_new; i_range++) { /* Handles indices stop_step skips over (no thickness between them) */
                dt->depth_interval_index[i_range] = i - 1;
            }
            i_range_old = i_range_new;
        }
    }
    dt->depth_interval_index[n_ranges] = dt->n - 1;
}

void des_table_print(FILE *f, const des_table *dt) {
    if(!dt) {
        return;
    }
    for(size_t i = 0; i < dt->n; i++) {
        const des *des = &(dt->t[i]);
        fprintf(f, "DES %4zu %3zu %10.3lf %10.3lf %8.3lf\n", i, des->d.i, des->d.x / C_TFU, des->E / C_KEV, sqrt(des->S) / C_KEV);
    }
    for(size_t i = 0; i < dt->n_ranges; i++) {
        fprintf(f, "DES INDEX %3zu [%3zu,%3zu]\n", i, dt->depth_interval_index[i],  dt->depth_interval_index[i+1]);
    }
}

depth des_table_find_depth(const des_table *dt, size_t *i_des, depth depth_prev, ion *incident) {
    size_t i;
    assert(dt);
    assert(dt->n > 0);
    double E = incident->E;
    for(i = *i_des; i < dt->n - 1; i++) {
        if(dt->t[i].d.i != dt->t[i + 1].d.i) { /* This means last point of this layer */
            E = dt->t[i].E;
            *i_des = i + 1; /* The +1 prevents stopping at this same layer boundary on the next call */
#ifdef DEBUG
            fprintf(stderr, "Layer boundary at %g tfu. Setting E = %g keV and i_des = %zu. index = %zu\n", dt->t[i].d.x / C_TFU, E / C_KEV, *i_des, dt->depth_interval_index[dt->t[i + 1].d.i]);
#endif
            break;
        }
        if(dt->t[i].E < E) {/* i is the index of the first element in dt->t that has energy below E. So i-1 should have energy above E. */
            *i_des = i;
            break;
        }
    }
    if(i == dt->n - 1) {
        if(E < dt->t[i].E) {
#ifdef DEBUG
            fprintf(stderr, "Energy %g keV is below last point in table (%g keV).\n", E / C_KEV, dt->t[i].E / C_KEV);
#endif
            E = dt->t[i].E;
        }
    }
    assert(i > 0);
    const des *des_low = &dt->t[i - 1]; /* Closer to surface, higher energy */
    const des *des_high = &dt->t[i]; /* Deeper, lower energy */
#ifdef DEBUG
    fprintf(stderr, "des_low = %g tfu (i = %zu), des_high = %g tfu (i = %zu), depth_prev = %g tfu (i = %zu)\n", des_low->d.x / C_TFU, des_low->d.i, des_high->d.x / C_TFU, des_high->d.i, depth_prev.x / C_TFU, depth_prev.i);
#endif

    double E_diff = E - des_low->E;
    double E_interval = des_high->E - des_low->E;
    double S_interval = des_high->S - des_low->S;
    double d_interval = depth_diff(des_low->d, des_high->d);
    double frac = (E_diff/E_interval); /* zero if close to low bin */
    assert(frac >= 0.0 && frac <= 1.0);
    depth d_out;
    if(depth_prev.i != des_high->d.i) { /* First point after crossing layer */
        d_out.i = des_high->d.i;
    } else {
        d_out.i = depth_prev.i;
    }
    d_out.x = des_low->d.x + frac * (d_interval); /* Linear interpolation */
#if 0
    incident->E = des_low->E + frac * E_interval;
#else
    incident->E = E;
#endif
    incident->S = des_low->S + frac * S_interval;
#ifdef DEBUG
    fprintf(stderr, "d_out = %g tfu (i = %zu), E = %g keV, S = %g keV\n", d_out.x / C_TFU, d_out.i, incident->E / C_KEV, sqrt(incident->S) * C_FWHM / C_KEV);
#endif
    return d_out;
}

void des_table_set_ion_depth(const des_table *dt, ion *ion, depth d) {
    size_t lo = dt->depth_interval_index[d.i];
    size_t hi = dt->depth_interval_index[d.i + 1];
    size_t mi;
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "Depth %g = tfu, index = %zu. Probably in DES table between %zu and %zu.\n", d.x / C_TFU, d.i, lo, hi);
#endif
    while (hi - lo > 1) {
        mi = (hi + lo) / 2;
        if (d.x >= dt->t[mi].d.x) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    double depth_diff = d.x - dt->t[lo].d.x; /* Positive */
    double depth_interval = dt->t[hi].d.x - dt->t[lo].d.x;
    double frac = depth_diff/depth_interval; /* zero at low bin, 1.0 at high bin */
    double frac_inv = 1.0 - frac;
    double E = dt->t[lo].E * frac_inv + dt->t[hi].E * frac;
    double S = dt->t[lo].S * frac_inv + dt->t[hi].S * frac;
    assert(depth_diff >= 0.0);
    ion->E = E;
    ion->S = S;
    assert(S >= 0.0);
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "lo = %zu, mi = %zu, hi = %zu, depth_diff = %g tfu, depth_interval = %g tfu, frac = %lf, E = %g keV, S = %g keV\n",
            lo, mi, hi, depth_diff / C_TFU, depth_interval / C_TFU, frac, E / C_KEV, sqrt(S) / C_KEV);
#endif
}

des_table *des_table_compute(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample) {
    ion ion = *incident;
    des_table *dt = des_table_init(DES_TABLE_INITIAL_ALLOC);
    size_t i = 0;
    depth d_before;
    depth d_after = depth_start;
    do {
        d_before = d_after;
        des *des = &dt->t[i];
        des->E = ion.E;
        des->S = ion.S;
        des->d = d_before;
        if(d_before.x >= sample->thickness) {
            i++; /* TODO: Safe? */
            break;
        }

        d_after = stop_step(ws, &ion, sample, d_before, stop_step_calc_incident(ws, &ion));
        if(i == dt->n) {
#ifdef DEBUG
            fprintf(stderr, "DES table reallocation, size %zu reached.\n", dt->n);
#endif
            des_table_realloc(dt, dt->n * 2);
        }
        i++;
    } while(ion.E > ws->emin);
    des_table_realloc(dt, i); /* Shrinks to size */
    des_table_rebuild_index(dt);
    return dt;
}


double stop_sample(const sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, const depth depth, double E) {
    const double em = E * incident->mass_inverse;
    double S1 = 0.0;
#ifdef NEUTRONS_EXIST
    if(incident->Z == 0)
        return 0.0; /* Return something neutron-specific (or whatever) */
#endif
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
        if(type == GSTO_STO_TOT) {
            S1 += c * (
                    jibal_gsto_get_em(ws->gsto, GSTO_STO_ELE, incident->Z, target->Z, em)
                    #ifdef NUCLEAR_STOPPING_FROM_JIBAL
                    + jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, target->Z, target->mass)
                    #else
                    + ion_nuclear_stop(incident, target, ws->params->nucl_stop_accurate)
#endif
            );
        } else {
            S1 += c * (jibal_gsto_get_em(ws->gsto, type, incident->Z, target->Z, em));
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
    if(incident->Z == 0) {
        k1 = 0.0;
    } else {
        k1 = stop_sample(ws, incident, sample, ws->stopping_type, depth, E);
    }
    double h_max = h_max_perp * incident->inverse_cosine_theta; /*  we can take bigger steps since we are going sideways. Note that inverse_cosine_theta can be negative and in this case h_max should also be negative so h_max is always positive! */
    double h = (step / k1); /* (energy) step should always be positive, as well as k1, so depth step h (not perpendicular, but "real" depth) is always positive  */
    if(k1 < STOP_STEP_MINIMUM_STOPPING) { /* We are (almost) dividing by zero if there is no or very little stopping. Assume stopping is zero and make a jump. */
#ifdef DEBUG_STOP_STEP
        fprintf(stderr, "stop_step returns no progress, because k1 = %g eV/tfu (x = %.3lf tfu, E = %.3lg keV)\n", k1/C_EV_TFU, depth.x/C_TFU, E/C_KEV);
#endif
        h = STOP_STEP_DEPTH_FALLBACK;
        if(h > h_max_perp) {
            return depth_next;
        }
        depth.x += h;
        return depth;
    }
    assert(h_max >= 0.0);
    assert(h > 0.0);
    struct depth halfdepth;
    struct depth fulldepth;
    halfdepth.i = depth.i;
    if(h >= h_max) { /* Depth step would take us beyond a depth range. We stop exactly on the boundary */
        h = h_max;
        halfdepth.x = depth.x + h_max_perp / 2.0;
        fulldepth = depth_next;
        if(h < 0.001 * C_TFU) {
            return fulldepth;
        }
    } else {
        double h_perp = h * incident->cosine_theta; /* x + h_perp is the actual perpendicular depth */
        halfdepth.x = depth.x + h_perp / 2.0;
        fulldepth.i = depth.i;
        fulldepth.x = depth.x + h_perp;
    }


    if(ws->params->rk4) {
        k2 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k1);
        k3 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k2);
        k4 = stop_sample(ws, incident, sample, ws->stopping_type, fulldepth, E - h * k3);
        stop = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    } else {
        stop = k1;
    }
    assert(stop > 0.0);
    dE = -1.0 * h * stop; /* Energy change in thickness "h". It is always negative! */
#ifndef NO_STATISTICAL_STRAGGLING
    double s_ratio = stop_sample(ws, incident, sample, ws->stopping_type, fulldepth, E + dE) / k1; /* Ratio of stopping for non-statistical broadening. TODO: at x? */
    incident->S *= pow2(s_ratio);
#endif
    incident->S += h * stop_sample(ws, incident, sample, GSTO_STO_STRAGG, halfdepth, E + (0.5 * dE)); /* Straggling, calculate at mid-energy */
    incident->E += dE;
    return fulldepth; /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}

double cross_section_straggling_fixed(const sim_reaction *sim_r, const prob_dist *pd, double E, double S) {
    const double std_dev = sqrt(S);
    double cs_sum = 0.0;
    for(size_t i = 0; i < pd->n; i++) {
        prob_point *pp = &(pd->points[i]);
        double E_stragg = E + pp->x * std_dev;
        double cs = pp->p * sim_r->cross_section(sim_r, E_stragg);
        cs_sum += cs;
#ifdef DEBUG_CS_STRAGG
        fprintf(stderr, "%zu %10g %10g %10g %10g\n", i, pp->x, E_stragg/C_KEV, pp->p, cs/C_MB_SR);
#endif
    }
#ifdef DEBUG_CS_STRAGG
    double unweighted = sim_r->cross_section(sim_r, E);
        double diff = (cs_sum-unweighted)/unweighted;
        fprintf(stderr, "Got cs %.7lf mb/sr (unweighted by straggling %.7lf mb/sr) diff %.7lf%%.\n", cs_sum/C_MB_SR, unweighted/C_MB_SR, 100.0*diff);
#endif
    return cs_sum;
}


double cross_section_straggling(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, const prob_dist *pd, double E, double S) {
#ifdef DEBUG_CS_STRAGG
    fprintf(stderr, "CS stragg, E = %g keV, S = %g keV (std. dev)\n", E/C_KEV, std_dev/C_KEV);
#endif
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
#ifdef DEBUG_CS_VERBOSE
    fprintf(stderr, "cs_stragg_function(), E = %g keV, S = %g keV FWHM, E_mean = %g keV, a = %g. Result %g mb/sr.\n", x/C_KEV, p->sigma*C_FWHM/C_KEV, p->E_mean/C_KEV, a, (0.398942280401432703/p->sigma) * result/C_MB_SR );
#endif
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
#ifdef DEBUG_CS_VERBOSE
    fprintf(stderr, "Integrated from %g keV to %g keV in %zu steps (limit %zu), got (%g +- %g) mb/sr\n", E_low/C_KEV, E_high/C_KEV, w->size, w->limit, result/C_MB_SR, error/C_MB_SR);
#endif
    return result;
}

double cs_function(double x, void * params) {
    struct cs_int_params *p = (struct cs_int_params *) params;
    p->d.x = p->d_before->x + p->stop_slope * (x - p->E_front); /* Depth assuming constant stopping inside brick */
    double c = get_conc(p->sample, p->d, p->sim_r->i_isotope);
    double sigma;
    double S = p->S_front + p->stragg_slope * (x - p->E_front);
#ifdef DEBUG_CS_VERBOSE
        fprintf(stderr, "S = %g keV, E = %g keV\n", sqrt(S)/C_KEV, x/C_KEV);
#endif
        sigma = cross_section_straggling(p->sim_r, p->w, p->stragg_int_accuracy, p->cs_stragg_pd, x, S);

#ifdef DEBUG
    fprintf(stderr, "Depth %g tfu, energy %g keV, sigma  %g mb/sr, c %g %%\n", p->d.x/C_TFU, x/C_KEV, sigma/C_MB_SR, c/C_PERCENT);
#endif
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
#ifdef DEBUG
    fprintf(stderr, "E from %g keV to %g keV, diff %g keV\n", E_front/C_KEV, E_back/C_KEV, (E_front - E_back) / C_KEV);
        fprintf(stderr, "stop avg %g eV/tfu\n", params.stop_slope/(C_EV/C_TFU));
        fprintf(stderr, "integration result          = % 18g\n", result);
        fprintf(stderr, "integration estimated error = % 18g\n", error);
        fprintf(stderr, "integration intervals       = %zu\n", ws->w_int_cs->size);
        fprintf(stderr, "final result                = %g mb/sr\n", final/C_MB_SR);
#endif
    return final;
}

double cross_section_concentration_product_fixed(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after,
                                           double S_front, double S_back) {
    depth d;
    d.i = d_after->i;
    const double x_step = (d_after->x - d_before->x) * ws->params->cs_frac;
    const double E_step = (E_back - E_front) * ws->params->cs_frac;
    const double S_step = (S_back - S_front) * ws->params->cs_frac;
    double sum = 0.0;
    for(size_t i = 1; i <= ws->params->cs_n_steps; i++) { /* Compute cross section and concentration product in several "sub-steps" */
        d.x = d_before->x + x_step * i;
        double E = E_front + E_step * i;
#ifdef DEBUG_CS_VERBOSE
        fprintf(stderr, "i = %zu, E = %g keV, (E_front = %g keV, E_back = %g keV)\n", i, E/C_KEV, E_front/C_KEV, E_back/C_KEV);
#endif
        double c = get_conc(sample, d, sim_r->i_isotope);
        double S = S_front + S_step * i;
        double sigma = cross_section_straggling(sim_r, ws->w_int_cs_stragg, ws->params->int_cs_stragg_accuracy, ws->params->cs_stragg_pd, E, S);
#ifdef DEBUG_CS_VERBOSE
        fprintf(stderr, "S = %g, E = %g keV\n", S, E/C_KEV);
#endif
        sum += sigma * c;
    }
    return sample->ranges[d.i].yield * sum / (ws->params->cs_n_steps * 1.0);
}

double cross_section_concentration_product_new(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back) {
    const int type = 1; /* TODO: pick either adaptive or something else */
    const double E_step_max = 10.0 * C_KEV;
    const double depth_step_max = 100.0 * C_TFU;
    const double S_avg_FWHM = C_FWHM * sqrt((S_front + S_back)/2.0); /* units of energy... instread of C_FWHM we should have a configurable "fudge factor" here */
    const double E_step_nominal = -1.0 * GSL_MIN_DBL(E_step_max, S_avg_FWHM); /* Actual step is negative */
    const double d_diff = d_after->x - d_before->x;
    double sigma;
    double c; /* Concentration (but only if constant!) */
    if(sample->no_conc_gradients) {
        c = *sample_conc_bin(sample, d_before->i, sim_r->i_isotope); /* Concentration weighting was not done inside the loop */
        if(c < CONC_TOLERANCE) {
            return 0.0;
        }
    }
    depth d = *d_before;
    if(type == 0) {
        sigma = cross_section_concentration_product_adaptive(ws, sample, sim_r, E_front, E_back, d_before, d_after, S_front, S_back);
    } else {
        double E_diff = E_back - E_front;
        size_t n_steps = GSL_MAX(ceil((E_diff)/E_step_nominal), ceil(d_diff/depth_step_max));
        double frac = 1.0/(1.0*(n_steps+1));
        const double x_step = (d_after->x - d_before->x) * frac;
        const double E_step = (E_back - E_front) * frac;
        const double S_step = (S_back - S_front) * frac;
#ifdef DEBUG
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
#ifdef DEBUG
            fprintf(stderr, "i = %zu, E = %g keV, S = %g keV, (E_front = %g keV, E_back = %g keV), sigma = %g mb/sr\n", i, E/C_KEV, C_FWHM * sqrt(S)/C_KEV, E_front/C_KEV, E_back/C_KEV, sigma_partial / C_MB_SR);
#endif
        }
        sigma = sum / n_steps;
        if(sample->no_conc_gradients) {
            sigma *= c; /* Concentration weighting was not done inside the loop */
        }
    }
    return sample->ranges[d.i].yield * sigma;
}

double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after,
                                           double S_front, double S_back) {
    if(ws->params->mean_conc_and_energy) { /* This if-branch is fast and also serves as a testing branch, since it is a lot easier to understand... */
        const depth d_halfdepth = {.x = (d_before->x + d_after->x) /
                                        2.0, .i = d_after->i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double c = get_conc(sample, d_halfdepth, sim_r->i_isotope);
        if(c < ABUNDANCE_THRESHOLD)
            return 0.0;
        const double E_mean = (E_front + E_back) / 2.0;
        assert(sim_r->cross_section);
        double sigma = sim_r->cross_section(sim_r, E_mean);
        return sigma * c * sample->ranges[d_halfdepth.i].yield;
    } else if(ws->params->cs_n_steps == 0) {
        return cross_section_concentration_product_adaptive(ws, sample, sim_r, E_front, E_back, d_before, d_after, S_front, S_back);
    } else {
        return cross_section_concentration_product_fixed(ws, sample, sim_r, E_front, E_back, d_before, d_after, S_front, S_back);
    }
    return 0.0;
}

void exit_from_sample(ion *p, const depth depth_start, const sim_workspace *ws, const sample *sample) {
    depth d = depth_start;
    while(1) { /* Exit from sample (hopefully) */
#ifdef DEBUG_REACTION
        fprintf(stderr, "  Exiting... depth = %g tfu (i = %zu)\n", d.x, d.i);
#endif
        if(p->inverse_cosine_theta > 0.0 && d.x >= (sample->thickness - DEPTH_TOLERANCE)) { /* Exit through back (transmission) */
            break;
        }
        if(p->inverse_cosine_theta < 0.0 && d.x <= DEPTH_TOLERANCE) { /* Exit (surface, front of sample) */
            break;
        }
        depth d_after = stop_step(ws, p, sample, d, ws->params->stop_step_exiting == 0.0 ? ws->params->stop_step_fudge_factor * (p->E * 0.07 + sqrt(p->S) + 10.0 * C_KEV)
                                                                                        : ws->params->stop_step_exiting); /* TODO: 7% of energy plus straggling plus 10 keV is a weird rule. Automatic stop size should be based more on required accuracy in stopping. */
        if(p->E < ws->emin) {
#ifdef DEBUG_REACTION
            fprintf(stderr,
                            "  Reaction %zu with %s: Energy below EMIN when surfacing from %.3lf tfu, break break.\n",
                            i, r->r->target->name, d_after.x / C_TFU);
#endif
            return;
        }
        d = d_after;
    }
}

double stop_step_calc_incident(const sim_workspace *ws, const ion *ion) {
    if(ws->params->stop_step_incident > 0.0) {
        return ws->params->stop_step_incident;
    }
    double step = sqrt(ion->S) + ws->params->stop_step_add;
    if(step < STOP_STEP_MIN_FALLBACK) {
        return STOP_STEP_MIN_FALLBACK;
    }
    if(ws->params->stop_step_max > 0.0 && step > ws->params->stop_step_max) {
        return ws->params->stop_step_max;
    }
    return step;
}

double stop_step_calculate(const sim_workspace *ws, const ion *ion) { /* Calculate stop step to take */
    if(ws->params->stop_step_incident > 0) {
        return ws->params->stop_step_incident;
    }
    //double E_step = ws->params.stop_step_incident == 0.0?ws->params.stop_step_fudge_factor*sqrt(detector_resolution(ws->det, ion1.isotope, ion1.E)+ion1.S):ws->params.stop_step_incident;
    double broad = sqrt(ion->S) + ws->params->stop_step_add;
    if(broad < ws->params->stop_step_min) {
        return ws->params->stop_step_min;
    }
    return ws->params->stop_step_fudge_factor * broad; /* Fudge factor also affects the minimum stop step */
}

void simulate_reaction_new_routine(const ion *incident, const depth depth_start, sim_workspace *ws, const sample *sample, const des_table *dt, const geostragg_vars *g, sim_reaction *sim_r) {
    ion ion1 = *incident; /* Shallow copy */
    des_table_set_ion_depth(dt, &ion1, depth_start); /* TODO: is this necessary? */
    simulate_init_reaction(sim_r, sample, g, ws->emin, ion1.E);
    depth d_before, d_after = depth_start;
    size_t i_des = 0;
    brick *b = NULL, *b_prev = NULL;
    int skipped, crossed;
    for(size_t i_brick = 0; i_brick < sim_r->n_bricks; i_brick++) {
        assert(ion1.S >= 0.0);
        skipped = FALSE;
        crossed = FALSE;
        b_prev = b;
        b = &sim_r->bricks[i_brick];
        d_before = d_after;
        d_after = des_table_find_depth(dt, &i_des, d_before, &ion1);
        if(ion1.E < ws->emin) {
#ifdef DEBUG
            fprintf(stderr, "Energy minimum.\n");
#endif
            sim_r->last_brick = i_brick - 1;
            break;
        }
        if(i_brick == 0 || d_after.i > d_before.i) { /* There was a layer (depth range) crossing. If stop_step() took this into account when making DES table the only issue is the .i index. depth (.x) is not changed. */
#ifndef NO_SKIP_EMPTY_RANGES
            double conc_start = *sample_conc_bin(sample, d_after.i, sim_r->i_isotope);
            double conc_stop = *sample_conc_bin(sample, d_after.i + 1, sim_r->i_isotope);
#ifdef DEBUG
            fprintf(stderr, "Brick %zu crosses into range %zu, d_before = %g tfu.\n", i_brick, d_after.i, d_before.x / C_TFU);
            fprintf(stderr, "Concentration varies between %g%% and %g%%\n", conc_start * 100.0, conc_stop * 100.0);
#endif
            if(conc_start < CONC_TOLERANCE && conc_stop < CONC_TOLERANCE) { /* This isotope concentration is zero in this layer, skip to next one */
                size_t i_des_skip = dt->depth_interval_index[d_after.i + 1];
                des *des_skip = &(dt->t[i_des_skip]);
                ion1.E = des_skip->E; /* TODO: long skips may make E_deriv calculation below inaccurate */
                ion1.S = des_skip->S;
                d_after = des_skip->d;
                skipped = TRUE;
#ifdef DEBUG
                fprintf(stderr, "Skipped to %g.\n", d_after.x / C_TFU);
#endif
            }
#endif // NO_SKIP_EMPTY_RANGES
            d_before.i = d_after.i;
            crossed = TRUE;
        } else if(d_after.i < d_before.i){
            fprintf(stderr, "This shouldn't happen. But it did. E = %g keV\n", ion1.E / C_KEV);
        }
        assert(ion1.S >= 0.0);
        b->d = d_after;
        b->E_0 = ion1.E;
        b->S_0 = ion1.S;
        if(ws->params->geostragg) {
            b->S_geo_x = geostragg(ws, sample, sim_r, &(g->x), d_after, b->E_0);
            b->S_geo_y = geostragg(ws, sample, sim_r, &(g->y), d_after, b->E_0);
        }
        sim_reaction_product_energy_and_straggling(sim_r, &ion1); /* sets sim_r->p */
        assert(sim_r->p.S >= 0.0);
        b->E_r = sim_r->p.E;
        b->S_r = sim_r->p.S;

        exit_from_sample(&sim_r->p, d_after, ws, sample);

        if(ws->det->foil) { /* Energy loss in detector foil */
            depth d_foil = {.i = 0, .x = 0.0};
            ion ion_foil = *&sim_r->p;
            ion_set_angle(&ion_foil, 0.0, 0.0); /* Foils are not tilted. We use a temporary copy of "p" to do this step. */
            exit_from_sample(&ion_foil, d_foil, ws, ws->det->foil);
            b->E = ion_foil.E;
            b->S = ion_foil.S;
        } else {
            b->E = sim_r->p.E;
            b->S = sim_r->p.S;
        }

        if(b->E < ws->emin) {
            sim_r->last_brick = i_brick - 1;
            break;
        }

        double E_deriv;
        double sigma_conc;
        double d_diff;

        if(b_prev) {
            d_diff = depth_diff(d_before, d_after);
            if(d_diff == 0) { /* Zero thickness brick, so no energy change either */
                E_deriv = b_prev->deriv * 1.5; /* TODO: Can't calculate, assume. The 1.5 is a safety factor. TODO: may lead to exponential growth...  */
#ifdef DEBUG
                fprintf(stderr, "Zero thickness brick (after skip?)\n");
#endif
                sigma_conc = 0.0;
            } else {
                sigma_conc = cross_section_concentration_product_new(ws, sample, sim_r, b_prev->E_0, b->E_0, &d_before, &d_after, b_prev->S_0, b->S_0); /* Product of concentration and sigma for isotope i_isotope target and this reaction. */
                assert(b_prev->E_0 > b->E_0);
                double E0_diff = b_prev->E_0 - b->E_0;
                double E_diff = b_prev->E - b->E;
                E_deriv = fabs(E_diff / E0_diff); /* how many keVs does the reaction product energy change for each keV of incident ion1 energy change */
#ifdef DEBUG
                fprintf(stderr, "should deriv be %g? (based on difference to previous brick in incident energy %g keV and detector energy %g keV)\n", E_deriv, E0_diff / C_KEV, E_diff / C_KEV); /* TODO: this is accurate only near surface */
#endif
                assert(E_deriv > 0.0);
            }
        } else { /* First brick */
            sigma_conc = 0.0;
            d_diff = 0.0;
            double k_incident = stop_sample(ws, &ion1, sample, ws->stopping_type, d_after, b->E_0);
            double k_exiting = stop_sample(ws, &sim_r->p, sample, ws->stopping_type, d_after, b->E_r);
            double K = reaction_product_energy(sim_r->r, sim_r->theta, b->E_0) / b->E_0;
            double k_incident_eff = k_incident * ion1.inverse_cosine_theta;
            double k_exiting_eff = k_exiting * sim_r->p.inverse_cosine_theta;

            E_deriv = fabs((K * k_incident * ion1.inverse_cosine_theta - k_exiting * sim_r->p.inverse_cosine_theta) / (k_incident * ion1.inverse_cosine_theta)); /* TODO: NaN if reaction is not possible */
#ifdef DEBUG
            fprintf(stderr, "should deriv be %g? (calculated using stopping powers and kinematics: %g eV/tfu and %g eV/tfu)\n", E_deriv, k_incident_eff / C_EV_TFU, k_exiting_eff / C_EV_TFU); /* TODO: this is accurate only near surface */
#endif

        }
        if(E_deriv < ENERGY_DERIVATIVE_MIN) {
            E_deriv = ENERGY_DERIVATIVE_MIN;
        }
        b->deriv = E_deriv;
        b->thick = d_diff;
        b->Q = ion1.inverse_cosine_theta * sigma_conc * b->thick;
        b->sc = sigma_conc;
#ifdef DEBUG
        fprintf(stderr, "%12s %3zu %3zu:%10.3lf %3zu:%10.3lf %3zu %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %.3lf %10.3e %8.3lf\n",
                sim_r->r->target->name, i_brick, d_before.i, d_before.x / C_TFU, d_after.i, d_after.x / C_TFU, i_des,
                b->E_0 / C_KEV, sqrt(b->S_0) / C_KEV,
                b->E_r / C_KEV, sqrt(b->S_r) / C_KEV,
                b->E / C_KEV, sqrt(b->S) / C_KEV,
                E_deriv, get_conc(sample, d_after, sim_r->i_isotope) * 100.0, sigma_conc / C_MB_SR, b->Q, b->deriv);
#endif

        assert(!isnan(ion1.E));
        if(d_after.x >= sim_r->max_depth) {
            sim_r->last_brick = i_brick;
            break;
        }
        if(!skipped) {
            double S_sigma = sqrt(detector_resolution(ws->det, sim_r->p.isotope, b->E) + b->S);
            assert(S_sigma > 0.0);
            double E_change = -ws->params->brick_width_sigmas * S_sigma / E_deriv;
            assert(E_change < 0.0);
            ion1.E += E_change;
            if(E_deriv > 9.999) {
#ifdef DEBUG
                fprintf(stderr, "Using a high deriv energy changes by %g keV.\n", E_change/C_KEV);
#endif
            }
        }
    }
}

int simulate(const ion *incident, const depth depth_start, sim_workspace *ws, const sample *sample) {
    geostragg_vars g = geostragg_vars_calculate(ws, incident);
    des_table *dt = des_table_compute(incident, depth_start, ws, sample); /* Depth, energy and straggling of incident ion */
#ifdef DEBUG
    des_table_print(stderr, dt);
#endif
    for(size_t i_reaction = 0; i_reaction < ws->n_reactions; i_reaction++) {
        simulate_reaction_new_routine(incident, depth_start, ws, sample, dt, &g, ws->reactions[i_reaction]);
    }
    des_table_free(dt);
    sim_workspace_histograms_calculate(ws);
    return EXIT_SUCCESS;
}

int simulate_old(const ion *incident, const depth depth_start, sim_workspace *ws, const sample *sample) { /* Ion is expected to be in the sample coordinate system at starting depth. Also note that sample may be slightly different (e.g. due to roughness) to ws->sim->sample */
    assert(sample->n_ranges);
    int warnings = 0;
#ifdef DEBUG
    fprintf(stderr, "Sample thickness is %g\n", sample->thickness/C_TFU);
#endif
    size_t i_depth;
    ion ion = *incident; /* Shallow copy of the incident ion */
    geostragg_vars g = geostragg_vars_calculate(ws, incident);
#ifdef DEBUG
    fprintf(stderr, "Simulate from depth %g tfu (index %zu), detector theta = %g deg, calculated theta = %g deg. %zu reactions.\n", depth_start.x/C_TFU, depth_start.i, ws->det->theta/C_DEG, g.scatter_theta/C_DEG, ws->n_reactions);
    fprintf(stderr, "Ion energy at start %g keV, straggling %g keV FWHM.\n", incident->E/C_KEV, C_FWHM * sqrt(incident->S) / C_KEV);
#endif
    depth d_before = depth_start;
    for(size_t i = 0; i < ws->n_reactions; i++) {
#ifdef DEBUG
        fprintf(stderr, "Initializing reaction %zu\n", i);
#endif
        simulate_init_reaction(ws->reactions[i], sample, &g, ws->emin, ion.E);
    }
    if(fabs(ion.cosine_theta) < 1e-6) {
#ifdef DEBUG
        fprintf(stderr, "Ion was going sideways in the sample, we nudged it a little.\n");
#endif
        ion.theta += 0.01 * C_DEG;
        ion.cosine_theta = cos(ion.theta);
    }
    for(i_depth = 0; i_depth < ws->n_bricks; i_depth++) {
#ifdef DEBUG
        fprintf(stderr, "\r %04zu %8.3lf keV %8.3lf tfu. ", i_depth, ion.E/C_KEV, d_before.x / C_TFU);
#endif
        if(warnings > SIMULATE_WARNING_LIMIT) {
            jabs_message(MSG_ERROR, stderr, "Warning limit reached. Won't calculate anything.\n");
            break;
        }
        if(d_before.x >= sample->thickness) /* We're in too deep. */
            break;
        if(ion.inverse_cosine_theta < 0.0 && d_before.x < 0.001 * C_TFU) /* We're coming out of the sample and about to exit */
            break;
        if(ion.E < ws->emin) {
#ifdef DEBUG
            fprintf(stderr, "Break due to low energy (%.3lf keV < %.3lf keV), x = %.3lf, i_range = %zu.\n", ion.E/C_KEV, ws->sim->emin/C_KEV, d_before.x/C_TFU, d_before.i);
#endif
            break;
        }
        const double E_front = ion.E;
        const double S_front = ion.S;
        double E_step = stop_step_calculate(ws, &ion);
        depth d_after;
        if(i_depth) {
            d_after = stop_step(ws, &ion, sample, d_before, E_step);
        } else {
            d_after = d_before;
        }
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "After:  %g tfu in range %zu\n", d_after.x/C_TFU, d_after.i);
#endif
        const double d_diff = depth_diff(d_before, d_after);
        /* DEPTH BIN [x, x+d_diff) */
        const double E_back = ion.E;
        const double S_back = ion.S;
        if(i_depth && fabs(d_diff) < 0.001 * C_TFU && E_front - E_back < 0.001 * C_KEV) {
            jabs_message(MSG_WARNING, stderr,
                         "Warning: no or very little progress was made (E step (goal) %g keV, E from %g keV to E = %g keV, depth = %g tfu, d_diff = %g tfu), check stopping or step size.\n",
                         E_step / C_KEV, E_front / C_KEV, E_back / C_KEV, d_before.x / C_TFU, d_diff / C_TFU);
            //sample_print(stderr, sample, FALSE);
            //d_before.x += incident->inverse_cosine_theta*0.0001*C_TFU;
            d_before = d_after;
            d_before.x += 0.002 * C_TFU;
            warnings++;
            continue;
        }

#ifdef DEBUG_VERBOSE
        double E_diff = E_front-E_back;
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3lf keV to  %8.3lf keV (diff %6.4lf keV)\n", x/C_TFU, (x+h)/C_TFU, E_front/C_KEV, ws->ion.E/C_KEV, E_diff/C_KEV);
#endif
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(ion.S) / C_KEV);
#endif
        int alive = 0;
        for(size_t i = 0; i < ws->n_reactions; i++) {
            if(ws->reactions[i]->stop) {
                continue;
            }
            simulate_reaction(ws->reactions[i], ws, sample, &g, i_depth, d_before, d_after, &ion, E_front, S_front, E_back, S_back, d_diff);
            if(!ws->reactions[i]->stop) {
                alive++;
            }
        }
        d_before = d_after;
        if(!alive) {
#ifdef DEBUG
            fprintf(stderr, "\nAll reactions have ceased by depth %g.\n", d_before.x/C_TFU);
#endif
            break;
        }
    }
    if(i_depth == ws->n_bricks - 1) {
        jabs_message(MSG_WARNING, stderr, "Warning: used all allocated bricks (%zu), reached depth of %g tfu. Simulated spectra may be partial.\n", ws->n_bricks, d_before.x / C_TFU);
    }
    sim_workspace_histograms_calculate(ws);
    return EXIT_SUCCESS;
}

void simulate_reaction(sim_reaction *sim_r, const sim_workspace *ws, const sample *sample, const geostragg_vars *g, size_t i_depth, const depth d_before, const depth d_after, const ion *incident, double E_front, double S_front, double E_back, double S_back, double d_diff) {
#ifdef N_BRICKS_REDUNDANT_CHECK
    if(i_depth >= r->n_bricks) { /* This check is currently not necessary, since r->n_bricks == ws->n_bricks */
                fprintf(stderr, "Too many bricks (%zu max). Data partial.\n", r->n_bricks);
                r->stop = TRUE;
                continue;
            }
#endif
    brick *b = &sim_r->bricks[i_depth];
    if(!ws->params->ds && d_before.x >= sim_r->max_depth) { /* Reactions stop when we are too deep in the sample, unless, of course, if DS is enabled. TODO: check optimizations for DS */
#ifdef DEBUG
        fprintf(stderr, "\nReaction with %s stops, because maximum depth is reached at x = %.3lf tfu.\n", sim_r->r->target->name, d_before.x / C_TFU); /* TODO: give reactions a name */
#endif
        b->Q = 0.0;
        sim_r->stop = TRUE;
        return;
    }
#if 0
    if(incident->E + 3.0*sqrt(incident->S) < sim_r->r->E_min) { /* Beam average energy is three sigmas below reaction minimum, we can stop calculation (due to straggling weighted cross sections we can't stop immediately below E_min). */
        b->Q = 0.0;
        sim_r->stop = TRUE;
        return;
    }
#endif
    b->d = d_before;
    b->E_0 = incident->E; /* Sort of energy just before the reaction. */
    b->S_0 = incident->S;

#ifdef DEBUG_REACTION
    fprintf(stderr, "Reaction %s (%zu): %s\n", reaction_name(r->r), i, r->r->target->name);
#endif
    if(ws->params->geostragg) {
        b->S_geo_x = geostragg(ws, sample, sim_r, &(g->x), d_after, incident->E);
        b->S_geo_y = geostragg(ws, sample, sim_r, &(g->y), d_after, incident->E);
    }
    sim_reaction_product_energy_and_straggling(sim_r, incident);
    assert(sim_r->p.E > 0.0);
    b->E_r = sim_r->p.E;
    b->S_r = sim_r->p.S;
    exit_from_sample(&sim_r->p, d_after, ws, sample);

    if(ws->det->foil) { /* Energy loss in detector foil */
        depth d_foil = {.i = 0, .x = 0.0};
        ion ion_foil = *&sim_r->p;
        ion_set_angle(&ion_foil, 0.0, 0.0); /* Foils are not tilted. We use a temporary copy of "p" to do this step. */
        exit_from_sample(&ion_foil, d_foil, ws, ws->det->foil);
        b->E = ion_foil.E;
        b->S = ion_foil.S;
    } else {
        b->E = sim_r->p.E;
        b->S = sim_r->p.S;
    }
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "Reaction %2zu depth from %8.3lf tfu to %8.3lf tfu, E = %8.3lf keV, Straggling: eloss %7.3lf keV, geo %7.3lf keV\n", i, d_before.x/C_TFU, d_after.x/C_TFU, b->E/C_KEV, sqrt(b->S)/C_KEV, sqrt(b->S_geo_x+b->S_geo_y)/C_KEV);
#endif
    if(sim_r->p.E < ws->emin) {
        sim_r->stop = TRUE;
        b->Q = 0.0;
        return;
    }
    b->thick = d_diff;
    double sigma_conc;
    if(i_depth) {
        sigma_conc = cross_section_concentration_product(ws, sample, sim_r, E_front, E_back, &d_before, &d_after, S_front, S_back); /* Product of concentration and sigma for isotope i_isotope target and this reaction. */
    } else {
        sigma_conc = 0.0;
    }
    if(sigma_conc > 0.0) {
        if(d_after.i == sample->n_ranges - 2) {
            sigma_conc *= ws->sim->channeling_offset + ws->sim->channeling_slope * (E_front + E_back) / 2.0;
        }
        b->Q = incident->inverse_cosine_theta * sigma_conc * b->thick;
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
    sim_r->last_brick = i_depth;
}

void simulate_init_reaction(sim_reaction *sim_r, const sample *sample, const geostragg_vars *g, double E_min, double E_max) {
    if(!sim_r) {
        fprintf(stderr, "Simulation reaction is NULL\n");
        return;
    }
    sim_r->last_brick = 0;
    sim_r->stop = FALSE;
    ion_set_angle(&sim_r->p, g->theta_product, g->phi_product);
    sim_reaction_recalculate_internal_variables(sim_r, g->scatter_theta, E_min, E_max);
    if(sim_r->stop) {
        sim_r->max_depth = 0.0;
        return;
    }
    sim_r->max_depth = sample_isotope_max_depth(sample, sim_r->i_isotope);
    sim_reaction_reset_bricks(sim_r);
    if(sim_r->i_isotope >= sample->n_isotopes) { /* No target isotope for reaction. */
        sim_r->stop = TRUE;
    }
#ifdef DEBUG
    fprintf(stderr, "Simulation reaction: %s %s. Max depth %g tfu. i_isotope=%zu, stop = %i.\n\n"
            , reaction_name(sim_r->r), sim_r->r->target->name, sim_r->max_depth / C_TFU, sim_r->i_isotope, sim_r->stop);
    ion_print(stderr, &sim_r->p);
#endif
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
        if(l > 4 && strncmp(filename + l - 4, ".csv", 4) == 0) { /* For CSV: print header line */
            sep = ','; /* and set the separator! */
            fprintf(f, "\"Channel\",\"Energy (keV)\",\"Simulated\"");
            if(exp) {
                fprintf(f, ",\"Experimental\"");
            }
            for(size_t j = 0; j < ws->n_reactions; j++) {
                const reaction *r = ws->reactions[j]->r;
                fprintf(f, ",\"%s (%s)\"", r->target->name, reaction_name(r));
            }
            fprintf(f, "\n");
        }
    }
    for(size_t i = 0; i < ws->n_channels; i++) {
        fprintf(f, "%zu%c%.3lf%c", i, sep, detector_calibrated(ws->det, JIBAL_ANY_Z, i) / C_KEV,
                sep); /* Channel, energy. TODO: Z-specific calibration can have different energy (e.g. for a particular reaction). */
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
        for(size_t j = 0; j < ws->n_reactions; j++) {
            gsl_histogram *histo = ws->reactions[j]->histo;
            if(i >= histo->n || histo->bin[i] == 0.0) {
                fprintf(f, "%c0", sep);
            } else {
                fprintf(f, "%c%e", sep, histo->bin[i]);
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
    fit_params_add_parameter(params, &sim->fluence, "fluence", "", 1.0); /* This must be the first parameter always, as there is a speedup in the fit routine */
    fit_params_add_parameter(params, &sim->channeling_offset, "channeling", "", 1.0);
    fit_params_add_parameter(params, &sim->channeling_slope, "channeling_slope", "1/keV", 1.0 / C_KEV);
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
                         (Z == JIBAL_ANY_Z) ? "" : "_",
                         (Z == JIBAL_ANY_Z) ? "" : jibal_element_name(fit->jibal->elements, Z),
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
            size_t range_index = i_range + 1; /* Human readable indexing */
            if(r->x > 0.0) {
                snprintf(param_name, param_name_max_len, "thick%zu", range_index);
                fit_params_add_parameter(params, &(r->x), param_name, "tfu", C_TFU);
            }
            if(r->yield != 1.0) {
                snprintf(param_name, param_name_max_len, "yield%zu", range_index);
                fit_params_add_parameter(params, &(r->yield), param_name, "", 1.0);
            }
            if(r->bragg != 1.0) {
                snprintf(param_name, param_name_max_len, "bragg%zu", range_index);
                fit_params_add_parameter(params, &(r->bragg), param_name, "", 1.0);
            }
            if(r->stragg != 1.0) {
                snprintf(param_name, param_name_max_len, "stragg%zu", range_index);
                fit_params_add_parameter(params, &(r->stragg), param_name, "", 1.0);
            }
            if(r->rough.model != ROUGHNESS_NONE && r->rough.x > 0.0) {
                snprintf(param_name, param_name_max_len, "rough%zu", range_index);
                fit_params_add_parameter(params, &(r->rough.x), param_name, "tfu", C_TFU);
            }
            for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
                if(*sample_model_conc_bin(sm, i_range, i_mat) < CONC_TOLERANCE) /* Don't add fit variables for negative, zero or very low concentrations. */
                    continue;
                snprintf(param_name, param_name_max_len, "conc%zu_%s", range_index, sm->materials[i_mat]->name);
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
        const sim_reaction *r = ws->reactions[i];
        fprintf(f, "#Reaction %zu: %s %s\n", i+1, reaction_name(r->r), r->r->target->name);
        if(r->r->filename) {
            fprintf(f, "#Filename: %s\n", r->r->filename);
        }
        fprintf(f, "#i  brick    depth    thick      E_0  S_0(el)      E_r  S_r(el)   E(det)    S(el)    S(geo)    S(sum) sigma*conc            Q  dE(det)/dE_0\n");
        for(size_t j = 0; j <= r->last_brick; j++) {
            brick *b = &r->bricks[j];
            fprintf(f, "%2zu %4zu %10.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %10.1lf %10.3lf %12.6g %8.3lf\n",
                    i, j, b->d.x / C_TFU, b->thick / C_TFU,
                    b->E_0 / C_KEV, C_FWHM * sqrt(b->S_0) / C_KEV,
                    b->E_r / C_KEV, C_FWHM * sqrt(b->S_r) / C_KEV,
                    b->E / C_KEV, C_FWHM * sqrt(b->S) / C_KEV,
                    C_FWHM * sqrt(b->S_geo_x + b->S_geo_y) / C_KEV, C_FWHM * b->S_sum / C_KEV,
                    b->sc / C_MB_SR, b->Q * ws->fluence * ws->det->solid,
                    b->deriv
                    );
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
#ifdef DEBUG
    fprintf(stderr, "%zu rough layers\n", n_rl);
#endif
    if(!n_rl) {
        ion_set_angle(&ws->ion, 0.0 * C_DEG, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        return simulate(&ws->ion, depth_seek(ws->sample, 0.0 * C_TFU), ws, ws->sample);
    }
    struct sample *sample_rough = sample_copy(ws->sample);
    size_t i_rl = 0;
    thick_prob_dist **tpds = malloc(sizeof(thick_prob_dist *) * n_rl);
    for(size_t i = 0; i < ws->sample->n_ranges; i++) {
        sample_range *r = &(ws->sample->ranges[i]);
        if(r->rough.model == ROUGHNESS_NONE)
            continue;
        tpds[i_rl] = NULL;
        if(r->rough.model == ROUGHNESS_GAMMA) {
#ifdef DEBUG
            fprintf(stderr, "Range %zu is rough (gamma), amount %g tfu, n = %zu spectra\n", i, ws->sample->ranges[i].rough.x/C_TFU, ws->sample->ranges[i].rough.n);
#endif
            tpds[i_rl] = thickness_probability_table_gamma(r->x, r->rough.x, r->rough.n);
        } else if(r->rough.model == ROUGHNESS_FILE) {
#ifdef DEBUG
            fprintf(stderr, "Range %zu is rough (FILE), filename = \"%s\"\n", i, r->rough.file->filename);
#endif
            tpds[i_rl] = thickness_probability_table_copy(r->rough.file->tpd);
        }
        if(tpds[i_rl]) {
            tpds[i_rl]->i_range = i;
#ifdef DEBUG
            fprintf(stderr, "TPD (i_range %zu) for depth %zu (%.3lf tfu nominal), roughness %.3lf tfu:\n", tpds[i_rl]->i_range, i, ws->sample->ranges[i].x/C_TFU, ws->sample->ranges[i].rough.x/C_TFU);
            thickness_probability_table_print(stderr, tpds[i_rl]);
#endif
            i_rl++;
        }
    }
    size_t iter_total = 1;
    for(i_rl = 0; i_rl < n_rl; i_rl++) { /* Calculate cumulative product (number of subspectra) */
        thick_prob_dist *tpd = tpds[i_rl];
        tpd->modulo = iter_total;
        iter_total *= tpd->n;
    }
    for(size_t i_iter = 0; i_iter < iter_total; i_iter++) {
#ifdef DEBUG
        fprintf(stderr, "Roughness step %zu/%zu.\n", i_iter+1, iter_total);
#endif
        double p = 1.0;
        for(size_t i_range = 0; i_range < ws->sample->n_ranges; i_range++) { /* Reset ranges for every iter */
            sample_rough->ranges[i_range].x = ws->sample->ranges[i_range].x;
        }
        for(i_rl = 0; i_rl < n_rl; i_rl++) {
            thick_prob_dist *tpd = tpds[i_rl]; /* One particular thickness probability distribution ("i"th one) */
            size_t j = (i_iter / tpd->modulo) % tpd->n; /* "j"th roughness element */
            thick_prob *pj = &tpds[i_rl]->p[j]; /* ..is this one */
            p *= pj->prob; /* Probability is multiplied by the "i"th roughness, element "j" to get the subspectra weight */
            double x_diff = pj->x - ws->sample->ranges[tpd->i_range].x; /* Amount to change thickness of this and all subsequent layers */
#ifdef DEBUG
            fprintf(stderr, "Modifying ranges from %zu to %zu by %g tfu.\n", tpd->i_range, ws->sample->n_ranges, x_diff/C_TFU);
#endif
            for(size_t i_range = tpd->i_range; i_range < ws->sample->n_ranges; i_range++) {
                sample_rough->ranges[i_range].x += x_diff;
            }
        }
#ifdef DEBUG
        fprintf(stderr, "Weight %.6lf.\n", p);
        sample_print(NULL, sample_rough, FALSE);
#endif
        ws->fluence = p * p_sr;
        ion_set_angle(&ws->ion, 0.0, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        sample_thickness_recalculate(sample_rough);
        status = simulate(&ws->ion, depth_seek(ws->sample, 0.0), ws, sample_rough);
        if(status != EXIT_SUCCESS)
            break;
    }
    for(size_t i = 0; i < n_rl; i++) {
        thickness_probability_table_free(tpds[i]);
    }
    sample_free(sample_rough);
    free(tpds);
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
    sim_calc_params_update(ws->params);
    jabs_message(MSG_ERROR, stderr, "\n");
    const jibal_isotope *incident = ws->sim->beam_isotope;
    while(1) {
        double E_front = ion1.E;
        if(E_front < ws->emin)
            break;
        depth d_after = stop_step(ws, &ion1, ws->sample, d_before, ws->params->stop_step_fudge_factor * sqrt(detector_resolution(ws->det, ion1.isotope, ion1.E) + ion1.S)); /* TODO: step? */
        double thick_step = depth_diff(d_before, d_after);
        const depth d_halfdepth = {.x = (d_before.x + d_after.x) /
                                        2.0, .i = d_after.i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double E_back = ion1.E;
        const double E_mean = (E_front + E_back) / 2.0;

        jabs_message(MSG_VERBOSE, stderr, "DS depth from %9.3lf tfu to %9.3lf tfu, E from %6.1lf keV to %6.1lf keV. p*sr = %g\n", d_before.x / C_TFU, d_after.x / C_TFU, E_front / C_KEV,
                     E_back / C_KEV, fluence);
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
#ifdef DEBUG
                    if(d_before.x == 0.0) {
                    fprintf(stderr, "DS polar %.3lf, azi %.3lf, scatter %.3lf\n", ds_polar/C_DEG, ds_azi/C_DEG, scattering_angle(&ion2, ws)/C_DEG);
                }
#endif
                    if(scattering_angle(&ion2, ws) > 19.99999 * C_DEG) {
                        if(simulate(&ion2, d_halfdepth, ws, ws->sample)) {
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        if(ws->sample->ranges[ws->sample->n_ranges - 1].x - d_after.x < 0.01 * C_TFU)
            break;
        d_before = d_after;
    }
    sim_workspace_calculate_sum_spectra(ws);
    return EXIT_SUCCESS;
}
