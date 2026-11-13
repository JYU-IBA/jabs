/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

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
#include "message.h"
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
            jabs_message(MSG_WARNING, "Warning: probably too deep! This is a bug!\n");
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


void stop_ms_step(const ms_model *ms, ion *ion, const sample *sample, const depth d, double E, double h) {
    if(!ms || ion->Z == 0 || h <= 0.0 || E <= 0.0) {
        return;
    }
    int Z2_prev = 0;
    double c_sum = 0.0;
    for(size_t i_isotope = 0; i_isotope <= sample->n_isotopes; i_isotope++) { /* Isotopes of the same element are summed if they are consecutive. One extra iteration to handle the last element. */
        int Z2 = 0;
        double c = 0.0;
        if(i_isotope < sample->n_isotopes) {
            Z2 = sample->isotopes[i_isotope]->Z;
            c = sample->no_conc_gradients ? *sample_conc_bin(sample, d.i, i_isotope) : get_conc(sample, d, i_isotope);
        }
        if(Z2 != Z2_prev) {
            if(Z2_prev > 0 && c_sum > CONC_TOLERANCE) {
                ms_step_element(ms, &ion->ms, ion->Z, Z2_prev, c_sum, E, h);
            }
            c_sum = 0.0;
            Z2_prev = Z2;
        }
        c_sum += c;
    }
}

static double stop_gsto(const ion *incident, gsto_stopping_type type, int Z2, double em, double ln_em) { /* Electronic stopping or straggling of incident ion in element Z2, see ion_gsto_table_init(). ln_em = log(em). */
    const jabs_ion_gsto_data *gd = &incident->ion_gsto->gsto_data[Z2];
    const jabs_gsto_table *t = (type == GSTO_STO_STRAGG) ? &gd->stragg : &gd->stop;
    assert(t->file);
    return ion_gsto_table_eval_ln(t, em, ln_em);
}

static void stop_sample_sw(const jabs_stop *stop, const jabs_stop *stragg, const ion *incident, const sample *sample, const depth depth, double E, stop_sw *sw) { /* Stopping, straggling and stopping gradient in one go */
    sw->S = 0.0;
    sw->W = 0.0;
    sw->S_x = 0.0;
    if(incident->Z == 0) {
        return;
    }
    const double em = E * incident->mass_inverse;
    const double ln_em = log(em); /* Once for all tables */
    double inverse_width = 0.0; /* Concentrations are linear within a range. Zero if concentrations are constant. */
    if(!sample->no_conc_gradients) {
        double width = sample->ranges[depth.i + 1].x - sample->ranges[depth.i].x;
        if(width > 0.0) {
            inverse_width = 1.0 / width;
        }
    }
    int Z2_old = -1; /* Not valid */
    double S_Z2 = 0.0, W_Z2 = 0.0;
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        const jibal_isotope *target = sample->isotopes[i_isotope];
        double c = get_conc(sample, depth, i_isotope);
        double dc = 0.0;
        if(inverse_width > 0.0) {
            const double *cbin = sample_conc_bin(sample, depth.i, i_isotope);
            dc = (cbin[sample->n_isotopes] - cbin[0]) * inverse_width;
        }
        if(c < CONC_TOLERANCE && dc == 0.0)
            continue;
        int Z2 = target->Z;
        assert(Z2 <= incident->ion_gsto->Z2_max);
        if(Z2 != Z2_old) { /* Electronic stopping and straggling are the same for all isotopes of the same element */
            S_Z2 = stop_gsto(incident, stop->type, Z2, em, ln_em);
            W_Z2 = stop_gsto(incident, stragg->type, Z2, em, ln_em);
            Z2_old = Z2;
        }
        double S = S_Z2;
        if(stop->type == GSTO_STO_TOT) {
            S += ion_nuclear_stop(incident, target, E);
        }
        sw->S += c * S;
        sw->S_x += dc * S;
        sw->W += c * W_Z2;
    }
    const sample_range *r = &sample->ranges[depth.i];
    sw->S *= r->bragg;
    sw->S_x *= r->bragg;
    sw->W *= r->stragg;
}

static inline double stop_u_derivative(const stop_sw *sw, double u, double cosine_theta) { /* du/dl, see stop_step() */
    if(sw->S <= 0.0) {
        return 0.0;
    }
    return sw->W / pow2(sw->S) - 2.0 * u * sw->S_x * cosine_theta / sw->S;
}

depth stop_step(const jabs_stop *stop, const jabs_stop *stragg, ion *incident, const sample *sample, depth depth_before, double step) {
    double k1, k2, k3, k4, stopping, dE, E;
    depth depth_next = stop_next_crossing(incident, sample, &depth_before);
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
    /* Energy E and u = S_variance / S^2 are integrated along the path (l), where S_variance is the energy loss variance (incident->S). Two ions at the same depth with
     * energy difference delta obey d(delta)/dl = -(dS/dE) delta, therefore d(S_variance)/dl = W - 2 (dS/dE) S_variance. Derivative of stopping with respect to energy
     * cancels out from du/dl = W / S^2 - 2 u (dS/dl)_E / S, where (dS/dl)_E is the change of stopping due to concentration gradients (zero in layered samples).
     * k1...k4 are slopes of energy loss (stopping) at various x (depth) and E. Note convention: positive values, i.e. -dE/dx! */
    E = incident->E;
    stop_sw sw1;
    stop_sample_sw(stop, stragg, incident, sample, depth_before, E, &sw1);
    k1 = sw1.S;
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


    const double u = incident->S / pow2(k1);
    double du; /* Mean du/dl over the step */
    double g1 = stop_u_derivative(&sw1, u, incident->cosine_theta);
    if(stop->rk4) {
        stop_sw sw2, sw3, sw4;
        stop_sample_sw(stop, stragg, incident, sample, halfdepth, E - (h / 2.0) * k1, &sw2);
        double g2 = stop_u_derivative(&sw2, u + (h / 2.0) * g1, incident->cosine_theta);
        k2 = sw2.S;
        stop_sample_sw(stop, stragg, incident, sample, halfdepth, E - (h / 2.0) * k2, &sw3);
        double g3 = stop_u_derivative(&sw3, u + (h / 2.0) * g2, incident->cosine_theta);
        k3 = sw3.S;
        stop_sample_sw(stop, stragg, incident, sample, fulldepth, E - h * k3, &sw4);
        double g4 = stop_u_derivative(&sw4, u + h * g3, incident->cosine_theta);
        k4 = sw4.S;
        stopping = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        du = (g1 + 2.0 * g2 + 2.0 * g3 + g4) / 6.0;
#ifdef DEBUG
        if(stopping <= 0.0) {
            DEBUGMSG("Stopping is negative: RK4 (k1 = %g, k2 = %g, k3 = %g, k4 = %g) eV/tfu\n",
                     k1 / C_EV_TFU, k2 / C_EV_TFU, k3 / C_EV_TFU, k4 / C_EV_TFU);
        }
#endif
    } else { /* Energy: Euler. u: midpoint. */
        stopping = k1;
        stop_sw sw_mid;
        stop_sample_sw(stop, stragg, incident, sample, halfdepth, E - (h / 2.0) * k1, &sw_mid);
        du = stop_u_derivative(&sw_mid, u + (h / 2.0) * g1, incident->cosine_theta);
    }
#ifdef DEBUG
    if(stopping <= 0.0) {
        DEBUGMSG("Stopping is %g eV/tfu at E = %g keV, k1 = %g eV/tfu, h = %g tfu, step = %g keV\n", stopping / C_EV_TFU, E / C_KEV, k1 / C_EV_TFU, h / C_TFU, step / C_KEV);
    }
    assert(stopping > 0.0);
#endif
    dE = -1.0 * h * stopping; /* Energy change in thickness "h". It is always negative! */
    double stop_end = stop_sample(stop, incident, sample, fulldepth, E + dE);
    incident->S = GSL_MAX_DBL(0.0, u + h * du) * pow2(stop_end);
    if(stop->ms) {
        double p_ratio = stop_end / k1; /* Energy differences propagate over the step like straggling above */
        incident->ms.P *= p_ratio;
        incident->ms.I = incident->ms.I * p_ratio + h * stop_end; /* Exact for constant composition: I = l * S(E) */
        stop_ms_step(stop->ms, incident, sample, halfdepth, E + (0.5 * dE), h);
    }
    incident->E += dE;
    return fulldepth; /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}

/* Dormand-Prince 5(4): nodes, stages (the last row is the fifth order solution) and error estimate coefficients (difference of fifth and fourth order solutions).
 * The last stage is evaluated at the fifth order solution and it is reused as the first stage of the next step. */
static const double dp_c[7] = {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};
static const double dp_a[7][6] = {
        {0.0},
        {1.0 / 5.0},
        {3.0 / 40.0, 9.0 / 40.0},
        {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0},
        {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0},
        {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0},
        {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}
};
#define STOP_ADAPTIVE_SAFETY (0.9)
#define STOP_ADAPTIVE_FACTOR_MIN (0.2) /* Smallest and largest change in step from one step to the next */
#define STOP_ADAPTIVE_FACTOR_MAX (2.0)
static const double dp_e[7] = {71.0 / 57600.0, 0.0, -71.0 / 16695.0, 71.0 / 1920.0, -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0};

depth stop_step_adaptive(const jabs_stop *stop, const jabs_stop *stragg, ion *incident, const sample *sample, depth depth_before, double step, stop_adaptive *state) {
    /* Same differential equations for E and u as in stop_step(), but integrated with Dormand-Prince 5(4) and the step (path length) is chosen so that the estimated
     * error of energy is about stop->accuracy and relative error of u is about stop->stragg_accuracy. The step is also limited by energy loss "step" and by range
     * boundaries. The error estimate of a step only determines the next step, steps are not rejected (and retried with a shorter step). This keeps the result a
     * continuous function of the initial conditions and the sample, which is important for fitting (numerical derivatives). */
    depth depth_next = stop_next_crossing(incident, sample, &depth_before);
    double h_max_perp = depth_next.x - depth_before.x;
    depth_before.i = depth_next.i;
    const double E = incident->E;
    stop_sw sw[7];
    if(state->valid && state->d_end.i == depth_before.i && state->d_end.x == depth_before.x && state->E_end == E) {
        sw[0] = state->sw_end;
    } else {
        stop_sample_sw(stop, stragg, incident, sample, depth_before, E, &sw[0]);
        state->n_evals++;
    }
    state->valid = FALSE;
    const double k1 = sw[0].S;
    if(k1 < STOP_STEP_MINIMUM_STOPPING) { /* See stop_step() */
        if(STOP_STEP_DEPTH_FALLBACK > h_max_perp) {
            return depth_next;
        }
        depth_before.x += STOP_STEP_DEPTH_FALLBACK;
        return depth_before;
    }
    const double h_max = h_max_perp * incident->inverse_cosine_theta; /* Path length to range boundary, positive */
    double h = step / k1;
    if(state->h > 0.0) {
        h = GSL_MIN_DBL(h, state->h);
    }
    h = GSL_MAX_DBL(h, STOP_STEP_ABSOLUTE_MINIMUM_STEP);
    const double u = incident->S / pow2(k1);
    double k[7], g[7];
    k[0] = k1;
    g[0] = stop_u_derivative(&sw[0], u, incident->cosine_theta);
    double E_new = E, u_new = u;
    depth fulldepth;
    const double h_ctrl = h; /* Step suggested by error control (or initial guess) */
    int boundary = FALSE;
    if(h >= h_max) {
        h = h_max;
        boundary = TRUE;
        if(h < 0.001 * C_TFU) {
            return depth_next;
        }
    }
    double factor = 1.0;
    while(1) {
        fulldepth.i = depth_before.i;
        fulldepth.x = boundary ? depth_next.x : depth_before.x + h * incident->cosine_theta;
        int ok = TRUE;
        for(int s = 1; s < 7; s++) {
            double sk = 0.0, sg = 0.0;
            for(int j = 0; j < s; j++) {
                sk += dp_a[s][j] * k[j];
                sg += dp_a[s][j] * g[j];
            }
            E_new = E - h * sk;
            u_new = u + h * sg;
            if(E_new <= 0.0) {
                ok = FALSE;
                break;
            }
            depth d = fulldepth;
            if(dp_c[s] < 1.0) {
                d.x = depth_before.x + dp_c[s] * (fulldepth.x - depth_before.x);
            }
            stop_sample_sw(stop, stragg, incident, sample, d, E_new, &sw[s]);
            state->n_evals++;
            k[s] = sw[s].S;
            g[s] = stop_u_derivative(&sw[s], u_new, incident->cosine_theta);
        } /* After the last stage E_new and u_new are the fifth order solution */
        double err = 0.0; /* Error relative to tolerance */
        factor = STOP_ADAPTIVE_FACTOR_MIN;
        if(ok) {
            double err_E = 0.0, err_u = 0.0;
            for(int j = 0; j < 7; j++) {
                err_E += dp_e[j] * k[j];
                err_u += dp_e[j] * g[j];
            }
            err = h * fabs(err_E) / stop->accuracy;
            double u_scale = stop->stragg_accuracy * GSL_MAX_DBL(fabs(u), fabs(u_new));
            if(u_scale > 0.0) {
                err = GSL_MAX_DBL(err, h * fabs(err_u) / u_scale);
            }
            factor = err > 0.0 ? STOP_ADAPTIVE_SAFETY * pow(err, -0.2) : 1e10;
            factor = GSL_MIN_DBL(GSL_MAX_DBL(factor, STOP_ADAPTIVE_FACTOR_MIN), STOP_ADAPTIVE_FACTOR_MAX);
            break;
        }
        state->n_rejected++; /* Energy went below zero, should not happen often since steps end at emin */
        h *= STOP_ADAPTIVE_FACTOR_MIN;
        boundary = FALSE;
        if(h <= STOP_STEP_ABSOLUTE_MINIMUM_STEP) {
            incident->E = GSL_MAX_DBL(0.0, E - h * k1);
            fulldepth.x = depth_before.x + h * incident->cosine_theta;
            return fulldepth;
        }
    }
    state->h = h * factor;
    if(boundary) { /* Step truncated at a boundary does not tell much about the next step. This is continuous with respect to h_ctrl when h_ctrl = h_max. */
        state->h = GSL_MAX_DBL(state->h, h_ctrl * GSL_MIN_DBL(factor, 1.0));
    }
    const double stop_end = k[6];
    incident->S = GSL_MAX_DBL(0.0, u_new) * pow2(stop_end);
    if(stop->ms) { /* See stop_step() */
        double p_ratio = stop_end / k1;
        incident->ms.P *= p_ratio;
        incident->ms.I = incident->ms.I * p_ratio + h * stop_end;
        depth halfdepth = {.i = depth_before.i, .x = (depth_before.x + fulldepth.x) / 2.0};
        stop_ms_step(stop->ms, incident, sample, halfdepth, (E + E_new) / 2.0, h);
    }
    incident->E = GSL_MAX_DBL(0.0, E_new);
    state->sw_end = sw[6];
    state->d_end = fulldepth;
    state->E_end = incident->E;
    state->valid = (E_new > 0.0);
    return fulldepth;
}

double stop_sample(const jabs_stop *stop, const ion *incident, const sample *sample, const depth depth, double E) {
    const double em = E * incident->mass_inverse;
    int Z2_old = -1; /* Not valid */
    double S_Z2 = 0.0;
    double out = 0.0;
    if(incident->Z == 0) {
        return 0.0;
    }
    const double ln_em = log(em); /* Once for all tables */
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        const jibal_isotope *target = sample->isotopes[i_isotope];
        double c = get_conc(sample, depth, i_isotope);
        if(c < CONC_TOLERANCE)
            continue;
        int Z2 = target->Z;
        assert(Z2 <= incident->ion_gsto->Z2_max);
        if(Z2 != Z2_old) { /* Electronic stopping (or straggling) is the same for all isotopes of the same element */
            S_Z2 = stop_gsto(incident, stop->type, Z2, em, ln_em);
            Z2_old = Z2;
        }
        double S = S_Z2;
        if(stop->type == GSTO_STO_TOT) {
            S += ion_nuclear_stop(incident, target, E);
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

double stop_step_max(const jabs_stop_step_params *params, const ion *ion) { /* Largest step, used as a limit for adaptive steps */
    if(params->step > 0.0) {
        return params->step;
    }
    if(params->max_rel > 0.0) {
        return GSL_MIN_DBL(params->max, params->max_rel * ion->E);
    }
    return params->max;
}

double stop_step_calc(const jabs_stop_step_params *params, const ion *ion) {
    if(params->step > 0.0) {
        return params->step;
    }
    double step = params->sigmas * sqrt(ion->S);
    step = GSL_MIN_DBL(step, params->max);
    step = GSL_MAX_DBL(step, params->min);
    if(params->max_rel > 0.0) {
        step = GSL_MIN_DBL(step, params->max_rel * ion->E);
    }
    return step;
}

int stop_sample_exit(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params_exiting, ion *p, const depth depth_start, const sample *sample) {
    return stop_sample_exit_remaining(stop, stragg, params_exiting, p, depth_start, sample, NULL);
}

static double stop_remaining_path(const ion *p, const sample *sample, depth d) { /* Path length from depth d to the surface the ion is heading to */
    double x = (p->inverse_cosine_theta > 0.0) ? sample->thickness - d.x : d.x;
    return GSL_MAX_DBL(x, 0.0) * fabs(p->inverse_cosine_theta);
}

stop_range_tables *stop_range_tables_new(void) {
    return calloc(1, sizeof(stop_range_tables));
}

static void stop_range_table_free(stop_range_table *t) {
    if(!t) {
        return;
    }
    free(t->isotopes);
    free(t->conc);
    free(t->E);
    free(t->R);
    free(t->S);
    free(t->Q);
    free(t->dQ);
    free(t);
}

void stop_range_tables_free(stop_range_tables *rt) {
    if(!rt) {
        return;
    }
    for(size_t i = 0; i < rt->n; i++) {
        stop_range_table_free(rt->t[i]);
    }
    free(rt->t);
    free(rt);
}

static int sample_range_conc_constant(const sample *sample, size_t i_range) {
    if(i_range + 1 >= sample->n_ranges) {
        return FALSE;
    }
    if(sample->no_conc_gradients) {
        return TRUE;
    }
    const double *c = sample_conc_bin(sample, i_range, 0);
    for(size_t i = 0; i < sample->n_isotopes; i++) {
        if(c[i] != c[i + sample->n_isotopes]) {
            return FALSE;
        }
    }
    return TRUE;
}

static int stop_range_table_matches(const stop_range_table *t, const jabs_stop *stop, const jabs_stop *stragg, const ion *p, const sample *sample, size_t i_range, double emin) {
    if(t->isotope != p->isotope || t->ion_gsto != p->ion_gsto || t->stop_type != stop->type || t->stragg_type != stragg->type || t->emin != emin || t->n_isotopes != sample->n_isotopes) {
        return FALSE;
    }
    const sample_range *r = &sample->ranges[i_range];
    if(t->bragg != r->bragg || t->stragg != r->stragg) {
        return FALSE;
    }
    const double *c = sample_conc_bin(sample, i_range, 0);
    for(size_t i = 0; i < sample->n_isotopes; i++) {
        if(t->isotopes[i] != sample->isotopes[i] || t->conc[i] != c[i]) {
            return FALSE;
        }
    }
    return TRUE;
}

static stop_range_table *stop_range_table_make(const jabs_stop *stop, const jabs_stop *stragg, const ion *p, const sample *sample, size_t i_range, double emin) {
    /* Nodes are equally spaced in ln(E) from emin to the maximum energy of stopping data. Path length R and straggling integral Q are integrated between nodes
     * (Gauss-Legendre, three points). Between nodes they are interpolated with cubic Hermite polynomials using the derivatives (1/S and W/S^3) at the nodes,
     * likewise energy as a function of path length (derivative S). In constant composition E and u = (energy loss variance) / S^2 obey dE/dl = -S and
     * du/dl = W / S^2 (see stop_step()), therefore crossing path length l from E_in is E_out = R^-1(R(E_in) - l) and u_out = u_in + Q(E_in) - Q(E_out). */
    stop_range_table *t = calloc(1, sizeof(stop_range_table));
    if(!t) {
        return NULL;
    }
    t->isotope = p->isotope;
    t->ion_gsto = p->ion_gsto;
    t->stop_type = stop->type;
    t->stragg_type = stragg->type;
    t->emin = emin;
    t->n_isotopes = sample->n_isotopes;
    t->bragg = sample->ranges[i_range].bragg;
    t->stragg = sample->ranges[i_range].stragg;
    t->isotopes = malloc(t->n_isotopes * sizeof(jibal_isotope *));
    t->conc = malloc(t->n_isotopes * sizeof(double));
    if(!t->isotopes || !t->conc) {
        stop_range_table_free(t);
        return NULL;
    }
    const double *c = sample_conc_bin(sample, i_range, 0);
    for(size_t i = 0; i < t->n_isotopes; i++) {
        t->isotopes[i] = sample->isotopes[i];
        t->conc[i] = c[i];
    }
    t->valid = FALSE;
    const double emax = GSL_MIN_DBL(p->ion_gsto->emax, E_MAX);
    if(p->Z == 0 || emin <= 0.0 || emax <= emin) {
        return t;
    }
    const double ln_range = log(emax / emin);
    const size_t n = (size_t)ceil(ln_range / STOP_RANGE_TABLE_DLNE) + 1;
    const double dlne = ln_range / (double)(n - 1);
    t->E = malloc(n * sizeof(double));
    t->R = malloc(n * sizeof(double));
    t->S = malloc(n * sizeof(double));
    t->Q = malloc(n * sizeof(double));
    t->dQ = malloc(n * sizeof(double));
    if(!t->E || !t->R || !t->S || !t->Q || !t->dQ) {
        stop_range_table_free(t);
        return NULL;
    }
    const depth d_mid = {.x = (sample->ranges[i_range].x + sample->ranges[i_range + 1].x) / 2.0, .i = i_range};
    stop_sw sw;
    for(size_t i = 0; i < n; i++) {
        t->E[i] = (i == n - 1) ? emax : emin * exp(dlne * (double)i);
        stop_sample_sw(stop, stragg, p, sample, d_mid, t->E[i], &sw);
        if(sw.S < STOP_STEP_MINIMUM_STOPPING) {
            return t; /* Not valid */
        }
        t->S[i] = sw.S;
        t->dQ[i] = sw.W / (sw.S * sw.S * sw.S);
    }
    static const double gl_x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834}; /* sqrt(3/5) */
    static const double gl_w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    t->R[0] = 0.0;
    t->Q[0] = 0.0;
    for(size_t i = 0; i < n - 1; i++) {
        const double E_mid = (t->E[i] + t->E[i + 1]) / 2.0, half = (t->E[i + 1] - t->E[i]) / 2.0;
        double sum_R = 0.0, sum_Q = 0.0;
        for(int k = 0; k < 3; k++) {
            stop_sample_sw(stop, stragg, p, sample, d_mid, E_mid + half * gl_x[k], &sw);
            if(sw.S < STOP_STEP_MINIMUM_STOPPING) {
                return t;
            }
            sum_R += gl_w[k] / sw.S;
            sum_Q += gl_w[k] * sw.W / (sw.S * sw.S * sw.S);
        }
        t->R[i + 1] = t->R[i] + half * sum_R;
        t->Q[i + 1] = t->Q[i] + half * sum_Q;
    }
    t->n = n;
    t->emax = emax;
    t->ln_emin = log(emin);
    t->inv_dlne = 1.0 / dlne;
    t->valid = TRUE;
    DEBUGMSG("Range table of %s in range %zu made, %zu nodes from %g keV to %g keV.", p->isotope->name, i_range, n, emin / C_KEV, emax / C_KEV);
    return t;
}

static const stop_range_table *stop_range_table_get(const jabs_stop *stop, const jabs_stop *stragg, const ion *p, const sample *sample, size_t i_range, double emin) {
    stop_range_tables *rt = stop->range_tables;
    for(size_t i = 0; i < rt->n; i++) {
        if(stop_range_table_matches(rt->t[i], stop, stragg, p, sample, i_range, emin)) {
            return rt->t[i];
        }
    }
    stop_range_table *t = stop_range_table_make(stop, stragg, p, sample, i_range, emin);
    if(!t) {
        return NULL;
    }
    stop_range_table **tables = realloc(rt->t, (rt->n + 1) * sizeof(stop_range_table *));
    if(!tables) {
        stop_range_table_free(t);
        return NULL;
    }
    rt->t = tables;
    rt->t[rt->n] = t;
    rt->n++;
    return t;
}

static inline double stop_hermite(double y0, double y1, double d0, double d1, double h, double u) { /* Cubic Hermite polynomial, values y and derivatives d at the ends of an interval of width h, u in [0, 1] */
    const double u2 = u * u, u3 = u2 * u;
    return y0 * (2.0 * u3 - 3.0 * u2 + 1.0) + y1 * (3.0 * u2 - 2.0 * u3) + h * (d0 * (u3 - 2.0 * u2 + u) + d1 * (u3 - u2));
}

static size_t stop_range_table_index(const stop_range_table *t, double E) { /* Index of node at or below E, E in [emin, emax] */
    double x = (log(E) - t->ln_emin) * t->inv_dlne;
    size_t i = x > 0.0 ? (size_t)x : 0;
    if(i > t->n - 2) {
        i = t->n - 2;
    }
    return i;
}

static void stop_range_table_eval(const stop_range_table *t, double E, double *R, double *Q) {
    const size_t i = stop_range_table_index(t, E);
    const double h = t->E[i + 1] - t->E[i];
    const double u = (E - t->E[i]) / h;
    *R = stop_hermite(t->R[i], t->R[i + 1], 1.0 / t->S[i], 1.0 / t->S[i + 1], h, u);
    *Q = stop_hermite(t->Q[i], t->Q[i + 1], t->dQ[i], t->dQ[i + 1], h, u);
}

static double stop_range_table_energy(const stop_range_table *t, double R) { /* Inverse of path length R(E), R in [0, R[n - 1]] */
    size_t lo = 0, hi = t->n - 1;
    while(hi - lo > 1) {
        size_t mi = (lo + hi) / 2;
        if(R >= t->R[mi]) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    const double h = t->R[lo + 1] - t->R[lo];
    return stop_hermite(t->E[lo], t->E[lo + 1], t->S[lo], t->S[lo + 1], h, (R - t->R[lo]) / h);
}

static int stop_range_table_cross(const jabs_stop *stop, const jabs_stop *stragg, ion *p, const sample *sample, depth *d, double emin, double *remaining) {
    /* Moves exiting ion p from depth d across the range to the next range boundary (towards the surface it is heading to) using range tables. Returns 1 if the ion
     * crossed the range (d is updated), -1 if it stopped in the range (path length it still had to go is stored in remaining, if not NULL) and 0 if tables
     * can not be used (composition is not constant, energy is outside table), then the ion must be stepped. */
    if(p->E <= emin) {
        return 0;
    }
    const depth d_next = stop_next_crossing(p, sample, d);
    const size_t i_range = d_next.i;
    if(!sample_range_conc_constant(sample, i_range)) {
        return 0;
    }
    const stop_range_table *t = stop_range_table_get(stop, stragg, p, sample, i_range, emin);
    if(!t || !t->valid || p->E > t->emax) {
        return 0;
    }
    const depth d_mid = {.x = (sample->ranges[i_range].x + sample->ranges[i_range + 1].x) / 2.0, .i = i_range};
    const double l = (d_next.x - d->x) * p->inverse_cosine_theta; /* Path length to range boundary, positive */
    double R_in, Q_in;
    stop_range_table_eval(t, p->E, &R_in, &Q_in);
    const double u_in = p->S / pow2(stop_sample(stop, p, sample, d_mid, p->E));
    if(R_in <= l) { /* Stops in this range */
        if(remaining) {
            const depth d_stop = {.x = d->x + R_in * p->cosine_theta, .i = i_range};
            *remaining = stop_remaining_path(p, sample, d_stop);
        }
        p->S = (u_in + Q_in) * pow2(t->S[0]);
        p->E = emin;
        return -1;
    }
    const double E_out = stop_range_table_energy(t, R_in - l);
    double R_out, Q_out;
    stop_range_table_eval(t, E_out, &R_out, &Q_out);
    p->S = GSL_MAX_DBL(0.0, u_in + Q_in - Q_out) * pow2(stop_sample(stop, p, sample, d_mid, E_out));
    p->E = E_out;
    *d = d_next;
    return 1;
}

int stop_sample_exit_remaining(const jabs_stop *stop, const jabs_stop *stragg, const jabs_stop_step_params *params_exiting, ion *p, const depth depth_start, const sample *sample, double *remaining) {
    depth d = depth_start;
    if(remaining) {
        *remaining = 0.0;
    }
    int last = FALSE;
    double emin = GSL_MAX_DBL(stop->emin, p->ion_gsto->emin);
    const int adaptive = stop->rk4 && stop->accuracy > 0.0;
    stop_adaptive state = {0};
    while(1) { /* Exit from sample (hopefully) */
        if(p->inverse_cosine_theta > 0.0 && d.x >= (sample->thickness - DEPTH_TOLERANCE)) { /* Exit through back (transmission) */
            return 0;
        }
        if(p->inverse_cosine_theta < 0.0 && d.x <= DEPTH_TOLERANCE) { /* Exit (surface, front of sample) */
            return 0;
        }
        if(last) {
            DEBUGMSG("Last step taken, E = %g keV, depth still %.3lf tfu, break break.", p->E / C_KEV, d.x / C_TFU);
            break;
        }
        if(stop->range_tables && !stop->ms) {
            int crossed = stop_range_table_cross(stop, stragg, p, sample, &d, emin, remaining);
            if(crossed > 0) {
                continue;
            } else if(crossed < 0) {
                return -1;
            }
        }
        double E_step = stop_step_calc(params_exiting, p);
        if(adaptive && state.h > 0.0) { /* First adaptive step is the fixed step, then step is limited by error control and maximum step */
            E_step = stop_step_max(params_exiting, p);
        }

        if(p->E - E_step <= emin) {
            if(p->E <= emin) {
                DEBUGMSG("Energy %g keV is below %g keV. Assuming ion stops.\n", p->E / C_KEV, emin / C_KEV);
                if(remaining) {
                    *remaining = stop_remaining_path(p, sample, d);
                }
                return -1;
            }
            E_step = p->E - emin;
            DEBUGMSG("Energy %g keV will go below emin = %g keV if we take the suggested step. Taking one final step of (approximately) %g keV at depth %g tfu and seeing what happens then.",
                            p->E / C_KEV, emin / C_KEV, E_step / C_KEV, d.x / C_TFU);
            last = TRUE;
        }
        depth d_after = adaptive ? stop_step_adaptive(stop, stragg, p, sample, d, E_step, &state) : stop_step(stop, stragg, p, sample, d, E_step);
        if(last && p->E - emin > STOP_LAST_STEP_REMAINDER * E_step) { /* Last step was truncated at range boundary or energy loss was smaller than expected (stopping decreases with energy) */
            DEBUGMSG("Last step ended at depth %g tfu with E = %g keV, which is still above emin = %g keV, continuing.", d_after.x / C_TFU, p->E / C_KEV, emin / C_KEV);
            last = FALSE;
        }
        d = d_after;
    }
    if(remaining) {
        *remaining = stop_remaining_path(p, sample, d);
    }
    return -1; /* Stopped inside the sample */
}
