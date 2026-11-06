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

static double stop_gsto(const ion *incident, gsto_stopping_type type, int Z2, double em) { /* Electronic stopping or straggling of incident ion in element Z2, see ion_gsto_table_init() */
    const jabs_ion_gsto_data *gd = &incident->ion_gsto->gsto_data[Z2];
    const jabs_gsto_table *t = (type == GSTO_STO_STRAGG) ? &gd->stragg : &gd->stop;
    assert(t->file);
    return ion_gsto_table_eval(t, em);
}

static void stop_sample_sw(const jabs_stop *stop, const jabs_stop *stragg, const ion *incident, const sample *sample, const depth depth, double E, stop_sw *sw) { /* Stopping, straggling and stopping gradient in one go */
    sw->S = 0.0;
    sw->W = 0.0;
    sw->S_x = 0.0;
    if(incident->Z == 0) {
        return;
    }
    const double em = E * incident->mass_inverse;
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
            S_Z2 = stop_gsto(incident, stop->type, Z2, em);
            W_Z2 = stop_gsto(incident, stragg->type, Z2, em);
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
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        const jibal_isotope *target = sample->isotopes[i_isotope];
        double c = get_conc(sample, depth, i_isotope);
        if(c < CONC_TOLERANCE)
            continue;
        int Z2 = target->Z;
        assert(Z2 <= incident->ion_gsto->Z2_max);
        if(Z2 != Z2_old) { /* Electronic stopping (or straggling) is the same for all isotopes of the same element */
            S_Z2 = stop_gsto(incident, stop->type, Z2, em);
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
    depth d = depth_start;
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
        double E_step = stop_step_calc(params_exiting, p);
        if(adaptive && state.h > 0.0) { /* First adaptive step is the fixed step, then step is limited by error control and maximum step */
            E_step = stop_step_max(params_exiting, p);
        }

        if(p->E - E_step <= emin) {
            if(p->E <= emin) {
                DEBUGMSG("Energy %g keV is below %g keV. Assuming ion stops.\n", p->E / C_KEV, emin / C_KEV);
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
    return -1; /* Stopped inside the sample */
}
