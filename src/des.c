/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <math.h>
#include <jibal_cross_section.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include "jabs_debug.h"
#include "defaults.h"
#include "des.h"
#include "stop.h"

extern inline const des *des_table_min_energy_bin(const des_table *dt);

des_table *des_table_init(size_t n) {
    des_table *dt = malloc(sizeof(des_table));
    dt->t = NULL;
    dt->depth_interval_index = NULL;
    dt->n = 0;
    dt->n_ranges = 0;
    dt->depth_increases = TRUE;
    if(des_table_realloc(dt, n)) {
        free(dt);
        return NULL;
    }
    return dt;
}
int des_table_realloc(des_table *dt, size_t n) {
    if(n > DES_TABLE_MAX_SIZE) {
        DEBUGMSG("DES table requested size %zu larger than allowed %i", n, DES_TABLE_MAX_SIZE);
        free(dt->t);
        dt->t = NULL;
    } else {
        DEBUGMSG("DES table %p realloc to %zu elements", (void *) dt, n);
        dt->t = realloc(dt->t, n * sizeof(des));
    }
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
    if(dt->depth_interval_index) {
        free(dt->depth_interval_index);
    }
    dt->depth_interval_index = calloc(dt->n_ranges + 1, sizeof(size_t));

    size_t i_range_old = dt->t[0].d.i;
    dt->depth_interval_index[i_range_old] = 0;
    if(dt->depth_increases) {
        for(size_t i = 1; i < dt->n; i++) {
            const des *des = &(dt->t[i]);
            if(des->d.i > i_range_old) { /* index increases (des table has increasing depth) */
                for(size_t i_range = i_range_old + 1; i_range <= des->d.i && i_range < dt->n_ranges; i_range++) { /* Handles indices step skips over (no thickness between them) */
                    dt->depth_interval_index[i_range] = i - 1;
                }
                i_range_old = des->d.i;
            }
        }
        dt->depth_interval_index[dt->n_ranges] = dt->n - 1; /* Last point, last index */
    } else {
        for(size_t i = 1; i < dt->n; i++) {
            const des *des = &(dt->t[i]);
            if(des->d.i < i_range_old) {
                for(size_t i_range = des->d.i; i_range < i_range_old; i_range++) {
                    dt->depth_interval_index[i_range + 1] = i;
                }
                i_range_old = des->d.i;
            }
        }
        dt->depth_interval_index[i_range_old] = dt->n - 1;
    }
}

void des_table_print(FILE *f, const des_table *dt) {
    if(!dt) {
        return;
    }
    fprintf(f, "DES    i d.i        d.x        d.E      d.S\n");
    for(size_t i = 0; i < dt->n; i++) {
        const des *des = &(dt->t[i]);
        fprintf(f, "DES %4zu %3zu %10.3lf %10.3lf %8.3lf\n", i, des->d.i, des->d.x / C_TFU, des->E / C_KEV, sqrt(des->S) / C_KEV);
    }
    fprintf(f, "DES DEPTH INCREASES: %s\n", dt->depth_increases?"TRUE":"FALSE");
    if(dt->depth_interval_index) {
        for(size_t i = 0; i < dt->n_ranges; i++) {
            fprintf(f, "DES INDEX %3zu [%4zu, %4zu]\n", i, dt->depth_interval_index[i], dt->depth_interval_index[i + 1]);
        }
    } else {
        fprintf(f, "DES INDEX DOES NOT EXIST\n");
    }
}

depth des_table_find_depth(const des_table *dt, size_t *i_des, depth depth_prev, ion *incident) {
    size_t i;
    assert(dt);
    assert(dt->n > 0);
    double E = incident->E;
    DEBUGVERBOSEMSG("Where is E = %g keV in DES table? Start index %zu. Previous depth %g tfu (range %zu)",
            E / C_KEV, *i_des, depth_prev.x / C_TFU, depth_prev.i);
    for(i = *i_des; i < dt->n - 1; i++) {
        if(dt->t[i].d.i != dt->t[i + 1].d.i) { /* This means last point of this layer */
            E = dt->t[i].E;
            *i_des = i + 1; /* The +1 prevents stopping at this same layer boundary on the next call */
            DEBUGVERBOSEMSG("Layer boundary at %g tfu. Setting E = %g keV and i_des = %zu. index = %zu", dt->t[i].d.x / C_TFU, E / C_KEV, *i_des, dt->depth_interval_index[dt->t[i + 1].d.i]);
            break;
        }
        if(dt->t[i].E < E) {/* i is the index of the first element in dt->t that has energy below E. So i-1 should have energy above E. */
            DEBUGVERBOSEMSG("Breaking, i_des = i = %zu, because first element %.12e J (%g keV) <  E = %.12e J (%g keV)", i, dt->t[i].E, dt->t[i].E / C_KEV,  E, E / C_KEV);
            *i_des = i;
            break;
        }
    }
    if(i == dt->n - 1) {
        if(E < dt->t[i].E) {
            DEBUGVERBOSEMSG("Energy %g keV is below last point in table (%g keV). Changing energy.", E / C_KEV, dt->t[i].E / C_KEV);
            E = dt->t[i].E;
        }
    } else if(i == 0) {
        if(E > dt->t[i].E) {
            DEBUGVERBOSEMSG("Energy %g keV is above first point in table (%g keV). Changing energy.", E / C_KEV, dt->t[i].E / C_KEV);
            E = dt->t[i].E;
        }
        i = 1;
    }
    assert(i > 0);
    const des *des_low = &dt->t[i - 1]; /* Closer to surface, higher energy */
    const des *des_high = &dt->t[i]; /* Deeper, lower energy */
    DEBUGVERBOSEMSG("dt->t[%zu] = %g tfu (i = %zu), dt->t[%zu] = %g tfu (i = %zu), depth_prev = %g tfu (i = %zu)",
            i - 1,  des_low->d.x / C_TFU, des_low->d.i, i, des_high->d.x / C_TFU, des_high->d.i, depth_prev.x / C_TFU, depth_prev.i);
    double E_diff = E - des_low->E;
    double E_interval = des_high->E - des_low->E;
    double S_interval = des_high->S - des_low->S;
    double d_interval = depth_diff(des_low->d, des_high->d);
    double frac;
    if(fabs(E_interval) < 0.1 * C_EV) { /* Prevent div by zero */
        frac = 0.0;
    } else {
        frac = (E_diff / E_interval); /* zero if close to low bin */
    }
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
    incident->ms.w = des_low->ms.w + frac * (des_high->ms.w - des_low->ms.w);
    incident->ms.sum_a = des_low->ms.sum_a + frac * (des_high->ms.sum_a - des_low->ms.sum_a);
    incident->ms.sum_b = des_low->ms.sum_b + frac * (des_high->ms.sum_b - des_low->ms.sum_b);
    incident->ms.I = des_low->ms.I + frac * (des_high->ms.I - des_low->ms.I);
    incident->ms.P = des_low->ms.P + frac * (des_high->ms.P - des_low->ms.P);
    incident->flux = des_low->flux + frac * (des_high->flux - des_low->flux);
    DEBUGVERBOSEMSG("d_out = %g tfu (i = %zu), E = %g keV, S = %g keV", d_out.x / C_TFU, d_out.i, incident->E / C_KEV, sqrt(incident->S) * C_FWHM / C_KEV);
    return d_out;
}

double des_total_cross_section(const jibal_isotope *incident, const jibal_isotope *target, double theta_min, double E) {
    const int n = DES_BEAM_LOSS_ANGLE_STEPS;
    const double u_min = log(theta_min), u_max = log(C_PI);
    const double du = (u_max - u_min) / n;
    double sum = 0.0;
    for(int i = 0; i < n; i++) { /* Midpoint rule in u = ln(theta), d(Omega) = 2 pi sin(theta) theta du */
        double theta = exp(u_min + (i + 0.5) * du);
        sum += jibal_cross_section_rbs(incident, target, theta, E, JIBAL_CS_ANDERSEN) * 2.0 * C_PI * sin(theta) * theta * du;
    }
    return sum;
}

static double des_beam_loss(const sample *sample, const ion *ion, const depth d, double E, double path, double theta_min) { /* Probability of scattering to counted dual scattering angles over path length */
    double sum = 0.0;
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        double c = get_conc(sample, d, i_isotope);
        if(c < CONC_TOLERANCE) {
            continue;
        }
        sum += c * des_total_cross_section(ion->isotope, sample->isotopes[i_isotope], theta_min, E);
    }
    return sum * path;
}

des_table *des_table_compute(const jabs_stop *stop, const jabs_stop *stragg, const sim_calc_params *scp, const sample *sample, const ion *incident, depth depth_start, double emin, double beam_loss_angle_min) {
    ion ion = *incident;
    ion_ms_reset(&ion.ms); /* Path of incident ion starts here */
    ion.flux = 1.0;
    des_table *dt = des_table_init(DES_TABLE_INITIAL_ALLOC);
    if(!dt) {
        return NULL;
    }
    size_t i = 0;
    depth d_before;
    depth d_after = depth_start;
    dt->depth_increases = (ion.inverse_cosine_theta > 0.0);
    assert(ion.inverse_cosine_theta <= -1.0 || ion.inverse_cosine_theta >= 1.0);
    DEBUGMSG("Computing DES table, incident ion is %s. start_depth = %g tfu, E = %g keV, angle in sample %g deg (1/cos = %.6lf).", incident->isotope->name, depth_start.x / C_TFU, ion.E / C_KEV,
            ion.theta / C_DEG, ion.inverse_cosine_theta);
    if(incident->ion_gsto->emin < emin) {
        emin = incident->ion_gsto->emin;
        DEBUGMSG("DES table emin set from ion emin = %g keV", emin / C_KEV);
    }
    do {
        if(i == dt->n) {
            DEBUGMSG("DES table reallocation, size %zu reached when E = %g keV.", dt->n, ion.E / C_KEV);
            if(des_table_realloc(dt, dt->n * 2)) {
                break;
            }
        }
        d_before = d_after;
        des *des = &dt->t[i];
        des->E = ion.E;
        des->S = ion.S;
        des->ms = ion.ms;
        des->flux = ion.flux;
        des->d = d_before;
        if((ion.inverse_cosine_theta > 0.0 && d_before.x >= sample->thickness - DEPTH_TOLERANCE) || (ion.inverse_cosine_theta < 0.0 && d_before.x < DEPTH_TOLERANCE)) {
            DEBUGMSG("DES table calculation stops at %g tfu (i = %zu).", d_before.x / C_TFU, i);
            i++;
            break;
        }
        double E_before = ion.E;
        d_after = stop_step(stop, stragg, &ion, sample, d_before, stop_step_calc(&scp->incident_stop_params, &ion));
        if(beam_loss_angle_min > 0.0) {
            depth d_half = {.x = (d_before.x + d_after.x) / 2.0, .i = d_after.i};
            double path = fabs(d_after.x - d_before.x) * fabs(ion.inverse_cosine_theta);
            ion.flux *= exp(-des_beam_loss(sample, &ion, d_half, (E_before + ion.E) / 2.0, path, beam_loss_angle_min));
        }
        i++;
    } while(ion.E > emin);
    if(dt->n) {
        des_table_realloc(dt, i); /* Shrinks to size */
        if(i > 0 ) {
            dt->n_ranges = GSL_MAX(depth_start.i, dt->t[i - 1].d.i) + 1;
        } else {
            dt->n_ranges = 0;
        }
        des_table_rebuild_index(dt);
        return dt;
    } else {
        DEBUGSTR("DES table unused, freeing.");
        des_table_free(dt);
        return NULL;
    }
}

void des_set_ion(const des *des, ion *ion) {
    ion->E = des->E;
    ion->S = des->S;
    ion->ms = des->ms;
    ion->flux = des->flux;
}

depth des_next_range(const des_table *dt, ion *incident, depth d) {
    assert(d.i < dt->n_ranges);
    size_t i_range_next;
    if(incident->inverse_cosine_theta > 0.0) { /* Going deeper */
        i_range_next = d.i + 1;
    } else {
        if(d.i == 0) {
            i_range_next = 0;
        } else {
            i_range_next = d.i - 1;
        }
    }
    size_t i_des_skip = dt->depth_interval_index[i_range_next];
    assert(i_des_skip < dt->n);
    des *des_skip = &(dt->t[i_des_skip]);
    incident->E = des_skip->E; /* TODO: long skips may make E_deriv calculation inaccurate */
    incident->S = des_skip->S;
    incident->ms = des_skip->ms;
    incident->flux = des_skip->flux;
    return des_skip->d;
}
