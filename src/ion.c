/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_minmax.h>
#include <jibal_stragg.h>
#include "defaults.h"
#include "ion.h"
#include "rotate.h"
#include "jabs_debug.h"

void ion_reset(ion *ion) {
    ion->isotope = NULL;
    ion->E = 0.0;
    ion->S = 0.0;
    ion->Z = 0;
    ion->mass = 0.0;
    ion->theta = 0.0;
    ion->phi = 0.0;
    ion->cosine_theta = 1.0;
    ion->cosine_phi = 1.0;
    ion->inverse_cosine_theta = 1.0; /* These need to be set matching to the angles */
    ion->inverse_cosine_phi = 1.0;
    ion->nucl_stop = NULL; /* Beware of memory leaks! */
    ion->ion_gsto = NULL; /* Beware of memory leaks! */
    ion_ms_reset(&ion->ms);
    ion->flux = 1.0;
}

void ion_ms_reset(ion_ms *ms) {
    ms->w = 0.0;
    ms->sum_a = 0.0;
    ms->sum_b = 0.0;
    ms->I = 0.0;
    ms->P = 1.0;
}

void ion_set_isotope(ion *ion, const jibal_isotope *isotope) {
    if(!isotope)
        return;
    ion->isotope = isotope;
    ion->mass = isotope->mass;
    ion->mass_inverse = 1.0/isotope->mass;
    ion->Z = isotope->Z;
}

void ion_set_angle(ion *ion, double theta, double phi) {
    if(!ion)
        return;
    if(theta == ion->theta && phi == ion->phi) {
        return;
    }
    ion->theta = fmod(theta, C_2PI);
    if(ion->theta > C_PI) { /* Limit theta to [0, pi] */
        ion->theta = C_2PI - ion->theta;
        phi += C_PI;
    }
    ion->phi = fmod(phi, C_2PI);
    ion->cosine_theta = cos(ion->theta);
    ion->cosine_phi = cos(ion->phi);
    ion->inverse_cosine_theta = 1.0/ion->cosine_theta;
    ion->inverse_cosine_phi = 1.0/ion->cosine_phi;
}

double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, double E) {
#ifdef DEBUG
    assert(isotope->i < ion->nucl_stop->n_isotopes);
    assert(ion->nucl_stop->t[isotope->i].target == isotope);
#endif
    nucl_stop_pair x = ion->nucl_stop->t[isotope->i];
    return x.k * nuclear_stopping_universal(ion->nucl_stop, x.eps0 * E);
}

void ion_rotate(ion *ion, double theta2, double phi2) { /* Wrapper for rotate() */
    double theta, phi;
    rotate(theta2, phi2, ion->theta, ion->phi, &theta, &phi);
    ion_set_angle(ion, theta, phi);
}

void ion_print(FILE *f, const ion *ion) {
    fprintf(f, "ion %s (Z=%i, mass=%.3lf u), E = %.3lf keV, theta = %.3lf deg (cos_theta = %.3lf, 1/cos_theta = %.3lf), phi = %.3lf deg. Nuclear stopping isotopes %zu calculated for %s.\n",
            ion->isotope->name, ion->Z, ion->mass/C_U, ion->E/C_KEV,  ion->theta/C_DEG, ion->cosine_theta, ion->inverse_cosine_theta, ion->phi/C_DEG, ion->nucl_stop->n_isotopes, ion->nucl_stop->incident->name);
}

void ion_gsto_table_init(jabs_gsto_table *t, const gsto_file_t *file, const double *data, double unit_factor) {
    /* Interpolation of stopping or straggling tables uses a natural cubic spline in ln(em), evaluated as a cubic Hermite polynomial with derivatives at nodes from
     * the spline. Unlike linear interpolation, this has a continuous first and second derivative with respect to energy, which is important for the accuracy of
     * Runge-Kutta integration (see stop_step()) when steps are longer than table node spacing. Near steps in the data (e.g. SRIM2013 for some Z1, Z2) and where
     * data is flat, derivatives are limited to prevent overshoot, then the second derivative is not continuous at those nodes. Tables with only two nodes or
     * nodes at zero energy are interpolated linearly in em. */
    t->file = file;
    t->data = data;
    t->unit_factor = unit_factor;
    t->lnem = NULL;
    t->deriv = NULL;
    t->lnem0 = 0.0;
    t->inv_dlnem = 0.0;
    size_t n = file->xpoints;
    if(n < 3 || file->em[0] <= 0.0) {
        return;
    }
    t->lnem = malloc(n * sizeof(double));
    t->deriv = malloc(n * sizeof(double));
    for(size_t i = 0; i < n; i++) {
        t->lnem[i] = log(file->em[i]);
    }
    t->lnem0 = t->lnem[0];
    double dlnem = (t->lnem[n - 1] - t->lnem[0]) / (double)(n - 1);
    int uniform = TRUE;
    for(size_t i = 0; i < n; i++) {
        if(fabs(t->lnem[i] - (t->lnem0 + dlnem * (double)i)) > 1e-6 * dlnem) {
            uniform = FALSE;
            break;
        }
    }
    if(uniform) {
        t->inv_dlnem = 1.0 / dlnem;
    }
    double *M = calloc(n, sizeof(double)); /* Second derivatives at nodes, natural spline: M[0] = M[n-1] = 0 */
    double *c = malloc(n * sizeof(double)); /* Tridiagonal solver workspace */
    for(size_t i = 1; i < n - 1; i++) { /* Thomas algorithm: h_(i-1) M_(i-1) + 2 (h_(i-1) + h_i) M_i + h_i M_(i+1) = 6 (slope_i - slope_(i-1)) */
        double h0 = t->lnem[i] - t->lnem[i - 1];
        double h1 = t->lnem[i + 1] - t->lnem[i];
        double r = 6.0 * ((data[i + 1] - data[i]) / h1 - (data[i] - data[i - 1]) / h0);
        double b = 2.0 * (h0 + h1) - h0 * (i > 1 ? c[i - 1] : 0.0);
        c[i] = h1 / b;
        M[i] = (r - h0 * (i > 1 ? M[i - 1] : 0.0)) / b;
    }
    for(size_t i = n - 2; i >= 1; i--) {
        M[i] -= c[i] * M[i + 1];
    }
    for(size_t i = 0; i < n - 1; i++) {
        double h = t->lnem[i + 1] - t->lnem[i];
        t->deriv[i] = (data[i + 1] - data[i]) / h - h * (2.0 * M[i] + M[i + 1]) / 6.0;
    }
    double h = t->lnem[n - 1] - t->lnem[n - 2];
    t->deriv[n - 1] = (data[n - 1] - data[n - 2]) / h + h * (M[n - 2] + 2.0 * M[n - 1]) / 6.0;
    free(M);
    free(c);
    for(size_t i = 0; i < n; i++) { /* Hyman filter: where data is monotonic, derivatives are limited so that interpolation is monotonic too. Prevents overshoot near steps in data. */
        double d_left = i > 0 ? (data[i] - data[i - 1]) / (t->lnem[i] - t->lnem[i - 1]) : 0.0;
        double d_right = i < n - 1 ? (data[i + 1] - data[i]) / (t->lnem[i + 1] - t->lnem[i]) : 0.0;
        double d_min;
        if(i == 0) {
            d_min = d_right;
        } else if(i == n - 1) {
            d_min = d_left;
        } else if(d_left * d_right > 0.0) {
            d_min = fabs(d_left) < fabs(d_right) ? d_left : d_right;
        } else if(d_left * d_right == 0.0) { /* Flat on one side */
            d_min = 0.0;
        } else {
            continue; /* Extremum, spline derivative is not changed */
        }
        if(t->deriv[i] * d_min < 0.0) {
            t->deriv[i] = 0.0;
        } else if(fabs(t->deriv[i]) > 3.0 * fabs(d_min)) {
            t->deriv[i] = 3.0 * d_min;
        }
    }
}

void ion_gsto_table_free(jabs_gsto_table *t) {
    free(t->lnem);
    free(t->deriv);
}

double ion_gsto_table_eval(const jabs_gsto_table *t, double em) { /* Interpolated value (SI units) at em (energy per mass), see ion_gsto_table_init() */
    const gsto_file_t *file = t->file;
    const double *data = t->data;
    const size_t n = file->xpoints;
    if(em <= file->em[0]) {
        DEBUGMSG("Energy %g keV/u below GSTO file minimum %g keV/u, returning value corresponding to minimum", em / (C_KEV / C_U), file->em[0] / (C_KEV / C_U));
        return t->unit_factor * data[0];
    }
    if(em >= file->em[n - 1]) {
        DEBUGMSG("Energy %g keV/u above GSTO file maximum %g keV/u, returning value corresponding to maximum", em / (C_KEV / C_U), file->em[n - 1] / (C_KEV / C_U));
        return t->unit_factor * data[n - 1];
    }
    if(!t->deriv) {
        int lo = jibal_gsto_em_to_index(file, em);
        lo = GSL_MIN_INT(GSL_MAX_INT(lo, 0), (int)n - 2);
        return t->unit_factor * jibal_linear_interpolation(file->em[lo], file->em[lo + 1], data[lo], data[lo + 1], em);
    }
    const double x = log(em);
    size_t lo;
    if(t->inv_dlnem > 0.0) {
        lo = (size_t)((x - t->lnem0) * t->inv_dlnem);
    } else { /* Binary search */
        size_t hi = n - 1;
        lo = 0;
        while(hi - lo > 1) {
            size_t mi = (lo + hi) / 2;
            if(x >= t->lnem[mi]) {
                lo = mi;
            } else {
                hi = mi;
            }
        }
    }
    if(lo > n - 2) { /* Rounding */
        lo = n - 2;
    }
    const double h = t->lnem[lo + 1] - t->lnem[lo];
    const double u = (x - t->lnem[lo]) / h;
    const double u2 = u * u, u3 = u2 * u;
    return t->unit_factor * (data[lo] * (2.0 * u3 - 3.0 * u2 + 1.0) + data[lo + 1] * (3.0 * u2 - 2.0 * u3) + h * (t->deriv[lo] * (u3 - 2.0 * u2 + u) + t->deriv[lo + 1] * (u3 - u2)));
}

jabs_ion_gsto *ion_gsto_new(const jibal_isotope *incident, const jibal_gsto *gsto) {
    jabs_ion_gsto *ig = malloc(sizeof(jabs_ion_gsto));
    ig->incident = incident;
    ig->gsto_data = calloc(gsto->Z2_max + 1, sizeof(jabs_ion_gsto_data));
    ig->emin = 0.0;
    ig->emax = E_MAX;
    ig->refcount = 1;
    ig->Z2_max = gsto->Z2_max;
    if(incident->Z == 0) { /* Leave everything empty for neutrons */
        return ig;
    }
    DEBUGMSG("New gsto data for incident %s. Z2 max is %i", incident->name, gsto->Z2_max);
    for(int Z2 = 1; Z2 <= gsto->Z2_max; Z2++) {
        jabs_ion_gsto_data *gd = &ig->gsto_data[Z2];
        const gsto_file_t *stopfile = jibal_gsto_get_assigned_file(gsto, GSTO_STO_ELE, incident->Z, Z2);
        if(stopfile) {
            ion_gsto_table_init(&gd->stop, stopfile, jibal_gsto_file_get_data(stopfile, incident->Z, Z2), stopfile->stounit == GSTO_STO_UNIT_EV15CM2 ? C_EV_TFU : 1.0);
            double emin = jibal_gsto_xunit_to_energy(stopfile->xunit, stopfile->xmin, incident->mass);
            double emax = jibal_gsto_xunit_to_energy(stopfile->xunit, stopfile->xmax, incident->mass);
            DEBUGMSG("Stop data Z1 = %i, Z2 = %i, file: %s, emin = %g keV, emax %g keV", incident->Z, Z2, stopfile->name, emin / C_KEV, emax / C_KEV);
            ig->emin = GSL_MAX_DBL(ig->emin, emin);
            if(emax > 0.0) {
                ig->emax = GSL_MIN_DBL(ig->emax, emax);
            }
        } else {
            DEBUGVERBOSEMSG("No stopfile, Z1 = %i, Z2 = %i", incident->Z, Z2);
        }
        const gsto_file_t *straggfile = jibal_gsto_get_assigned_file(gsto, GSTO_STO_STRAGG, incident->Z, Z2);
        if(straggfile) {
            ion_gsto_table_init(&gd->stragg, straggfile, jibal_gsto_file_get_data(straggfile, incident->Z, Z2), straggfile->straggunit == GSTO_STRAGG_UNIT_BOHR ? jibal_stragg_bohr(incident->Z, Z2) : 1.0);
            double emin = jibal_gsto_xunit_to_energy(straggfile->xunit, straggfile->xmin, incident->mass);
            double emax = jibal_gsto_xunit_to_energy(straggfile->xunit, straggfile->xmax, incident->mass);
            DEBUGMSG("Stragg data Z1 = %i, Z2 = %i, file: %s, emin = %g keV, emax %g keV", incident->Z, Z2, straggfile->name, emin / C_KEV, emax / C_KEV);
            ig->emin = GSL_MAX_DBL(ig->emin, emin);
            if(emax > 0.0) {
                ig->emax = GSL_MIN_DBL(ig->emax, emax);
            }
        }  else {
            DEBUGVERBOSEMSG("No straggfile, Z1 = %i, Z2 = %i", incident->Z, Z2);
        }
    }
    DEBUGMSG("Ion %s stop E range [%g keV, %g keV]", incident->name, ig->emin / C_KEV, ig->emax / C_KEV);
    return ig;
}

jabs_ion_gsto *ion_gsto_shared(jabs_ion_gsto *ion_gsto) {
    if(!ion_gsto) {
        return NULL;
    }
    assert(ion_gsto->refcount > 0);
    ion_gsto->refcount++;
    DEBUGMSG("Shallow copy of ion_gsto = %p made (incident: %s), refcount is %i", (void *)ion_gsto, ion_gsto->incident->name, ion_gsto->refcount);
    return ion_gsto;
}

void ion_gsto_free(jabs_ion_gsto *ion_gsto) {
    if(!ion_gsto) {
        return;
    }
    ion_gsto->refcount--;
    DEBUGVERBOSEMSG("ion_gsto_free(ion_gsto = %p), incident %s, called, refcount is now %i", (void *)ion_gsto, ion_gsto->incident->name, ns->refcount);
    assert(ion_gsto->refcount >= 0);
    if(ion_gsto->refcount == 0) {
        for(int Z2 = 0; Z2 <= ion_gsto->Z2_max; Z2++) {
            ion_gsto_table_free(&ion_gsto->gsto_data[Z2].stop);
            ion_gsto_table_free(&ion_gsto->gsto_data[Z2].stragg);
        }
        free(ion_gsto->gsto_data);
        free(ion_gsto);
    }
}
