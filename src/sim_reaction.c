/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <string.h>
#include <assert.h>
#include "spectrum.h"
#include "sim_reaction.h"

sim_reaction *sim_reaction_init(const ion *incident_ion, const jibal_isotope *isotopes, const sample *sample, const detector *det, const reaction *r, size_t n_channels, size_t n_bricks) {
    if(!r) {
        return NULL;
    }
    assert(r->product);
    sim_reaction *sim_r = malloc(sizeof(sim_reaction));
    sim_r->r = r;
    ion *p = &sim_r->p;
    ion_reset(p);
    sim_r->max_depth = 0.0;
    sim_r->i_isotope = sample->n_isotopes; /* Intentionally not valid */

    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        if(sample->isotopes[i_isotope] == r->target) {
#ifdef DEBUG
            fprintf(stderr, "Reaction target isotope %s is isotope number %zu in sample.\n", r->target->name, i_isotope);
#endif
            sim_r->i_isotope = i_isotope;
        }
    }
    sim_r->histo = gsl_histogram_alloc(n_channels); /* free'd by sim_workspace_free */
    spectrum_set_calibration(sim_r->histo, det, r->product->Z); /* Setting histogram with Z-specific (or as fallback, default) calibration. */
    gsl_histogram_reset(sim_r->histo);
    sim_r->n_bricks = n_bricks;
    sim_r->bricks = calloc(sim_r->n_bricks, sizeof(brick));
    ion_set_isotope(p, r->product);
    if(p->isotope == incident_ion->isotope) {
        p->nucl_stop = nuclear_stopping_shared_copy(incident_ion->nucl_stop);
    } else {
        p->nucl_stop = nuclear_stopping_new(p->isotope, isotopes);
    }
    sim_reaction_set_cross_section_by_type(sim_r);
    return sim_r;
}

void sim_reaction_free(sim_reaction *sim_r) {
    if(!sim_r) {
        return;
    }
    if(sim_r->histo) {
        gsl_histogram_free(sim_r->histo);
        sim_r->histo = NULL;
    }
    if(sim_r->bricks) {
        free(sim_r->bricks);
        sim_r->bricks = NULL;
    }
    nuclear_stopping_free(sim_r->p.nucl_stop);
    free(sim_r);
}

void sim_reaction_recalculate_internal_variables(sim_reaction *sim_r, const sim_calc_params *params, double theta, double E_min, double E_max) {
    /* Calculate variables for Rutherford (and Andersen) cross sections. This is done for all reactions, even if they are not RBS or ERD reactions! */
    (void) E_min; /* Energy range could be used to set something (in the future) */
    (void) E_max;
    if(!sim_r || !sim_r->r)
        return;
    const jibal_isotope *incident = sim_r->r->incident;
    const jibal_isotope *target = sim_r->r->target;
    sim_r->E_cm_ratio = target->mass / (incident->mass + target->mass);
    sim_r->mass_ratio = incident->mass / target->mass;
    sim_r->theta = theta;
    if(sim_r->r->Q == 0.0) {
        sim_r->K = reaction_product_energy(sim_r->r, sim_r->theta, 1.0);
    } else {
        sim_r->K = 0.0;
    }
    sim_r->cs_constant = 0.0;
    sim_r->theta_cm = 0.0; /* Will be recalculated, if possible */
    reaction_type type = sim_r->r->type;

    if(!reaction_is_possible(sim_r->r, params, theta)) {
        sim_r->stop = TRUE;
#ifdef DEBUG
        fprintf(stderr, "Reaction not possible, returning.\n");
#endif
        return;
    }
    if(type == REACTION_RBS) {
        sim_r->theta_cm = sim_r->theta + asin(sim_r->mass_ratio * sin(sim_r->theta));
    }
    if(type == REACTION_RBS_ALT) {
        sim_r->theta_cm = C_PI - (asin(sim_r->mass_ratio * sin(sim_r->theta)) - sim_r->theta);
    }
    if(type == REACTION_RBS || type ==REACTION_RBS_ALT) {
        sim_r->cs_constant = fabs((pow2(sin(sim_r->theta_cm))) / (pow2(sin(sim_r->theta)) * cos(sim_r->theta_cm - sim_r->theta)) *
                                  pow2((incident->Z * C_E * target->Z * C_E) / (4.0 * C_PI * C_EPSILON0)) *
                                  pow4(1.0 / sin(sim_r->theta_cm / 2.0)) * (1.0 / 16.0));
    }
    if(type == REACTION_ERD) { /* ERD */
        sim_r->theta_cm = C_PI - 2.0 * sim_r->theta;
        sim_r->cs_constant = pow2(incident->Z * C_E * target->Z * C_E / (8 * C_PI * C_EPSILON0)) * pow2(1.0 + incident->mass / target->mass) * pow(cos(sim_r->theta), -3.0) * pow2(sim_r->E_cm_ratio);
    }
    if(sim_r->r->cs == JABS_CS_ANDERSEN) {
        sim_r->r_VE_factor = 48.73 * C_EV * incident->Z * target->Z * sqrt(pow(incident->Z, 2.0 / 3.0) + pow(target->Z, 2.0 / 3.0)); /* Factors for Andersen correction */
        sim_r->r_VE_factor2 = pow2(0.5 / sin(sim_r->theta_cm / 2.0));
    }
#ifdef DEBUG
    fprintf(stderr, "Reaction recalculated, theta = %g deg, theta_cm = %g deg, K = %g (valid for RBS and ERD). Q = %g MeV.\n", sim_r->theta/C_DEG, sim_r->theta_cm/C_DEG, sim_r->K, sim_r->r->Q / C_MEV);
#endif
}

void sim_reaction_reset_bricks(sim_reaction *sim_r) {
    memset(sim_r->bricks, 0, sizeof(brick) * sim_r->n_bricks);
}

void sim_reaction_set_cross_section_by_type(sim_reaction *sim_r) {
    switch(sim_r->r->type) {
        case REACTION_RBS:
            /* Falls through */
        case REACTION_RBS_ALT:
            /* Falls through */
        case REACTION_ERD:
            sim_r->cross_section = sim_reaction_cross_section_rutherford;
            break;
        case REACTION_FILE:
            sim_r->cross_section = sim_reaction_cross_section_tabulated;
            break;
#ifdef JABS_PLUGINS
        case REACTION_PLUGIN:
            sim_r->cross_section = sim_reaction_cross_section_plugin;
            break;
#endif
        default:
#ifdef DEBUG
            fprintf(stderr, "Unknown reaction type %i, cross_section() function pointer is NULL!\n", sim_r->r->type);
#endif
            sim_r->cross_section = NULL;
    }
}

double sim_reaction_andersen(const sim_reaction *sim_r, double E_cm) {
    const double r_VE = sim_r->r_VE_factor / E_cm;
    return pow2(1 + 0.5 * r_VE) / pow2(1 + r_VE + sim_r->r_VE_factor2 * pow2(r_VE));
}

double sim_reaction_cross_section_rutherford(const sim_reaction *sim_r, double E) {
#ifdef CROSS_SECTIONS_FROM_JIBAL
    return jibal_cross_section_erd(sim_r->r->incident, sim_r->r->target, sim_r->theta, E, sim_r->r->cs);
#else
    const reaction *r = sim_r->r;
    if(E > r->E_max || E < r->E_min)
        return 0.0;
    const double E_cm = sim_r->E_cm_ratio * E;
    double sigma_r = sim_r->cs_constant / pow2(E_cm) ;
    switch(r->cs) {
        case JABS_CS_RUTHERFORD:
            return sigma_r;
        case JABS_CS_ANDERSEN:
            return sigma_r * sim_reaction_andersen(sim_r, E_cm);
        default:
            return 0.0;
    }
#endif
}

double sim_reaction_cross_section_tabulated(const sim_reaction *sim_r, double E) {
    size_t lo, mi, hi;
    const reaction *r = sim_r->r;
    const struct reaction_point *t = r->cs_table;
    hi = r->n_cs_table - 1;
    lo = 0;
    if(E < t[lo].E || E > t[hi].E) {
#ifdef REACTIONS_FALL_BACK
        return sim_reaction_cross_section_rutherford(sim_r, E); /* Fall back quietly to analytical formulae outside tabulated values */
#else
        return 0.0;
#endif
    }
    while (hi - lo > 1) {
        mi = (hi + lo) / 2;
        if (E >= t[mi].E) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    return t[lo].sigma+((t[lo+1].sigma-t[lo].sigma)/(t[lo+1].E-t[lo].E))*(E-t[lo].E);
}

#ifdef JABS_PLUGINS
double sim_reaction_cross_section_plugin(const sim_reaction *sim_r, double E) {
    jabs_plugin_reaction *r = sim_r->r->plugin_r;
    return r->cs(r, sim_r->theta, E);
}
#endif

void sim_reaction_product_energy_and_straggling(sim_reaction *r, const ion *incident) {
    if(r->r->Q == 0.0) {
        assert(r->K > 0 && r->K <= 1.0);
        assert(incident->E > 1.0 * C_EV && incident->E < 1000.0 * C_MEV);
        r->p.E = incident->E * r->K;
        r->p.S = incident->S * pow2(r->K);
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Product energy %g keV, eloss straggling %g keV FWHM. Calculated using K = %g\n", r->p.E/C_KEV, C_FWHM * sqrt(r->p.S) / C_KEV, r->K);
#endif
        return;
    }
    r->p.E = reaction_product_energy(r->r, r->theta, incident->E);
    double epsilon = 0.001*C_KEV;
    double deriv = (reaction_product_energy(r->r, r->theta, incident->E+epsilon) - r->p.E)/(epsilon); /* TODO: this derivative could be solved analytically */
    r->p.S = incident->S * pow2(deriv) * incident->E;
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "deriv %g, E_out/E %g, E_out = %g keV, E = %g keV\n", deriv, r->p.E / incident->E, r->p.E/C_KEV, incident->E/C_KEV);
#endif
}

void sim_reaction_print_bricks(FILE *f, const sim_reaction *r, double psr) {
    fprintf(f, "#Reaction %s %s\n", reaction_name(r->r), r->r->target->name);
    if(r->r->filename) {
        fprintf(f, "#Filename: %s\n", r->r->filename);
    }
    fprintf(f, "#brick    depth    thick      E_0  S_0(el)      E_r  S_r(el)   E(det)    S(el)    S(geo)    S(sum) sigma*conc              Q  dE(det)/dE_0\n");
    for(size_t j = 0; j <= r->last_brick; j++) {
        brick *b = &r->bricks[j];
        fprintf(f, "%4zu %10.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %10.1lf %10.3lf %14.6e %8.3lf\n",
                j, b->d.x / C_TFU, b->thick / C_TFU,
                b->E_0 / C_KEV, C_FWHM * sqrt(b->S_0) / C_KEV,
                b->E_r / C_KEV, C_FWHM * sqrt(b->S_r) / C_KEV,
                b->E / C_KEV, C_FWHM * sqrt(b->S) / C_KEV,
                C_FWHM * sqrt(b->S_geo_x + b->S_geo_y) / C_KEV, C_FWHM * b->S_sum / C_KEV,
                b->sc / C_MB_SR, b->Q * psr,
                b->deriv
        );
    }
}
