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
#include "jabs_debug.h"
#include "spectrum.h"
#include "sim_reaction.h"
#include "scatint.h"
#include "defaults.h"

sim_reaction *sim_reaction_init(const sample *sample, const detector *det, const reaction *r, size_t n_channels, size_t n_bricks) {
    if(!r) {
        return NULL;
    }
    assert(r->product);
    sim_reaction *sim_r = calloc(1, sizeof(sim_reaction));
    sim_r->r = r;
    ion *p = &sim_r->p;
    ion_reset(p);
    sim_r->i_isotope = sample->n_isotopes; /* Intentionally not valid */

    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        if(sample->isotopes[i_isotope] == r->target) {
            DEBUGMSG("Reaction target isotope %s is isotope number %zu in sample.", r->target->name, i_isotope);
            sim_r->i_isotope = i_isotope;
        }
    }
    sim_r->histo = jabs_histogram_alloc(n_channels); /* free'd by sim_workspace_free */
    calibration_apply_to_histogram(detector_get_calibration(det, r->product->Z), sim_r->histo); /* Setting histogram with Z-specific (or as fallback, default) calibration. */
    jabs_histogram_reset(sim_r->histo);
    sim_r->n_bricks = n_bricks;
    sim_r->bricks = calloc(sim_r->n_bricks, sizeof(brick));
    ion_set_isotope(p, r->product);
    p->nucl_stop = r->nucl_stop; /* We just borrow this */
    assert(p->nucl_stop);
    assert(p->nucl_stop->incident = p->isotope);
    p->ion_gsto = r->ion_gsto; /* We just borrow this */
    assert(p->ion_gsto);
    assert(p->ion_gsto->incident = p->isotope);
    sim_reaction_set_cross_section_by_type(sim_r);
    return sim_r;
}

void sim_reaction_free(sim_reaction *sim_r) {
    if(!sim_r) {
        return;
    }
    if(sim_r->histo) {
        jabs_histogram_free(sim_r->histo);
        sim_r->histo = NULL;
    }
    if(sim_r->bricks) {
        free(sim_r->bricks);
        sim_r->bricks = NULL;
    }
    free(sim_r->cs_table);
    free(sim_r);
}

int sim_reaction_recalculate_internal_variables(sim_reaction *sim_r, const sim_calc_params *params, double theta, double emin, double emin_incident, double emax_incident) {
    /* Calculate variables for Rutherford (and Andersen) cross sections. This is done for all reactions, even if they are not RBS or ERD reactions! */
    /* E_min and E_max are assumed to be minimum and maximum (respectively) for the incident ion energy (given by simulation, stopping data availability etc) */
    if(!sim_r || !sim_r->r) {
        return EXIT_FAILURE;
    }
    const jibal_isotope *incident = sim_r->r->incident;
    const jibal_isotope *target = sim_r->r->target;
    sim_r->E_cm_ratio = target->mass / (incident->mass + target->mass);
    sim_r->mass_ratio = incident->mass / target->mass;
    sim_r->theta = theta;
    sim_r->K = 0.0;
    sim_r->cs_constant = 0.0;
    sim_r->theta_cm = 0.0; /* Will be recalculated, if possible */
    sim_r->emin_incident = GSL_MAX_DBL(emin, emin_incident);
    sim_r->emin_incident = GSL_MAX_DBL(sim_r->emin_incident, sim_r->r->E_min);
    sim_r->emax_incident = GSL_MIN_DBL(emax_incident, sim_r->r->E_max); /* Note that sometimes going below this energy is required (straggling) */

    sim_r->emin_product = GSL_MAX_DBL(reaction_product_energy(sim_r->r, sim_r->theta, sim_r->emin_incident), sim_r->p.ion_gsto ? sim_r->p.ion_gsto->emin : 0.0);
    sim_r->emax_product = GSL_MIN_DBL(reaction_product_energy(sim_r->r, sim_r->theta, sim_r->emax_incident), sim_r->p.ion_gsto ? sim_r->p.ion_gsto->emax : 0.0);
    reaction_type type = sim_r->r->type;

    if(!reaction_is_possible(sim_r->r, params, theta)) {
        sim_r->stop = TRUE;
        DEBUGSTR("Reaction not possible, returning.");
        return EXIT_SUCCESS; /* Not a failure as such */
    }
    if(sim_r->r->Q == 0.0) {
        sim_r->K = reaction_product_energy(sim_r->r, sim_r->theta, 1.0);
    } else {
        sim_r->K = 0.0;
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
    } else if(sim_r->r->cs == JABS_CS_LECUYER) {
        sim_r->lecuyer_factor = 48.73 * C_EV * incident->Z * pow(target->Z, 4.0/3.0); /* Factor 0.049 given in NIM 160 (1979) 337-346 */
    }
    if(sim_r->r->cs != JABS_CS_NONE && (type == REACTION_RBS || type == REACTION_RBS_ALT || type == REACTION_ERD)) {
        if(params->screening_tables || sim_r->r->cs == JABS_CS_UNIVERSAL) { /* Universal needs to be precalculated to be remotely useful. Forcing it! */
            if(sim_reaction_recalculate_screening_table(sim_r)) {
                jabs_message(MSG_ERROR, "Recalculating screening table for reaction %s failed.\n", reaction_name(sim_r->r));
                DEBUGMSG("Recalculating screening table for reaction %s failed.", reaction_name(sim_r->r));
                return EXIT_FAILURE;
            }
        }
    } else {
        sim_reaction_reset_screening_table(sim_r);
    }

    DEBUGMSG("Reaction %s recalculated, theta = %g deg, theta_cm = %g deg, K = %g (valid for RBS and ERD). Q = %g MeV. E_incident = [%g keV, %g keV], E_product = [%g keV, %g keV]",
             sim_r->r->name, sim_r->theta/C_DEG, sim_r->theta_cm/C_DEG, sim_r->K, sim_r->r->Q / C_MEV,
             sim_r->emin_incident / C_KEV, sim_r->emax_incident / C_KEV,
             sim_r->emin_product / C_KEV, sim_r->emax_product / C_KEV
    );
    return EXIT_SUCCESS;
}

int sim_reaction_recalculate_screening_table(sim_reaction *sim_r) {
    jabs_reaction_cs cs = sim_r->r->cs;
    size_t n;
    sim_reaction_reset_screening_table(sim_r);
    scatint_params *sp = NULL;
    n = SCREENING_TABLE_ELEMENTS; /* TODO: make more clever */
    potential_type pt = POTENTIAL_NONE;
    switch(cs) { /* Set potential type for those cross sections that use scatint (= numerical scattering integral solver). If analytical solution is used (e.g. Andersen), keep as POTENTIAL_NONE */
        case JABS_CS_ANDERSEN:
            break;
        case JABS_CS_UNIVERSAL:
            pt = POTENTIAL_UNIVERSAL;
            break;
        case JABS_CS_TEST:
            pt = POTENTIAL_TEST;
            break;
        case JABS_CS_THOMASFERMI:
            pt = POTENTIAL_TF_SOMMERFELD;
            break;
        case JABS_CS_LECUYER:
            break;
        default:
            return EXIT_SUCCESS;
    }
    if(pt != POTENTIAL_NONE) {
        sp = scatint_init(sim_r->r->type, pt, sim_r->r->incident, sim_r->r->target);
        if(!sp) {
            return EXIT_FAILURE;
        }
        scatint_set_theta(sp, sim_r->theta);
    }
    sim_r->cs_table = malloc(n * sizeof(struct reaction_point));
    if(!sim_r->cs_table) {
        return EXIT_FAILURE;
    }
    sim_r->n_cs_table = n;
    double emin = sim_r->emin_incident * 0.9; /* Safety factor included, note that screening table is *just* a screening table, cross section below sim_r->r->E_min should be zero, but screening is not! */
    double emax = sim_r->emax_incident * 1.1;
    sim_r->cs_estep = (emax - emin)/(n - 1);
    DEBUGMSG("Computing screening, reaction %s, from %g keV to %g keV with %zu steps of %g keV", sim_r->r->name, emin / C_KEV, emax / C_KEV, n, sim_r->cs_estep / C_KEV);
    for(size_t i = 0; i < n; i++) {
        double E = emin + (sim_r->cs_estep) * i;
        const double E_cm = sim_r->E_cm_ratio * E;
        double sigma_r = sim_r->cs_constant / pow2(E_cm) ;
        struct reaction_point *rp = &sim_r->cs_table[i];
        rp->E = E;
        switch(cs) {
            case JABS_CS_ANDERSEN:
                rp->sigma = sim_reaction_andersen(sim_r, E_cm);
                break;
            case JABS_CS_LECUYER:
                rp->sigma = sim_reaction_lecuyer(sim_r, E_cm);
                break;
            case JABS_CS_TEST:
            case JABS_CS_THOMASFERMI:
            case JABS_CS_UNIVERSAL:
                scatint_set_energy(sp, E);
                rp->sigma  = scatint_sigma(sp) / sigma_r; /* Note: only screening correction, not cross section! */
                break;
            default: /* Not reached */
                rp->sigma  = 0.0;
        }
        if(rp->sigma == 0.0) {
            jabs_message(MSG_ERROR, "Could not compute screening correction for E = %g keV (emin = %g keV, emax = %g keV), Rutherford = %g mb/sr, reaction %s\n", E / C_KEV, emin / C_KEV, emax / C_KEV, sigma_r / C_MB_SR, reaction_name(sim_r->r));
            DEBUGMSG("Could not compute screening correction for E = %g keV, reaction %s\n", E / C_KEV, sim_r->r->name)
            scatint_params_free(sp);
            return EXIT_FAILURE;
        }
    }
    scatint_params_free(sp);
    DEBUGVERBOSEMSG("Screening table with %zu points recalculated, energy in range [%g keV, %g keV]", sim_r->n_cs_table, emin / C_KEV, emax / C_KEV);
    return EXIT_SUCCESS;
}

void sim_reaction_reset_screening_table(sim_reaction *sim_r) {
    free(sim_r->cs_table);
    sim_r->n_cs_table = 0;
}

double sim_reaction_evaluate_screening_table(const sim_reaction *sim_r, double E) {
    size_t i = (E - sim_r->cs_table[0].E)/sim_r->cs_estep;
    if(i >= sim_r->n_cs_table) { /* We have to extrapolate, either E is too high or too low (unsigned i has underflowed). This obviously should not happen (frequently). */
        DEBUGVERBOSEMSG("Screening table is being extrapolated, because E = %g keV!", E / C_KEV);
        if(E < sim_r->emin_incident) {
            return sim_r->cs_table[0].sigma;
        } else {
            return sim_r->cs_table[sim_r->n_cs_table - 1].sigma;
        }
    }
    return sim_r->cs_table[i].sigma + (sim_r->cs_table[i+1].sigma - sim_r->cs_table[i].sigma)*((E - sim_r->cs_table[i].E)/sim_r->cs_estep);
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

double sim_reaction_lecuyer(const sim_reaction *sim_r, double E_cm) {
    return (1.0 - sim_r->lecuyer_factor/E_cm);
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
    if(sim_r->cs_table) {
        return sigma_r * sim_reaction_evaluate_screening_table(sim_r, E);
    }
    switch(r->cs) {
        case JABS_CS_RUTHERFORD:
            return sigma_r;
        case JABS_CS_ANDERSEN:
            return sigma_r * sim_reaction_andersen(sim_r, E_cm);
        case JABS_CS_UNIVERSAL:
            return 0.0; /* Universal screening should be combined with screening table calculation */
        case JABS_CS_LECUYER:
            return sigma_r * sim_reaction_lecuyer(sim_r, E_cm);
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
    fprintf(f, "#Reaction %s\n", reaction_name(r->r));
    if(r->r->filename) {
        fprintf(f, "#Filename: %s\n", r->r->filename);
    }
    fprintf(f, "#  i v      depth    thick      E_0   S_0(el)      E_r   S_r(el)      E_s  S_s(el)   E(det)    S(el)    S(geo)    S(sum) sigma*conc     Q (counts)  dE(det)/dE_0\n");
    for(size_t j = 0; j <= r->last_brick; j++) {
        brick *b = &r->bricks[j];
        fprintf(f, "%4zu %1i %10.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf  %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %10.1lf %10.3lf %14.6e %8.3lf\n",
                j, b->valid, b->d.x / C_TFU, b->thick / C_TFU,
                b->E_0 / C_KEV, C_FWHM * sqrt(b->S_0) / C_KEV,
                b->E_r / C_KEV, C_FWHM * sqrt(b->S_r) / C_KEV,
                b->E_s / C_KEV, C_FWHM * sqrt(b->S_s) / C_KEV,
                b->E / C_KEV, C_FWHM * sqrt(b->S) / C_KEV,
                C_FWHM * sqrt(b->S_geo_x + b->S_geo_y) / C_KEV, C_FWHM * b->S_sum / C_KEV,
                b->sc / C_MB_SR, b->Q * psr,
                b->deriv
        );
    }
}
