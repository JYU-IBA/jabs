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
#include "defaults.h"
#include "generic.h"
#include "message.h"
#include "simulation_workspace.h"
#include "spectrum.h"


void sim_workspace_init_reactions(sim_workspace *ws) {
    const simulation *sim = ws->sim;
    ws->reactions = calloc(sim->n_reactions, sizeof (sim_reaction *));
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        ws->reactions[i_reaction] = sim_reaction_init(sim->sample, ws->det, sim->reactions[i_reaction], ws->n_channels, ws->n_bricks);
        ws->n_reactions++;
    }
}

void sim_workspace_calculate_number_of_bricks(sim_workspace *ws) {
    size_t n_bricks = 0;
    const detector *det = ws->det;
    const simulation *sim = ws->sim;
    if(ws->params->n_bricks_max) {
        n_bricks = ws->params->n_bricks_max;
    } else {
        if(det->type == DETECTOR_ENERGY) {
            if(ws->params->incident_stop_params.step == 0.0) { /* Automatic incident step size */
                double E_min = GSL_MAX_DBL(ws->params->incident_stop_params.min, ws->params->brick_width_sigmas * sqrt(detector_resolution(ws->det, sim->beam_isotope, sim->beam_E))); /* Trying to guesstimate a lower limit for incident ion energy steps */
                n_bricks = (int) ceil(sim->beam_E / E_min) + ws->sample->n_ranges;
            } else {
                n_bricks = (int) ceil(sim->beam_E / ws->params->incident_stop_params.step) + ws->sample->n_ranges; /* This is conservative */
            }
        } else { /* TODO: maybe something more clever is needed here :) */
            n_bricks = BRICKS_DEFAULT;
        }
        if(n_bricks > BRICKS_MAX) {
            n_bricks = BRICKS_MAX;
            jabs_message(MSG_WARNING, stderr,  "Number of bricks limited to %zu.\n", n_bricks);
        }
    }
    ws->n_bricks = n_bricks;
    DEBUGMSG("Number of bricks in workspace: %zu", n_bricks);
}

sim_workspace *sim_workspace_init(const jibal *jibal, const simulation *sim, const detector *det) {
    if(!jibal || !sim || !det) {
        jabs_message(MSG_ERROR, stderr,  "No JIBAL, sim or det. Guru thinks: %p, %p %p.\n", jibal, sim, det);
        return NULL;
    }
    if(!sim->sample) {
        jabs_message(MSG_ERROR, stderr,  "No sample has been set. Will not initialize workspace.\n");
        return NULL;
    }
    sim_workspace *ws = calloc(1, sizeof(sim_workspace));
    ws->sim = sim;
    ws->emin = sim->emin;
    ws->fluence = sim->fluence;
    ws->det = det;
    ws->sample = sim->sample;
    ws->params = sim_calc_params_defaults(NULL);
    sim_calc_params_copy(sim->params,  ws->params);
    sim_calc_params_update(ws->params);
    ws->gsto = jibal->gsto;
    ws->isotopes = jibal->isotopes;
    ws->n_reactions = 0; /* Will be incremented later */

    if(sim->params->beta_manual && sim->params->ds) {
        jabs_message(MSG_WARNING, stderr,  "Manual exit angle is enabled in addition to dual scattering. This is an unsupported combination. Manual exit angle calculation will be disabled.\n");
        ws->params->beta_manual  = FALSE;
    }

    if(sim->params->beta_manual && sim->params->geostragg) {
        jabs_message(MSG_WARNING, stderr,  "Manual exit angle is enabled in addition to geometric scattering. This is an unsupported combination. Geometric straggling calculation will be disabled.\n");
        ws->params->geostragg = FALSE;
    }
    ws->ion = sim->ion; /* Shallow copy, but that is ok */
    ws->ion.E = sim->beam_E;
    ws->ion.S = pow2(sim->beam_E_broad / C_FWHM);

    sim_workspace_recalculate_n_channels(ws, ws->sim);

    if(ws->n_channels == 0) {
        jabs_message(MSG_ERROR, stderr,  "Number of channels in workspace is zero. Aborting initialization.\n");
        free(ws);
        return NULL;
    }

    if(detector_sanity_check(det, ws->n_channels)) {
        jabs_message(MSG_ERROR, stderr, "Detector %zu failed sanity check!", det->name);
        free(ws);
        return NULL;
    }

    ws->histo_sum = jabs_histogram_alloc(ws->n_channels);
    if(!ws->histo_sum) {
        return NULL;
    }
    spectrum_set_calibration(ws->histo_sum, ws->det->calibration); /* Calibration (assuming default calibration) can be set now. */
    jabs_histogram_reset(ws->histo_sum); /* This is not necessary, since contents should be set after simulation is over (successfully). */

    sim_workspace_calculate_number_of_bricks(ws);
    sim_workspace_init_reactions(ws);

    if(ws->params->cs_adaptive) { /* Actually integrate, allocate workspace for this */
        DEBUGMSG("cs_n_steps = 0, allocating integration workspace w_int_cs with %zu max intervals.", ws->params->int_cs_max_intervals);
        ws->w_int_cs = gsl_integration_workspace_alloc(ws->params->int_cs_max_intervals);
    } else {
        ws->w_int_cs = NULL;
    }
    if(ws->params->cs_n_stragg_steps == 0) {
        DEBUGMSG("cs_n_stragg_steps = 0, allocating integration workspace w_int_cs_stragg with %zu max intervals.", ws->params->int_cs_stragg_max_intervals);
        ws->w_int_cs_stragg = gsl_integration_workspace_alloc(ws->params->int_cs_stragg_max_intervals);
    } else {
        ws->w_int_cs_stragg = NULL;
    }
    ws->stop.type = GSTO_STO_TOT;
    ws->stop.gsto = jibal->gsto;
    ws->stop.rk4 = ws->params->rk4;
    ws->stop.nuclear_stopping_accurate = ws->params->nuclear_stopping_accurate;
    ws->stop.emin = ws->emin;
    ws->stragg = ws->stop; /* Copy */
    ws->stragg.type = GSTO_STO_STRAGG;
    return ws;
}

void sim_workspace_free(sim_workspace *ws) {
    if(!ws)
        return;
    sim_calc_params_free(ws->params);
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction_free(ws->reactions[i]);
    }
    free(ws->reactions);
    gsl_integration_workspace_free(ws->w_int_cs);
    gsl_integration_workspace_free(ws->w_int_cs_stragg);
    jabs_histogram_free(ws->histo_sum);
    free(ws);
}

void sim_workspace_recalculate_n_channels(sim_workspace *ws, const simulation *sim) { /* TODO: assumes calibration function is increasing */
    size_t n_max = CHANNELS_ABSOLUTE_MIN; /* Always simulate at least CHANNELS_ABSOLUTE_MIN channels */
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        const reaction *r = sim->reactions[i_reaction];
        if(!reaction_is_possible(r, ws->params, ws->det->theta)) {
            DEBUGMSG("Reaction %zu (target %s, type %s) is not possible when theta = %g deg. Skipping.", i_reaction + 1, r->target->name,
                    reaction_type_to_string(r->type), ws->det->theta / C_DEG);
            continue;
        }
        double E = reaction_product_energy(r, ws->det->theta, sim->beam_E + ws->params->sigmas_cutoff * sim->beam_E_broad / C_FWHM); /* Highest energy we could possibly reach, maybe */
        double E_safer = E + ws->params->sigmas_cutoff * sqrt(detector_resolution(ws->det, r->product, E)); /* Add resolution sigma to max energy */
        E_safer *= 1.1; /* and 10% for good measure! */
        const calibration *cal = detector_get_calibration(ws->det, r->product->Z);
        double E_old = calibration_eval(cal, n_max - 1);
        while(n_max < CHANNELS_ABSOLUTE_MAX) {
            double E_det = calibration_eval(cal, n_max);
            if(E_det < E_old) {
                DEBUGMSG("Monotonicity fail at %zu, %g keV < %g keV\n", n_max, E_det / C_KEV, E_old / C_KEV); /* This monotonicity check does not guarantee the monotonicity of all calibrations every time, so further checks are necessary */
                jabs_message(MSG_ERROR, stderr, "Could not determine number of channels required in simulation while processing reaction %zu, since calibrations must be monotonous, but channel %zu gives %.14g keV, which is less than %.14g keV on the previous channel.\n", i_reaction + 1, n_max, E_det / C_KEV, E_old / C_KEV);
                char *cal_str = calibration_to_string(cal);
                jabs_message(MSG_ERROR, stderr, "The offending calibration is currently \"%s\"\n", cal_str);
                free(cal_str);
                ws->n_channels = 0;
                return;
            }
            if(E_det > E_safer) {
                break;
            }
            E_old = E_det;
            n_max++;
        } /* Increase number of channels until we hit goal energy. */
        DEBUGMSG("Reaction %zu: %s, max E = %g keV -> %g keV after resolution and safety factor have been added, current n_max %zu (channels)",
                 i_reaction + 1, r->name, E / C_KEV, E_safer / C_KEV, n_max);
    }
    if(n_max >= CHANNELS_ABSOLUTE_MAX) {
        DEBUGMSG("Number of channels in workspace would exceed the built-in limit %i!\n", CHANNELS_ABSOLUTE_MAX);
        n_max = 0;
    }
    ws->n_channels = n_max;
    DEBUGMSG("Number of channels in workspace set to %zu\n", ws->n_channels);
}

void sim_workspace_calculate_sum_spectra(sim_workspace *ws) {
    for(size_t j = 0; j < ws->n_reactions; j++) {
        sim_reaction *sim_r = ws->reactions[j];
        if(sim_r->n_convolution_calls == 0 || !sim_r->histo) {
            continue;
        }
        size_t n = GSL_MIN(ws->histo_sum->n, sim_r->histo->n);
        for(size_t i = 0; i < n; i++) {
            ws->histo_sum->bin[i] += sim_r->histo->bin[i];
        }
    }
}

void sim_workspace_histograms_reset(sim_workspace *ws) {
    jabs_histogram_reset(ws->histo_sum);
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = ws->reactions[i];
        if(!r)
            continue;
        jabs_histogram_reset(r->histo);
    }
}

size_t sim_workspace_histograms_calculate(sim_workspace *ws) {
    size_t n_meaningful = 0;
#if 0
#pragma omp parallel default(none) shared(ws) shared(n_meaningful)
#pragma omp for nowait
#endif
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = ws->reactions[i];
        if(!r)
            continue;
        assert(r->last_brick < r->n_bricks);
        if(r->last_brick == 0) {
            continue;
        }
        bricks_calculate_sigma(ws->det, r->p.isotope, r->bricks, r->last_brick);
        const calibration *c = detector_get_calibration(ws->det, r->p.Z);
#ifdef DEBUG_BRICK_OUTPUT
        fprintf(stdout, "BRICK #Reaction %zu: %s, last brick %zu\n", i + 1, r->r->name, r->last_brick);
#endif
        bricks_convolute(r->histo, c, r->bricks, r->last_brick, ws->fluence * ws->det->solid, ws->params->sigmas_cutoff, r->emin, ws->params->gaussian_accurate);
#ifdef DEBUG_BRICK_OUTPUT
        fprintf(stdout, "BRICK \nBRICK \n");
#endif
        r->n_convolution_calls++;
#if 0
#pragma omp critical
#endif
        n_meaningful++;
    }
    return n_meaningful;
}

void sim_workspace_histograms_scale(sim_workspace *ws, double scale) {
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = ws->reactions[i];
        if(!r || !r->histo)
            continue;
        jabs_histogram_scale(r->histo, scale);
    }
}


int sim_workspace_print_spectra(const result_spectra *spectra, const char *filename) {
    char sep = ' ';
    if(!spectra || spectra->n_spectra == 0) {
        jabs_message(MSG_ERROR, stderr, "No spectra!\n");
        return EXIT_FAILURE;
    }
    jabs_histogram *h_sum = result_spectra_simulated_histo(spectra);
    jabs_histogram *h_exp = result_spectra_experimental_histo(spectra);
    if(!h_sum) {
        jabs_message(MSG_ERROR, stderr, "No simulated spectrum!\n");
        return EXIT_FAILURE;
    }
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f) {
        return EXIT_FAILURE;
    }
    if(filename) {
        if(strncmp(jabs_file_extension_const(filename), ".csv", 4) == 0) { /* For CSV: print header line */
            sep = ','; /* and set the separator! */
            fprintf(f, "\"Channel\",\"Energy (keV)\"");
            for(size_t i_spectrum = 0; i_spectrum < spectra->n_spectra; i_spectrum++) {
                fprintf(f, ",\"%s\"", spectra->s[i_spectrum].name);
            }
            fprintf(f, "\n");
        }
    }
    size_t n_ch;
    double *range;
    if(h_exp && h_exp->n > h_sum->n) { /* Energy for output is taken from experimental spectrum if it exists and extends to higher energy than simulated spectrum */
        range = h_exp->range;
        n_ch = h_exp->n;
    } else {
        range = h_sum->range;
        n_ch = h_sum->n;
    }
    for(size_t ch = 0; ch < n_ch; ch++) {
        fprintf(f, "%zu%c%.3lf", ch, sep, range[ch] / C_KEV); /* Channel, energy. TODO: Z-specific calibration can have different energy (e.g. for a particular reaction). */
        for(size_t j = 0; j < spectra->n_spectra; j++) {
            jabs_histogram *h = spectra->s[j].histo;
            if(!h || ch >= h->n) {
                fprintf(f, "%c0", sep);
            } else {
                fprintf(f, "%c%.8g", sep, h->bin[ch]);
            }
        }
        fprintf(f, "\n");
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

int sim_workspace_print_bricks(const sim_workspace *ws, const char *filename) {
    double psr = ws->fluence * ws->det->solid;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f) {
        return EXIT_FAILURE;
    }
    fprintf(f, "#This is a bricks file containing up to %zu reactions.\n", ws->n_reactions);
    fprintf(f, "#Each brick is convoluted with a gaussian to make an energy spectrum.\n");
    fprintf(f, "#In case of roughness, this represents the bricks of the last simulated roughness subspectrum.\n");
    fprintf(f, "#Turn roughness and DS off if you want to understand the contents of this file.\n");
    fprintf(f, "#Each brick has two edges, the high energy edge and the low energy edge.\n");
    fprintf(f, "#\"depth\" is the depth of the low energy edge in tfu (1e15 at./cm2), thick is thickness of the brick (difference in depth to high energy edge) in the same units.\n");
    fprintf(f, "#\"E_0\" and \"S_0(el)\" are the energy and electronic energy loss straggling in keVs FWHM of the incident ion at the particular depth (low energy edge).\n");
    fprintf(f, "#\"E_r\" and \"S_r(el)\" are the same, but for the reaction product at the particular depth. For elastic scattering E_r = K * E_0.\n");
    fprintf(f, "#\"E(det)\" and \"S(el)\" are the same, but for the reaction product at the detector (after detector foil, if one is set).\n");
    fprintf(f, "#\"S(geo)\" is the broadening due to finite beam spot and detector size (keV FWHM) at the surface if calculated, zero otherwise.\n");
    fprintf(f, "#\"S(sum)\" is the sum of all broadening (electronic, geometric and detector resolution) in keV FWHM. This is the actual broadening used for convolution.\n");
    fprintf(f, "#When convoluting, JaBS convolutes a the brick (a box from low energy to high energy edge with area of Q counts) by a gaussian with varying sigma (different S(sum) for either edge)\n");
    fprintf(f, "#\"sigma * conc\" is the concentration of the target isotope multiplied by the (straggling weighted) cross section, in mb/sr. This is representative of the entire brick.\n");
    fprintf(f, "#\"Q\" is the number of counts in the brick. It is the sigma*conc product multiplied by effective thickness of the brick (thick/cos(a)), number of incident particles and solid angle.\n");
    fprintf(f, "#\"dE(det)/dE_0\" is the (estimate of) derivative dE(det)/dE_0. This is used to figure out the length of energy steps (incident energy, E_0) to take to get comfortable brick spacing.\n");
    fprintf(f, "#Thanks for reading! If you have any questions, please email Jaakko...\n\n");
    for(size_t i = 0; i < ws->n_reactions; i++) {
        const sim_reaction *r = ws->reactions[i];
        fprintf(f, "#Reaction %zu, %zu bricks.\n", i + 1, r->last_brick);
        if(r->last_brick == 0) {
            continue;
        }
        sim_reaction_print_bricks(f, r, psr);
        fprintf(f, "\n\n");
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}
