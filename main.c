/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See LICENSE.txt for the full license.

 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <jibal.h>
#include <jibal_cs.h>
#include <jibal_defaults.h>
#include <gsl/gsl_histogram.h>

#include "options.h"
#include "sample.h"
#include "simulation.h"
#include "spectrum.h"
#include "fit.h"
#include "jabs.h"
#include "script.h"
#include "defaults.h"

int main(int argc, char **argv) {
    fprintf(stderr, "JaBS version %s. Copyright (C) 2021 Jaakko Julin.\n", jabs_version()); /* These are printed when running non-interactively with just command line parameters */
    fprintf(stderr, "Compiled using JIBAL %s, current library version %s.\n\n", JIBAL_VERSION, jibal_version());
    global_options *global = global_options_alloc();
    simulation *sim = sim_init();
    clock_t start, end;
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error,
                jibal_error_string(jibal->error));
        return 1;
    }
    global->jibal = jibal;
    sim->beam_isotope = jibal_isotope_find(jibal->isotopes, NULL, 2, 4); /* Default: 4He */
    read_options(global, sim, &argc, &argv);
    sim_sanity_check(sim);

    gsl_histogram *exp = NULL;
    if(global->exp_filename) {
        exp = read_experimental_spectrum(global->exp_filename, sim->det);
        if(!exp) {
            fprintf(stderr, "Error! Can not open file \"%s\".\n", global->exp_filename);
            return EXIT_FAILURE;
        }
    }
    sample_model *sm;
    if(global->sample_filename) {
        sm = sample_model_from_file(jibal, global->sample_filename);
        if(!sm) {
            fprintf(stderr, "Could not load a sample model from file \"%s\".\n", global->sample_filename);
            return EXIT_FAILURE;
        }
        if(argc != 0) {
            fprintf(stderr, "Unexpected command line parameters, total %i, starting with %s.\n", argc, *argv);
            return EXIT_FAILURE;
        }
    } else if(argc > 0) {
        sample_model *sm_raw  = sample_model_from_argv(jibal, argc, argv);
#ifdef DEBUG
        fprintf(stderr, "Sample model, as read from command line.\n");
        sample_model_print(stderr, sm_raw);
#endif
        sm = sample_model_split_elements(sm_raw);
        sample_model_free(sm_raw);
        if(!sm) {
            fprintf(stderr, "Error in reading sample model from command line.\n");
            return EXIT_FAILURE;
        }
    } else {
        return script_process(jibal, stdin);
    }
    fputs(COPYRIGHT_STRING, stderr);
    sample *sample = sample_from_sample_model(sm);
    if(global->verbose) {
        sample_model_print(stderr, sm);
        fprintf(stderr, "\nSimplified sample model for simulation:\n");
        sample_print(stderr, sample, global->print_isotopes);
    }
    if(sm->n_ranges == 0 || sm->n_materials == 0) {
        fprintf(stderr, "Can not simulate nothing.\n");
    }
    reaction **reactions = NULL;
    if(global->verbose) {
        fprintf(stderr, "Default RBS cross section model used: %s\n", jibal_cross_section_name(jibal->config->cs_rbs));
        fprintf(stderr, "Default ERD cross section model used: %s\n", jibal_cross_section_name(jibal->config->cs_erd));
        fprintf(stderr, "\n");
    }
    reactions = make_reactions(sample, sim, global->rbs?jibal->config->cs_rbs:JIBAL_CS_NONE, global->erd?jibal->config->cs_erd:JIBAL_CS_NONE);
    if(global->reaction_filenames) {
        if(process_reaction_files(jibal->isotopes, reactions, global->reaction_filenames, global->n_reaction_filenames)) {
            fprintf(stderr, "Could not process all reaction files. Aborting.\n");
            return EXIT_FAILURE;
        }
    }

    if(!reactions || reactions[0] == NULL ) {
        fprintf(stderr, "No reactions, nothing to do.\n");
        return EXIT_FAILURE;
    } else {
        if(global->verbose) {
            fprintf(stderr, "%zu reactions.\n", reaction_count(reactions));
            reactions_print(stderr, reactions);
        }
    }

    if(assign_stopping(jibal->gsto, sim, sample, reactions)) {
        return EXIT_FAILURE;
    }
    if(global->verbose) {
        jibal_gsto_print_assignments(jibal->gsto);
        jibal_gsto_print_files(jibal->gsto, 1);
    }
    jibal_gsto_load_all(jibal->gsto);
    simulation_print(stderr, sim);
    fprintf(stderr, "\nSTARTING SIMULATION... NOW! Hold your breath!\n");
    fflush(stderr);
    start = clock();
    sim_workspace *ws = NULL;
    struct fit_stats fit_stats;
    if(global->fit) {
        fit_data *f_data = fit_data_new(jibal, sim, exp, sm, reactions, global->fit_vars, global->fit_low, global->fit_high, global->print_iters);
        if(!f_data) {
            fprintf(stderr, "No parameters to fit!\n");
            return EXIT_FAILURE;
        }
        fprintf(stderr, "Fit range [%lu, %lu]\n", f_data->low_ch, f_data->high_ch);
        fit_stats = fit(exp, f_data);
        fprintf(stderr, "\nFinal parameters:\n");
        simulation_print(stderr, sim);
        fprintf(stderr, "\nFinal composition:\n");
        sample_print(stderr, f_data->sample, global->print_isotopes);
        sample_areal_densities_print(stderr, f_data->sample, global->print_isotopes);
        fprintf(stderr, "\nFinal sample model:\n");
        sample_model_print(stderr, sm);
        ws = f_data->ws;
        fprintf(stderr,"CPU time used for actual simulation: %.3lf s.\n", fit_stats.cputime_actual);
        fprintf(stderr,"Per spectrum simulation: %.3lf ms.\n", 1000.0*fit_stats.cputime_actual/fit_stats.n_evals);
        fit_data_free(f_data);
    } else {
        ws = sim_workspace_init(sim, reactions, sample, jibal);
        simulate_with_ds(ws);
    }
    if(!ws) {
        fprintf(stderr, "ERROR! Simulation workspace does not exist after the simulation. This implies something failed spectacularly.\nNo output generated.\n");
        return EXIT_FAILURE;
    }
    if(exp) {
        set_spectrum_calibration(exp, ws->sim.det); /* Update the experimental spectra to final calibration */
    }
    print_spectra(global->out_filename, ws, exp);

    output_bricks(global->bricks_filename, ws);

    if(global->detector_out_filename) {
        FILE *f_det;
        if((f_det = fopen(global->detector_out_filename, "w"))) {
            detector_print(f_det, ws->sim.det);
        } else {
            fprintf(stderr, "Could not write detector to file \"%s\".\n", global->detector_out_filename);
        }
    }
    if(global->sample_out_filename) {
        FILE *f_sout;
        if((f_sout = fopen(global->sample_out_filename, "w"))) {
            sample_model_print(f_sout, sm);
        } else {
            fprintf(stderr, "Could not write sample to file \"%s\".\n", global->sample_out_filename);
        }
    }

    sim_workspace_free(ws);
    end = clock();
    double cputime_total =(((double) (end - start)) / CLOCKS_PER_SEC);
    fprintf(stderr, "...finished!\n\n");
    fprintf(stderr, "Total CPU time: %.3lf s.\n", cputime_total);
    for(reaction **r = reactions; *r != NULL; r++) {
        reaction_free(*r);
    }
    free(reactions);
    sample_model_free(sm);
    sim_free(sim);
    sample_free(sample);
    jibal_free(jibal);
    free(exp);
    global_options_free(global);
    return EXIT_SUCCESS;
}
