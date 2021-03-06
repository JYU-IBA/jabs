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
#include <string.h>

#include <jibal.h>
#include <jibal_cs.h>

#include "defaults.h"
#include "generic.h"
#include "options.h"
#include "sample.h"
#include "simulation.h"
#include "spectrum.h"
#include "fit.h"
#include "script.h"
#include "script_session.h"
#include "script_command.h"

#include "prob_dist.h"
int main(int argc, char * const *argv) {
#ifdef DEBUG
    fprintf(stderr, "Note: this is a debug build!\n");
#endif
#if 0
    for(size_t n = 1; n < 10; n++) {
        fprintf(stderr, "n = %zu\n", n);
        prob_dist *pd = prob_dist_gaussian(n);
        prob_dist_free(pd);
        fprintf(stderr, "\n");
    }
    exit(EXIT_SUCCESS);
#endif


    int script_files = FALSE;
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return EXIT_FAILURE;
    }
    script_session *session = script_session_init(jibal, NULL);
    if(!session) {
        fprintf(stderr, "Can not initialize session.\n");
        return EXIT_FAILURE;
    }
    struct fit_data *fit_data = session->fit;
    simulation *sim = fit_data->sim;
    cmdline_options *cmd_opt = cmdline_options_init();
    read_options(jibal, sim, cmd_opt, &argc, &argv);
    if(!sim->beam_isotope) {
        sim->beam_isotope = jibal_isotope_find(jibal->isotopes, NULL, 2, 4); /* 4He default */
    }
    if(cmd_opt->fast) {
        sim_calc_params_defaults_fast(sim->params);
    }
    sim_calc_params_ds(sim->params, cmd_opt->ds);
    sim->params->stop_step_incident = cmd_opt->stop_step_incident;
    sim->params->stop_step_exiting = cmd_opt->stop_step_exiting;
    sim->params->depthsteps_max = cmd_opt->depthsteps_max;
    sim->rbs = cmd_opt->rbs;
    sim->erd = cmd_opt->erd;
    sim_calc_params_update(sim->params);
    session->output_filename = strdup_non_null(cmd_opt->output_filename);
    roi range = {.i_det = 0, .low = cmd_opt->fit_low, .high = cmd_opt->fit_high};
    fit_data_fit_range_add(fit_data, &range); /* We add just this one range (or none) */
    sim_sanity_check(sim);
    if(cmd_opt->exp_filename) {
        session->fit->exp[0] = spectrum_read_detector(cmd_opt->exp_filename, sim->det[0]);
        if(!session->fit->exp[0]) {
            fprintf(stderr, "Error! Can not open file \"%s\".\n", cmd_opt->exp_filename);
            return EXIT_FAILURE;
        }
    }
    if(cmd_opt->sample_filename) {
        fit_data->sm = sample_model_from_file(jibal, cmd_opt->sample_filename);
        if(!fit_data->sm) {
            fprintf(stderr, "Could not load a sample model from file \"%s\".\n", cmd_opt->sample_filename);
            return EXIT_FAILURE;
        }
    } else if(argc > 0 && strcmp(argv[0], "sample") == 0) {
        argc--;
        argv++;
        fit_data->sm  = sample_model_from_argv(jibal, &argc, &argv);
        if(argc != 0) {
            fprintf(stderr, "Error in reading sample model from command line (%i args remain)\n", argc);
            return EXIT_FAILURE;
        }
        argc = 0;
    }
    if(argc > 0) {
        script_files = TRUE;
    } else if(!fit_data->sm) { /* No sample file, no sample and no files given on command line, fallback to interactive (or script) mode */
        cmd_opt->interactive = TRUE;
    }
    if(!cmd_opt->interactive) {
        greeting(FALSE);
    }
    int status = 0;
    if(cmd_opt->interactive || script_files) {
        if(script_files) {
            for(int i = 0; i < argc; i++) {
                if(script_session_load_script(session, argv[i])) {
                    return EXIT_FAILURE;
                }
                status = script_process(session);
#ifdef DEBUG
                fprintf(stderr, "Script %i/%i given from command line has been processed. Status: %s\n", i+1, argc, script_command_status_to_string(status));
#endif
                if(status != SCRIPT_COMMAND_SUCCESS) {

                    return status;
                }
            }
        }
#ifdef DEBUG
        fprintf(stderr, "All scripts given from command line have been processed.\n");
#endif
        if(cmd_opt->interactive) {
            greeting(TRUE);
            if(script_session_load_script(session, NULL)) {
                return EXIT_FAILURE;
            }
            status = script_process(session);
        }
    } else { /* Non-interactive, pure command line mode. Run a single sim or fit. */
        if(!session->output_filename) {
            session->output_filename = strdup("-"); /* If no output filename given, set it to "-" (interpreted as stdout) */
        }
        if(fit_data->sim->rbs) {
            sim_reactions_add_auto(fit_data->sim, fit_data->sm, REACTION_RBS, sim_cs(fit_data->sim, REACTION_RBS)); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
        }
        if(fit_data->sim->erd) {
            sim_reactions_add_auto(fit_data->sim, fit_data->sm, REACTION_ERD, sim_cs(fit_data->sim, REACTION_ERD));
        }
        for(size_t i = 0; i < cmd_opt->n_reaction_filenames; i++) {
            sim_reactions_add_r33(fit_data->sim, jibal->isotopes, cmd_opt->reaction_filenames[i]);
        }
        if(cmd_opt->fit) {
            fprintf(stderr, "Running a fit in non-interactive mode.\n");
            status = script_fit(session, 1, &cmd_opt->fit_vars);
        } else {
            fprintf(stderr, "Running a simulation in non-interactive mode.\n");
            status = script_simulate(session, 0, NULL);
        }
    }
    script_session_free(session);
    cmdline_options_free(cmd_opt);
    jibal_free(jibal);
    return status;
}
