/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

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

#include "jabs_debug.h"
#include "defaults.h"
#include "generic.h"
#include "git.h"
#include "options.h"
#include "sample.h"
#include "simulation.h"
#include "spectrum.h"
#include "fit.h"
#include "script.h"
#include "script_session.h"
#include "script_command.h"
#include "idf2jbs.h"

int idf2jbs(int argc, char * const *argv) {
    char *filename_out = NULL;
    if(argc == 1) {
        idf_error idferr = idf_parse(argv[0], &filename_out);
        if(idferr == IDF2JBS_SUCCESS) {
            fprintf(stderr, "Success. Wrote script to file \"%s\"\n", filename_out);
            free(filename_out);
            return EXIT_SUCCESS;
        } else {
            fprintf(stderr, "IDF2JBS failed with error code %i (%s).\n", idferr, idf_error_code_to_str(idferr));
            return EXIT_FAILURE;
        }
    } else {
        fprintf(stderr, "Usage (idf2jbs): jabs idf2jbs <idf file>\nNote that idf file must have suffix .xml or .idf.\nOn successful run .jbs (simulation script) and .dat (spectrum) files will be created.\n");
        return EXIT_FAILURE;
    }
}

int main(int argc, char * const *argv) {
    DEBUGMSG("This is a debug build, version %s", jabs_version());
    if(git_populated()) {
        DEBUGMSG("%sGit build. branch %s, commit %s, date %s",
                 git_dirty() ? "Dirty " : "", git_branch(), git_commit_sha1(), git_commit_date());
    }
    if(argc >= 2 && strcmp(argv[1], "idf2jbs") == 0) {
        argc -= 2;
        argv += 2;
        return idf2jbs(argc, argv);
    }
    int script_files = FALSE;
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return EXIT_FAILURE;
    }
    script_session *session = script_session_init(jibal, NULL);
    if(!session) {
        fprintf(stderr, "Can not initialize session.\n");
        jibal_free(jibal);
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
    sim->params->incident_stop_params.step = cmd_opt->incident_stop_step;
    sim->params->exiting_stop_params.step = cmd_opt->exiting_stop_step;
    sim->params->n_bricks_max = cmd_opt->depthsteps_max;
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
                DEBUGMSG("Script %i/%i given from command line has been processed. Status: %s", i+1, argc, script_command_status_to_string(status));
                if(status != SCRIPT_COMMAND_SUCCESS) {
                    return status;
                }
            }
        }
        DEBUGSTR("All scripts given from command line have been processed.");
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
            sim_reactions_add_auto(fit_data->sim, fit_data->sm, REACTION_RBS, sim_cs(fit_data->sim, REACTION_RBS), TRUE); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
            sim_reactions_add_auto(fit_data->sim, fit_data->sm, REACTION_RBS_ALT, sim_cs(fit_data->sim, REACTION_RBS_ALT), TRUE);
        }
        if(sim_do_we_need_erd(fit_data->sim)) {
            sim_reactions_add_auto(fit_data->sim, fit_data->sm, REACTION_ERD, sim_cs(fit_data->sim, REACTION_ERD), TRUE);
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
    return status == SCRIPT_COMMAND_FAILURE ? EXIT_FAILURE : EXIT_SUCCESS;
}
