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

#include "options.h"
#include "sample.h"
#include "simulation.h"
#include "spectrum.h"
#include "fit.h"
#include "script.h"
#include "generic.h"
#include "defaults.h"

int main(int argc, char **argv) {
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
    sample_model *sm = fit_data->sm;
    cmdline_options *cmd_opt = cmdline_options_init();
    read_options(jibal, sim, cmd_opt, &argc, &argv);
    if(!sim->beam_isotope) {
        sim->beam_isotope = jibal_isotope_find(jibal->isotopes, NULL, 2, 4); /* 4He default */
    }
    sim->params = sim_calc_params_defaults(cmd_opt->ds, cmd_opt->fast);
    sim->params.stop_step_incident = cmd_opt->stop_step_incident;
    sim->params.stop_step_exiting = cmd_opt->stop_step_exiting;
    sim->params.depthsteps_max = cmd_opt->depthsteps_max;
    session->output_filename = strdup_non_null(cmd_opt->output_filename);
    session->bricks_out_filename = strdup_non_null(cmd_opt->bricks_filename);
    session->sample_out_filename = strdup_non_null(cmd_opt->sample_filename);
    session->detector_out_filename = strdup_non_null(cmd_opt->detector_out_filename);
    fit_range range = {.low = cmd_opt->fit_low, .high = cmd_opt->fit_high};
    fit_range_add(fit_data, &range); /* We add just this one range (or none) */
    fit_data->print_iters = cmd_opt->print_iters;
    sim_sanity_check(sim);
    if(cmd_opt->exp_filename) {
        session->fit->exp = read_experimental_spectrum(cmd_opt->exp_filename, sim->det);
        if(!session->fit->exp) {
            fprintf(stderr, "Error! Can not open file \"%s\".\n", cmd_opt->exp_filename);
            return EXIT_FAILURE;
        }
    }
    if(cmd_opt->sample_filename) {
        sm = sample_model_from_file(jibal, cmd_opt->sample_filename);
        if(!sm) {
            fprintf(stderr, "Could not load a sample model from file \"%s\".\n", cmd_opt->sample_filename);
            return EXIT_FAILURE;
        }
    } else if(argc > 0 && strcmp(argv[0], "sample") == 0) {
        argc--;
        argv++;
        sm  = sample_model_from_argv(jibal, argc, argv);
        if(!sm) {
            fprintf(stderr, "Error in reading sample model from command line.\n");
            return EXIT_FAILURE;
        }
        argc = 0;
    }
    if(argc > 0) {
        script_files = TRUE;
    } else if(!sm) { /* No sample file, no sample and no files given on command line, fallback to interactive (or script) mode */
        cmd_opt->interactive = TRUE;
    }
    if(sm) {
        sim->sample = sample_from_sample_model(sm);
    }
    if(!cmd_opt->interactive) {
        greeting(FALSE);
    }
    int status = 0;
    if(cmd_opt->interactive || script_files) {
        if(script_files) {
            for(int i = 0; i < argc; i++) {
                FILE *f_script = fopen(argv[i], "r");
                if(!f_script) {
                    fprintf(stderr, "Can not open script from file \"%s\".\n", argv[i]);
                    return EXIT_FAILURE;
                }
                fprintf(stderr, "\nRunning script \"%s\"\n\n", argv[i]);
                status = script_process(session, f_script);
                fclose(f_script);
                if(status) {
                    fprintf(stderr, "Error running script \"%s\"\n", argv[i]);
                    return status;
                }
                fprintf(stderr, "Finished running script \"%s\"\n", argv[i]);
            }
        }
        if(cmd_opt->interactive) {
            greeting(TRUE);
            status = script_process(session, stdin);
        }
    } else { /* Non-interactive, pure command line mode. Run a single sim or fit. */
        if(cmd_opt->fit) {
            status = script_fit(session, 1, &cmd_opt->fit_vars);
        } else {
            status = script_simulate(session, 0, NULL);
        }
    }
    script_session_free(session);
    cmdline_options_free(cmd_opt);
    jibal_free(jibal);
    return status;
}
