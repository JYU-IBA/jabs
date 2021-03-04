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
#include "layers.h"
#include "spectrum.h"
#include "fit.h"
#include "jabs.h"

int main(int argc, char **argv) {
    global_options global = {.jibal = NULL, .out_filename = NULL, .verbose = 0, .exp_filename = NULL,
                             .bricks_filename = NULL, .fit = 0, .fit_low = 0, .fit_high = 0, .fit_vars = NULL};
    simulation *sim = sim_init();
    clock_t start, end;
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return 1;
    }
    global.jibal = jibal;
    sim->beam_isotope = jibal_isotope_find(jibal->isotopes, NULL, 2, 4); /* Default: 4He */
    read_options(&global, sim, &argc, &argv);
    if(argc < 2) {
        usage();
        return EXIT_FAILURE;
    }
    fprintf(stderr, "JaBS version %s. Copyright (C) 2021 Jaakko Julin.\n", jabs_version());
    fprintf(stderr, "Compiled using JIBAL %s, current library version %s.\n", jibal_VERSION, jibal_version());
    fprintf(stderr, "JaBS comes with ABSOLUTELY NO WARRANTY.\n"
                    "This is free software, and you are welcome to redistribute it under certain conditions.\n"
                    "Run 'jabs -h' for more information.\n\n");
    sim_sanity_check(sim);

    gsl_histogram *exp = NULL;
    if(global.exp_filename) {
        exp = read_experimental_spectrum(global.exp_filename, 16384); /* TODO: number of channels? */
        set_experimental_spectrum_calibration(exp, sim); /* TODO: this is not really used anywhere. If it is used with fits etc it needs to be updated. */
    }

    FILE *f;
    if(global.out_filename) {
        f = fopen(global.out_filename, "w");
        if (!f) {
            fprintf(stderr, "Can't open file \"%s\" for output.\n", global.out_filename);
            return EXIT_FAILURE;
        }
    } else {
        f = stdout;
    }

    size_t n_layers = 0;
    jibal_layer **layers = read_layers(jibal, argc, argv, &n_layers);
    if(!layers)
        return EXIT_FAILURE;

    sample *sample = sample_from_layers(layers, n_layers);
    if(!sample || sample->n_isotopes == 0)
        return EXIT_FAILURE;
    sample_print(stderr, sample);

    reaction *reactions = make_rbs_reactions(sample, sim);
    fprintf(stderr, "\n");
    reactions_print(stderr, reactions);
    if(reactions[0].type == REACTION_NONE) {
        fprintf(stderr, "No reactions, nothing to do.\n");
        return EXIT_FAILURE;
    }

    if(assign_stopping(jibal->gsto, sim, sample)) {
        return EXIT_FAILURE;
    }
    jibal_gsto_print_assignments(jibal->gsto);
    jibal_gsto_load_all(jibal->gsto);
    fprintf(stderr, "Default RBS cross section model used: %s\n", jibal_cs_rbs_name(jibal->config));

    sim_calculate_geometry(sim);
    simulation_print(stderr, sim);
    fprintf(stderr, "\nSTARTING SIMULATION... NOW! Hold your breath!\n");
    fflush(stderr);
    start = clock();
    sim_workspace *ws = NULL;
    struct fit_stats fit_stats;
    if(global.fit) {
        struct fit_data fit_data;
        fit_data.n_iters_max = 150;
        fit_data.low_ch = global.fit_low;
        if(fit_data.low_ch <= 0)
            fit_data.low_ch = (int)(exp->n*0.1);
        fit_data.high_ch = global.fit_high;
        if(fit_data.high_ch <= 0)
            fit_data.high_ch = exp->n - 1;
        fprintf(stderr, "Fit range [%i, %i]\n", fit_data.low_ch, fit_data.high_ch);
        fit_data.jibal = jibal;
        fit_data.sim = sim;
        fit_data.exp = exp;
        fit_data.sample = NULL;
        fit_data.layers = layers;
        fit_data.n_layers = n_layers;
        fit_data.reactions = reactions;
        fit_data.ws = NULL;
        fit_data.fit_params = fit_params_new();
        add_fit_params(&global, sim, layers, n_layers, fit_data.fit_params);
        if(fit_data.fit_params->n == 0) {
            fprintf(stderr, "No parameters to fit!\n");
            return EXIT_FAILURE;
        }
        fit_stats = fit(exp, &fit_data);
        fprintf(stderr, "\nFinal parameters:\n");
        simulation_print(stderr, sim);
        fprintf(stderr, "\nFinal composition:\n");
        sample_print(stderr, fit_data.sample);
        ws = fit_data.ws;
        fprintf(stderr,"CPU time used for actual simulation: %.3lf s.\n", fit_data.cputime_actual);
        fprintf(stderr,"Per spectrum simulation: %.3lf ms.\n", 1000.0*fit_data.cputime_actual/fit_stats.n_evals);
        fit_params_free(fit_data.fit_params);
    } else {
        ws = sim_workspace_init(sim, reactions, sample, jibal);
        no_ds(ws, sample);
    }
    if(!ws) {
        fprintf(stderr, "Unexpected error.\n");
        return EXIT_FAILURE;
    }
    print_spectra(f, &global, ws, sample, exp);

    output_bricks(global.bricks_filename, ws);

    sim_workspace_free(ws);
    end = clock();
    double cputime_total =(((double) (end - start)) / CLOCKS_PER_SEC);
    fprintf(stderr, "...finished!\n\n");
    fprintf(stderr, "Total CPU time: %.3lf s.\n", cputime_total);
    if(f != stdout) {
        fclose(f);
    }
    layers_free(layers, n_layers);
    sim_free(sim);
    sample_free(sample);
    jibal_free(jibal);
    free(exp);
    return EXIT_SUCCESS;
}
