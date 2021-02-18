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
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <jibal.h>
#include <jibal_gsto.h>
#include <jibal_stop.h>
#include <jibal_kin.h>
#include <jibal_cs.h>
#include <gsl/gsl_histogram.h>

#include "options.h"
#include "sample.h"
#include "ion.h"
#include "brick.h"
#include "simulation.h"
#include "layers.h"
#include "spectrum.h"
#include "fit.h"

#define CONCENTRATION_CUTOFF 1e-8

#define NUMBER_OF_SIMULATIONS 1

double stop_sample(sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, double x, double E) {
    double em=E/incident->mass;
    int i_isotope;

    double S1 = 0.0;

    get_concs(ws, sample, x, ws->c);
    for(i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        if(ws->c[i_isotope] < CONCENTRATION_CUTOFF)
             continue;
        if(type == GSTO_STO_TOT) {
            S1 += ws->c[i_isotope] * (
                    jibal_gsto_get_em(ws->gsto, GSTO_STO_ELE, incident->Z, sample->isotopes[i_isotope]->Z, em)
                    +jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, sample->isotopes[i_isotope]->Z, sample->isotopes[i_isotope]->mass)
            );
        } else {
            S1 += ws->c[i_isotope] * (
                    jibal_gsto_get_em(ws->gsto, type, incident->Z, sample->isotopes[i_isotope]->Z, em)
            );
        }
    }
    //assert(S1 > 0.0);
    return S1;
}

double stop_step(sim_workspace *ws, ion *incident, const sample *sample, double x, double h_max, double step) {
    /* positive max depth step (h_max) also gives the direction. Energy step (step) should be negative if regular stopping is done. */
    double k1, k2, k3, k4, stop, dE, E, h_max_orig = h_max;
    int maxstep = 0;
    /* k1...k4 are slopes of energy loss (stopping) at various x (depth) and E. Note convention: positive values, i.e. -dE/dx! */
    E = incident->E;
    k1 = stop_sample(ws, incident, sample, ws->stopping_type, x, E);
    if(k1 < 0.001*C_EV_TFU) { /* Fail on positive values, zeroes (e.g. due to zero concentrations) and too small negative values */
#ifdef DEBUG
        fprintf(stderr, "stop_step returns 0.0, because k1 = %g eV/tfu (x = %.3lf tfu, E = %.3lg keV)\n", k1/C_EV_TFU, x/C_TFU, E/C_KEV);
#endif
        return 0.0;
    }
    h_max *=  incident->inverse_cosine; /* h_max is the perpendicular distance, but we can take bigger steps (note scaling elsewhere) */
    double h = (step / k1); /* step should always be positive, as well as k1  */
    if(h >= fabs(h_max)) {
        maxstep = 1;
        h = fabs(h_max) * 0.999;
    }
    double h_abs = h;
    assert(h_abs > 0.0);
    h = copysign(h, h_max); /* h is positive when going deeper, i.e. h + x is deeper than x . */

    double h_perp = h*incident->cosine; /* x + h_perp is the actual perpendicular depth */
    if(ws->rk4) {
        k2 = stop_sample(ws, incident, sample, ws->stopping_type, x + (h_perp / 2.0), E - (h_abs / 2.0) * k1);
        k3 = stop_sample(ws, incident, sample, ws->stopping_type, x + (h_perp / 2.0), E - (h_abs / 2.0) * k2);
        k4 = stop_sample(ws, incident, sample, ws->stopping_type, x + h_perp, E - h_abs * k3);
        stop = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    } else {
        stop = k1;
    }
    assert(stop > 0.0);
    dE =  -1.0*h_abs * stop; /* Energy change in thickness "h". It is always negative! */
#if 0
    fprintf(stderr, "%s stop = %.3lf eV/tfu ( x = %.3lf, h = %.3lf, h_max = %.3lf), E = %.3lf keV, h = %6.3lf  dE = %.5lf keV\n", incident->isotope->name, stop/C_EV_TFU, x/C_TFU, h/C_TFU, h_max/C_TFU, E/C_KEV, h/C_TFU, dE/C_KEV);
#endif
#ifdef DEBUG
    if(fabs(stop) < 0.1*C_EV_TFU) {
        fprintf(stderr, "Not good!\n");
        return 0.0;
    }
#endif
#ifndef STATISTICAL_STRAGGLING
   double s_ratio = stop_sample(ws, incident, sample, ws->stopping_type, x, E+dE)/k1; /* Ratio of stopping for non-statistical broadening. TODO: at x? */
#ifdef DEBUG
    //if((s_ratio)*(s_ratio) < 0.9 || (s_ratio)*(s_ratio) > 1.1) { /* Non-statistical broadening. */
    //   fprintf(stderr, "YIKES, s_ratio = %g, sq= %g\n", s_ratio, (s_ratio)*(s_ratio));
    //}
#endif
   incident->S *= (s_ratio)*(s_ratio);
#endif
   incident->S += h_abs*stop_sample(ws, incident, sample, GSTO_STO_STRAGG, x + (h_perp/2.0), (E+dE/2)); /* Straggling, calculate at mid-energy */

   assert(isnormal(incident->S));

    if(maxstep) {
        incident->E += dE/0.999;
        return h_max_orig;
    } else {
        incident->E += dE;
        return h*incident->cosine;
    }
      /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}

void simulate(sim_workspace *ws, const simulation *sim, reaction *reactions, const sample *sample) {
    double x;
    ion ion = sim->ion;
    assert(sample->n_ranges);
    double thickness = sample->cranges[sample->n_ranges-1];
    double next_crossing = sample->cranges[1];
    double h_max;
    int i_range = 0;
    int i_reaction;
    for(i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        reaction *r = &reactions[i_reaction];
        r->p.E = sim->ion.E * r->K;
        r->E = r->p.E;
        r->S = 0.0;
        r->p.S = 0.0;
        r->stop = 0;
    }
    for (x = 0.0; x < thickness;) {
        while (i_range < sample->n_ranges-1 && x >= sample->cranges[i_range+1]) {
            i_range++;
#ifdef DEBUG
            fprintf(stderr, "Crossing to range %i = [%g, %g)\n", i_range, sample->cranges[i_range]/C_TFU, sample->cranges[i_range+1]/C_TFU);
#endif
            next_crossing = sample->cranges[i_range+1];
        }
        if(ion.E < sim->emin) {
            fprintf(stderr, "Return due to low energy.\n");
            break;
        }
        h_max = next_crossing - x;
        if(h_max < 0.0001*C_TFU) {
            x += 0.0001*C_TFU;
            fprintf(stderr, "Step too small. Let's nudge forward! x=%lf tfu\n", x/C_TFU);
            continue;
        }

        double E_front = ion.E;
        double h = stop_step(ws, &ion, sample, x, h_max, sim->stop_step_incident);
        assert(h > 0.0);
        /* DEPTH BIN [x, x+h) */
        double E_back = ion.E;
#ifdef DEBUG
        double E_diff = E_front-E_back;
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3lf keV to  %8.3lf keV (diff %6.4lf keV)\n", x/C_TFU, (x+h)/C_TFU, E_front/C_KEV, ion.E/C_KEV, E_diff/C_KEV);
#endif
        double S_back = ion.S;
        double E_mean = (E_front+E_back)/2.0;
#ifdef DEBUG
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(ion.S) / C_KEV);
#endif
        for (i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
            reaction *r = &reactions[i_reaction];
            if(r->stop)
                continue;
            r->p.E = ion.E * r->K;
            r->p.S = ion.S * r->K;
            assert(r->p.E > 0.0);

            if(x >= r->max_depth) {
#ifdef DEBUG
                fprintf(stderr, "Reaction %i with %s stops, because maximum depth is reached.\n", i_reaction, reactions[i_reaction].isotope->name);
#endif
                r->stop = 1;
                continue;
            }

            double x_out;
            int i_range_out = i_range;
            for (x_out = x + h; x_out > 0.0;) { /* Calculate energy and straggling of backside of slab */
                double h_out_max = sample->cranges[i_range_out] - x_out;
                double h_out = stop_step(ws, &r->p, sample, x_out-0.00001*C_TFU, h_out_max, sim->stop_step_exiting); /* FIXME: 0.0001*C_TFU IS A STUPID HACK */
                x_out += h_out;
                if( r->p.E < sim->emin) {
                    r->stop = 1;
#ifdef DEBUG
                    fprintf(stderr, "Reaction %i with %s: Energy below EMIN when surfacing from %.3lf tfu, break break. Last above was %.3lf keV\n",i_reaction, reactions[i_reaction].isotope->name, (x+h)/C_TFU, r->E/C_KEV);
#endif
                    break;
                }
                assert(h_out < 0.0);
                while (i_range_out > 0 && x_out <= sample->cranges[i_range_out]) {
                    i_range_out--;
#ifdef DEBUG_VERBOSE
                    fprintf(stderr, "Outgoing from reaction %i crossing to range %i = [%g, %g) when x_out = %g tfu.\n", i_reaction, i_range_out, sample->cranges[i_range_out]/C_TFU, sample->cranges[i_range_out+1]/C_TFU, x_out/C_TFU);
#endif
                }
            }
            double c = get_conc(ws, sample, x+h/2.0, r->i_isotope); /* TODO: x+h/2.0 is actually exact for linearly varying concentration profiles. State this clearly somewhere. */
            if(c > CONCENTRATION_CUTOFF && r->p.E > sim->emin) {/* TODO: concentration cutoff? TODO: it->E should be updated when we start calculating it again?*/
                double sigma = jibal_cross_section_rbs(ion.isotope, r->isotope, sim->theta, E_mean, JIBAL_CS_ANDERSEN);
                double Q = c * ws->p_sr_cos_alpha * sigma * h; /* TODO: worst possible approximation... */
#ifdef dfDEBUG
                fprintf(stderr, "    %s: E_scatt = %.3lf, E_out = %.3lf (prev %.3lf, sigma = %g mb/sr, Q = %g (c = %.4lf%%)\n",
                        r->isotope->name, E_back * r->K/C_KEV, r->p.E/C_KEV, r->E/C_KEV, sigma/C_MB_SR, Q, c*100.0);
#endif
                assert(sigma > 0.0);
                assert(r->p.S > 0.0);
                assert(Q < 1.0e7 || Q > 0.0);
                brick_int(sqrt(r->p.S + sim->energy_resolution), sqrt(r->S + sim->energy_resolution), r->p.E, r->E, ws->histos[i_reaction], Q);
            }
            r->E = r->p.E;
            r->S = r->p.S;
        }
#if 0
        fprintf(stderr, "Surf: from x_out = %g tfu, gives E_out = %g MeV (prev was %g MeV, diff %g keV), stragg = %g keV\n", (x-h)/C_TFU, E_out/C_MEV, E_out_prev/C_MEV, (E_diff)/C_KEV, sqrt(S_out)/C_KEV);
        fprintf(stderr, "\n");
#endif
        if(!isnormal(ion.E)) {
            fprintf(stderr, "SOMETHING DOESN'T LOOK RIGHT HERE.\n");
            return;
        }
        x += h;
        ion.S = S_back;
        ion.E = E_back;
    }
}

reaction *make_rbs_reactions(const sample *sample, const simulation *sim, int *n_reactions) { /* Note that sim->ion needs to be set! */
    int i;
    *n_reactions = 0; /* we calculate this */
    reaction * reactions = malloc(sample->n_isotopes*sizeof(reaction)); /* TODO: possible memory leak */
    for(i = 0; i < sample->n_isotopes; i++) {
        reaction *r = &reactions[sim->n_reactions];
        r->isotope = sample->isotopes[i];
        r->i_isotope = i;
        ion_set_isotope(&r->p, sim->ion.isotope); /* TODO: for RBS, yes... */
        ion_set_angle(&r->p, sim->beta);
        double theta_max=asin(r->isotope->mass/sim->ion.isotope->mass);
        if(sim->ion.isotope->mass >= r->isotope->mass && sim->theta > theta_max) {
#ifdef DEBUG
            fprintf(stderr, "RBS with %s is not possible (theta max %g deg)\n", r->isotope->name, theta_max);
#endif
            continue;
        }
        r->K = jibal_kin_rbs(sim->ion.mass, r->isotope->mass, sim->theta, '+'); /* TODO: this is too hard coded for RBS right now */
        r->p.E = sim->ion.E * r->K;
        r->E = r->p.E;
        r->p.S = 0.0;
        r->max_depth = sample_isotope_max_depth(sample, r->i_isotope);
        r->stop = 0;
        (*n_reactions)++;
    };
    return reactions;
}

int assign_stopping(jibal_gsto *gsto, simulation *sim, sample *sample) {
    int i;
    for(i = 0; i < sample->n_isotopes; i++) {
        if (!jibal_gsto_auto_assign(gsto, sim->ion.Z, sample->isotopes[i]->Z)) { /* TODO: this is only for RBS */
            fprintf(stderr, "Can not assign stopping.\n");
            return 1;
        }
    }
    return 0;
}

void print_spectra(FILE *f, const global_options *global,  const simulation *sim, const sim_workspace *ws, const sample *sample, const reaction *reactions, const gsl_histogram *exp) {
    int i, j;
    int csv = 0;
    char sep = ' ';
    if(global->out_filename) {
        size_t l = strlen(global->out_filename);
        if(l > 4 && strncmp(global->out_filename+l-4, ".csv", 4) == 0) { /* For CSV: print header line */
            csv = 1;
            sep = ','; /* and set the separator! */
            fprintf(f, "\"Channel\",\"Simulated\"");
            if(exp) {
                fprintf(f, ",\"Experimental\"");
            }
            for(j = 0; j < sim->n_reactions; j++) {
                const reaction *r = &reactions[j];
                fprintf(f, ",\"%s\"", sample->isotopes[r->i_isotope]->name);
            }
            fprintf(f, "\n");
        }
    }
    for(i = 0; i < ws->n_channels; i++) {
        double sum = 0.0;
        for (j = 0; j < sim->n_reactions; j++) { /* Sum comes always first, which means we have to compute it first. */
            sum += ws->histos[j]->bin[i];
        }
        if(sum == 0.0) {
            fprintf(f, "%i%c0", i, sep); /* Tidier output with a clean zero */
        } else {
            fprintf(f, "%i%c%e", i, sep, sum);
        }
        if(exp) {
            if(i < exp->n) {
                fprintf(f, "%c%g", sep, exp->bin[i]);
            } else {
                fprintf(f, "%c0", sep);
            }
        }
        for (j = 0; j < sim->n_reactions; j++) {
            if(ws->histos[j]->bin[i] == 0.0) {
                fprintf(f,"%c0", sep);
            } else {
                fprintf(f, "%c%e", sep, ws->histos[j]->bin[i]);
            }
        }
        fprintf(f, "\n");
    }
}

void add_fit_params(global_options *global, simulation *sim, jibal_layer **layers, const int n_layers, fit_params *params) {
#ifdef DEBUG
    fprintf(stderr, "fitvars = %s\n", global->fit_vars);
#endif
    if(!global->fit_vars)
        return;
    char *token, *s, *s_orig;
    s_orig = s = strdup(global->fit_vars);
    assert(s != NULL);
    while ((token = strsep(&s, ",")) != NULL) { /* parse comma separated list of parameters to fit */
#ifdef DEBUG
        fprintf(stderr, "Thing to fit: \"%s\"\n", token);
#endif
        if(strcmp(token, "calib") == 0) {
            fit_params_add_parameter(params, &sim->energy_slope); /* TODO: prevent adding already added things */
            fit_params_add_parameter(params, &sim->energy_offset);
            fit_params_add_parameter(params, &sim->energy_resolution);
        }
        if(strcmp(token, "slope") == 0) {
            fit_params_add_parameter(params, &sim->energy_slope);
        }
        if(strcmp(token, "offset") == 0) {
            fit_params_add_parameter(params, &sim->energy_offset);
        }
        if(strcmp(token, "reso") == 0) {
            fit_params_add_parameter(params, &sim->energy_resolution);
        }
        if(strcmp(token, "fluence") == 0) {
            fit_params_add_parameter(params, &sim->p_sr);
        }
        if(strncmp(token, "thickness", 9) == 0 && strlen(token) > 9) {
            int i_layer = atoi(token+9);
            if(i_layer >= 1 && i_layer <= n_layers) {
                fit_params_add_parameter(params, &layers[i_layer-1]->thickness);
            } else {
                fprintf(stderr, "No layer %i (parsed from \"%s\")\n", i_layer, token);
            }
        }
    }
    free(s_orig);
}

int main(int argc, char **argv) {
    fprintf(stderr, "JaBS version %s. Copyright (C) 2021 Jaakko Julin.\n", jabs_version());
    fprintf(stderr, "JaBS comes with ABSOLUTELY NO WARRANTY.\n"
                    "This is free software, and you are welcome to redistribute it under certain conditions.\n"
                    "Run 'jabs -h' for more information.\n\n");
    global_options global = {.jibal = NULL, .out_filename = NULL, .verbose = 0, .exp_filename = NULL,
                             .fit = 0, .fit_low = 0, .fit_high = 0, .fit_vars = NULL};
    simulation *sim = sim_init();
    clock_t start, end;
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return 1;
    }
    global.jibal = jibal;
    ion_set_isotope(&sim->ion, jibal_isotope_find(jibal->isotopes, NULL, 2, 4)); /* Default: 4He */
    read_options(&global, sim, &argc, &argv);
    if(argc < 2) {
        usage();
        return EXIT_FAILURE;
    }
    sim_sanity_check(sim);

    gsl_histogram *exp = NULL;
    if(global.exp_filename) {
        exp = read_experimental_spectrum(global.exp_filename, 16384); /* TODO: number of channels? */
        set_experimental_spectrum_calibration(exp, sim);
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

    reaction *reactions = make_rbs_reactions(sample, sim, &sim->n_reactions);
    fprintf(stderr, "\n");
    reactions_print(stderr, reactions, sim->n_reactions);

    if(assign_stopping(jibal->gsto, sim, sample)) {
        return EXIT_FAILURE;
    }
    jibal_gsto_print_assignments(jibal->gsto);
    jibal_gsto_load_all(jibal->gsto);


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
            fit_data.high_ch = exp->n-1;
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
        if (fit_data.ws) {
            print_spectra(f, &global, sim, fit_data.ws, sample, reactions, exp);
        }
        ws = fit_data.ws;
        fprintf(stderr,"CPU time used for actual simulation: %.3lf s.\n", fit_data.cputime_actual);
        fprintf(stderr,"Per spectrum simulation: %.3lf ms.\n", 1000.0*fit_data.cputime_actual/fit_stats.n_evals);
        fit_params_free(fit_data.fit_params);
    } else {
        ws = sim_workspace_init(sim, sample, jibal->gsto);
        simulate(ws, sim, reactions, sample);
    }
    if(!ws) {
        fprintf(stderr, "Unexpected error.\n");
        return EXIT_FAILURE;
    }
    print_spectra(f, &global, sim, ws, sample, reactions, exp);
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
