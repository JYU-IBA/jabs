#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <jibal.h>
#include <jibal_gsto.h>
#include <jibal_stop.h>
#include <jibal_stragg.h>
#include <jibal_kin.h>
#include <jibal_cs.h>
#include <gsl/gsl_histogram.h>

#include "sample.h"
#include "ion.h"
#include "brick.h"
#include "simulation.h"

#define ALPHA (0.0*C_DEG)
#define BETA (0.0*C_DEG)
#define THETA (170*C_DEG)
#define DETECTOR_RESOLUTION (15.0*C_KEV/C_FWHM)
#define PARTICLES_SR (1.0e12)

#define E_MIN (100.0*C_KEV)
#define N_LAYERS_MAX 100
#define HISTOGRAM_BIN (1.0*C_KEV)
#define STOP_STEP_INCIDENT (5.0*C_KEV)
#define STOP_STEP_EXITING (25.0*C_KEV)

#define CONCENTRATION_CUTOFF 1e-8

#define NUMBER_OF_SIMULATIONS 20


double stop_sample(sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, double x, double E) {
    double em=E/incident->mass;
    int i_isotope;

    double S1 = 0.0;

    get_concs(ws, sample, x, ws->c);
    for(i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        //double c = get_conc(sample, x, i_isotope); /* TODO: this could be sped up */
        double c = ws->c[i_isotope];
        if(c < CONCENTRATION_CUTOFF)
             continue;
        if(type == GSTO_STO_TOT) {
            S1 += c * (
                    jibal_gsto_get_em(ws->gsto, GSTO_STO_ELE, incident->Z, sample->isotopes[i_isotope]->Z, em)
                    +jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, sample->isotopes[i_isotope]->Z, sample->isotopes[i_isotope]->mass)
            );
        } else {
            S1 += c * (
                    jibal_gsto_get_em(ws->gsto, type, incident->Z, sample->isotopes[i_isotope]->Z, em)
            );
        }
    }
    //assert(S1 > 0.0);
    return S1;
}

#if 0
void recalculate_concs(sim_isotope *its, double x) {
    int i_isotope;
    double sum = 0.0;
    int n_isotopes = 0;
    for(i_isotope = 0; its[i_isotope].isotope != NULL; i_isotope++) {
        n_isotopes++;
        its[i_isotope].c = get_conc(x, &its[i_isotope]);
        sum += its[i_isotope].c;
    }
#ifdef DEBUG
    if(sum < 0.99 || sum > 1.01) {
        fprintf(stderr, "ISSUE AT x=%g with %i isotopes, x=%12.8lf, sum of concs is %g\n", x, n_isotopes, x, sum);
    }
#endif
}
#endif

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

void rbs(sim_workspace *ws, const simulation *sim, reaction *reactions, const sample *sample) {
    double x;
    ion ion = sim->ion;
    assert(sample->n_ranges);
    double thickness = sample->cranges[sample->n_ranges-1];
    double next_crossing = sample->cranges[1];
    double h_max;
    int i_range = 0;
    int i_reaction;
#ifdef DEBUG
    fprintf(stderr, "Thickness %g tfu, E = %g MeV\n", thickness/C_TFU, ion.E/C_MEV);
#endif
    for(i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        reaction *r = &reactions[i_reaction];
        r->p.E = sim->ion.E * r->K;
        r->E = r->p.E;
        r->S = 0.0;
        r->p.S = 0.0;
        r->stop = 0;
        ws->histos[i_reaction] = gsl_histogram_calloc_uniform(sim->n_channels, 0 * C_KEV, sim->histogram_bin * sim->n_channels); /* free'd by sim_workspace_free */

    }
    for (x = 0.0; x < thickness;) {
        while (i_range < sample->n_ranges-1 && x >= sample->cranges[i_range+1]) {
            i_range++;
#ifdef DEBUG
            fprintf(stderr, "Crossing to range %i = [%g, %g)\n", i_range, sample->cranges[i_range]/C_TFU, sample->cranges[i_range+1]/C_TFU);
#endif
            next_crossing = sample->cranges[i_range+1];
        }
        if(ion.E < E_MIN) {
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
        double h = stop_step(ws, &ion, sample, x, h_max, STOP_STEP_INCIDENT);
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
                double h_out = stop_step(ws, &r->p, sample, x_out-0.00001*C_TFU, h_out_max, STOP_STEP_EXITING); /* FIXME: 0.0001*C_TFU IS A STUPID HACK */
                x_out += h_out;
                if( r->p.E < E_MIN) {
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
            if(c > CONCENTRATION_CUTOFF && r->p.E > E_MIN) {/* TODO: concentration cutoff? TODO: it->E should be updated when we start calculating it again?*/
                double sigma = jibal_cross_section_rbs(ion.isotope, r->isotope, THETA, E_mean, JIBAL_CS_ANDERSEN);
                double Q = c * sim->p_sr_cos_alpha * sigma * h; /* TODO: worst possible approximation... */
#ifdef dfDEBUG
                fprintf(stderr, "    %s: E_scatt = %.3lf, E_out = %.3lf (prev %.3lf, sigma = %g mb/sr, Q = %g (c = %.4lf%%)\n",
                        r->isotope->name, E_back * r->K/C_KEV, r->p.E/C_KEV, r->E/C_KEV, sigma/C_MB_SR, Q, c*100.0);
#endif
                assert(sigma > 0.0);
                assert(r->p.S > 0.0);
                assert(Q < 1.0e7 || Q > 0.0);
                brick_int(sqrt(r->p.S + DETECTOR_RESOLUTION * DETECTOR_RESOLUTION), sqrt(r->S + DETECTOR_RESOLUTION * DETECTOR_RESOLUTION), r->p.E, r->E, ws->histos[i_reaction], Q);
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
        r->K = jibal_kin_rbs(sim->ion.mass, r->isotope->mass, sim->theta, '+'); /* TODO: this is too hard coded for RBS right now */
        r->p.E = sim->ion.E * r->K;
        r->E = r->p.E;
        r->p.S = 0.0;
        r->max_depth = sample->cranges[sample->n_ranges-1]-1000.0*C_TFU;

        r->stop = 0;
#ifdef DEBUG
        fprintf(stderr, "Reaction i=%i, with %s: K = %.5lf, max depth : %.3lf tfu\n", sim->n_reactions, r->isotope->name, r->K, r->max_depth/C_TFU);
#endif
        (*n_reactions)++;
    };
    return reactions;
}

jibal_layer **read_layers(jibal *jibal, int argc, char **argv, int *n_layers) {
    jibal_layer **layers = calloc(N_LAYERS_MAX, sizeof(jibal_layer *));
    *n_layers=0;
    while (argc >= 2 && *n_layers < N_LAYERS_MAX) {
        jibal_layer *layer = jibal_layer_new(jibal_material_create(jibal->elements, argv[0]),
                                             jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, argv[1]));
        if (!layer) {
            fprintf(stderr, "Not a valid layer: %s!\n", argv[0]);
            free(layers);
            return NULL;
        }

        layers[*n_layers] = layer;
        argc -= 2;
        argv += 2;
        (*n_layers)++;
    }
    return layers;
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

void print_spectra(FILE *f, const simulation *sim, const sim_workspace *ws) {
    int i, j;
    for(i = 0; i < sim->n_channels; i++) {
        double sum = 0.0;
        for (j = 0; j < sim->n_reactions; j++) {
            sum += ws->histos[j]->bin[i];
        }
        fprintf(f, "%4i %e", i, sum);
        for (j = 0; j < sim->n_reactions; j++) {
            fprintf(f, " %e", ws->histos[j]->bin[i]);
        }
        fprintf(f, "\n");
    }
}

int main(int argc, char **argv) {
    simulation sim;
    int i,j,k;
    clock_t start, end;
    double cpu_time_used;
#if 0
    erf_Q_test();
    return 0;
#endif
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return 1;
    }
    if(argc < 5) {
        fprintf(stderr, "Not enough arguments! Usage: %s: isotope energy file material thickness material2 thickness2...\n", argv[0]);
        return 1;
    }
    sim.incident = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    if(!sim.incident) {
        fprintf(stderr, "No isotope %s found.\n", argv[1]);
    }
    sim.ion.E = jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, argv[2]);
#ifdef DEBUG
    fprintf(stderr, "Beam E = %g MeV, incident ion Z = %i, mass = %.4lf u\n", sim.ion.E/C_MEV, sim.incident->Z, sim.incident->mass/C_U);
#endif
    if (sim.ion.E > 1000.0*C_MEV || sim.ion.E < 10*C_KEV) {
        fprintf(stderr, "Hmm...? Check your numbers.\n");
        return -1;
    }
    sim.histogram_bin = HISTOGRAM_BIN;
    sim.n_channels= ceil(1.1 * sim.ion.E / sim.histogram_bin);
    const char *filename = argv[3];
    FILE *f = fopen(filename, "w");
    if(!f) {
        fprintf(stderr, "Can't open file \"%s\" for output.\n", filename);
        return EXIT_FAILURE;
    }
    argc -= 4;
    argv += 4;

    int n_layers = 0;
    jibal_layer **layers = read_layers(jibal, argc, argv, &n_layers);
    if(!layers)
        return EXIT_FAILURE;

    sample sample = sample_from_layers( layers, n_layers);
    if(sample.n_isotopes == 0)
        return EXIT_FAILURE;

    sample_print(stderr, &sample);

    sim.alpha = ALPHA;
    sim.beta = BETA;
    sim.theta = THETA;
    sim.p_sr = PARTICLES_SR; /* TODO: particles * sr / cos(alpha) */
    sim.p_sr_cos_alpha = sim.p_sr / cos(sim.alpha);
    sim.ion.S = 0.0; /* TODO: initial beam broadening goes here */
    ion_set_isotope(&sim.ion, sim.incident);
    ion_set_angle(&sim.ion, ALPHA);

    reaction *reactions = make_rbs_reactions(&sample, &sim, &sim.n_reactions);
    if(assign_stopping(jibal->gsto, &sim, &sample)) {
        return EXIT_FAILURE;
    }
    jibal_gsto_print_assignments(jibal->gsto);
    jibal_gsto_load_all(jibal->gsto);

    start = clock();
    for (int n = 0; n < NUMBER_OF_SIMULATIONS; n++) {
        reaction *r = malloc(sim.n_reactions * sizeof(reaction));
        memcpy(r, reactions, sim.n_reactions * sizeof(reaction)); /* TODO: when we stop mutilating the reactions we can stop doing this */
        sim_workspace *ws = sim_workspace_init(&sim, &sample, jibal->gsto);
        rbs(ws, &sim, reactions, &sample);
#if 0
        if(n == NUMBER_OF_SIMULATIONS-1)
#else
        if(n == 0)
#endif
            print_spectra(f, &sim, ws);
        sim_workspace_free(ws);
        free(r);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    fprintf(stderr, "CPU time: %g ms per spectrum, for actual simulation of %i spectra.\n", cpu_time_used*1000.0/NUMBER_OF_SIMULATIONS, NUMBER_OF_SIMULATIONS);
    fclose(f);
    for(i = 0; i < n_layers; i++) {
        jibal_layer_free(layers[i]);
    }
    sample_free(&sample);
    free(layers);
    jibal_free(jibal);
    return 0;
}
