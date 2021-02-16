#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include <jibal.h>
#include <jibal_gsto.h>
#include <jibal_stop.h>
#include <jibal_kin.h>
#include <jibal_cs.h>
#include <gsl/gsl_histogram.h>

#include "defaults.h"
#include "sample.h"
#include "ion.h"
#include "brick.h"
#include "simulation.h"

#define ENERGY (2.0*C_MEV)
#define ALPHA (15.0*C_DEG)
#define BETA (0.0*C_DEG)
#define THETA (165.0*C_DEG)
#define DETECTOR_RESOLUTION (15.0*C_KEV/C_FWHM)
#define PARTICLES_SR (1.0e12)

#define E_MIN (100.0*C_KEV)
#define N_LAYERS_MAX 100
#define HISTOGRAM_BIN (1.0*C_KEV)
#define STOP_STEP_INCIDENT (5.0*C_KEV)
#define STOP_STEP_EXITING (25.0*C_KEV)

#define CONCENTRATION_CUTOFF 1e-8

#define NUMBER_OF_SIMULATIONS 1

#define USAGE_STRING "Usage simu [-E <energy>] <material1> <thickness1> [<material2> <thickness2> ...]\n\nExample: simu -E 2MeV --alpha 10deg --beta 0deg -theta 170deg Au 500tfu SiO2 1000tfu Si 10000tfu\n"

const char *simu_version() {
        return simu_VERSION;
}

void usage() {
    fprintf(stderr, USAGE_STRING);
}

typedef struct {
    jibal *jibal;
    int verbose;
    char *out_filename;
} global_options;

void read_options(global_options *global, simulation *sim, int *argc, char ***argv) {
    static struct option long_options[] = {
            {"help",      no_argument,       NULL, 'h'},
            {"version",   no_argument,       NULL, 'V'},
            {"verbose",   optional_argument, NULL, 'v'},
            {"out",       required_argument, NULL, 'o'},
            {"ion",       required_argument, NULL, 'I'},
            {"energy",    required_argument, NULL, 'E'},
            {"alpha",     required_argument, NULL, 'a'},
            {"beta",      required_argument, NULL, 'b'},
            {"theta",     required_argument, NULL, 't'},
            {"fluence",   required_argument, NULL, 'F'},
            {"resolution",required_argument, NULL, 'R'},
            {"step_incident",required_argument, NULL, 'S'},
            {"step_exiting",required_argument, NULL, '0'},
            {"fast", optional_argument, NULL, 'f'},
            {NULL, 0,                NULL,   0}
    };
    while (1) {
        int option_index = 0;
        char c = getopt_long(*argc, *argv, "hvVE:o:a:b:t:I:F:R:S:f", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case 'f':
                if (optarg)
                    sim->fast = atoi(optarg);
                else
                    sim->fast++;
                break;
            case '0':
                sim->stop_step_exiting = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'S':
                sim->stop_step_incident = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'a':
                sim->alpha = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 'b':
                sim->beta = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 't':
                sim->theta = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
                break;
            case 'V':
                printf("%s\n", simu_VERSION);
                exit(EXIT_SUCCESS);
                break; /* Unnecessary */
            case 'v':
                if (optarg)
                    global->verbose = atoi(optarg);
                else
                    global->verbose++;
                break;
            case 'o':
                global->out_filename = optarg;
                break;
            case 'E':
                sim->ion.E = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'I':
                ion_set_isotope(&sim->ion, jibal_isotope_find(global->jibal->isotopes, optarg, 0, 0));
                if(!sim->ion.isotope) {
                    fprintf(stderr, "%s is not a valid isotope.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'F':
                sim->p_sr = strtod(optarg, NULL);
                break;
            case 'R':
                sim->energy_resolution = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg)/C_FWHM;
                sim->energy_resolution *= sim->energy_resolution; /* square */
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
                break;
        }
    }
    *argc -= optind;
    *argv += optind;
}

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

void rbs(sim_workspace *ws, const simulation *sim, reaction *reactions, const sample *sample) {
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
        ws->histos[i_reaction] = gsl_histogram_calloc_uniform(sim->n_channels, sim->energy_offset, sim->energy_offset+sim->energy_slope * sim->n_channels); /* free'd by sim_workspace_free */

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

void print_spectra(FILE *f, const global_options *global,  const simulation *sim, const sim_workspace *ws, const sample *sample, const reaction *reactions) {
    int i, j;
    int csv = 0;
    char sep = ' ';
    if(global->out_filename) {
        size_t l = strlen(global->out_filename);
        if(l > 4 && strncmp(global->out_filename+l-4, ".csv", 4) == 0) {
            csv = 1;
            sep = ',';
            fprintf(f, "\"Channel\",\"Sum\"");
            for(j = 0; j < sim->n_reactions; j++) {
                const reaction *r = &reactions[j];
                fprintf(f, ",\"%s\"", sample->isotopes[r->i_isotope]->name);
            }
            fprintf(f, "\n");
        }
    }
    for(i = 0; i < sim->n_channels; i++) {
        double sum = 0.0;
        for (j = 0; j < sim->n_reactions; j++) {
            sum += ws->histos[j]->bin[i];
        }
        if(sum == 0.0) {
            fprintf(f, "%i%c0", i, sep);
        } else {
            fprintf(f, "%i%c%e", i, sep, sum);
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

void layers_free(jibal_layer **layers, int n_layers) {
    int i;
    for(i = 0; i < n_layers; i++) {
        jibal_layer_free(layers[i]);
    }
    free(layers);
}

void reactions_print(FILE *f, reaction *reactions, int n_reactions) {
    int i;
    for (int i = 0; i < n_reactions; i++) {
        reaction *r = &reactions[i];
        fprintf(stderr, "Reaction %3i/%i: RBS with %5s (isotope id %3i): K = %7.5lf, max depth = %9.3lf tfu\n", i+1, n_reactions,
                r->isotope->name, r->i_isotope,
                r->K, r->max_depth / C_TFU);
    }

}

void simulation_print(FILE *f, const simulation *sim) {
    fprintf(stderr, "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->ion.isotope->name, sim->ion.isotope->Z, sim->ion.isotope->A, sim->ion.isotope->mass/C_U);
    fprintf(stderr, "E = %.3lf\n", sim->ion.E/C_MEV);
    fprintf(stderr, "alpha = %.3lf deg\n", sim->alpha/C_DEG);
    fprintf(stderr, "beta = %.3lf deg\n", sim->beta/C_DEG);
    fprintf(stderr, "theta = %.3lf deg\n", sim->theta/C_DEG);
    fprintf(stderr, "particles * sr = %e\n", sim->p_sr);
    fprintf(stderr, "calibration offset = %.3lf keV\n", sim->energy_offset/C_KEV);
    fprintf(stderr, "calibration slope = %.5lf keV\n", sim->energy_slope/C_KEV);
    fprintf(stderr, "detector resolution = %.3lf keV FWHM\n", sqrt(sim->energy_resolution)*C_FWHM/C_KEV);
    fprintf(stderr, "step for incident ions = %.3lf keV\n", sim->stop_step_incident/C_KEV);
    fprintf(stderr, "step for exiting ions = %.3lf keV\n", sim->stop_step_exiting/C_KEV);
    fprintf(stderr, "fast level = %i\n", sim->fast);
}

int main(int argc, char **argv) {
    global_options global = {.jibal = NULL, .out_filename = NULL, .verbose = 0};
    simulation sim = {.alpha = ALPHA, .beta = BETA, .theta = THETA,
                      .p_sr = PARTICLES_SR, .energy_resolution = DETECTOR_RESOLUTION*DETECTOR_RESOLUTION,
                      .stop_step_incident = STOP_STEP_INCIDENT, .stop_step_exiting = STOP_STEP_EXITING,
                      .fast = 0,
                      .ion = {.E = ENERGY }};
    sample sample;
    clock_t start, end;
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return 1;
    }
    global.jibal = jibal;
    sim.ion.E = 2.0 * C_MEV;
    ion_set_isotope(&sim.ion, jibal_isotope_find(jibal->isotopes, NULL, 2, 4)); /* Default: 4He */
    ion_set_angle(&sim.ion, 0.0 * C_DEG);
    read_options(&global, &sim, &argc, &argv);
    if(argc < 2) {
        usage();
        return EXIT_FAILURE;
    }
    if(!sim.ion.isotope) {
        fprintf(stderr, "No valid isotope given for the beam.\n");
    }
    if (sim.ion.E > 1000.0*C_MEV || sim.ion.E < 10*C_KEV) {
        fprintf(stderr, "Hmm...? Check your numbers.\n");
        return -1;
    }
    sim.energy_slope = HISTOGRAM_BIN;
    sim.energy_offset = 0.0*C_KEV;
    sim.n_channels= ceil((1.1 * sim.ion.E - sim.energy_offset)/ sim.energy_slope);
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

    int n_layers = 0;
    jibal_layer **layers = read_layers(jibal, argc, argv, &n_layers);
    if(!layers)
        return EXIT_FAILURE;

    sample = sample_from_layers(layers, n_layers);
    if(sample.n_isotopes == 0)
        return EXIT_FAILURE;
    sample_print(stderr, &sample);

    sim.p_sr_cos_alpha = sim.p_sr / cos(sim.alpha);
    sim.ion.S = 0.0; /* TODO: initial beam broadening goes here */
    ion_set_angle(&sim.ion, sim.alpha);

    reaction *reactions = make_rbs_reactions(&sample, &sim, &sim.n_reactions);
    fprintf(stderr, "\n");
    reactions_print(stderr, reactions, sim.n_reactions);

    if(assign_stopping(jibal->gsto, &sim, &sample)) {
        return EXIT_FAILURE;
    }
    jibal_gsto_print_assignments(jibal->gsto);
    jibal_gsto_load_all(jibal->gsto);


    simulation_print(stderr, &sim);
    fprintf(stderr, "\nSTARTING SIMULATION... NOW! Hold your breath!\n");
    fflush(stderr);
    start = clock();
    for (int n = 0; n < NUMBER_OF_SIMULATIONS; n++) {
        reaction *r = malloc(sim.n_reactions * sizeof(reaction));
        memcpy(r, reactions, sim.n_reactions *
                             sizeof(reaction)); /* TODO: when we stop mutilating the reactions we can stop doing this */
        sim_workspace *ws = sim_workspace_init(&sim, &sample, jibal->gsto);
        if (sim.fast) {
            ws->stopping_type = GSTO_STO_ELE;
            ws->rk4 = 0;
        }
        rbs(ws, &sim, reactions, &sample);
#if 0
        if(n == NUMBER_OF_SIMULATIONS-1)
#else
        if(n == 0)
#endif
            print_spectra(f, &global, &sim, ws, &sample,     reactions);
        sim_workspace_free(ws);
        free(r);
    }
    end = clock();
    double cputime_spectrum_ms = (((double) (end - start)) / CLOCKS_PER_SEC)*1000.0/NUMBER_OF_SIMULATIONS;
    fprintf(stderr, "...finished!\n\n");
    fprintf(stderr, "CPU time per spectrum: %.3lf ms.%s\n", cputime_spectrum_ms, cputime_spectrum_ms<10.0?" That was fast!":"");
    if(f != stdout) {
        fclose(f);
    }
    layers_free(layers, n_layers);
    sample_free(&sample);
    jibal_free(jibal);
    return EXIT_SUCCESS;
}
