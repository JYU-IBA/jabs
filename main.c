#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <jibal.h>
#include <jibal_gsto.h>
#include <jibal_stop.h>
#include <jibal_stragg.h>
#include <jibal_kin.h>
#include <jibal_cs.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_histogram.h>

#define ALPHA (0.0*C_DEG)
#define BETA (10.0*C_DEG)
#define THETA (170.0*C_DEG)
#define DETECTOR_RESOLUTION (15.0*C_KEV/C_FWHM)
#define PARTICLES_SR (1.0e12)

#define E_MIN (100.0*C_KEV)
#define N_LAYERS_MAX 100
#define HISTOGRAM_BIN (1.0*C_KEV)
#define STOP_STEP_INCIDENT (5.0*C_KEV)
#define STOP_STEP_EXITING (5.0*C_KEV)

#define CONCENTRATION_CUTOFF 1e-8

#define NUMBER_OF_SIMULATIONS 1

typedef struct {
    size_t n;
    double *ranges;
    int i; /* accelerator */
} conc_range;

typedef struct {
    const jibal_isotope *isotope;
    double E;
    double K; /* TODO: this should be part of a "reaction" and not the "isotope" */
    double c; /* accelerator of concentration, value assigned by recalculate_concs() and used by stopping related things etc */
    double step;
    gsl_histogram *h; /* energy histogram */
    double *conc; /* concentration values corresponding to ranges given in an associated (shared!) conc_range */
    double max_depth;
    conc_range *r;
} sim_isotope;

typedef struct {
    const jibal_isotope *isotope;
    double E;
    double mass;
    int Z;
    double angle; /* Traversing matter "straight on" when angle = 0, getting stuck sideways if angle = 90.0*C_DEG */
    double inverse_cosine; /* Inverse cosine of angle. Traversing matter "straight on" means 1.0 and going sidewways approaches infinity. */
} sim_ion;

typedef struct {
    int n_isotopes;
    sim_isotope *its;
    conc_range crange;
} sample;

typedef struct {
    int n_channels;
    double histogram_bin;
    const jibal_isotope *incident;
    double E;
    double p_sr;
    double p_sr_cos_alpha; /* particles * sr / cos(alpha) */
    double alpha;
    double beta;
    double theta;
} simulation;



double get_conc(double x, const sim_isotope *it) {
    int lo, mi, hi;
    conc_range *r = it->r;
    if(x >= r->ranges[r->i] && x < r->ranges[r->i+1]) { /* Use saved bin. We'll probably receive a lot of repeated calls to the same bin, so this avoid cost of other range checking and binary searches. */
        if(it->conc[r->i] == it->conc[r->i+1]) /* Constant */
            return it->conc[r->i];
        else /* Linear interpolation */
            return it->conc[r->i] + ((it->conc[r->i + 1] - it->conc[r->i]) / (r->ranges[r->i + 1] - r->ranges[r->i])) *
                                (x - r->ranges[r->i]);
    }
    if(x < r->ranges[0] || x >= r->ranges[r->n-1]) { /* Out of bounds concentration is zero. Maybe the execution doesn't go here if all goes as planned, so this could be changed to an assert. */
        return 0.0;
    }
    hi = r->n;
    lo = 0;
    while (hi - lo > 1) {
        mi = (hi + lo) / 2;
        if (x >= r->ranges[mi]) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    r->i = lo; /* Save value for next time */
    if(r->ranges[lo+1]-r->ranges[lo] == 0) /* Zero width. Return value of left side. */
        return it->conc[lo];
    return it->conc[lo] + ((it->conc[lo+1] - it->conc[lo])/(r->ranges[lo+1] - r->ranges[lo])) * (x - r->ranges[lo]);
}

conc_range conc_rance_alloc(size_t n) {
    //conc_range *r = malloc(sizeof(conc_range));
    conc_range r = {.n = n, .i = 0};
    r.ranges = calloc(n, sizeof(double ));
    return r;
}

sim_isotope *sim_isotopes_alloc(size_t n_isotopes, size_t n_channels, conc_range *r) {
    sim_isotope *isotopes = calloc(n_isotopes+1, sizeof(sim_isotope));
    int i;
    for (i = 0; i < n_isotopes; i++) {
        sim_isotope *isotope = &isotopes[i];
        isotope->isotope = NULL;
        isotope->E = 0.0;
        isotope->h = gsl_histogram_alloc(n_channels);
        isotope->r = r;
        isotope->conc = calloc(r->n, sizeof (double));
    }
    isotopes[n_isotopes].isotope = NULL;
    return isotopes;
}

inline double erfc(double x) {
    return exp(-1.0950081470333*x-0.75651138383854*x*x);
    /* Tsay, WJ., Huang, C.J., Fu, TT. et al. J Prod Anal 39, 259–269 (2013). https://doi.org/10.1007/s11123-012-0283-1 */
}
inline double erf_Q(double x);

double erf_Q(double x) {
    return x < 0.0 ? 1.0-0.5*erfc(-1.0*x/sqrt(2.0)) : 0.5*erfc(x/sqrt(2.0));
}

void erf_Q_test() {
    double x;
    int i;
    double sigma = 1.0;
    double x_0 = 100.0;
    double x_L = 10.0;
    double x_H = 15.0;
    double sum = 0.0;
    double step = 0.03;
    for(i=0; i <10000; i++) {
        x = -50.0 + i*step;
        double fx = gsl_sf_erf_Q((x_H-x)/sigma);
        double fx_low = gsl_sf_erf_Q((x_L-x)/sigma);
        fx = erf_Q((x_H-x)/sigma);
        fx_low = erf_Q((x_L-x)/sigma);
        double out = (fx_low-fx)/(x_H-x_L); /* One count, convolution of gaussian rectangle */
        sum += out*step; /* This is the integral! */
        fprintf(stdout, "%g %g %g %12.8lf\n", x, fx, out, sum);
    }
}

void brick_int(double sigma, double E_low, double E_high, gsl_histogram *h, double Q) { /* Energy spectrum h, integrate over rectangular brick convoluted with a gaussian */
    int i;
#ifndef NO_OPTIMIZE_BRICK
    size_t lo;
    size_t hi;
    gsl_histogram_find(h, E_low, &lo);
    gsl_histogram_find(h, E_high, &hi);
    double foo = (E_high-E_low)/(HISTOGRAM_BIN);
    int n = ceil(foo*5);
    if(n < 50)
        n = 50;
    if(lo > n)
        lo -= n;
    else
        lo = 0;
    if(hi > h->n-n)
        hi = h->n;
    else
        hi += n;
    for(i = lo; i < hi; i++) {
#else
    for(i = 0; i < h->n; i++) {
#endif
        double w = h->range[i+1] - h->range[i];
        double E = (h->range[i] + h->range[i+1])/2.0; /* Approximate gaussian at center */
#ifdef ERF_Q_FROM_GSL
        double y = (gsl_sf_erf_Q((E_low-E)/sigma)-gsl_sf_erf_Q((E_high-E)/sigma))/(E_high-E_low);
#else
        double y = (erf_Q((E_low-E)/sigma)-erf_Q((E_high-E)/sigma))/(E_high-E_low);
#endif
        h->bin[i] += y*w*Q;
    }
}

double stop_target(jibal_gsto *workspace, const jibal_isotope *incident, sim_isotope *target, double E) { /* Call recalculate_concs() before this */
    double em=E/incident->mass;
    int i, i_isotope;

    double S1 = 0.0;
    for(i_isotope = 0; target[i_isotope].isotope != NULL; i_isotope++) {
        sim_isotope *it = &target[i_isotope];
        if(it->c < CONCENTRATION_CUTOFF)
             continue;
        S1 += it->c * (
                jibal_gsto_get_em(workspace, GSTO_STO_ELE, incident->Z, it->isotope->Z, em)
                +jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, it->isotope->Z, it->isotope->mass)
                );
    }
    assert(S1 > 0.0);
    return S1;
}

double stragg_target(jibal_gsto *workspace, const jibal_isotope *incident, sim_isotope *target, double E) { /* Call recalculate_concs() before this */
    double em=E/incident->mass;
    int i_isotope;

    double S2 = 0.0;
    for(i_isotope = 0; target[i_isotope].isotope != NULL; i_isotope++) {
        sim_isotope *it = &target[i_isotope];
        if(it->c <= CONCENTRATION_CUTOFF)
            continue;
        S2 +=  it->c*jibal_gsto_get_em(workspace, GSTO_STO_STRAGG, incident->Z, it->isotope->Z, em);
    }
    return S2;
}

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

double stop_step(jibal_gsto *workspace, const jibal_isotope *incident, sim_isotope *target, double *h, double h_max, double *E, double *S, double step) {
    double k1, k2, k3, k4, stop, dE;
    k1 = -1.0*stop_target(workspace, incident, target, *E);
    if(k1 > -0.01*C_EV_TFU) { /* Fail on positive values, zeroes (e.g. due to zero concentrations) and too small negative values */
        *h = 0.0;
        return 0.0;
    }
    *h = (-1.0*step / k1); /* we try to aim for some energy loss */
    if(*h > h_max)
        *h = h_max; /* but we have some other limitations too */
#ifdef NO_ULTRASMALL_STEPS
    if(*h < 0.1*C_TFU) {
        *h = 0.1*C_TFU;
    }
#endif
    k2 = -1.0*stop_target(workspace, incident, target, *E + (*h / 2) * k1);
    k3 = -1.0*stop_target(workspace, incident, target, *E + (*h / 2) * k2);
    k4 = -1.0*stop_target(workspace, incident, target, *E + (*h) * k3);
    stop = (k1 + 2 * k2 + 2 * k3 + k4)/6;
    //stop = k1;
    dE = (*h)*stop; /* Energy change in thickness "h" */
    //fprintf(stderr, "stop = %g eV/tfu, E = %g keV, h = %6.3lf  dE = %g keV\n", stop/C_EV_TFU, *E/C_KEV, *h/C_TFU, dE/C_KEV);
#ifdef DEBUG
    if(fabs(stop) < 0.1*C_EV_TFU) {
        fprintf(stderr, "Not good!\n");
        return 0.0;
    }
#endif
    double s_ratio = stop_target(workspace, incident, target, *E+dE)/(k1); /* Ratio of stopping */
#ifdef DEBUG
    if((s_ratio)*(s_ratio) < 0.9 || (s_ratio)*(s_ratio) > 1.1) {
        fprintf(stderr, "YIKES, s_ratio = %g, sq= %g\n", s_ratio, (s_ratio)*(s_ratio));
    }
#endif
    *S *= (s_ratio)*(s_ratio);
    *S += fabs(*h)*stragg_target(workspace, incident, target, (*E+dE/2)); /* Straggling, calculate at mid-energy */
    *E += dE;
    assert(isnormal(*S));
    return dE; /* TODO: return something useful */
}

void rbs(jibal_gsto *workspace, simulation *sim, sample *sample, sim_ion *ion) {
    double E = sim->E;
    double S=0.0;
    double x;
    double h = workspace->stop_step;
    conc_range *crange = &sample->crange;
    sim_isotope *its = sample->its;
    int i_isotope;
    int n_isotopes = sample->n_isotopes;
    for(i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
        sim_isotope *it = &its[i_isotope];
        it->K = jibal_kin_rbs(ion->mass, it->isotope->mass, sim->theta, '+'); /* TODO: this is too hard coded for RBS right now */
        it->E = E * it->K;
        it->step = h;
    };
    double thickness = crange->ranges[crange->n-1];
    double next_crossing = crange->ranges[1];
    double h_max;
    int i_range = 0;
#ifdef DEBUG
    fprintf(stderr, "Thickness %g tfu, stop step %g tfu, E = %g MeV\n", thickness/C_TFU, h/C_TFU, E/C_MEV);
#endif
    for (x = 0.0; x <= thickness;) {
#if 1
        while (x >= crange->ranges[i_range+1]) {
            i_range++;
            if (i_range >= crange->n-1) {
#ifdef DEBUG
                fprintf(stderr, "return due to last range, last (uncalculated) x = %g tfu\n", x/C_TFU);
#endif
                return;
            }
#ifdef DEBUG
            fprintf(stderr, "Crossing to range %i = [%g, %g)\n", i_range, crange->ranges[i_range]/C_TFU, crange->ranges[i_range+1]/C_TFU);
#endif
            next_crossing = crange->ranges[i_range+1];
        }
#endif
        if(E < E_MIN) {
            fprintf(stderr, "Return due to low energy.\n");
            return;
        }
        h_max = next_crossing - x;
        if(h_max < 0.0001*C_TFU) {
            x += 0.0001*C_TFU;
            fprintf(stderr, "Step too small. Let's nudge forward! x=%lf tfu\n", x/C_TFU);
            continue;
        }

        double E_front = E;
        recalculate_concs(its, x);
        stop_step(workspace, ion->isotope, its, &h, h_max, &E, &S, STOP_STEP_INCIDENT);
        assert(h > 0.0);
        /* DEPTH BIN [x, x+h) */
        double E_back = E;
        double E_diff = E_front-E_back;
#if 0
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3lf keV to  %8.3lf keV (diff %6.4lf keV)\n", x/C_TFU, (x+h)/C_TFU, E_front/C_KEV, E/C_KEV, E_diff/C_KEV);
#endif
        double S_back = S;
        double E_mean = (E_front+E_back)/2.0;
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(S) / C_KEV);
#endif
        for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it = &its[i_isotope];
            if( x > it->max_depth)
                continue;
            recalculate_concs(its, x);
            double c = it->c;
            double S_out = S_back * it->K;
            double E_out = E_back * it->K;
            double x_out;
            int i_range_out = i_range;
            for (x_out = x + h; x_out > 0.0;) { /* Calculate energy and straggling of backside of slab */
                //fprintf(stderr, "Surfacing... x_out = %g tfu... ", x_out/C_TFU);
                recalculate_concs(its, x_out-0.00001*C_TFU); /* TODO: x_out may be close to an area where concentrations are zero, therefore no stopping... */
                double remaining = x_out - crange->ranges[i_range_out];
                if(remaining < 0.1*C_TFU) {
                    x_out = crange->ranges[i_range_out];
                    i_range_out--;
                    continue;
                } else {
                    stop_step(workspace, ion->isotope, its, &(it->step), remaining, &E_out, &S_out, STOP_STEP_EXITING);
                }
                x_out -= it->step;
                if(it->step == 0.0)
                    break;
                if(E_out < E_MIN)
                    break;
            }
            if(c > CONCENTRATION_CUTOFF && E_out > E_MIN) {/* TODO: concentration cutoff? TODO: it->E should be updated when we start calculating it again?*/
                const jibal_isotope *isotope = it->isotope;
                double sigma = jibal_cross_section_rbs(ion->isotope, isotope, THETA, E_mean, JIBAL_CS_ANDERSEN);
                double Q = c * sim->p_sr_cos_alpha * sigma * h; /* TODO: worst possible approximation... */
#ifdef DEBUG
                fprintf(stderr, "    %s: E_scatt = %.3lf, E_out = %.3lf, sigma = %g mb/sr, Q = %g\n", isotope->name, E_back * it->K/C_KEV, E_out/C_KEV, sigma/C_MB_SR, Q);
#endif
                assert(sigma > 0.0);
                assert(S_out > 0.0);
                assert(Q < 1.0e7 || Q > 0.0);
                brick_int(sqrt(S_out + DETECTOR_RESOLUTION * DETECTOR_RESOLUTION), E_out, it->E, it->h, Q);
            }
            it->E = E_out;
        }
#if 0
        fprintf(stderr, "Surf: from x_out = %g tfu, gives E_out = %g MeV (prev was %g MeV, diff %g keV), stragg = %g keV\n", (x-h)/C_TFU, E_out/C_MEV, E_out_prev/C_MEV, (E_diff)/C_KEV, sqrt(S_out)/C_KEV);
        fprintf(stderr, "\n");
#endif
        if(!isnormal(E)) {
            fprintf(stderr, "SOMETHING DOESN'T LOOK RIGHT HERE.\n");
            return;
        }
        x += h;
        S = S_back;
        E = E_back;
    }
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
sample sample_from_layers(jibal *jibal, simulation *sim, jibal_layer **layers, int n_layers) {
        sample s = {.n_isotopes = 0, .its = NULL};
        int i, j, k;
    s.crange = conc_rance_alloc(2*n_layers);
    conc_range *crange = &s.crange;
    int i_isotope, n_isotopes=0;
    for(i = 0; i < n_layers; i++) {
        jibal_layer *layer = layers[i];
        if(!jibal_gsto_auto_assign_material(jibal->gsto, sim->incident, layer->material)) {
            fprintf(stderr, "Couldn't assign stopping.\n");
            return s;
        }
#ifdef DEBUG
        fprintf(stderr, "Layer %i/%i. Thickness %g tfu\n", i+1, n_layers, layer->thickness/C_TFU);
        jibal_material_print(stderr, layer->material);
#endif
        crange->ranges[2*i] = i?crange->ranges[2*i-1]:0.0;
        crange->ranges[2*i+1] = crange->ranges[2*i] + layer->thickness;
        for (j = 0; j < layer->material->n_elements; ++j) {
            n_isotopes += layer->material->elements[j].n_isotopes;
        }
    }
#ifdef DEBUG
    for (i = 0; i < crange->n; i++) {
        fprintf(stderr, "ranges[%i]  = %g\n", i, crange->ranges[i]);
    }
    fprintf(stderr, "Total %i isotopes (each isotope in each layer treated differently...)\n", n_isotopes);
#endif
    sim_isotope *its = calloc(n_isotopes + 1, sizeof(sim_isotope));
    its[n_isotopes].isotope = NULL; /* Last isotope of isotopes is NULL. For-loops without knowledge of n_isotopes are possible. */
    i_isotope = 0;
    for (i = 0; i < n_layers; i++) {
        jibal_layer *layer = layers[i];
        for (j = 0; j < layer->material->n_elements; ++j) {
            jibal_element *element = &layer->material->elements[j];
            for (k = 0; k < element->n_isotopes; k++) {
                assert(i_isotope < n_isotopes);
                sim_isotope *it = &its[i_isotope];
                it->isotope = element->isotopes[k];
                it->conc = calloc(crange->n, sizeof(double));
                it->conc[2 * i] = element->concs[k] * layer->material->concs[j];
                it->conc[2 * i + 1] = element->concs[k] * layer->material->concs[j];
                it->max_depth = crange->ranges[2 * i + 1];
                it->h = gsl_histogram_calloc_uniform(sim->n_channels, 0 * C_KEV, sim->histogram_bin * sim->n_channels);
                it->r = crange;
                it->E = 0.0;
                i_isotope++;
            }
        }
    }
    s.its = its;
    s.n_isotopes = n_isotopes;
    return s;
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
    if(sim.incident) {
#ifdef DEBUG
        fprintf(stderr, "Mass %.3lf u\n", sim.incident->mass/C_U);
#endif
    } else {
        fprintf(stderr, "No isotope %s found.\n", argv[1]);
    }
    sim.E = jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, argv[2]);
#ifdef DEBUG
    fprintf(stderr, "Beam E = %g MeV\n", sim.E/C_MEV);
#endif
    if (sim.E > 1000.0*C_MEV || sim.E < 10*C_KEV) {
        fprintf(stderr, "Hmm...? Check your numbers.\n");
        return -1;
    }
    sim.histogram_bin = HISTOGRAM_BIN;
    sim.n_channels= ceil(1.1 * sim.E / sim.histogram_bin);
    sim.p_sr = PARTICLES_SR; /* TODO: particles * sr / cos(alpha) */
    sim.alpha = ALPHA;
    sim.beta = BETA;
    sim.theta = THETA;
    sim.p_sr_cos_alpha = sim.p_sr / cos(sim.alpha);
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

    sample sample = sample_from_layers(jibal, &sim, layers, n_layers);
    if(sample.n_isotopes == 0)
        return EXIT_FAILURE;
    int i_isotope;
    int n_isotopes = sample.n_isotopes;
    conc_range *crange = &sample.crange;
    sim_isotope *its = sample.its;
    fprintf(stderr, "DEPTH(tfu) ");
    for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
        fprintf(stderr, "%8s ", its[i_isotope].isotope->name);
    }
    fprintf(stderr, "\n");
    for (i = 0; i < crange->n; i++) {
        fprintf(stderr, "%10.3lf", crange->ranges[i]/C_TFU);
        for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it =  &its[i_isotope];
            fprintf(stderr, " %8.4lf", it->conc[i]*100.0);
        }
        fprintf(stderr, "\n");
    }

#ifdef PRINT_SAMPLE_MODEL
    for (i = 0; i < 100; ++i) {
        double d = 10.0*C_TFU*i;
        fprintf(stderr, "%10.3lf", d/C_TFU);
        for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it =  &its[i_isotope];
            double x = get_conc(d, it);
            fprintf(stderr, " %8.4lf", x*100.0);
        }
        fprintf(stderr, "\n");
    }
#endif
    jibal_gsto_print_assignments(jibal->gsto);
    jibal_gsto_load_all(jibal->gsto);

    start = clock();
    for (int n = 0; n < NUMBER_OF_SIMULATIONS; n++) {
        sim_ion ion;
        ion.E = sim.E;
        ion.isotope = sim.incident;
        ion.angle = 0.0;
        ion.inverse_cosine = 1.0; /* TODO: calc these */
        ion.mass = ion.isotope->mass;
        ion.Z = ion.isotope->Z;
        rbs(jibal->gsto, &sim, &sample, &ion);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    for(i = 0; i < sim.n_channels; i++) {
        double sum = 0.0;
        for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it = &its[i_isotope];
            sum += it->h->bin[i];
        }
        fprintf(f, "%4i %e", i, sum);
        for(i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it = &its[i_isotope];
            fprintf(f, " %e", it->h->bin[i]);

        }
        fprintf(f, "\n");
    }
    fprintf(stderr, "CPU time: %g ms per spectrum, for actual simulation of %i spectra.\n", cpu_time_used*1000.0/NUMBER_OF_SIMULATIONS, NUMBER_OF_SIMULATIONS);
    fclose(f);
    for(i = 0; i < n_layers; i++) {
        jibal_layer_free(layers[i]);
    }
    free(layers);
    jibal_free(jibal);
    for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
        sim_isotope *it = &its[i_isotope];
        free(it->conc);
        gsl_histogram_free(it->h);
    }
    free(its);
    free(crange->ranges);
    return 0;
}
