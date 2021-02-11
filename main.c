#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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

#define E_MIN (100.0*C_KEV)
#define N_LAYERS_MAX 100
#define HISTOGRAM_CHANNELS 4000

typedef struct {
    size_t n;
    double *ranges;
    int i; /* accelerator */
} conc_range;

typedef struct {
    const jibal_isotope *isotope;
    double E;
    double K; /* TODO: this should be part of a "reaction" and not the "isotope" */
    double c; /* accelerator of concentration */
    double step;
    gsl_histogram *h;
    double *conc; /* concentration values corresponding to ranges given in an associated (shared!) conc_range */
    conc_range *r;
} sim_isotope;


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

conc_range *conc_rance_alloc(size_t n) {
    conc_range *r = malloc(sizeof(conc_range));
    r->n = n;
    r->ranges = calloc(n, sizeof(double ));
    r->i = 0;
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

/* TODO:
 *  - make conc profiles for each isotope to be simulated.
 *    Ranges from layer model? Linear varying (like MCERD) or histogram-like (layer-model)?
 *
 * - Allocate energy histograms.
 *   These should have ranges independently from other energy histograms, since the same bin can correspond to a different energy.
 *   Energy histograms can be summed this way, too.
 * - Provide stopping calls, like jibal_stop(), but using the new material description.
 * - List of reactions? By default lets start with RBS from every target isotope (if RBS possible)
 * - Track "bricking" process in some helper varible (like E_out_prev)
 * */


void foo() {
    double x;
    int i;
    double sigma = 1.0;
    double x_0 = 100.0;
    double x_L = 50.0;
    double x_H = 55.0;
    double sum = 0.0;
    double step = 0.2;
    for(i=0; i <10000; i++) {
        x = -1000.0 + i*step;
        double fx = gsl_sf_erf_Q((x_H-x)/sigma);
        double fx_low = gsl_sf_erf_Q((x_L-x)/sigma);
        double out = (fx_low-fx)/(x_H-x_L); /* One count, convolution of gaussian rectangle */
        sum += out*step; /* This is the integral! */
        fprintf(stdout, "%g %g %g %g\n", x, fx, out, sum);
    }
}

void brick_int(double sigma, double E_low, double E_high, gsl_histogram *h, double Q) { /* Energy spectrum h, integrate over rectangular brick convoluted with a gaussian */
    int i;
    size_t lo;
    size_t hi;
#ifndef NO_OPTIMIZE_BRICK
    gsl_histogram_find(h, E_low, &lo);
    gsl_histogram_find(h, E_low, &hi);
    if(lo > 20)
        lo -= 20;
    else
        lo = 0;
    if(hi > h->n-20)
        hi = h->n;
    else
        hi += 20;
#else
    lo = 0;
    hi = h->n;
#endif
    for(i = lo; i < hi; i++) {
        double w = h->range[i+1] - h->range[i];
        double E = (h->range[i] + h->range[i+1])/2.0; /* Approximate gaussian at center */
        double y = (gsl_sf_erf_Q((E_low-E)/sigma)-gsl_sf_erf_Q((E_high-E)/sigma))/(E_high-E_low);
        h->bin[i] += y*w*Q;
    }
}

double stop_target(jibal_gsto *workspace, const jibal_isotope *incident, sim_isotope *target, double E) { /* Call recalculate_concs() before this */
    double em=E/incident->mass;
    int i, i_isotope;

    double S1 = 0.0;
    for(i_isotope = 0; target[i_isotope].isotope != NULL; i_isotope++) {
        sim_isotope *it = &target[i_isotope];
        S1 += it->c * (
                jibal_gsto_get_em(workspace, GSTO_STO_ELE, incident->Z, it->isotope->Z, em)
                +jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, it->isotope->Z, it->isotope->mass)
                );
    }
    return S1;
}

double stragg_target(jibal_gsto *workspace, const jibal_isotope *incident, sim_isotope *target, double E) { /* Call recalculate_concs() before this */
    double em=E/incident->mass;
    int i, i_isotope;

    double S2 = 0.0;
    for(i_isotope = 0; target[i_isotope].isotope != NULL; i_isotope++) {
        sim_isotope *it = &target[i_isotope];
        S2 +=  it->c*jibal_gsto_get_em(workspace, GSTO_STO_STRAGG, incident->Z, it->isotope->Z, em);
    }
    return S2;
}

void recalculate_concs(sim_isotope *its, double x) {
    int i_isotope;
    for(i_isotope = 0; its[i_isotope].isotope != NULL; i_isotope++) {
        its[i_isotope].c = get_conc(x, &its[i_isotope]);
    }
}

double stop_step(jibal_gsto *workspace, const jibal_isotope *incident, sim_isotope *target, double h, double *E, double *S) {
    double k1, k2, k3, k4, stop, dE;
    k1 = -1.0*stop_target(workspace, incident, target, *E); /* */
    k2 = -1.0*stop_target(workspace, incident, target, *E + (h / 2) * k1);
    k3 = -1.0*stop_target(workspace, incident, target, *E + (h / 2) * k2);
    k4 = -1.0*stop_target(workspace, incident, target, *E + h * k3);
    stop = (k1 + 2 * k2 + 2 * k3 + k4)/6;
    //stop = k1;
    dE = h*stop; /* Energy change in thickness "h" */
    //fprintf(stderr, "stop = %g eV/tfu, E = %g keV,  dE = %g keV\n", stop/C_EV_TFU, *E/C_KEV, dE/C_KEV);
    double s_ratio = stop_target(workspace, incident, target, *E+dE)/(k1); /* Ratio of stopping */
    *S *= (s_ratio)*(s_ratio);
    *S += fabs(h)*stragg_target(workspace, incident, target, (*E+dE/2)); /* Straggling, calculate at mid-energy */
    *E += dE;
    return dE; /* TODO: return something useful */
}

void scatter(const jibal_isotope *incident, const jibal_isotope *target, double E) {
    double sigma = jibal_cross_section_rbs(incident, target, THETA, E, JIBAL_CS_ANDERSEN);
}




double overlayer(jibal_gsto *workspace, const jibal_isotope *incident, sim_isotope *target, double p_sr, double E_0, double *S, conc_range *crange) {
    double E = E_0;
    double E_prev = E;
    double dE;
    double x;
    double h = workspace->stop_step;
    int n_isotopes;
    for(n_isotopes = 0; target[n_isotopes].isotope != NULL; n_isotopes++) {
        sim_isotope *it = &target[n_isotopes];
        it->K = jibal_kin_rbs(incident->mass, it->isotope->mass, THETA, '+'); /* TODO: this is too hard coded for RBS right now */
        it->E = E * it->K;
        it->step = h;
    };
    double thickness = crange->ranges[crange->n-1];
#ifdef DEBUG
    fprintf(stderr, "Thickness %g tfu, stop step %g tfu, E_0 = %g MeV\n", thickness/C_TFU, h/C_TFU, E_0/C_MEV);
#endif
    for (x = 0.0; x <= thickness;) {
        if(x+h > thickness || E < E_MIN) { /* Last step may be partial */ /* TODO: sharp transitions from one range bin to other?! */
            break; /* TODO: BRUTAL */
        }
        recalculate_concs(target, x);
        /* DEPTH BIN [x, x-h) */
#if 1
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3g keV\n", x/C_TFU, (x+h)/C_TFU, E/C_KEV);
#endif
        double E_front = E;
        stop_step(workspace, incident, target, h, &E, S);
        double S_back = *S;
        double E_back = E;
        double E_mean = (E_front+E_back)/2.0;
#ifndef NO_ADAPTIVE_STEPPING_FOR_INCIDENT
        double E_diff = E_front-E_back;
        if(h > 1.0*C_TFU && E_diff > 2.0*C_KEV || E_diff < 0.5*C_KEV) {
            fprintf(stderr, "E_diff too large or too small: %g keV with step %g tfu!\n", E_diff/C_KEV, h/C_TFU);
            h *= 1.0*C_KEV/E_diff;
            E = E_front;
            fprintf(stderr, "New step %g tfu, back to E = %g MeV!\n\n", h/C_TFU, E/C_MEV);
            continue;
        }
#endif
#ifdef DEBUG
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(*S) / C_KEV);
#endif
        int i_isotope;
        for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it = &target[i_isotope];
            if(it->c < 1e-12) /* TODO: concentration cutoff? TODO: it->E should be updated when we start calculating it again?*/
                continue;
            const jibal_isotope *isotope = it->isotope;
            double sigma = jibal_cross_section_rbs(incident, isotope, THETA, E_mean,JIBAL_CS_ANDERSEN);
            double Q = it->c * p_sr * sigma * h; /* TODO: worst possible approximation... */
            double S_out = S_back * it->K;
            double E_out = E_back * it->K;
            double x_out;
            for (x_out = x + h; x_out >= 0.0;) { /* Calculate energy and straggling of backside of slab */
                    //fprintf(stderr, "Surfacing... x_out = %g tfu... ", x_out/C_TFU);
                    if(x_out < it->step) {
                        stop_step(workspace, incident, target, x_out, &E_out, &S_out);
                        break;
                    } else {
                        stop_step(workspace, incident, target, it->step, &E_out, &S_out);
                        x_out -= it->step;
                    }
            }
#ifdef DEBUG
            fprintf(stderr, "    %s: E_scatt = %.3lf, E_out = %.3lf, sigma = %g mb/sr\n", isotope->name, E_back * it->K/C_KEV, E_out/C_KEV, sigma/C_MB_SR);
#endif
            brick_int(sqrt(*S+DETECTOR_RESOLUTION*DETECTOR_RESOLUTION), E_out, it->E, it->h, Q);
            it->E = E_out;
        }
#if 0

        fprintf(stderr, "Surf: from x_out = %g tfu, gives E_out = %g MeV (prev was %g MeV, diff %g keV), stragg = %g keV\n", (x-h)/C_TFU, E_out/C_MEV, E_out_prev/C_MEV, (E_diff)/C_KEV, sqrt(S_out)/C_KEV);
        fprintf(stderr, "\n");
#endif
        if(!isnormal(E)) {
            fprintf(stderr, "SOMETHING DOESN'T LOOK RIGHT HERE.\n");
            return 0.0;
        }
        x += h;
        *S = S_back;
        E = E_back;
    }
    return E;
}

int main(int argc, char **argv) {
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return 1;
    }
    if(argc < 5) {
        fprintf(stderr, "Not enough arguments! Usage: %s: isotope energy material thickness material2 thickness2\n", argv[0]);
        return 1;
    }
    const jibal_isotope *incident = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    if(incident) {
        fprintf(stderr, "Mass %.3lf u\n", incident->mass/C_U);
    } else {
        fprintf(stderr, "No isotope %s found.\n", argv[1]);
    }
    double E = jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, argv[2]);
    fprintf(stderr, "Beam E = %g MeV\n", E/C_MEV);
    if (E > 1000.0*C_MEV || E < 10*C_KEV) {
        fprintf(stderr, "Hmm...? Check your numbers.\n");
        return -1;
    }
    argc -= 3;
    argv += 3;
    jibal_layer **layers = calloc(N_LAYERS_MAX, sizeof(jibal_layer *));
    int n_layers=0;
    while (argc >= 2 && n_layers < N_LAYERS_MAX) {
        jibal_layer *layer = jibal_layer_new(jibal_material_create(jibal->elements, argv[0]),
                                             jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, argv[1]));
        if (!layer) {
            fprintf(stderr, "NO LAYER %s!\n", argv[0]);
            return (EXIT_FAILURE);
        }
        if (!jibal_gsto_auto_assign_material(jibal->gsto, incident, layer->material)) {
            fprintf(stderr, "Couldn't assign stopping.\n");
            return (EXIT_FAILURE);
        }
        layers[n_layers] = layer;
        argc -= 2;
        argv += 2;
        n_layers++;
    }
    int i,j,k;
    conc_range *crange = conc_rance_alloc(2*n_layers);
    int i_isotope, n_isotopes=0;
    for(i = 0; i < n_layers; i++) {
        jibal_layer *layer = layers[i];
        fprintf(stderr, "Layer %i/%i. Thickness %g tfu\n", i+1, n_layers, layer->thickness/C_TFU);
        jibal_material_print(stderr, layer->material);
        crange->ranges[2*i] = i?crange->ranges[2*i-1]:0.0;
        crange->ranges[2*i+1] = crange->ranges[2*i] + layer->thickness;
        for (j = 0; j < layer->material->n_elements; ++j) {
            n_isotopes += layer->material->elements[j].n_isotopes;
        }
    }
    for (i = 0; i < crange->n; i++) {
        fprintf(stderr, "ranges[%i]  = %g\n", i, crange->ranges[i]);
    }
    fprintf(stderr, "Total %i isotopes (each isotope in each layer treated differently...)\n", n_isotopes);

    sim_isotope *its = calloc(n_isotopes+1, sizeof(sim_isotope));
    its[n_isotopes].isotope = NULL; /* Last isotope of isotopes is NULL. For-loops without knowledge of n_isotopes are possible. */
    i_isotope=0;
    for(i = 0; i < n_layers; i++) {
        jibal_layer *layer = layers[i];
        for (j = 0; j < layer->material->n_elements; ++j) {
            jibal_element *element = &layer->material->elements[j];
            for(k = 0; k < element->n_isotopes; k++) {
                assert(i_isotope < n_isotopes);
                sim_isotope *it = &its[i_isotope];
                it->isotope = element->isotopes[k];
                it->conc = calloc(crange->n, sizeof (double));
                it->conc[2*i] = element->concs[k] * layer->material->concs[j];
                it->conc[2*i+1] = element->concs[k] * layer->material->concs[j];
                it->h = gsl_histogram_calloc_uniform(HISTOGRAM_CHANNELS, 0*C_KEV, 4000*C_KEV);
                it->r = crange;
                it->E = 0.0;
                //fprintf(stderr, "i: %i, j: %i, k: %i, i_isotope: %i, name: %s\n", i, j, k, i_isotope, it->isotope->name);
                i_isotope++;
            }
        }
    }
    fprintf(stderr, "\nDEPTH(tfu) ");
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
        fprintf(stdout, "%10.3lf", d/C_TFU);
        for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it =  &its[i_isotope];
            double x = get_conc(d, it);
            fprintf(stdout, " %8.4lf", x*100.0);
        }
        fprintf(stdout, "\n");
    }
#endif

    jibal_gsto_print_assignments(jibal->gsto);
    jibal_gsto_load_all(jibal->gsto);

    //gsl_histogram *histo = gsl_histogram_calloc_uniform(2000, 100*C_KEV, 2100*C_KEV);

    double S = 0.0;
    double p_sr = 1.0e12; /* TODO: particles * sr / cos(alpha) */
    overlayer(jibal->gsto, incident, its, p_sr, E, &S, crange);

    for(i = 0; i < HISTOGRAM_CHANNELS; i++) {
        fprintf(stdout, "%4i", i);
        double sum = 0.0;
        for (i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
            sim_isotope *it = &its[i_isotope];
            fprintf(stdout, " %10.3lf", it->h->bin[i]);
            sum += it->h->bin[i];
        }
        fprintf(stdout, " %10.3lf\n", sum);
    }

    jibal_free(jibal);
    //gsl_histogram_fprintf(stdout, histo, "%g", "%g");
    //gsl_histogram_free(histo);
    return 0;
}
