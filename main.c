#include <stdlib.h>
#include <stdio.h>
#include <jibal.h>
#include <jibal_gsto.h>
#include <jibal_stop.h>
#include <jibal_stragg.h>
#include <jibal_kin.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_histogram.h>

#define ALPHA (0.0*C_DEG)
#define BETA (10.0*C_DEG)
#define THETA (170.0*C_DEG)
#define DETECTOR_RESOLUTION (15.0*C_KEV/C_FWHM)

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
    for(i=0; i < h->n; i++) {
        double w = h->range[i+1] - h->range[i];
        double E = (h->range[i] + h->range[i+1])/2.0; /* Approximate gaussian at center */
        double y = (gsl_sf_erf_Q((E_low-E)/sigma)-gsl_sf_erf_Q((E_high-E)/sigma))/(E_high-E_low);
        h->bin[i] += y*w*Q;
    }
}

double stop_step(jibal_gsto *workspace, const jibal_isotope *incident, const jibal_material *material, double h, double *E, double *S) {
    double k1, k2, k3, k4, stop, dE;
    k1 = -1.0*jibal_stop(workspace, incident, material, *E); /* */
    k2 = -1.0*jibal_stop(workspace, incident, material, *E + (h / 2) * k1);
    k3 = -1.0*jibal_stop(workspace, incident, material, *E + (h / 2) * k2);
    k4 = -1.0*jibal_stop(workspace, incident, material, *E + h * k3);
    stop = (k1 + 2 * k2 + 2 * k3 + k4)/6;
    dE = h*stop; /* Energy change in thickness "h" */
    //fprintf(stderr, "stop = %g eV/tfu, E = %g keV,  dE = %g keV\n", stop/C_EV_TFU, *E/C_KEV, dE/C_KEV);
    double s_ratio = jibal_stop(workspace, incident, material, *E+dE)/(k1); /* Ratio of stopping */
    *S *= (s_ratio)*(s_ratio);
    *S += fabs(h)*jibal_stragg(workspace, incident, material, (*E+dE/2)); /* Straggling, calculate at mid-energy */
    *E += dE;
    return 0.0; /* TODO: return something useful */
}


double overlayer(jibal_gsto *workspace, const jibal_isotope *incident, const jibal_layer *layer, double p_sr, double E_0, double *S, gsl_histogram *histo) {
    double E = E_0;
    double dE;
    double x;
    double h = workspace->stop_step;
#ifdef DEBUG
    fprintf(stderr, "Thickness %g tfu, stop step %g tfu, E_0 = %g MeV\n", layer->thickness/C_TFU, h/C_TFU, E_0/C_MEV);
#endif
    double K = jibal_kin_rbs(incident->mass, layer->material->elements[0].isotopes[0]->mass, THETA, '+');
    double E_out_prev = K*E;
    for (x = layer->thickness; x >= 0.0;) {
        if(x-h < 0.0) { /* Last step may be partial */
            h=x;
            if(h < workspace->stop_step/1e6) {
                break;
            }
        }


        /* DEPTH BIN [x, x-h) */
        fprintf(stderr, "x = %g, x-h = %g\n", x/C_TFU, (x-h)/C_TFU);
        double x_front = x;
        double x_back = x-h;
#if 0
        double k1, k2, k3, k4, stop;
        k1 = factor*jibal_stop(workspace, incident, layer->material, E);
        k2 = factor*jibal_stop(workspace, incident, layer->material, E + (h / 2) * k1);
        k3 = factor*jibal_stop(workspace, incident, layer->material, E + (h / 2) * k2);
        k4 = factor*jibal_stop(workspace, incident, layer->material, E + h * k3);
        stop = (k1 + 2 * k2 + 2 * k3 + k4)/6;
        dE = h*stop; /* Energy change in thickness "h" */
        double s_ratio = jibal_stop(workspace, incident, layer->material, E+dE)/(k1/factor); /* Ratio of stopping */
        *S *= (s_ratio)*(s_ratio); /* Non-statistical broadening due to energy dependent stopping */
        *S += h*jibal_stragg(workspace, incident, layer->material, (E+dE/2)); /* Straggling, calculate at mid-energy */
        E += dE;
#endif
        double S_front = *S;
        double E_front = E;
        stop_step(workspace, incident, layer->material, h, &E, S);
        double S_back = *S;
        double E_back = E;
        double E_mean = (E_front+E_back)/2.0;
        double sigma = jibal_cross_section_rbs(incident, layer->material->elements[0].isotopes[0], THETA, E_mean, JIBAL_CS_ANDERSEN);
        double Q = p_sr * sigma * h; /* TODO: worst possible approximation... */
#ifdef DEBUG
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sigma = %g mb/sr, sqrt(S) = %g keV, Q = %g\n", E_front/C_MEV, E_back/C_MEV, E_mean/C_MEV, sigma/C_MB_SR, sqrt(*S)/C_KEV, Q);
#endif
        double S_out = S_back * K;
        double E_out = E_back * K;
        double x_out;
        for (x_out = x-h; x_out <= layer->thickness;) { /* Calculate energy and straggling of backside of slab */
            //fprintf(stderr, "Surfacing... x_out = %g tfu... ", x_out/C_TFU);
            stop_step(workspace, incident, layer->material, h, &E_out, &S_out);
            x_out += h;
        }
        double E_diff = E_out_prev-E_out;
        if(h > 1.0*C_TFU && E_diff > 5.0*C_KEV || E_diff < 2.0*C_KEV) {
            fprintf(stderr, "E_diff too large or too small: %g keV with step %g tfu!\n", E_diff/C_KEV, h/C_TFU);
            h *= 3.5*C_KEV/E_diff;
            E = E_front;
            fprintf(stderr, "New step %g tfu, back to E = %g MeV!\n\n", h/C_TFU, E/C_MEV);
            continue;
        }
        fprintf(stderr, "Surf: from x_out = %g tfu, gives E_out = %g MeV (prev was %g MeV, diff %g keV), stragg = %g keV\n", (x-h)/C_TFU, E_out/C_MEV, E_out_prev/C_MEV, (E_diff)/C_KEV, sqrt(S_out)/C_KEV);
        fprintf(stderr, "\n");
        //fprintf(stdout, "%.5lf %.5lf\n", E_out/C_KEV, Q);
        brick_int(sqrt(*S+DETECTOR_RESOLUTION*DETECTOR_RESOLUTION), E_out, E_out_prev, histo, Q);
        if(!isnormal(E)) {
            fprintf(stderr, "SOMETHING DOESN'T LOOK RIGHT HERE.\n");
            return 0.0;
        }
        x -= h;
        *S = S_back;
        E = E_back;
        E_out_prev = E_out;
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
        fprintf(stderr, "Not enough arguments! Usage: %s: isotope material thickness energy\n", argv[0]);
        return 1;
    }
    const jibal_isotope *incident = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    jibal_layer *layer = jibal_layer_new(jibal_material_create(jibal->elements, argv[2]), jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, argv[3]));
    if(incident) {
        fprintf(stderr, "Mass %.3lf u\n", incident->mass/C_U);
    }
    if(!layer) {
        fprintf(stderr, "NO LAYER!\n");
        return (EXIT_FAILURE);
    }
    if(!jibal_gsto_auto_assign_material(jibal->gsto, incident, layer->material)) {
        fprintf(stderr, "Couldn't assign stopping.\n");
        return (EXIT_FAILURE);
    }
    jibal_gsto_print_assignments(jibal->gsto);
    jibal_gsto_load_all(jibal->gsto);
    gsl_histogram *histo = gsl_histogram_calloc_uniform(2000, 100*C_KEV, 2100*C_KEV);
    double E = jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, argv[4]);
    double S = 0.0;
    double p_sr = 1.0e12; /* TODO: particles * sr / cos(alpha) */
    overlayer(jibal->gsto, incident, layer, p_sr, E, &S, histo);
    jibal_free(jibal);
    gsl_histogram_fprintf(stdout, histo, "%g", "%g");
    gsl_histogram_free(histo);
    return 0;
}
