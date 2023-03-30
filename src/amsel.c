#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <jibal.h>
#include <jibal_masses.h>
#include <jibal_stop.h>
#include <jibal_stragg.h>

typedef struct amsel_angular_spread {
    double tau;
    int n;
    double w;
    double v_w;
    double gamma_b;
} amsel_angular_spread_p;

typedef struct amsel_table {
    int n_min; /* Smallest p->n in table (n of first element) */
    int n_max; /* Largest p->n in table (n of last element), inclusive */
    amsel_angular_spread_p *p; /* Table of values, number of elements is n */
    int n;
} amsel_table;

int amsel_table_index_n(const amsel_table *t, int n) { /* Returns index of t->p table, -1 on error */
    int i = n - t->n_min;
    fprintf(stderr, "index = %i\n", i);
    if(i < 0 || i >= t->n) {
        return -1;
    }
    return i;
}

int amsel_table_index_tau(const amsel_table *t, double tau) { /* Returns index in t->p table with p.tau *closest* to tau in log2 */
    double logtau = log2(tau);
    int n = (int) round(logtau); /* Note: nearest index! */
    //fprintf(stderr, " n = %i (log2 of tau is %g)\n", n, logtau);
    return amsel_table_index_n(t, n);
}

int amsel_table_set_element(amsel_table *t, int i, double w, double v_w, double gamma_b) {
    if(i < 0 || i >= t->n) {
        return -1;
    }
    amsel_angular_spread_p *p = &t->p[i];
    p->w = w;
    p->v_w = v_w;
    p->gamma_b = gamma_b;
    return 0;
}

amsel_table *amsel_table_generate() {
    amsel_table *t = calloc(1, sizeof(amsel_table));
    t->n = 23;
    t->n_min = -10;
    t->p = calloc(t->n, sizeof(amsel_angular_spread_p));
    for(int i = 0; i < t->n; i++) {
        amsel_angular_spread_p *p = &t->p[i];
        p->n = t->n_min + i;
        p->tau = pow(2.0, p->n);
    }
    amsel_table_set_element(t, 0, 2.894e-05, 0.626, 2.176);
    amsel_table_set_element(t, 1, 8.797e-05, 0.626, 2.177);
    amsel_table_set_element(t, 2, 2.665e-04, 0.627, 2.184);
    amsel_table_set_element(t, 3, 8.005e-04, 0.634, 2.175);
    amsel_table_set_element(t, 4, 2.362e-03, 0.649, 2.177);
    amsel_table_set_element(t, 5, 6.723e-03, 0.680, 2.172);
    amsel_table_set_element(t, 6, 1.799e-02, 0.734, 2.157);
    amsel_table_set_element(t, 7, 4.414e-02, 0.817, 2.126);
    amsel_table_set_element(t, 8, 9.811e-02, 0.923, 2.081);
    amsel_table_set_element(t, 9, 1.988e-01, 1.042, 2.030);
    amsel_table_set_element(t, 10, 3.730e-01, 1.160, 1.982);
    amsel_table_set_element(t, 11, 6.599e-01, 1.269, 1.940);
    amsel_table_set_element(t, 12, 1.117e+00, 1.363, 1.905);
    amsel_table_set_element(t, 13, 1.830e+00, 1.445, 1.877);
    amsel_table_set_element(t, 14, 2.923e+00, 1.513, 1.855);
    amsel_table_set_element(t, 15, 4.581e+00, 1.570, 1.828);
    amsel_table_set_element(t, 16, 7.705e+00, 1.618, 1.824);
    amsel_table_set_element(t, 17, 1.080e+01, 1.658, 1.812);
    amsel_table_set_element(t, 18, 1.633e+01, 1.692, 1.803);
    amsel_table_set_element(t, 19, 2.451e+01, 1.721, 1.795);
    amsel_table_set_element(t, 20, 3.657e+01, 1.745, 1.789);
    amsel_table_set_element(t, 21, 5.427e+01, 1.766, 1.784);
    amsel_table_set_element(t, 22, 8.020e+01, 1.784, 1.779);
    t->n_max = t->n_min + t->n - 1;
    return t;
}

void amsel_table_free(amsel_table *t) {
    if(!t) {
        return;
    }
    free(t->p);
    free(t);
}

double wa_hwhm(const amsel_table *t, double tau) {
    const double a = 0.5387;
    const double b = 0.3434;
    const double c = 1.3209;
    const double d = 1.0579;
    const double e = 13.185;
    int index = amsel_table_index_tau(t, tau);
    if(index < 0) {
        return 0.0;
    }
    const amsel_angular_spread_p *p = &t->p[index];
    double h1 = - a * tanh(b * (- e * e / fabs(e + log(tau) + c/d) + e)) + d;
    fprintf(stderr, "1/v = h1 = %g\n", h1);
    fprintf(stderr, "1.0/v_w = %g\n", 1.0/p->v_w);
    return p->w * pow(tau / p->tau, 1.0/p->v_w);
    //return p->w * pow(tau / p->tau, h1);
}

double screening_length_TF(int Z1) {
    return 0.8853 * C_BOHR_RADIUS * pow(Z1, -1.0/3.0);
}

double reduced_thickness_Z(double t, int Z) {
    double a = screening_length_TF(Z);
    return C_PI * pow2(a) * t;
}

double reduced_thickness_material(double t, jibal_material *material) {
    double result = 0.0;
    for(size_t i_element = 0; i_element < material->n_elements; i_element++) {
        result += material->concs[i_element] * reduced_thickness_Z(t, material->elements[i_element].Z);
    }
    return result;
}

double gamma_beta_TF(double tau) {
    const double a = 0.5584;
    const double b = 0.3591;
    const double c = 0.5982;
    const double d = 1.058;
    const double e = 6.7397;
    double h1 = - a * tanh(b * (- e * e / fabs(e + log(tau) + c/d) + e)) + d;
    return pow(1 + 1/h1, h1);
}

int amsel_calc(const jibal *jibal, const amsel_table *t, double depth, double phi, double E_lab, const jibal_isotope *incident, const jibal_material *material) {
    fprintf(stderr, "depth = %g tfu\n", depth / C_TFU);
    double thickness = depth / cos(phi);
    fprintf(stderr, "thickness = %g tfu\n", thickness / C_TFU);

    double angle_sum = 0.0;
    double angle_chord_sum = 0.0;
    jibal_layer *layer = jibal_layer_new(jibal_material_copy(material), thickness);
    double S = 0.0;
    double E_deep = jibal_layer_energy_loss_with_straggling(jibal->gsto, incident, layer, E_lab, -1.0, &S);
    fprintf(stderr, "E_lab (surface) = %g keV\n", E_lab / C_KEV);
    fprintf(stderr, "E_deep = %g keV\n", E_deep / C_KEV);
    double E_avg = (E_lab + E_deep) / 2.0;
    for(size_t i_element = 0; i_element < material->n_elements; i_element++) {
        int Z2 = material->elements[i_element].Z;
        double a = screening_length_TF(Z2);
        double tau =  C_PI * pow2(a) * thickness;
        double wa = wa_hwhm(t, tau);
        fprintf(stderr, "wa = %g (HWHM), angular spread\n", wa);
        double mu = E_avg * screening_length_TF(Z2) / (2.0 * incident->Z * C_E * Z2 * C_E);
        mu *= (4.0 * C_PI * C_EPSILON0); /* Damn it again! */
        double angle = wa / mu;
        fprintf(stderr, "angle = %g mrad (%g deg) HWHM\n", angle / 0.001, angle / C_DEG);
        double gamma_b = gamma_beta_TF(tau);
        fprintf(stderr, "Gamma_b = %g\n", gamma_b);
        double wb = wa / gamma_b;
        fprintf(stderr, "wb = %g (HWHM), chord\n", wb);
        double angle_chord = wb / mu;
        fprintf(stderr, "chord angle = %g mrad (%g deg) HWHM\n", angle_chord / 0.001, angle_chord / C_DEG); /* Lateral displacement is this angle times thickness */
        angle_sum += pow2(angle);
        angle_chord_sum += pow2(angle_chord);

    }
    double angle = 2.0 * sqrt(angle_sum); /* FWHM */
    double angle_chord = 2.0 * sqrt(angle_chord_sum);  /* FWHM */
    double dt = thickness * angle_chord * tan(phi);
    fprintf(stderr, "trajectory length difference %g tfu\n", dt / C_TFU);
    double stop = jibal_stop(jibal->gsto, incident, layer->material, E_deep);
    fprintf(stderr, "stopping force at depth = %g eV/tfu\n", stop / C_EV_TFU);
    double dE = dt * stop;
    fprintf(stdout, "%12.3lf %12.3lf %12.3lf %12.3lf %12.3lf\n", depth / C_TFU, thickness / C_TFU, angle / C_DEG, angle_chord / C_DEG, dE / C_KEV);
    jibal_layer_free(layer);
    return 0;
}

int main(int argc, char **argv) {
    if(argc < 5) {
        fprintf(stderr, "Not enough arguments! Usage: amsel <incident isotope> <target element> <max thickness> <incident E lab> <phi>\n");
        return EXIT_FAILURE;
    }
    jibal *jibal = jibal_init(NULL);
    const jibal_isotope *incident = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    if(!incident) {
        fprintf(stderr, "No such isotope: %s\n", argv[1]);
        return EXIT_FAILURE;
    }
    double thickness_max = 0.0;
    fprintf(stderr, "incident = %s\n", incident->name);
    jibal_unit_convert(jibal->units, JIBAL_UNIT_TYPE_LAYER_THICKNESS, argv[3], &thickness_max);
    double E_lab = 0.0;
    jibal_unit_convert(jibal->units, JIBAL_UNIT_TYPE_ENERGY, argv[4], &E_lab);
    fprintf(stderr, "E_lab = %g keV\n", E_lab / C_KEV);

    double phi = 0.0;
    jibal_unit_convert(jibal->units, JIBAL_UNIT_TYPE_ANGLE, argv[5], &phi);
    fprintf(stderr, "phi = %g deg (incidence angle, zero for beam perpendicular to surface)\n", phi / C_DEG);

    amsel_table *t = amsel_table_generate();
    size_t n_steps = 50;
    double thickness_step = thickness_max / (n_steps - 1);

    jibal_material *material = jibal_material_create(jibal->elements, argv[2]);
    if(!material) {
        fprintf(stderr, "Material could not be created from this: %s\n", argv[2]);
        return EXIT_FAILURE;
    }
    if(!jibal_gsto_auto_assign_material(jibal->gsto, incident, material)) {
        fprintf(stderr, "Couldn't assign stopping.\n");
        return EXIT_FAILURE;
    }
    jibal_gsto_load_all(jibal->gsto);
    for(size_t i_step = 0; i_step < n_steps; i_step++) {
        double thickness = thickness_step * i_step;
        amsel_calc(jibal, t, thickness, phi, E_lab, incident, material);
    }
    amsel_table_free(t);
    jibal_free(jibal);
    return EXIT_SUCCESS;
}
