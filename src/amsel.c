/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

/* Test tool for small angle multiple scattering (ms.c). Prints the angular spread, chord angle spread and energy spread
 * (path length differences only) of an ion beam as a function of depth in a homogeneous material. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <jibal.h>
#include <jibal_units.h>
#include <jibal_stop.h>
#include "ms.h"

int main(int argc, char **argv) {
    if(argc < 6) {
        fprintf(stderr, "Usage: amsel <incident isotope> <material> <max depth> <incident E lab> <angle of incidence> [TF|LJ] [n_out]\n");
        fprintf(stderr, "Angle of incidence is measured from the surface normal.\n");
        return EXIT_FAILURE;
    }
    jibal *jibal = jibal_init(NULL);
    if(!jibal || jibal->error) {
        fprintf(stderr, "Could not initialize JIBAL.\n");
        return EXIT_FAILURE;
    }
    const jibal_isotope *incident = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    if(!incident) {
        fprintf(stderr, "No such isotope: %s\n", argv[1]);
        return EXIT_FAILURE;
    }
    jibal_material *material = jibal_material_create(jibal->elements, argv[2]);
    if(!material) {
        fprintf(stderr, "Material could not be created from this: %s\n", argv[2]);
        return EXIT_FAILURE;
    }
    double depth_max = 0.0, E = 0.0, alpha = 0.0;
    jibal_unit_convert(jibal->units, JIBAL_UNIT_TYPE_LAYER_THICKNESS, argv[3], &depth_max);
    jibal_unit_convert(jibal->units, JIBAL_UNIT_TYPE_ENERGY, argv[4], &E);
    jibal_unit_convert(jibal->units, JIBAL_UNIT_TYPE_ANGLE, argv[5], &alpha);
    ms_model_type type = MS_MODEL_AMSEL_TF;
    if(argc > 6) {
        type = jibal_option_get_value(ms_model_option, argv[6]);
    }
    int n_out = 50;
    if(argc > 7) {
        n_out = atoi(argv[7]);
    }
    if(!jibal_gsto_auto_assign_material(jibal->gsto, incident, material)) {
        fprintf(stderr, "Couldn't assign stopping.\n");
        return EXIT_FAILURE;
    }
    jibal_gsto_load_all(jibal->gsto);
    ms_model *ms = ms_model_new(type);
    if(!ms) {
        fprintf(stderr, "Could not initialize multiple scattering model.\n");
        return EXIT_FAILURE;
    }
    fprintf(stderr, "%s in %s, E = %g keV, angle of incidence %g deg, model %s\n", incident->name, argv[2], E / C_KEV, alpha / C_DEG, ms_model_name(type));
    const int substeps = 50;
    const double inv_cos = 1.0 / cos(alpha);
    const double dl = depth_max * inv_cos / (n_out * substeps); /* Path length step */
    ion_ms s;
    s.w = 0.0;
    s.sum_a = 0.0;
    s.sum_b = 0.0;
    s.I = 0.0;
    s.P = 1.0;
    fprintf(stdout, "#   depth(tfu)   E(keV) angle_FWHM(deg) chord_FWHM(deg)      tau*         nu   dE_path_FWHM(keV)\n");
    for(int i_out = 0; i_out <= n_out; i_out++) {
        double C_b = 1.0;
        double tau_star = ms_tau_star(ms, &s);
        double nu = 0.0;
        ms_table_params(ms, tau_star, &C_b, NULL, &nu);
        double w_chord = s.w / C_b;
        double dE_path = s.I * tan(alpha) * w_chord; /* Path length l * (1 + beta * tan(alpha)), HWHM */
        fprintf(stdout, "%12.3lf %10.3lf %12.5lf %12.5lf %12.5g %10.4lf %12.4lf\n", i_out * dl * substeps / inv_cos / C_TFU, E / C_KEV,
                2.0 * s.w / C_DEG, 2.0 * w_chord / C_DEG, tau_star, nu, 2.0 * dE_path / C_KEV);
        if(i_out == n_out) {
            break;
        }
        for(int i = 0; i < substeps; i++) {
            double S = jibal_stop(jibal->gsto, incident, material, E);
            double E_mid = E - 0.5 * dl * S;
            double S_mid = jibal_stop(jibal->gsto, incident, material, E_mid);
            double E_new = E - dl * S_mid;
            if(E_new <= 0.0) {
                fprintf(stderr, "Ion stopped.\n");
                i_out = n_out;
                break;
            }
            double S_end = jibal_stop(jibal->gsto, incident, material, E_new);
            for(size_t i_elem = 0; i_elem < material->n_elements; i_elem++) {
                ms_step_element(ms, &s, incident->Z, material->elements[i_elem].Z, material->concs[i_elem], E_mid, dl);
            }
            s.I = s.I * S_end / S + dl * S_end;
            E = E_new;
        }
    }
    ms_model_free(ms);
    jibal_material_free(material);
    jibal_free(jibal);
    return EXIT_SUCCESS;
}
