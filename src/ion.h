/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

    Some parts of this source file under different license, see below!

 */
#ifndef JABS_ION_H
#define JABS_ION_H
#include <jibal_masses.h>
#include <jibal_gsto.h>
#include "nuclear_stopping.h"

#ifdef __cplusplus
extern "C" {
#endif
typedef struct jabs_gsto_assignment {
        const gsto_file_t *stopfile;
        const double *stopdata;
        const gsto_file_t *straggfile;
        const double *straggdata;
} jabs_ion_gsto_data;

typedef struct jabs_ion_gsto {
    jabs_ion_gsto_data *gsto_data; /* Array of gsto->Z2_max + 1 elements */
    double emin;
    double emax;
    const jibal_isotope *incident;
    int refcount;
    int Z2_max; /* Copy of gsto->Z2_max */
} jabs_ion_gsto;

typedef struct ion {
    const jibal_isotope *isotope;
    double E;
    double S;
    double mass;
    double mass_inverse;
    int Z;
    double theta; /* polar angle. theta = 0 is along the z axis (deeper into sample, perpendicular to surface) */
    double phi; /* azimuthal angle. phi = 0 is x-axis, phi = 90deg is y-axis. */
    double cosine_theta;
    double inverse_cosine_theta; /* Inverse cosine of theta. Traversing matter "straight on" means 1.0 and going sideways approaches +inf. Negative -1.0 to -inf for travelling "backwards". Values between (-1.0,1.0) are not allowed. */
    double cosine_phi;
    double inverse_cosine_phi;
    nuclear_stopping *nucl_stop;
    jabs_ion_gsto *ion_gsto;
} ion;

void ion_reset(ion *ion);
void ion_set_isotope(ion *ion, const jibal_isotope *isotope);
void ion_set_angle(ion *ion, double theta, double phi);
double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, int accurate);
void ion_rotate(ion *ion, double theta2, double phi2);
void ion_print(FILE *f, const ion *ion);
jabs_ion_gsto *ion_gsto_new(const jibal_isotope *incident, const jibal_gsto *gsto);
jabs_ion_gsto *ion_gsto_shared(jabs_ion_gsto *ion_gsto);
void ion_gsto_free(jabs_ion_gsto *ion_gsto);
#ifdef __cplusplus
}
#endif
#endif // JABS_ION_H
