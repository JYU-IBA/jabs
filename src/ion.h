/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

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
typedef struct jabs_gsto_table { /* Stopping or straggling of an incident ion in one element (Z2) */
    const gsto_file_t *file; /* NULL if not assigned */
    const double *data; /* Values at file->em */
    double unit_factor; /* Converts data to SI units */
    double *lnem; /* ln(em) at nodes, NULL if interpolation is linear in em */
    double *deriv; /* Derivative of data with respect to ln(em) at nodes, of a cubic spline through data in ln(em). NULL if interpolation is linear in em. */
    double lnem0; /* ln(em) of the first node */
    double inv_dlnem; /* Inverse of node spacing in ln(em) if nodes are equally spaced in ln(em) (speedup), otherwise zero */
} jabs_gsto_table;

typedef struct jabs_gsto_assignment {
    jabs_gsto_table stop;
    jabs_gsto_table stragg;
} jabs_ion_gsto_data;

typedef struct jabs_ion_gsto {
    jabs_ion_gsto_data *gsto_data; /* Array of gsto->Z2_max + 1 elements */
    double emin; /* Stopping and straggling data is available for all assigned target elements in energy range [emin, emax] */
    double emax;
    const gsto_file_t *emin_file; /* File (stopping or straggling) that determines emin, NULL if none */
    int emin_Z2; /* Target element of emin_file */
    const gsto_file_t *emax_file;
    int emax_Z2;
    const jibal_isotope *incident;
    int refcount;
    int Z2_max; /* Copy of gsto->Z2_max */
} jabs_ion_gsto;

typedef struct ion_ms { /* Small angle multiple scattering state of ion, accumulated along the path. See ms.h. */
    double w; /* HWHM of the projected angular distribution (radians) */
    double sum_a; /* Sums for weighted mean of scaling factor mu (Szilagyi et al. NIMB 100 (1995) 103, Eq. (50)), sum of tau * mu^(1-nu) */
    double sum_b; /* sum of tau * mu^(-nu) */
    double I; /* Energy change caused by relative path length change (dE/df, where path length is l * (1 + f)), accumulated along path */
    double P; /* Propagation factor of energy differences from start of path to current point (dE/dE_start) */
} ion_ms;

typedef struct ion {
    const jibal_isotope *isotope;
    double E;
    double S;
    ion_ms ms;
    double flux; /* Fraction of incident beam remaining, i.e. not lost due to (counted) dual scattering */
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
void ion_ms_reset(ion_ms *ms); /* Resets multiple scattering state, i.e. start of a new path */
void ion_set_isotope(ion *ion, const jibal_isotope *isotope);
void ion_set_angle(ion *ion, double theta, double phi);
double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, double E);
void ion_rotate(ion *ion, double theta2, double phi2);
void ion_print(FILE *f, const ion *ion);
void ion_gsto_table_init(jabs_gsto_table *t, const gsto_file_t *file, const double *data, double unit_factor);
void ion_gsto_table_free(jabs_gsto_table *t);
double ion_gsto_table_eval(const jabs_gsto_table *t, double em);
jabs_ion_gsto *ion_gsto_new(const jibal_isotope *incident, const jibal_gsto *gsto, const int *Z2_used); /* Z2_used is an array of gsto->Z2_max + 1 elements, data is used only for Z2 where Z2_used[Z2] is TRUE. NULL means all assigned Z2. */
jabs_ion_gsto *ion_gsto_shared(jabs_ion_gsto *ion_gsto);
int ion_gsto_warn_range(const jabs_ion_gsto *ion_gsto, double emin, double emax); /* Warns (jabs_message()) if stopping or straggling data does not cover energy range [emin, emax]. Simulation is cut off below ion_gsto->emin, values at ion_gsto->emax are used above it. Returns TRUE if warnings were issued. */
void ion_gsto_free(jabs_ion_gsto *ion_gsto);
#ifdef __cplusplus
}
#endif
#endif // JABS_ION_H
