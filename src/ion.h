#ifndef JABS_ION_H
#define JABS_ION_H
#include <jibal_masses.h>
#include "nuclear_stopping.h"

typedef struct {
    const jibal_isotope *isotope;
    double E;
    double S;
    double mass;
    double mass_inverse;
    int Z;
    double theta; /* polar angle. theta = 0 is along the z axis (deeper into sample, perpendicular to surface) */
    double phi; /* azimuthal angle. phi = 0 is x-axis, phi = 90deg is y-axis. */
    double cosine_theta;
    double inverse_cosine_theta; /* Inverse cosine of theta. Traversing matter "straight on" means 1.0 and going sideways approaches infinity. */
    double cosine_phi;
    double inverse_cosine_phi;
    nuclear_stopping *nucl_stop;
} ion;

void ion_reset(ion *ion);
void ion_set_isotope(ion *ion, const jibal_isotope *isotope);
void ion_set_angle(ion *ion, double theta, double phi);
double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, int accurate);
void ion_rotate(ion *ion, double theta2, double phi2);
void ion_print(FILE *f, const ion *ion);
#endif // JABS_ION_H
