#ifndef JABS_ION_H
#define JABS_ION_H
#include <jibal_masses.h>

typedef struct {
    double k;
    double eps0;
} nucl_stop_pair;

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
    nucl_stop_pair *nucl_stop;
    size_t nucl_stop_isotopes;
} ion;

void ion_reset(ion *ion);
void ion_set_isotope(ion *ion, const jibal_isotope *isotope);
void ion_set_angle(ion *ion, double theta, double phi);
double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, const jibal_isotope *isotopes, int accurate);
void ion_nuclear_stop_fill_params(ion *ion, const jibal_isotope *isotopes, int n_isotopes); /* allocates ion->nucl_stop, must be free'd */
void ion_rotate(ion *ion, double theta2, double phi2);
void ion_print(FILE *f, const ion *ion);
#endif // JABS_ION_H
