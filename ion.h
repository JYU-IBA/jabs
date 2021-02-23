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
    int Z;
    double angle; /* Traversing matter "straight on" when angle = 0, getting stuck sideways if angle = 90.0*C_DEG */
    double cosine;
    double inverse_cosine; /* Inverse cosine of angle. Traversing matter "straight on" means 1.0 and going sideways approaches infinity. */
    nucl_stop_pair *nucl_stop;
    int nucl_stop_isotopes;
} ion;

void ion_set_isotope(ion *ion, const jibal_isotope *isotope);
void ion_set_angle(ion *ion, double angle);
double ion_nuclear_stop(const ion *ion, const jibal_isotope *isotope, const jibal_isotope *isotopes);
void ion_nuclear_stop_fill_params(ion *ion, const jibal_isotope *isotopes, int n_isotopes); /* allocates ion->nucl_stop, must be free'd */
#endif // JABS_ION_H
