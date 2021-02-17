#ifndef JABS_ION_H
#define JABS_ION_H
#include <jibal_masses.h>

typedef struct {
    const jibal_isotope *isotope;
    double E;
    double S;
    double mass;
    int Z;
    double angle; /* Traversing matter "straight on" when angle = 0, getting stuck sideways if angle = 90.0*C_DEG */
    double cosine;
    double inverse_cosine; /* Inverse cosine of angle. Traversing matter "straight on" means 1.0 and going sideways approaches infinity. */
} ion;

void ion_set_isotope(ion *ion, const jibal_isotope *isotope);
void ion_set_angle(ion *ion, double angle);
#endif // JABS_ION_H
