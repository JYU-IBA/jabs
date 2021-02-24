#include <stdlib.h>
#include <jibal.h>
#include <jibal_masses.h>
#include <ion.h>
#include "rotate.h"

int main(int argc, char **argv) {
    if(argc < 3) {
        fprintf(stderr, "Not enough arguments. Usage: %s tilt detector_angle [detector_angle_phi]\n", argv[0]);
        return -1;
    }
    jibal *jibal = jibal_init(NULL);
    ion ion;
    double alpha, beta, theta_scatt; /* These will be calculated later */
    double sample_theta = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, argv[1]); /* Note that sample_theta can be negative! [-pi, pi] is ok. If the sample is tilted AWAY from the detector, use negative angles. */
    double sample_phi = 0.0 * C_DEG; /* 0 + n*PI is IBM, +/-90 deg + n*2*PI is Cornell */
    ion_reset(&ion);
    ion_set_isotope(&ion, jibal_isotope_find(jibal->isotopes, "4He", 0, 0));
    ion_print(stdout, &ion);
    ion_rotate(&ion, sample_theta, sample_phi); /* Ion in sample coordinate system (sample rotated to sample_theta along one axis. This is typical unless you do channeling.) */
    alpha = ion.theta;

    ion_print(stdout, &ion);
    double theta, phi;
    double detector_theta = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, argv[2]);
    double detector_phi = 0.0;
    if(argc >= 4) {
        detector_phi = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, argv[3]);
    }
    rotate(detector_theta, detector_phi, sample_theta, sample_phi, &theta, &phi); /* Detector in sample coordinate system, angles are detector in sample system. Note that for Cornell geometry phi = 90.0 deg! */
    fprintf(stdout, "Result: theta = %.3lf, phi = %.3lf\n", theta/C_DEG, phi/C_DEG);
    beta = (C_PI-theta);
    ion_rotate(&ion, -sample_theta, -sample_phi); /* Unrotate sample tilt */
    ion_rotate(&ion, detector_theta, detector_phi); /* Rotate to detector */
    theta_scatt = ion.theta;

    ion_print(stdout, &ion);
    fprintf(stdout, "Alpha = %.3lf deg.\n", alpha/C_DEG);
    fprintf(stdout, "Beta = %.3lf deg.\n", beta/C_DEG);
    fprintf(stdout, "Theta = %.3lf deg.\n", theta_scatt/C_DEG);
    jibal_free(jibal);
    return 0;
}
