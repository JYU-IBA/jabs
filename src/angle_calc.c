#include <stdlib.h>
#include <jibal.h>
#include <jibal_masses.h>
#include <ion.h>
#include "rotate.h"

void exit_calculator(double detector_theta, double detector_phi, double sample_theta, double sample_phi) {
    double alpha, beta, theta_scatt; /* These will be calculated later */
    ion ion;
    ion_reset(&ion);
    double theta, phi; /* temps */
    ion_rotate(&ion, sample_theta, sample_phi); /* Ion in sample coordinate system (sample rotated to sample_theta along one axis. This is typical unless you do channeling.) */
    alpha = ion.theta;
    rotate(detector_theta, detector_phi, sample_theta, sample_phi, &theta, &phi); /* Detector in sample coordinate system, angles are detector in sample system. Note that for Cornell geometry phi = 90.0 deg! */
    beta = (C_PI-theta);

    ion_rotate(&ion, -sample_theta, -sample_phi); /* Unrotate sample tilt */
    ion_rotate(&ion, detector_theta, detector_phi); /* Rotate to detector */
    theta_scatt = ion.theta; /* There is a more direct way to calculate this (hint detector_theta), but we'll double check it :) */
    fprintf(stdout, "alpha = %.3lf deg\n", alpha/C_DEG);
    fprintf(stdout, "beta = %.3lf deg\n", beta/C_DEG);
    fprintf(stdout, "theta = %.3lf deg\n", theta_scatt/C_DEG);
}

void ds_test(ion *ion, double detector_theta, double detector_phi, double sample_theta, double sample_phi) {
    //double theta, phi;
    ion_print(stdout, ion);

    ion_print(stdout, ion);
    fprintf(stdout, "\n");
    int ds_steps_polar = 10;
    int ds_steps_azi = 12;
    int n_ds = ds_steps_polar*ds_steps_azi;
    int i_ds;
    for(i_ds = 0; i_ds < n_ds; i_ds++) {
        int i_polar = i_ds / ds_steps_azi;
        double cosine = (ds_steps_polar - 2 * i_polar) / (1.0 * ds_steps_polar);
        double ds_polar = acos(cosine);
        int i_azi = i_ds % ds_steps_azi;
        double ds_azi = C_2_PI * (1.0 * i_azi) / (ds_steps_azi * 1.0);
        fprintf(stdout, "%3i %3i %3i %12g %12g %12g", i_ds, i_polar, i_azi, cosine, ds_polar/C_DEG, ds_azi/C_DEG);
        ion_set_angle(ion, 0.0, 0.0);
        ion_rotate(ion, sample_theta, sample_phi); /* Go to sample coordinates, sample tilted */



        ion_rotate(ion, ds_polar, ds_azi); /* Dual scattering: first scattering to some angle (scattering angle: ds_polar). Note that this does not follow SimNRA conventions. */
        //ion_print(stdout, ion);
        //fprintf(stdout, "DS first scattering: %.3lf deg\n", ds_polar / C_DEG);
        double theta_after_first = ion->theta;


        ion_rotate(ion, -ds_polar, -ds_azi); /* Undo scattering rotation */
        ion_rotate(ion, detector_theta, detector_phi); /* Rotate to detector (the angle is relative to lab so we still stay in sample coordinates. This rotation should actually be: 1) rotate to lab, 2) rotate to detector, 3) rotate to sample, but steps 1 and 3 cancel each other out, so we just do step 2.) */
        //ion_print(stdout, ion);
        double theta_after_second = ion->theta;
        double ds_polar_second = fabs(theta_after_second - theta_after_first); /* Second scattering angle */
        //fprintf(stdout, "DS second scattering: %.3lf deg\n", ds_polar_second / C_DEG);
        ion_rotate(ion, -sample_theta, -sample_phi); /* Rotate to lab */
        //ion_print(stdout, ion);


        /* Simulate a spectrum from depth x */
        //fprintf(stdout, "New ion alpha angle is %.3lf deg and second scattering angle is %.3lf\n", theta_after_first / C_DEG, ds_polar_second / C_DEG);
        //fprintf(stdout, "\n");
        fprintf(stdout, " %12g %12g\n", ds_polar_second/C_DEG, theta_after_first/C_DEG);
    }
}


int main(int argc, char **argv) {
    if(argc < 3) {
        fprintf(stderr, "Not enough arguments. Usage: %s tilt detector_angle [detector_angle_phi]\n", argv[0]);
        return -1;
    }
    jibal *jibal = jibal_init(NULL);

    double sample_theta = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, argv[1]); /* Note that sample_theta can be negative! [-pi, pi] is ok. If the sample is tilted AWAY from the detector, use negative angles. */
    double sample_phi = 0.0 * C_DEG; /* 0 + n*PI is IBM, +/-90 deg + n*2*PI is Cornell */



    double detector_theta = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, argv[2]);
    double detector_phi = 0.0;
    if(argc >= 4) {
        detector_phi = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, argv[3]);
    }
    fprintf(stderr, "sample_theta = %.3lf deg\ndetector_theta = %.3lf deg\ndetector_phi = %.3lf\n", sample_theta/C_DEG, detector_theta/C_DEG, detector_phi/C_DEG);

    exit_calculator(detector_theta, detector_phi, sample_theta, sample_phi);
    return 0;
    ion ion;
    ion_reset(&ion);
    ion.E = 2.0*C_MEV;
    ion_set_isotope(&ion, jibal_isotope_find(jibal->isotopes, "4He", 0, 0));
    ds_test(&ion, detector_theta, detector_phi, sample_theta, sample_phi);

    jibal_free(jibal);
    return 0;
}
