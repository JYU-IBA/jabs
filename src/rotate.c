/*
 * rotate() from rotate.c from MCERD, https://github.com/JYU-IBA/mcerd,
 * Some modifications made by Jaakko Julin.
 *
 * MCERD is Copyright (C) 1996-2020  Kai Arstila, 2015-2021 Jaakko Julin
 * and used here under terms of the GNU General Public License v2.
 *
 * For the rest
 *
    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.
 *
 *
 */

#include <math.h>
#include "jibal_units.h"

#define max(A, B)  ((A) > (B)) ? (A) : (B)
#define min(A, B)  ((A) < (B)) ? (A) : (B)

#include "rotate.h"

void rotate(double theta2, double phi2, double theta1, double phi1,
            double *theta, double *phi) {
/*
 *   Definition of the angles: theta2,phi2 is the angle of the
 *   second coordinate system in the first coordinate system.
 *   Theta1,phi1 is the specific direction in the second coordinate
 *   system, which direction in the first system is calculated
 *   in this routine (theta,phi).
 *
 *   This routine cannot be explained easily. Read eg. Goldstein
 *   about the Euler angles and coordinate transforms. Typical
 *   time for understanding this is three days ;-)
 */

    double cos_theta;
    double x, y, z, rx, ry, rz;
    double cosa1, cosa2, cosa3, sina1, sina2, sina3;

    double sin_theta, sin_fii, cos_fii;

    cos_theta = cos(theta1);
    sin_theta = sin(theta1);
    cos_fii = cos(phi1);
    sin_fii = sin(phi1);

    x = sin_theta * cos_fii;
    y = sin_theta * sin_fii;
    z = cos_theta;

    cosa1 = cos(theta2);
    sina1 = sin(theta2);

    sina2 = sin(phi2 + C_PI_2);
    cosa2 = cos(phi2 + C_PI_2);

    cosa3 = cosa2;
    sina3 = -sina2;

    rx = x * (cosa3 * cosa2 - cosa1 * sina2 * sina3) +
         y * (-sina3 * cosa2 - cosa1 * sina2 * cosa3) +
         z * sina1 * sina2;

    ry = x * (cosa3 * sina2 + cosa1 * cosa2 * sina3) +
         y * (-sina3 * sina2 + cosa1 * cosa2 * cosa3) -
         z * sina1 * cosa2;

    rz = x * sina1 * sina3 +
         y * sina1 * cosa3 +
         z * cosa1;

    rz = max(min(rz, 1.0), -1.0);

    *theta = acos(rz);
    if (rx != 0.0) {
        *phi = atan2(ry, rx);
    } else {
        *phi = 0.0;
    }
    *phi = fabs(fmod(*phi, C_2PI));
}

rot_vect rot_vect_from_angles(double theta, double phi) { /* Calculates 3D cartesian unit vector */
    double sintheta = sin(theta);
    double costheta = cos(theta);
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    rot_vect v = {sintheta * cosphi, sintheta * sinphi, costheta};
    return v;
}

double angle_tilt(double theta, double phi, char direction) { /* Directions 'x' and 'y' mean "horizontal" and "vertical", respectively. X and Y angles are measured from Z-axis and are orthogonal. */
    double angle = 0.0;
    rot_vect v = rot_vect_from_angles(theta, phi);
    if(direction == 'x') {
        angle = atan2(v.x, v.z);
    } else if(direction == 'y') {
        angle = atan2(v.y, v.z);
    }
    return angle;
}

double theta_tilt(double tilt_x, double tilt_y) {
    double theta = 0.0, phi = 0.0;
    rotate(tilt_x, 0.0, theta, phi, &theta, &phi);
    rotate(tilt_y, C_PI_2, theta, phi, &theta, &phi);
#ifdef DEBUG
    fprintf(stderr, "By rotating first by %.7lf deg (phi 0) and then %.7lf deg (phi 90 deg) we get angles %.7lf deg and %.7lf deg.\n",
            tilt_x/C_DEG, tilt_y/C_DEG,
            theta/C_DEG, phi/C_DEG);
#endif
    return theta;
}
