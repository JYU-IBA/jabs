/* rotate.c from MCERD, https://github.com/JYU-IBA/mcerd,
 * Some modifications made by Jaakko Julin.
 *
 * MCERD is Copyright (C) 1996-2020  Kai Arstila, 2015-2021 Jaakko Julin
 * and used here under terms of the GNU General Public License v2.
 *
 */

#include <math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

#ifndef M_PI_2
#define M_PI_2 (3.14159265358979323846264338327950288/2.0)
#endif

#ifndef TWO_PI
#define TWO_PI (2.0*M_PI)
#endif

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

    sina2 = sin(phi2 + M_PI_2);
    cosa2 = cos(phi2 + M_PI_2);

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
#if 0
    if (fabs(*phi) > TWO_PI) {
        *phi = fmod(*phi, 2.0 * PI);
    }
    if (*phi < 0.0) {
        *phi += 2.0 * PI;
    }
#endif
    *phi = fabs(fmod(*phi, TWO_PI));
}

