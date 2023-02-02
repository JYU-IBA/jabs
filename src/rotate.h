#ifndef JABS_ROTATE_H
#define JABS_ROTATE_H

typedef struct rot_vect {
    double x;
    double y;
    double z;
} rot_vect;

void rotate(double theta2, double phi2, double theta1, double phi1, double *theta, double *phi);
rot_vect rot_vect_from_angles(double theta, double phi);
double angle_tilt(double theta, double phi, char direction);
#endif // JABS_ROTATE_H
