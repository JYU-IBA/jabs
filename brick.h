#ifndef JABS_BRICK_H
#define JABS_BRICK_H

#include <gsl/gsl_histogram.h>

void brick_int(double sigma_low, double  sigma_high, double E_low, double E_high, gsl_histogram *h, double Q);

#endif // JABS_BRICK_H
