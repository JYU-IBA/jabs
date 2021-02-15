#ifndef _BRICK_H_
#define _BRICK_H_

#include <gsl/gsl_histogram.h>

void brick_int(double sigma_low, double  sigma_high, double E_low, double E_high, gsl_histogram *h, double Q);

#endif /* _BRICK_H_ */
