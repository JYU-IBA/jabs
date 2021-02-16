#include <math.h>
#ifdef ERF_Q_FROM_GSL
#include <gsl/gsl_sf.h>
#endif
#include "brick.h"

inline double erfc(double x) {
    return exp(-1.0950081470333*x-0.75651138383854*x*x);
    /* Tsay, WJ., Huang, C.J., Fu, TT. et al. J Prod Anal 39, 259–269 (2013). https://doi.org/10.1007/s11123-012-0283-1 */
}
extern inline double erf_Q(double x);

inline double erf_Q(double x) {
    return x < 0.0 ? 1.0-0.5*erfc(-1.0*x/M_SQRT2) : 0.5*erfc(x/M_SQRT2);
}

#define ERFQ_A (-1.0950081470333/M_SQRT2)
#define ERFQ_B (-0.75651138383854*0.5)
#define ERFQ_CUTOFF 5.0 // 5.0 is enough

double erf_Q_optim(double x);

inline double erf_Q_optim(double x) {
    return x < (-ERFQ_CUTOFF) ? 1.0 : x > (ERFQ_CUTOFF) ? 0.0 : x < 0.0? 1.0-0.5*exp(-1.0*ERFQ_A*x+ERFQ_B*x*x) : 0.5*exp(ERFQ_A*x+ERFQ_B*x*x);
}

void erf_Q_test() {
    double x;
    int i;
    double sigma = 1.0;
    double x_0 = 100.0;
    double x_L = 10.0;
    double x_H = 15.0;
    double sum = 0.0;
    double step = 0.03;
    for(i=0; i <10000; i++) {
        x = -50.0 + i*step;
#ifdef ERF_Q_FROM_GSL
        double fx = gsl_sf_erf_Q((x_H-x)/sigma);
        double fx_low = gsl_sf_erf_Q((x_L-x)/sigma);
#else
        double fx = erf_Q((x_H-x)/sigma);
        double fx_low = erf_Q((x_L-x)/sigma);
#endif
        double out = (fx_low-fx)/(x_H-x_L); /* One count, convolution of gaussian rectangle */
        sum += out*step; /* This is the integral! */
        fprintf(stdout, "%g %g %g %12.8lf\n", x, fx, out, sum);
    }
}

void brick_int(double sigma_low, double sigma_high, double E_low, double E_high, gsl_histogram *h, double Q) { /* Energy spectrum h, integrate over rectangular brick convoluted with a gaussian */
    int i;
#ifndef OPTIMIZE_BRICK
    size_t lo;
    size_t hi;
    gsl_histogram_find(h, E_low, &lo);
    gsl_histogram_find(h, E_high, &hi);

#ifdef SMART
    double foo = (E_high-E_low)/(3.0); /* TODO: smart choice */
    int n = ceil(foo*5);

    if(n < 50)
        n = 50;
#else
    int n = 50; /* TODO: stupid choice */
#endif

    if(lo > n)
        lo -= n;
    else
        lo = 0;
    if(hi > h->n-n-1)
        hi = h->n;
    else
        hi += n;
    for(i = lo; i <= hi; i++) {
#else
        for(i = 0; i < h->n; i++) {
#endif
        double w = h->range[i+1] - h->range[i];
        double E = (h->range[i] + h->range[i+1])/2.0; /* Approximate gaussian at center bin */
#ifdef ERF_Q_FROM_GSL
        double y = (gsl_sf_erf_Q((E_low-E)/sigma)-gsl_sf_erf_Q((E_high-E)/sigma))/(E_high-E_low);
#else
#ifdef SIGMA_MEAN
        double sigma = (sigma_low+sigma_high)/2.0;
        double y = (erf_Q((E_low-E)/sigma)-erf_Q((E_high-E)/sigma))/(E_high-E_low);
#else
        double y = (erf_Q_optim((E_low-E)/sigma_low)-erf_Q_optim((E_high-E)/sigma_high))/(E_high-E_low);
#endif
#endif
        h->bin[i] += y*w*Q;
    }
}