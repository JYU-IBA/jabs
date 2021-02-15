#include <math.h>
#ifdef ERF_Q_FROM_GSL
#include <gsl/gsl_sf.h>
#endif

#include "brick.h"

inline double erfc(double x) {
    return exp(-1.0950081470333*x-0.75651138383854*x*x);
    /* Tsay, WJ., Huang, C.J., Fu, TT. et al. J Prod Anal 39, 259–269 (2013). https://doi.org/10.1007/s11123-012-0283-1 */
}
inline double erf_Q(double x);

double erf_Q(double x) {
    return x < 0.0 ? 1.0-0.5*erfc(-1.0*x/sqrt(2.0)) : 0.5*erfc(x/sqrt(2.0));
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

void brick_int(double sigma, double E_low, double E_high, gsl_histogram *h, double Q) { /* Energy spectrum h, integrate over rectangular brick convoluted with a gaussian */
    int i;
#ifdef OPTIMIZE_BRICK
    size_t lo;
    size_t hi;
    gsl_histogram_find(h, E_low, &lo);
    gsl_histogram_find(h, E_high, &hi);
    double foo = (E_high-E_low)/(HISTOGRAM_BIN);
    int n = ceil(foo*5);
    if(n < 50)
        n = 50;
    if(lo > n)
        lo -= n;
    else
        lo = 0;
    if(hi > h->n-n)
        hi = h->n;
    else
        hi += n;
    for(i = lo; i < hi; i++) {
#else
        for(i = 0; i < h->n; i++) {
#endif
        double w = h->range[i+1] - h->range[i];
        double E = (h->range[i] + h->range[i+1])/2.0; /* Approximate gaussian at center */
#ifdef ERF_Q_FROM_GSL
        double y = (gsl_sf_erf_Q((E_low-E)/sigma)-gsl_sf_erf_Q((E_high-E)/sigma))/(E_high-E_low);
#else
        double y = (erf_Q((E_low-E)/sigma)-erf_Q((E_high-E)/sigma))/(E_high-E_low);
#endif
        h->bin[i] += y*w*Q;
    }
}
