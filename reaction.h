#ifndef _REACTION_H_
#define _REACTION_H_

#include <jibal_masses.h>
#include <gsl/gsl_histogram.h>
#include "ion.h"

typedef struct {
    const jibal_isotope *isotope; /* target isotope */
    int i_isotope; /* location of target isotope in concentration table */
    ion p; /* reaction product, e.g. in case of He-RBS this will be He. */
    double K;
    double max_depth;
    /* TODO: type? */
    /* TODO: cross section to use? */
    /* TODO: cross sections from files */
    double E; /* Previous energy bin, TODO: make a histogram instead! */
    double S; /* Previous straggling bin */
    int stop;
} reaction;

#endif //_REACTION_H_
