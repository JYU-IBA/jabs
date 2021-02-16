#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <gsl/gsl_histogram.h>
#include <jibal_masses.h>
#include <jibal_gsto.h>

#include "ion.h"
#include "sample.h"

typedef struct {
    double *c; /* Concentrations for n_isotopes at some arbitrary x */
    double c_x; /* at this x */
    jibal_gsto *gsto;
    int rk4;
    gsto_stopping_type stopping_type;
} sim_workspace;


typedef struct {
    const jibal_isotope *isotope; /* target isotope */
    int i_isotope; /* location of target isotope in concentration table */
    ion p; /* reaction product, e.g. in case of He-RBS this will be He. */
    double K;
    gsl_histogram *histo; /* energy histogram */
    double max_depth;
    /* TODO: type? */
    /* TODO: cross section to use? */
    /* TODO: cross sections from files */
    double E; /* Energy bin */
    double S;
    int stop;
} reaction;


typedef struct {
    int n_channels;
    int n_reactions;
    double histogram_bin;
    const jibal_isotope *incident;
    double p_sr;
    double p_sr_cos_alpha; /* particles * sr / cos(alpha) */
    double alpha;
    double beta;
    double theta;
    ion ion; /* TODO: this is altered by the simulation */
    reaction *reactions; /* size AT LEAST n_reactions */
} simulation;


sim_workspace *sim_workspace_init(sample *sample, jibal_gsto *gsto);
void sim_workspace_free(sim_workspace *ws);

#endif /* _SIMULATION_H_ */
