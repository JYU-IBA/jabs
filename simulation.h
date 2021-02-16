#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <gsl/gsl_histogram.h>
#include <jibal_masses.h>
#include <jibal_gsto.h>



#include "ion.h"
#include "reaction.h"

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
    ion ion;
} simulation;

typedef struct {
    int n_reactions; /* Same as sim->n_reactions, but we want to keep it here too, as it is used for allocations */
    double *c; /* Concentrations for n_isotopes at some arbitrary x */
    double c_x; /* at this x */
    int i_range_accel;
    jibal_gsto *gsto;
    int rk4;
    gsto_stopping_type stopping_type;
    gsl_histogram **histos; /* array of n_reactions */
} sim_workspace;

#include "sample.h"

sim_workspace *sim_workspace_init(const simulation *sim, sample *sample, jibal_gsto *gsto);
void sim_workspace_free(sim_workspace *ws);

#endif /* _SIMULATION_H_ */
