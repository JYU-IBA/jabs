#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <jibal_masses.h>
#include <jibal_gsto.h>

typedef struct {
    double *c; /* Concentrations for n_isotopes at some arbitrary x */
    double c_x; /* at this x */
    int i_range_accel;
    jibal_gsto *gsto;
    int rk4;
    gsto_stopping_type stopping_type;
} sim_workspace;

#include "ion.h"
#include "sample.h"
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
    const reaction *reactions; /* size AT LEAST n_reactions */
} simulation;

sim_workspace *sim_workspace_init(sample *sample, jibal_gsto *gsto);
void sim_workspace_free(sim_workspace *ws);

#endif /* _SIMULATION_H_ */
