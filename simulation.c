#include "simulation.h"

sim_workspace *sim_workspace_init(const simulation *sim, sample *sample, jibal_gsto *gsto) {
    sim_workspace *ws = malloc(sizeof(sim_workspace));
    ws->n_reactions = sim->n_reactions;
    ws->c = calloc(sample->n_isotopes, sizeof(double));
    ws->gsto = gsto;
    ws->rk4 = 1;
    ws->stopping_type = GSTO_STO_TOT;
    ws->i_range_accel = 0;
    ws->c_x = 0.0;
    get_concs(ws, sample, ws->c_x, ws->c);
    ws->histos = calloc(ws->n_reactions, sizeof(gsl_histogram *));
    return ws;
}

void sim_workspace_free(sim_workspace *ws) {
    int i;
    for(i = 0; i < ws->n_reactions; i++) {
        gsl_histogram_free(ws->histos[i]);
    }
    free(ws->c);
    free(ws->histos);
    free(ws);
}
