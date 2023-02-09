#include <jibal.h>
#include "ion.h"
#include "jabs.h"

int main(int argc, char **argv) {
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code %i (%s)\n", jibal->error, jibal_error_string(jibal->error));
        return EXIT_FAILURE;
    }
    ion testion;
    ion_reset(&testion);
    ion_set_isotope(&testion, jibal_isotope_find(jibal->isotopes, "4He", 0, 0));
    testion.E = 1.0 * C_MEV;
    testion.nucl_stop = nuclear_stopping_new(testion.isotope, jibal->isotopes);
#if 1
    ion_rotate(&testion, 180.0 * C_DEG, 0.0);
#endif
    sample_model *sm = sample_model_from_string(jibal, "Au 100tfu SiO2 1000tfu Si 300000tfu");
    sample *sample = sample_from_sample_model(sm);
#if 1
    depth depth_start = depth_seek(sample, 20000.0 * C_TFU);
#else
    depth depth_start = depth_seek(sample, 0.0 * C_TFU);
#endif
    simulation *sim = sim_init(jibal);
    sim->sample = sample;
    sim_workspace *ws = sim_workspace_init(jibal, sim, sim->det[0]);
    assign_stopping(jibal->gsto, sim);
    jibal_gsto_load_all(jibal->gsto);
    des_table *dt = des_table_compute(&ws->stop, &ws->stragg, ws->params, sample, &testion, depth_start, ws->emin);
    des_table_print(stderr, dt);
    reaction *r = reaction_make(testion.isotope, jibal_isotope_find(jibal->isotopes, "28Si", 0, 0), REACTION_RBS, JABS_CS_ANDERSEN);
    sim_reaction *sim_r = sim_reaction_init(sample, ws->det, r, ws->n_channels, ws->n_bricks);
    geostragg_vars g = geostragg_vars_calculate(&testion, 0.0, 0.0, ws->det, NULL, FALSE, FALSE);
    simulate_reaction(&testion, depth_start, ws, sample, dt, &g, sim_r);
    return EXIT_SUCCESS;
}
