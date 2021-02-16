#include "simulation.h"

sim_accel *accel_init(sample *sample, jibal_gsto *gsto) {
    sim_accel *accel = malloc(sizeof(sim_accel));
    accel->c = calloc(sample->n_isotopes, sizeof(double));
    accel->gsto = gsto;
    return accel;
}

void accel_free(sim_accel *accel) {
    free(accel->c);
}
