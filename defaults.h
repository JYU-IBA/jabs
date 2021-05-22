#ifndef JABS_DEFAULTS_H
#define JABS_DEFAULTS_H
#include <jibal_units.h>

/* Defaults for new simulations */
#define ENERGY (2.0*C_MEV)
#define ALPHA (0.0*C_DEG)
#define DETECTOR_THETA (165.0*C_DEG)
#define DETECTOR_PHI (0.0 * C_DEG) /* IBM geometry */
#define DETECTOR_RESOLUTION (15.0*C_KEV/C_FWHM)
#define PARTICLES_SR (1.0e12)
#define E_MIN (50.0*C_KEV)
#define ENERGY_SLOPE (1.0*C_KEV)
#define STOP_STEP_INCIDENT (0.0*C_KEV) /* Zero is automatic */
#define STOP_STEP_EXITING (0.0*C_KEV) /* Zero is automatic */

#define STOP_STEP_AUTO_FUDGE_FACTOR (1.0) /* Factor to automatic incident step size. TODO: make this runtime configurable */
#define SIMULATE_WARNING_LIMIT 10 /* Allowed number of non-critical warnings for each run of simulate() */

/* Other constants */
#define DEPTH_TOLERANCE (1.0e-6 * C_TFU)
#define GAMMA_ROUGHNESS_STEPS 21
#define CS_CONC_STEPS 3 /* Minimum 1, odd numbers preferred */
#define CS_STRAGG_HALF_N 3 /* Cross section weighting by straggling, number of steps is this times 2 + 1. Set to zero to disable (aka 1 step). */
#define DUAL_SCATTER_POLAR_STEPS 20
#define DUAL_SCATTER_AZI_STEPS 10
#define FIT_ITERS_MAX 100
#define FIT_XTOL (1e-7)
#define FIT_GTOL (1e-7)
#define FIT_FTOL (1e-7)
#endif // JABS_DEFAULTS_H
