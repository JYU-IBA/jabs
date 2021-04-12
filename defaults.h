#ifndef JABS_DEFAULTS_H
#define JABS_DEFAULTS_H
#include <jibal_units.h>

/* Defaults for new simulations */
#define ENERGY (2.0*C_MEV)
#define ALPHA (0.0*C_DEG)
#define BETA (15.0*C_DEG)
#define THETA (165.0*C_DEG)
#define DETECTOR_RESOLUTION (15.0*C_KEV/C_FWHM)
#define PARTICLES_SR (1.0e12)
#define E_MIN (100.0*C_KEV)
#define ENERGY_SLOPE (1.0*C_KEV)
#define STOP_STEP_INCIDENT (0.0*C_KEV) /* Zero is automatic */
#define STOP_STEP_EXITING (0.0*C_KEV) /* Zero is automatic */

/* Other constants */
#define DEPTH_TOLERANCE (1.0e-6 * C_TFU)
#define ROUGHNESS_STEPS 21
#endif // JABS_DEFAULTS_H
