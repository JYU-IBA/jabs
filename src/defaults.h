/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_DEFAULTS_H
#define JABS_DEFAULTS_H
#define COPYRIGHT_STRING "\n    This program is free software; you can redistribute it and/or modify \n    it under the terms of the GNU General Public License as published by\n    the Free Software Foundation; either version 2 of the License, or\n    (at your option) any later version.\n\n    See LICENSE.txt for the full license.\n\n"
#include <jibal_units.h>

#define SCRIPT_FILES_NESTED_MAX 8 /* How deep (number of levels) scripts can be nested, i.e. script file loads a script file loads a script file... */
#define SCRIPT_COMMANDS_NESTED_MAX 8 /* How deep (number of levels) script commands can be nested, i.e. script command has a subcommand, which has a subcommand... */
#define SCRIPT_COMMAND_MERGE_SORT_ARRAY_SIZE 16 /* Commands are sorted alphabetically using a merge sort, the algorithm uses a "Bottom-up" approach with a small fixed size array (can sort up to 2^SCRIPT_COMMAND_MERGE_SORT_ARRAY_SIZE elements) */
#define PROMPT "jabs> "

#define CHANNELS_ABSOLUTE_MIN 256
#define CHANNELS_MAX_DEFAULT 16384
#define CHANNELS_ABSOLUTE_MAX 131072

/* Defaults for new simulations */
#define ENERGY_DEFAULT (2.0*C_MEV)
#define ALPHA_DEFAULT (0.0*C_DEG)
#define DETECTOR_THETA_DEFAULT (165.0*C_DEG)
#define DETECTOR_PHI_DEFAULT (0.0 * C_DEG) /* IBM geometry */
#define DETECTOR_RESOLUTION_DEFAULT (15.0*C_KEV) /* FWHM */
#define DETECTOR_SOLID_DEFAULT (10.0 * C_MSR)
#define DETECTOR_DISTANCE_DEFAULT (100.0 * C_MM)
#define DETECTOR_LENGTH_DEFAULT (1000.0 * C_MM)
#define FLUENCE_DEFAULT (1.0e14)
#define E_MIN_DEFAULT (50.0*C_KEV)
#define ENERGY_SLOPE_DEFAULT (1.0*C_KEV)
#define DES_TABLE_INITIAL_ALLOC (1024)
#define DES_TABLE_MAX_SIZE (256 * DES_TABLE_INITIAL_ALLOC)
#define ENERGY_DERIVATIVE_MIN (0.1) /* dE(detected)/dE(incident) is forced to be higher or equal to this. */
#define ENERGY_DERIVATIVE_MAX (10.0)
#define BRICK_WIDTH_SIGMAS_DEFAULT (C_FWHM) /* Attempt to make bricks width energy difference equal to this many times the detector resolution and straggling sigmas combined (in quadrature (variances)) */
#define INCIDENT_STOP_STEP_DEFAULT (0.0*C_KEV) /* Zero is automatic */
#define INCIDENT_STOP_STEP_SIGMAS_DEFAULT (1.0)
#define INCIDENT_STOP_STEP_MIN_DEFAULT (0.5 * C_KEV) /* Default minimum stopping step. */
#define INCIDENT_STOP_STEP_MAX_DEFAULT (50.0 * C_KEV) /* Default maximum stopping step.  */
#define EXITING_STOP_STEP_DEFAULT (0.0*C_KEV) /* Zero is automatic */
#define EXITING_STOP_STEP_SIGMAS_DEFAULT (5.0)
#define EXITING_STOP_STEP_MIN_DEFAULT (20.0 * C_KEV) /* Default minimum stopping step. */
#define EXITING_STOP_STEP_MAX_DEFAULT (250.0 * C_KEV) /* Default maximum stopping step.  */
#define STOP_STEP_DEPTH_FALLBACK (100.0 * C_TFU) /* If stop step can't be determined using stopping reliable, do a maximum of this */
#define STOP_STEP_MINIMUM_STOPPING (0.1 * C_EV_TFU) /* Stopping below this is suspicious */
#define STOP_STEP_ABSOLUTE_MINIMUM_STEP (1.0 * C_TFU) /* In some cases smaller steps are taken, such as when approaching layer boundaries, but otherwise this should be smallest sanity check step. */
#define BRICKS_DEFAULT (1000)
#define BRICKS_MAX (10000)
#define SIMULATE_WARNING_LIMIT 10 /* Allowed number of non-critical warnings for each run of simulate() */
/* Other constants */
#define E_MIN (10.0 * C_KEV) /* Absolute minimum for everything */
#define E_MAX (1000.0*C_MEV) /* Absolute maximum for everything */
#define DEPTH_TOLERANCE (1.0e-6 * C_TFU)
#define ROUGH_TOLERANCE (0.1 * C_TFU) /* Roughness below this is equivalent to none */
#define CONC_TOLERANCE (1.0e-7)
#define GAMMA_ROUGHNESS_STEPS 21
#define ROUGHNESS_SUBSPECTRA_MAXIMUM 99
#define CS_STRAGG_STEPS 7 /* Number of steps used when weighting cross section by straggling. Odd numbers preferred. */
#define DUAL_SCATTER_POLAR_STEPS 21
#define DUAL_SCATTER_POLAR_SUBSTEPS 9
#define DUAL_SCATTER_AZI_STEPS 15
#define DUAL_SCATTER_INCIDENT_STOP_STEP_FACTOR_DEFAULT (2.0)
#define FIT_ITERS_MAX 100
#define FIT_XTOL (1e-7)
#define FIT_FAST_XTOL_MULTIPLIER (1.0)
#define FIT_CHISQ_TOL (1e-7) /* Relative change in chi squared to stop fitting */
#define FIT_FAST_CHISQ_TOL (1e-4)  /* Relative change in chi squared to stop fitting (fast fitting phase). This can be quite large, as turning on better physics changes the chisq. */
#define SIGMAS_CUTOFF (5.0)
#define SIGMAS_FAST_CUTOFF (3.5)
#define CS_CONC_MAX_INTEGRATION_INTERVALS 100
#define CS_CONC_INTEGRATION_ACCURACY (1e-5)
#define CS_STRAGG_MAX_INTEGRATION_INTERVALS 100
#define CS_STRAGG_INTEGRATION_ACCURACY (1e-5)
#define CS_ENERGY_STEP_MAX_DEFAULT (10.0 * C_KEV)
#define CS_DEPTH_STEP_MAX_DEFAULT (100.0 * C_TFU)
#define CS_STRAGG_STEP_FUDGE_FACTOR_DEFAULT (2.0)
#define SCREENING_TABLE_ELEMENTS (101)
#define REACTION_FILE_ANGLE_TOLERANCE_DEFAULT (1.0 * C_DEG)
#endif // JABS_DEFAULTS_H
