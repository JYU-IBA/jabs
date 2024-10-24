/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_FIT_H
#define JABS_FIT_H
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlinear.h>

#include <jibal.h>

#include "simulation.h"
#include "simulation_workspace.h"
#include "detector.h"
#include "sample.h"
#include "fit_params.h"
#include "simulation.h"
#include "reaction.h"
#include "sample.h"
#include "spectrum.h"
#include "histogram.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FIT_SUCCESS_CHISQ (2)
#define FIT_SUCCESS_DELTA (1)
#define FIT_SUCCESS (0)
#define FIT_ERROR_NONE (0)
#define FIT_ERROR_GENERIC (-1)
#define FIT_ERROR_MAXITER (-2)
#define FIT_ERROR_NO_PROGRESS (-3)
#define FIT_ERROR_SANITY (-4)
#define FIT_ERROR_WORKSPACE_INITIALIZATION (-5)
#define FIT_ERROR_IMPOSSIBLE (-6)
#define FIT_ERROR_ABORTED (-7)
#define FIT_ERROR_INIT (-8)

#define FIT_PHASE_FAST 1
#define FIT_PHASE_SLOW 2

struct fit_stats {
    int phase;
    size_t n_evals;
    size_t n_evals_iter; /* Number of function evaluations per iteration */
    size_t n_spectra;
    size_t n_spectra_iter; /* Number of spectra (full spectra with roughness etc) actually simulated per iteration call */
    double cputime_cumul;
    double cputime_iter;
    double chisq0;
    double chisq;
    double chisq_dof;
    double rcond;
    double norm;
    size_t iter;
    size_t iter_call;
    int error;
};

typedef struct roi {
    size_t i_det; /* detector sim->det[i_det] */
    size_t low;
    size_t high;
} roi;

typedef struct fit_data_det {
    size_t i_det;
    detector *det; /* Not stored here, same as sim->det[i_det], but we need a non-const pointer. TODO: make a copy */
    const jabs_histogram *exp; /* Not a copy! */
    roi *ranges; /* Same ranges as in fit_data, but only those relevant for "det". Full copies are made. */
    size_t n_ranges;
    size_t n_ch; /* in fit ranges */
    size_t f_offset;
    gsl_vector *f_iter; /* Stored residual vector values of f on first call of iter */
} fit_data_det;

typedef struct jacobian_space {
    simulation sim; /* Copy of simulation, partially shallow, overwritten on every Jacobian calculation, so don't store anything here! */
    detector **det; /* Array of detector pointers, detectors will be cloned here (if needed) */
    fit_variable *var; /* Active fit variable to be perturbed */
    gsl_vector *f_param; /* Residuals vector */
    double delta_inv; /* Inverse of perturbation 1/delta */
    size_t n_spectra_calculated; /* How many spectra (detectors) were actually computed for this var */
} jacobian_space; /* All the stuff needed to compute Jacobian */

typedef struct fit_data {
    fit_data_det *fdd; /* Detector specific stuff */
    result_spectra *spectra; /* all spectra (array of n_det_spectra), updates every iter at the start of iter. */
    size_t n_det_spectra; /* n_det (= n_fdd) when histograms were copied */
    jabs_histogram **exp; /* experimental data to be fitted, array of n_exp elements */
    size_t n_exp; /* same as sim->n_det, but this keeps track on how many spectra we have allocated in exp and ref */ /* TODO: remove */
    jabs_histogram *ref; /* reference spectra, exactly one */
    simulation *sim;
    const jibal *jibal; /* This shouldn't be here, but it is the only place I can think of */
    sample_model *sm;
    struct roi *fit_ranges; /* Array of fit_range, size n_fit_ranges, freed by fit_data_free() */
    fit_params *fit_params; /* Allocated with fit_data_new() and freed by fit_data_free() */
    size_t n_fit_ranges;
    size_t n_iters_max; /* Maximum number of iterations, in each fit phase */
    double xtol; /* Tolerance of step size */
    double chisq_tol; /* Chi squared relative change tolerance */
    double chisq_fast_tol; /* Chi squared relative change tolerance (fast phase) */
    gsl_multifit_nlinear_fdf *fdf;
    size_t dof; /* Degrees of freedom (calculated) */
    struct fit_stats stats; /* Fit statistics, updated as we iterate */
    int phase_start; /* Fit phase to start from (see FIT_PHASE -defines) */
    int phase_stop; /* Inclusive */
    int (*fit_iter_callback)(struct fit_stats stats);
    gsl_vector *f_iter;
    double h_df;
    jacobian_space *jspace;
    gsl_matrix *covar; /* Covariance matrix */
} fit_data;

void fit_data_det_residual_vector_set(const fit_data_det *fdd, const jabs_histogram *histo_sum, gsl_vector *f);
fit_data *fit_data_new(const jibal *jibal, simulation *sim);
int fit_data_jspace_init(fit_data *fit, size_t n_channels_in_fit);
void fit_data_jspace_free(fit_data *fit);
void fit_data_defaults(fit_data *f);
void fit_data_free(fit_data *fit); /* Doesn't free everything in fit, like sm, jibal, ... */
void fit_data_reset(fit_data *fit);
void fit_data_exp_reset(fit_data *fit);
void fit_data_print(const fit_data *fit, jabs_msg_level msg_level);
void fit_data_roi_print(const struct fit_data *fit_data, const struct roi *roi);
jabs_histogram *fit_data_exp(const fit_data *fit, size_t i_det);
jabs_histogram *fit_data_ref(const fit_data *fit_data);
fit_params *fit_params_all(fit_data *fit);
int fit_data_fdd_init(fit_data *fit);
void fit_data_fdd_free(fit_data *fit);
void fit_data_exp_alloc(fit_data *fit);
void fit_data_exp_free(fit_data *fit);
int fit_data_load_exp(struct fit_data *fit, size_t i_det, const char *filename);
int fit_data_set_sample_model(fit_data *fit, sample_model *sm_new); /* Sets sample model of the fit (copies pointer), freeing the old model. Resets reactions. */
jabs_histogram *fit_data_histo_sum(const fit_data *fit, size_t i_det);
void fit_data_spectra_copy_to_spectra_from_ws(result_spectra *s, const detector *det, const jabs_histogram *exp, const sim_workspace *ws); /* Makes deep copies of histograms */
int fit_data_spectra_alloc(fit_data *fit);
void fit_data_spectra_free(fit_data *fit);
int fit_data_add_det(struct fit_data *fit, detector *det);
fit_data_det *fit_data_fdd(const fit_data *fit, size_t i_det);
size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data);
struct fit_stats fit_stats_init(void);
int fit(fit_data *fit);
void fit_correlation_print(const gsl_matrix *covar, jabs_msg_level msg_level);
int fit_uncertainty_spectra(const fit_data *fit, const gsl_matrix *J, const gsl_matrix *covar, const gsl_vector *f, const gsl_vector *w, const char *filename); /* Calculates +/- uncertainty spectra and copies them to fit. If filename is not null, writes (debug) output in it. */
int fit_parameters_set_from_vector(struct fit_data *fit, const gsl_vector *x); /* Updates values in fit params as they are varied by the fit algorithm. */
int fit_function(const gsl_vector *x, void *params, gsl_vector *f);
int fit_determine_active_detectors(fit_data *fit);
int fit_sanity_check(const fit_data *fit);
void fit_iter_stats_update(struct fit_data *params, const gsl_multifit_nlinear_workspace *w);
void fit_iter_stats_print(const struct fit_stats *stats);
void fit_stats_print(const struct fit_stats *stats, jabs_msg_level msg_level);
int fit_data_fit_range_add(struct fit_data *fit_data, const struct roi *range); /* Makes a deep copy */
void fit_data_fit_ranges_free(struct fit_data *fit_data);
int fit_set_roi_from_string(roi *r, const char *str); /* Parses only low and high from "[low:high]". */
double fit_emin(struct fit_data *fit, size_t i_det); /* Returns lowest energy of fit ranges for detector i_det. Detectors, calibrations and fit ranges must be set before calling. */
const char *fit_error_str(int error);
int fit_range_compare(const void *a, const void *b);
#ifdef __cplusplus
}
#endif
#endif // JABS_FIT_H
