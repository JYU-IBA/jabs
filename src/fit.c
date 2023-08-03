/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

    Some parts of this source file under different license, see below!

 */

/*
 * This file is based partially based around the GSL examples for nonlinear least-squares fitting. See e.g.
 * https://www.gnu.org/software/gsl/doc/html/nls.html . Please note that GSL (https://www.gnu.org/software/gsl/) is
 * is distributed under the terms of the GNU General Public License (GPL).
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include "win_compat.h"
#include <string.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "jabs_debug.h"
#include "defaults.h"
#include "generic.h"
#include "jabs.h"
#include "spectrum.h"
#include "message.h"
#include "gsl_inline.h"
#include "histogram.h"
#include "fit.h"
#ifdef _OPENMP
#include <omp.h>
#endif


int fit_detector(const jibal *jibal, const fit_data_det *fdd, const simulation *sim, result_spectra *spectra, gsl_vector *f) {
    detector *det = sim_det(sim, fdd->i_det);
    if(!det) {
        jabs_message(MSG_ERROR, "No detector set.\n");
        return EXIT_FAILURE;
    }
    detector_update(det); /* non-const pointer inside const... not great */
    sim_workspace *ws = sim_workspace_init(jibal, sim, det);
    if(!ws) {
        jabs_message(MSG_ERROR, "Could not initialize workspace of detector %s.\n", det->name);
        return EXIT_FAILURE;
    }
    if(simulate_with_ds(ws)) {
        jabs_message(MSG_ERROR, "Simulation of detector %s spectrum failed!\n", det->name);
        return EXIT_FAILURE;
    }
    fit_data_det_residual_vector_set(fdd, ws->histo_sum, f);
#if 0
    if(iter_call == 1) { /* TODO: give spectra (and f_iter) as pointers, set these if pointer given */
        gsl_vector_memcpy(fdd->f_iter, f); /* Copy residual vector, this can be used on later calls if fdd is not active */
    }
#endif
    if(spectra) {
        result_spectra_free(spectra);
        fit_data_spectra_copy_to_spectra_from_ws(spectra, det, fdd->exp, ws);
    }
    sim_workspace_free(ws);
    return EXIT_SUCCESS;
}

void fit_deriv_prepare_jspaces(const gsl_vector *x, fit_data *fit) {
    jacobian_space *space = fit->jspace;
    for(size_t j = 0; j < fit->fit_params->n_active; j++) { /* Make copies of sim (shallow) with deep copies of detector and sample. This way we can parallelize. */
        struct jacobian_space *spc = &space[j];
        spc->n_spectra_calculated = 0;
        fit_variable *var = spc->var;
        double xj = gsl_vector_get(x, j);
        double delta = fit->h_df * fabs(xj);
        if(delta == 0.0) {
            delta = fit->h_df; /* TODO: this is what GSL does, but it doesn't always work */
        }
        spc->delta_inv = 1.0/delta;
        *(var->value) = (xj + delta) * var->value_orig; /* Perturb */
        spc->sim = *fit->sim; /* Shallow copy! */
        switch(var->type) { /* Deepen the copy based on the variable */
            case FIT_VARIABLE_SAMPLE:
                spc->sim.sample = sample_from_sample_model(fit->sm);
                break;
            case FIT_VARIABLE_DETECTOR:
                spc->sim.det = spc->det;
                spc->sim.det[var->i_det] = detector_clone(fit->sim->det[var->i_det]);
                break;
            default:
                break;
        }
        *(var->value) = (xj) * var->value_orig; /* Unperturb */
    }
}

void fit_deriv_cleanup_jspaces(fit_data *fit) {
    jacobian_space *space = fit->jspace;
    for(size_t j = 0; j <  fit->fit_params->n_active; j++) {
        struct jacobian_space *spc = &space[j];
        fit_variable *var = spc->var;
        fit->stats.n_spectra_iter += spc->n_spectra_calculated;
        switch(var->type) {
            case FIT_VARIABLE_SAMPLE:
                sample_free(spc->sim.sample);
                break;
            case FIT_VARIABLE_DETECTOR:
                detector_free(spc->sim.det[var->i_det]);
                break;
            default:
                break;
        }
    }
}

int fit_deriv_function(const gsl_vector *x, void *params, gsl_matrix *J) {
    struct fit_data *fit = (struct fit_data *) params;
    assert(fit->stats.iter_call >= 1); /* This function relies on data that was computed when iter_call == 1 */
    assert(fit->fdf->p == fit->fit_params->n_active);

    gsl_matrix_set_zero(J); /* We might only set parts of the Jacobian, the rest of the elements of the matrix should be zero. */

    fit_deriv_prepare_jspaces(x, fit);

    double start = jabs_clock();
    volatile int error = FALSE;
    const int n = (int) fit->fit_params->n_active;
    int j;
#pragma omp parallel default(none) shared(fit, J, error, n)
#pragma omp for schedule(dynamic)
    for(j = 0; j < n; j++) {
        //fprintf(stderr, "Thread id %i got %zu.\n", omp_get_thread_num(), j);
        struct jacobian_space *spc = &fit->jspace[j];
        fit_variable *var = spc->var;
        size_t i_start, i_stop;
        if(var->type == FIT_VARIABLE_DETECTOR) {
            i_start = var->i_det;
            i_stop = var->i_det;
        } else { /* All detectors */
            i_start = 0;
            i_stop = fit->sim->n_det - 1;
        }
        DEBUGMSG("Variable %s (i_v = %zu, j = %i) Jacobian. i_det = [%zu, %zu]", var->name, var->i_v, j, i_start, i_stop);
        for(size_t i_det = i_start; i_det <= i_stop; i_det++) {
            fit_data_det *fdd = &fit->fdd[i_det];
            gsl_vector_view v = gsl_vector_subvector(spc->f_param, fdd->f_offset, fdd->n_ch); /* The ROIs of this detector are a subvector of f_param (all channels in fit) */
            DEBUGMSG("Detector %zu (FDD %p), %zu ranges, offset: %zu, len: %zu", i_det, (void *) fdd, fdd->n_ranges, fdd->f_offset, fdd->n_ch);
            if(fdd->n_ranges == 0) {
                continue;
            }
            if(fit_detector(fit->jibal, fdd, &spc->sim, NULL, &v.vector)) {
                error = TRUE;
                break;
            }
            spc->n_spectra_calculated++;
            for(size_t i = 0; i < fdd->n_ch; i++) {
                double fnext = jabs_gsl_vector_get(&v.vector, i);
                double fi = jabs_gsl_vector_get(fit->f_iter, fdd->f_offset + i);
                jabs_gsl_matrix_set(J, i + fdd->f_offset, j, (fnext - fi) * spc->delta_inv);
            }
        }
        if(error) {
            continue;
        }
    }
    fit_deriv_cleanup_jspaces(fit);
    double end = jabs_clock();
    fit->stats.cputime_iter += (end - start);
    return error ? GSL_FAILURE : GSL_SUCCESS;
}

int fit_function(const gsl_vector *x, void *params, gsl_vector *f) {
    struct fit_data *fit = (struct fit_data *) params;
    fit->stats.iter_call++;
    if(fit->stats.iter_call == 1) { /* Clear stored spectra from previous iteration */
        for(size_t i = 0; i < fit->n_det_spectra; i++) {
            result_spectra_free(&fit->spectra[i]);
        }
    }
    if(fit_parameters_set_from_vector(fit, x)) { /* This also sets some helper numbers and arrays */
        return GSL_FAILURE;
    }
    DEBUGMSG("Fit iteration %zu call %zu. Size of vector x: %zu, f: %zu. %zu active fit parameters.",
             fit->stats.iter, fit->stats.iter_call, x->size, f->size, fit->fit_params->n_active);

    if(sim_sanity_check(fit->sim) == GSL_FAILURE) {
        fit->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }
    sample_free(fit->sim->sample); /* TODO: this is only necessary on first call of iter and when sample specific parameters are active */
    fit->sim->sample = sample_from_sample_model(fit->sm);
    if(!fit->sim->sample) {
        fit->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }
    double start = jabs_clock();
    volatile int error = FALSE;
    if(fit->sim->n_det == 1) { /* Skip OpenMP if only one workspace */
        fit_data_det *fdd = &fit->fdd[0];
        result_spectra *spectra = (fit->stats.iter_call == 1 ? &fit->spectra[fdd->i_det] : NULL);
        error = fit_detector(fit->jibal, fdd, fit->sim, spectra, f);
    } else {
        int i;
        const int n = (int) fit->sim->n_det;
#pragma omp parallel default(none) shared(fit, error, f, n)
#pragma omp for
        for(i = 0; i < n; i++) {
#ifdef _OPENMP
            DEBUGMSG("Thread id %i got %i.", omp_get_thread_num(), i);
#endif
            fit_data_det *fdd = &fit->fdd[i];
            if(fdd->n_ranges == 0) { /* This will NOT simulate those detectors that don't participate in fit! */
                continue;
            }
            result_spectra *spectra = (fit->stats.iter_call == 1 ? &fit->spectra[fdd->i_det] : NULL); /* Pass spectra pointer on first call */
            gsl_vector_view f_det = gsl_vector_subvector(f, fdd->f_offset, fdd->n_ch);
            if(fit_detector(fit->jibal, fdd, fit->sim, spectra, &f_det.vector)) {
                error = TRUE;
            }
        }
    }
    double end = jabs_clock();
    DEBUGMSG("Fit iteration %zu call %zu simulation done.", fit->stats.iter, fit->stats.iter_call);
    if(error) {
        jabs_message(MSG_ERROR, "Error in simulating (during fit).");
        return EXIT_FAILURE;
    }
    fit->stats.cputime_iter += (end - start);
    fit->stats.n_evals_iter++;
    fit->stats.n_spectra_iter += fit->sim->n_det;
    if(fit->stats.iter_call == 1) {
        gsl_vector_memcpy(fit->f_iter, f); /* Store residual vector on first call, this acts as baseline for everything we do later */
    }
    //gsl_vector_memcpy(f, fit->f);
    DEBUGSTR("Successful fit function call.");
    return GSL_SUCCESS;
}

int fit_sanity_check(const fit_data *fit) {
    if(sample_model_sanity_check(fit->sm)) {
        return GSL_FAILURE;
    }
    if(sim_sanity_check(fit->sim)) {
        return GSL_FAILURE;
    }
    return GSL_SUCCESS;
}

int fit_parameters_set_from_vector(struct fit_data *fit, const gsl_vector *x) {
    DEBUGMSG("Fit iteration has %zu active parameters from a total of %zu possible.", fit->fit_params->n_active, fit->fit_params->n)
    for(size_t i = 0; i < fit->fit_params->n; i++) {
        fit_variable *var = &fit->fit_params->vars[i];
        if(!var->active) {
            continue;
        }
        if(!isfinite(*(var->value))) {
            DEBUGMSG("Fit iteration %zu, call %zu, fit variable %zu (%s) is not finite.", fit->stats.iter, fit->stats.iter_call, var->i_v, var->name);
            fit->stats.error = FIT_ERROR_SANITY;
            return GSL_FAILURE;
        }
        *(var->value) = gsl_vector_get(x, var->i_v) * var->value_orig;
    }
    return GSL_SUCCESS;
}

void fit_iter_stats_update(struct fit_data *fit_data, const gsl_multifit_nlinear_workspace *w) {
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&fit_data->stats.rcond, w);
    gsl_blas_ddot(f, f, &fit_data->stats.chisq);
    fit_data->stats.norm = gsl_blas_dnrm2(f);
    fit_data->stats.chisq_dof = fit_data->stats.chisq / fit_data->dof;
    fit_data->stats.n_evals += fit_data->stats.n_evals_iter;
    fit_data->stats.n_spectra += fit_data->stats.n_spectra_iter;
    fit_data->stats.cputime_cumul += fit_data->stats.cputime_iter;
}

void fit_iter_stats_print(const struct fit_stats *stats) {
    jabs_message(MSG_INFO, "%4zu | %10.2e | %14.7e | %12.7lf | %7zu | %7zu | %10.3lf | %13.1lf |\n",
                 stats->iter, 1.0 / stats->rcond, stats->norm,
                 stats->chisq_dof, stats->n_evals, stats->n_spectra,
                 stats->cputime_cumul, 1000.0 * stats->cputime_iter / stats->n_spectra_iter);
}

void fit_stats_print(const struct fit_stats *stats, jabs_msg_level msg_level) {
    if(stats->chisq_dof > 0.0) {
        jabs_message(msg_level, "Final chisq/dof = %.7lf\n", stats->chisq_dof);
    }
}

int fit_data_fit_range_add(struct fit_data *fit_data, const struct roi *range) { /* Makes a deep copy */
    if(range->low == 0 && range->high == 0) {
        return EXIT_FAILURE;
    }
    if(range->high < range->low) {
        return EXIT_FAILURE;
    }
    fit_data->n_fit_ranges++;
    fit_data->fit_ranges = realloc(fit_data->fit_ranges, fit_data->n_fit_ranges * sizeof(roi));
    if(!fit_data->fit_ranges) {
        fit_data->n_fit_ranges = 0;
        return EXIT_FAILURE;
    }
    fit_data->fit_ranges[fit_data->n_fit_ranges - 1] = *range;
    qsort(fit_data->fit_ranges, fit_data->n_fit_ranges, sizeof(roi), fit_range_compare); /* Sorting by detector is important for fit residual vector to work */
    return EXIT_SUCCESS;
}

void fit_data_fit_ranges_free(struct fit_data *fit_data) {
    if(!fit_data)
        return;
    free(fit_data->fit_ranges);
    fit_data->fit_ranges = NULL;
    fit_data->n_fit_ranges = 0;
}

void fit_data_det_residual_vector_set(const fit_data_det *fdd, const jabs_histogram *histo_sum, gsl_vector *f) { /* Sets residual vector (f) based on histo_sum */
    size_t i_vec = 0;
    if(fdd->n_ranges == 0) {
        fprintf(stderr, "No ranges!\n");
    }
    for(size_t i_range = 0; i_range < fdd->n_ranges; i_range++) {
        const roi *range = &fdd->ranges[i_range];
        for(size_t i = range->low; i <= range->high; i++) {
            if(i >= histo_sum->n) { /* Outside range of simulated spectrum (simulated is zero) */
                jabs_gsl_vector_set(f, i_vec, fdd->exp->bin[i]);
            } else {
                jabs_gsl_vector_set(f, i_vec, fdd->exp->bin[i] - histo_sum->bin[i]);
            }
            i_vec++;
        }
    }
    DEBUGMSG("Set %zu elements of vector from %zu ranges.", i_vec, fdd->n_ranges);
}
fit_data *fit_data_new(const jibal *jibal, simulation *sim) {
    struct fit_data *fit = calloc(1, sizeof(struct fit_data));
    fit->jibal = jibal;
    fit->sim = sim;
    fit_data_exp_alloc(fit);
    fit_data_defaults(fit);
    return fit;
}

void fit_data_defaults(fit_data *f) {
    f->n_iters_max = FIT_ITERS_MAX;
    f->xtol = FIT_XTOL;
    f->chisq_tol = FIT_CHISQ_TOL;
    f->chisq_fast_tol = FIT_FAST_CHISQ_TOL;
    f->phase_start = FIT_PHASE_FAST;
    f->phase_stop = FIT_PHASE_SLOW;
}

int fit_data_jspace_init(fit_data *fit, size_t n_channels_in_fit) {
    if(!fit || !fit->fit_params || !fit->fit_params->n_active) {
        return EXIT_FAILURE;
    }
    fit->jspace =  calloc(fit->fit_params->n_active, sizeof(jacobian_space)); /* Space for each parameter for Jacobian calculation */
    for(size_t j = 0; j < fit->fit_params->n_active; j++) {
        jacobian_space *spc = &fit->jspace[j];
        spc->var = fit_params_find_active(fit->fit_params, j);
        assert(spc->var);
        if(spc->var->type == FIT_VARIABLE_DETECTOR) {
            spc->det = calloc(fit->sim->n_det, sizeof(detector *)); /* Array for detector pointers. This overrides array in sim, so it needs to be an array, although we only use one element of it. */
        } else {
            spc->det = NULL;
        }
        spc->f_param = gsl_vector_alloc(n_channels_in_fit); /* Residuals for this parameter */
    }
    return EXIT_SUCCESS;
}

void fit_data_jspace_free(fit_data *fit) {
    if(!fit) {
        return;
    }
    for(size_t j = 0; j < fit->fit_params->n_active; j++) {
        jacobian_space *spc = &fit->jspace[j];
        free(spc->det);
        gsl_vector_free(spc->f_param);
    }
    free(fit->jspace);
    fit->jspace = NULL;
}

void fit_data_free(fit_data *fit) {
    if(!fit)
        return;
    fit_data_reset(fit);
    fit_data_exp_free(fit);
    fit_data_fdd_free(fit);
    for(size_t i = 0; i < fit->n_det_spectra; i++) {
        result_spectra_free(&fit->spectra[i]);
    }
    free(fit->spectra);
    free(fit);
}

void fit_data_reset(fit_data *fit) {
    if(!fit) {
        return;
    }
    fit_data_fit_ranges_free(fit);
    fit_params_free(fit->fit_params);
    fit->fit_params = NULL;
    fit_data_spectra_free(fit);
}

void fit_data_exp_reset(fit_data *fit) {
    for(size_t i = 0; i < fit->n_exp; i++) {
        jabs_histogram_free(fit->exp[i]);
        fit->exp[i] = NULL;
    }
}

void fit_data_roi_print(const struct fit_data *fit_data, const struct roi *roi) {
    if(!fit_data) {
        return;
    }
    jabs_histogram *histo_sim = fit_data_histo_sum(fit_data, roi->i_det);
    jabs_histogram *histo_exp = fit_data_exp(fit_data, roi->i_det);
    jabs_histogram *histo_ref = fit_data->ref;
    size_t n_sim = jabs_histogram_channels_in_range(histo_sim, roi->low, roi->high);
    size_t n_exp = jabs_histogram_channels_in_range(histo_exp, roi->low, roi->high);
    size_t n_ref = jabs_histogram_channels_in_range(histo_ref, roi->low, roi->high);
    double sim_cts = jabs_histogram_roi(histo_sim, roi->low, roi->high);
    double exp_cts = jabs_histogram_roi(histo_exp, roi->low, roi->high);
    double ref_cts = jabs_histogram_roi(histo_ref, roi->low, roi->high);
    const detector *det = sim_det(fit_data->sim, roi->i_det);

    jabs_message(MSG_INFO, "          low = %12zu\n", roi->low);
    jabs_message(MSG_INFO, "         high = %12zu\n", roi->high);
    if(det) {
        jabs_message(MSG_INFO, "        E_low = %12.3lf keV (low energy edge of bin)\n", detector_calibrated(det, JIBAL_ANY_Z, roi->low) / C_KEV);
        jabs_message(MSG_INFO, "       E_high = %12.3lf keV (high energy edge of bin)\n", detector_calibrated(det, JIBAL_ANY_Z, roi->high + 1) / C_KEV);
    }
    if(histo_sim) {
        jabs_message(MSG_INFO, "        n_sim = %12zu (number of channels)\n", n_sim);
        jabs_message(MSG_INFO, "          sim = %12.8g (number of counts)\n", sim_cts);
    }
    if(histo_exp) {
        jabs_message(MSG_INFO, "        n_exp = %12zu (number of channels)\n", n_exp);
        jabs_message(MSG_INFO, "          exp = %12.8g (number of counts)\n", exp_cts);
        jabs_message(MSG_INFO, "    sqrt(exp) = %12.5lf\n", sqrt(exp_cts));
    }
    if(histo_ref) {
        jabs_message(MSG_INFO, "        n_ref = %12zu (number of channels)\n", n_ref);
        jabs_message(MSG_INFO, "          ref = %12.8g (number of counts)\n", ref_cts);
        jabs_message(MSG_INFO, "    sqrt(ref) = %12.5lf\n", sqrt(ref_cts));
    }
    if(histo_sim && histo_exp && sim_cts > 0 && exp_cts > 0) {
        jabs_message(MSG_INFO, "      exp-sim = %12.8g\n", exp_cts - sim_cts);
        jabs_message(MSG_INFO, "      sim/exp = %12.5lf\n", sim_cts / exp_cts);
        jabs_message(MSG_INFO, "      exp/sim = %12.5lf\n", exp_cts / sim_cts);
        jabs_message(MSG_INFO, "  1/sqrt(exp) = %12.5lf%%\n", 100.0 / sqrt(exp_cts));
        jabs_message(MSG_INFO, "(exp-sim)/exp = %12.5lf%%\n", 100.0 * (exp_cts - sim_cts) / exp_cts);
    }
    if(histo_sim && histo_ref && sim_cts > 0 && ref_cts > 0) {
        jabs_message(MSG_INFO, "      ref-sim = %12g\n", ref_cts - sim_cts);
        jabs_message(MSG_INFO, "      sim/ref = %12.5lf\n", sim_cts / ref_cts);
        jabs_message(MSG_INFO, "      ref/sim = %12.5lf\n", ref_cts / sim_cts);
        jabs_message(MSG_INFO, "  1/sqrt(ref) = %12.5lf%%\n", 100.0 / sqrt(ref_cts));
        jabs_message(MSG_INFO, "(ref-sim)/ref = %12.5lf%%\n", 100.0 * (ref_cts - sim_cts) / ref_cts);
    }
}

jabs_histogram *fit_data_exp(const fit_data *fit, size_t i_det) {
    if(!fit || !fit->exp)
        return NULL;
    if(i_det >= fit->sim->n_det)
        return NULL;
    assert(fit->n_exp == fit->sim->n_det);
    return fit->exp[i_det];
}

fit_params *fit_params_all(fit_data *fit) {
    simulation *sim = fit->sim;
    sample_model *sm = fit->sm;
    if(!sim)
        return NULL;
    size_t param_name_max_len = 256; /* Laziness. We use a fixed size temporary string. snprintf is used, so no overflows should occur, but very long names may be truncated. */
    char *param_name = malloc(sizeof(char) * param_name_max_len);
    fit_params *params = fit_params_new();
    fit_params_add_parameter(params, FIT_VARIABLE_BEAM, &sim->fluence, "fluence", "", 1.0, 0); /* This must be the first parameter always, as there is a speedup in the fit routine */
    fit_params_add_parameter(params, FIT_VARIABLE_GEOMETRY, &sim->sample_theta, "alpha", "deg", C_DEG, 0);
    fit_params_add_parameter(params, FIT_VARIABLE_BEAM, &sim->beam_E, "energy", "keV", C_KEV, 0);

    if(sm) {
        for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) {
            sample_range *r = &(sm->ranges[i_range]);
            size_t range_index = i_range + 1; /* Human readable indexing */
            if(r->x < 0.0) { /* Don't fit anything if thickness is zero (negative shouldn't be possible) */
                continue;
            }
            if(r->rough.model != ROUGHNESS_FILE) { /* Layer thickness is not a parameter with arbitrary roughness from a file */
                snprintf(param_name, param_name_max_len, "thick%zu", range_index);
                fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->x), param_name, "tfu", C_TFU, 0);
            }

            snprintf(param_name, param_name_max_len, "yield%zu", range_index);
            fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->yield), param_name, "", 1.0, 0);

            snprintf(param_name, param_name_max_len, "yield_slope%zu", range_index);
            fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->yield_slope), param_name, "", 1.0, 0);

            if(i_range == sm->n_ranges - 1) { /* Last range, add channeling "aliases" (=yield corrections) */
                snprintf(param_name, param_name_max_len, "channeling");
                fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->yield), param_name, "", 1.0, 0);

                snprintf(param_name, param_name_max_len, "channeling_slope");
                fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->yield_slope), param_name, "", 1.0, 0);
            }

            snprintf(param_name, param_name_max_len, "bragg%zu", range_index);
            fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->bragg), param_name, "", 1.0, 0);

            snprintf(param_name, param_name_max_len, "stragg%zu", range_index);
            fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->stragg), param_name, "", 1.0, 0);

            if(r->rough.model == ROUGHNESS_GAMMA && r->rough.x > 0.0) {
                snprintf(param_name, param_name_max_len, "rough%zu", range_index);
                fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, &(r->rough.x), param_name, "tfu", C_TFU, 0);
            }
            for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
                if(*sample_model_conc_bin(sm, i_range, i_mat) < CONC_TOLERANCE) /* Don't add fit variables for negative, zero or very low concentrations. */
                    continue;
                snprintf(param_name, param_name_max_len, "conc%zu_%s", range_index, sm->materials[i_mat]->name);
                fit_params_add_parameter(params, FIT_VARIABLE_SAMPLE, sample_model_conc_bin(sm, i_range, i_mat), param_name, "%", C_PERCENT, 0);
            }
        }
    }

    for(size_t i_det = 0; i_det < sim->n_det; i_det++) {
        detector *det = sim_det(sim, i_det);
        char *det_name = NULL;
        if(asprintf(&det_name, "det%zu_", i_det + 1) < 0) {
            return NULL;
        }
        snprintf(param_name, param_name_max_len, "%ssolid", det_name);
        fit_params_add_parameter(params, FIT_VARIABLE_DETECTOR, &det->solid, param_name, "msr", C_MSR, i_det);

        snprintf(param_name, param_name_max_len, "%stheta", det_name);
        fit_params_add_parameter(params, FIT_VARIABLE_DETECTOR, &det->theta, param_name, "deg", C_DEG, i_det);

        snprintf(param_name, param_name_max_len, "%sphi", det_name);
        fit_params_add_parameter(params, FIT_VARIABLE_DETECTOR, &det->phi, param_name, "deg", C_DEG, i_det);

        for(int Z = JIBAL_ANY_Z; Z <= det->cal_Z_max; Z++) {
            calibration *c = detector_get_calibration(det, Z);
            if(Z != JIBAL_ANY_Z && c == det->calibration) /* No Z-specific calibration */
                continue;
            assert(c);
            size_t n = calibration_get_number_of_params(c);
            for(int i = CALIBRATION_PARAM_RESOLUTION; i < (int) n; i++) {
                char *calib_param_name = calibration_param_name(c->type, i);
                snprintf(param_name, param_name_max_len, "%scalib%s%s_%s",
                         det_name,
                         (Z == JIBAL_ANY_Z) ? "" : "_",
                         (Z == JIBAL_ANY_Z) ? "" : jibal_element_name(fit->jibal->elements, Z),
                         calib_param_name);
                free(calib_param_name);
                fit_params_add_parameter(params, FIT_VARIABLE_DETECTOR, calibration_get_param_ref(c, i), param_name, "keV", C_KEV, i_det);
            }
        }
        free(det_name);
    }
    free(param_name);
    return params;
}

int fit_data_fdd_init(fit_data *fit) {
    if(!fit) {
        return EXIT_FAILURE;
    }
    if(fit->fdd) {
        free(fit->fdd);
    }
    fit->fdd = calloc(fit->sim->n_det, sizeof(fit_data_det));

    for(size_t i_roi = 0; i_roi < fit->n_fit_ranges; i_roi++) { /* Loop over rois, fdds, count rois and n of channels */
        const roi *roi = &fit->fit_ranges[i_roi];
        for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
            fit_data_det *fdd = &fit->fdd[i_det];
            if(i_det != roi->i_det) { /* This ROI is not for this fdd */
                continue;
            }
            fdd->n_ranges++;
            fdd->n_ch += (roi->high - roi->low + 1);
        }
    }

    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        fit_data_det *fdd = &fit->fdd[i_det];
        fdd->exp = fit_data_exp(fit, i_det);
        fdd->i_det = i_det;
        fdd->det = sim_det(fit->sim, i_det);
        fdd->f_iter = gsl_vector_alloc(fdd->n_ch);
        fdd->ranges = calloc(fdd->n_ranges, sizeof(roi));
        if(fdd->n_ranges == 0) {
            DEBUGMSG("Detector %zu no ranges, length: %zu", i_det + 1, fdd->n_ch);
        }
        size_t i = 0; /* Index of fdd roi */
        size_t i_vec = 0; /* fit residual vector index at beginning of roi */
        for(size_t i_roi = 0; i_roi < fit->n_fit_ranges; i_roi++) {
            const roi *roi = &fit->fit_ranges[i_roi];
            size_t l = (roi->high - roi->low + 1);
            if(i_det != roi->i_det) {
                i_vec += l;
                continue;
            }
            if(i == 0) {
                fdd->f_offset = i_vec;
                DEBUGMSG("FDD %zu, %zu ROIs, First ROI is %zu/%zu of all, subvector offset %zu, length %zu)", i_det, fdd->n_ranges, i_roi, fit->n_fit_ranges, fdd->f_offset, fdd->n_ch);
            }
            fdd->ranges[i] = *roi;
            i++;
            i_vec += l;
        }
    }
    return EXIT_SUCCESS;
}

void fit_data_fdd_free(fit_data *fit) {
    if(!fit || !fit->fdd) {
        return;
    }
    for(size_t i_fdd = 0; i_fdd < fit->sim->n_det; i_fdd++) {
        fit_data_det *fdd = &fit->fdd[i_fdd];
        if(!fdd) {
            continue;
        }
        gsl_vector_free(fdd->f_iter);
        free(fdd->ranges);
    }
    free(fit->fdd);
    fit->fdd = NULL;
}

void fit_data_exp_alloc(fit_data *fit) {
    if(!fit) {
        return;
    }
    size_t n_alloc = fit->sim->n_det;
    if(fit->n_exp == 0) { /* This could be handled by realloc too, but this is an easy way to get null pointers as contents. */
        fit->exp = calloc(n_alloc, sizeof(jabs_histogram *));
    } else {
        fit->exp = realloc(fit->exp, sizeof(jabs_histogram *) * n_alloc);
        for(size_t i = fit->n_exp; i < n_alloc; i++) { /* Reset newly allocated space */
            fit->exp[i] = NULL;
        }
    }
    fit->n_exp = n_alloc;
}

void fit_data_exp_free(fit_data *fit) {
    if(!fit->exp)
        return;
    fit_data_exp_reset(fit);
    free(fit->exp);
    fit->exp = NULL;
    fit->n_exp = 0;
}

int fit_data_load_exp(struct fit_data *fit, size_t i_det, const char *filename) {
    jabs_histogram *h = spectrum_read_detector(filename, sim_det(fit->sim, i_det));
    if(!h) {
        jabs_message(MSG_ERROR, "Reading spectrum from file \"%s\" was not successful.\n", filename);
        return EXIT_FAILURE;
    }
    if(fit->exp[i_det]) {
        jabs_histogram_free(fit->exp[i_det]);
    }
    fit->exp[i_det] = h;
    return EXIT_SUCCESS;
}

jabs_histogram *fit_data_histo_sum(const fit_data *fit, size_t i_det) {
    if(!fit)
        return NULL;
    if(i_det >= fit->n_det_spectra)
        return NULL;
    if(RESULT_SPECTRA_SIMULATED > fit->spectra[i_det].n_spectra)
        return NULL;
    return result_spectra_simulated_histo(&fit->spectra[i_det]);
}

void fit_data_spectra_copy_to_spectra_from_ws(result_spectra *spectra, const detector *det, const jabs_histogram *exp, const sim_workspace *ws) {
    spectra->n_spectra = ws->n_reactions + RESULT_SPECTRA_N_FIXED; /* Simulated, experimental + reaction spectra */
    spectra->s = calloc(spectra->n_spectra, sizeof(result_spectrum));
    result_spectrum_set(&spectra->s[RESULT_SPECTRA_SIMULATED], ws->histo_sum, "Simulated", NULL, REACTION_NONE);
    result_spectrum_set(&spectra->s[RESULT_SPECTRA_EXPERIMENTAL], exp, "Experimental", NULL, REACTION_NONE);
    calibration_apply_to_histogram(det->calibration, result_spectra_experimental_histo(spectra));
    DEBUGMSG("Copying %zu spectra (%zu reactions) from ws to spectra structure.", spectra->n_spectra, ws->n_reactions);
    for(size_t i = 0; i < ws->n_reactions; i++) {
        const reaction *r = ws->reactions[i]->r;
        result_spectrum *s = &spectra->s[RESULT_SPECTRA_REACTION_SPECTRUM(i)];
        s->histo = jabs_histogram_clone(ws->reactions[i]->histo);
        if(asprintf(&s->name, "%s (%s)", r->target->name, reaction_type_to_string(r->type)) < 0) {
            s->name = NULL;
        }
        s->target_isotope = r->target;
        s->type = r->type;
    }
}

int fit_data_spectra_alloc(fit_data *fit) {
    fit_data_spectra_free(fit);
    fit->spectra = calloc(fit->sim->n_det, sizeof(result_spectra));
    if(!fit->spectra) {
        fit->n_det_spectra = 0;
        return EXIT_FAILURE;
    }
    fit->n_det_spectra = fit->sim->n_det;
    return EXIT_SUCCESS;
}

void fit_data_spectra_free(fit_data *fit) {
    if(!fit) {
        return;
    }
    for(size_t i = 0; i < fit->n_det_spectra; i++) {
        result_spectra_free(&fit->spectra[i]);
    }
    free(fit->spectra);
    fit->spectra = NULL;
    fit->n_det_spectra = 0;
}

int fit_data_add_det(struct fit_data *fit, detector *det) {
    if(!fit || !det)
        return EXIT_FAILURE;
    if(sim_det_add(fit->sim, det)) {
        return EXIT_FAILURE;
    }
    fit_data_exp_alloc(fit); /* Number of detectors in sim changed, let these guys know it too */
    return EXIT_SUCCESS;
}

fit_data_det *fit_data_fdd(const fit_data *fit, size_t i_det) {
    if(!fit|| !fit->fdd)
        return NULL;
    if(i_det >= fit->sim->n_det)
        return NULL;
    return &fit->fdd[i_det];
}


size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data) {
    size_t sum = 0;
    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *r = &fit_data->fit_ranges[i];
        detector *det = sim_det(fit_data->sim, r->i_det);
        if(!det) {
            continue;
        }
        if(r->high >= det->channels) { /* Limited by detector. TODO: ideally ROIs should not be outside detector channels.. */
            sum += (det->channels - r->low);
        } else {
            sum += (r->high - r->low) + 1;
        }
    }
    return sum;
}

struct fit_stats fit_stats_init() {
    struct fit_stats s;
    s.n_evals = 0;
    s.n_evals_iter = 0;
    s.n_spectra = 0;
    s.n_spectra_iter = 0;
    s.cputime_cumul = 0.0;
    s.cputime_iter = 0.0;
    s.chisq0 = 0.0;
    s.chisq = 0.0;
    s.chisq_dof = 0.0;
    s.rcond = 0.0;
    s.iter = 0;
    s.iter_call = 0;
    s.error = FIT_ERROR_NONE;
    return s;
}

void fit_data_print(const fit_data *fit, jabs_msg_level msg_level) {
    if(!fit) {
        return;
    }
    if(fit->n_fit_ranges == 0) {
        jabs_message(MSG_ERROR, "No fit ranges.\n");
        return;
    }
    jabs_message(msg_level, "Fit ROI # | detector name |  low ch | high ch | exp counts | sim counts\n");
    for(size_t i = 0; i < fit->n_fit_ranges; i++) {
        roi *range = &fit->fit_ranges[i];
        double exp = jabs_histogram_roi(fit_data_exp(fit, range->i_det), range->low, range->high);
        double sim = jabs_histogram_roi(fit_data_histo_sum(fit, range->i_det), range->low, range->high);
        jabs_message(msg_level, " %8zu | %13s | %7zu | %7zu | %10g | %10g\n",
                     i + 1, detector_name(sim_det(fit->sim, range->i_det)), range->low, range->high,
                     exp, sim);

    }

    jabs_message(msg_level, "\nFit %zu ranges with %zu channels total.\n", fit->n_fit_ranges, fit_data_ranges_calculate_number_of_channels(fit));
}

int jabs_test_delta(const gsl_vector *dx, const gsl_vector *x, double epsabs, double epsrel) { /* test_delta() copied from GSL convergence.c and modified */
    int ok = TRUE;
    DEBUGVERBOSEMSG("Test deltas to x->size=%zu\n", x->size);
    for(size_t i = 0; i < x->size; i++) {
        double xi = gsl_vector_get(x, i);
        double dxi = gsl_vector_get(dx, i);
        double tolerance = epsabs + epsrel * fabs(xi);
        double rel = fabs(dxi) / tolerance; /* "How many times over the acceptable tolerance are we */
        DEBUGVERBOSEMSG("Test delta: i %zu, xi %g, dxi %g, tolerance %g, rel %g\n", i, xi, dxi, tolerance, rel);
        if(rel >= 1.0) {
            DEBUGVERBOSEMSG("Fails because %g > 1.0.\n", rel);
            ok = FALSE;
            break;
        }
    }
    if(ok)
        return GSL_SUCCESS;
    return GSL_CONTINUE;
}

int multifit_nlinear_print_jacobian(const gsl_multifit_nlinear_workspace *w, const char *filename) {
    gsl_matrix *J = gsl_multifit_nlinear_jac(w);
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f) {
        return EXIT_FAILURE;
    }
    for(size_t i = 0; i < J->size1; i++) {
        for(size_t j = 0; j < J->size2; j++) {
            double val = gsl_matrix_get(J, i, j);
            fprintf(f, " %12e", val);
        }
        fprintf(f, "\n");
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

int jabs_gsl_multifit_nlinear_driver(const size_t maxiter, const double xtol, const double chisq_tol, struct fit_data *fit_data, gsl_multifit_nlinear_workspace *w) {
    int status = 0;
    size_t iter;
    double chisq_dof_old;
    jabs_message(MSG_INFO, "iter |   cond(J)  |     |f(x)|     |   chisq/dof  |   evals | spectra | time cumul | time/spectrum |\n");
    jabs_message(MSG_INFO, "     |            |                |              |  cumul. |  cumul. |          s |            ms |\n");
    for(iter = 0; iter <= maxiter; iter++) {
        fit_data->stats.iter_call = 0;
        fit_data->stats.iter = iter;
        if(iter) {
            chisq_dof_old = fit_data->stats.chisq_dof;
            fit_data->stats.cputime_iter = 0.0;
            fit_data->stats.n_evals_iter = 0;
            fit_data->stats.n_spectra_iter = 0;
            status = gsl_multifit_nlinear_iterate(w);
            DEBUGMSG("Iteration status %i (%s)", status, gsl_strerror(status));
        }
        if(fit_data->stats.error) {
            return fit_data->stats.error;
        }
        if(status == GSL_ENOPROG && iter == 1) {
            return FIT_ERROR_NO_PROGRESS;
        }
        fit_iter_stats_update(fit_data, w);
        fit_iter_stats_print(&fit_data->stats);
#ifdef DEBUG
        if(fit_data->stats.phase > FIT_PHASE_FAST) {
            char *jacobian_filename;
            asprintf(&jacobian_filename, "jacobian_iter%zu.dat", fit_data->stats.iter);
            if(jacobian_filename) {
                multifit_nlinear_print_jacobian(w, jacobian_filename);
                free(jacobian_filename);
            }
        }
#endif
        if(fit_data->fit_iter_callback) {
            if(fit_data->fit_iter_callback(fit_data->stats)) {
                return FIT_ERROR_ABORTED;
            }
        }
        if(iter == 0)
            continue;
        /* test for convergence */
        status = jabs_test_delta(w->dx, w->x, xtol * xtol, xtol);
        if(status == GSL_SUCCESS) {
            return FIT_SUCCESS_DELTA;
        }
        double chisq_change = 1.0 - fit_data->stats.chisq_dof / chisq_dof_old;
        if(fit_data->stats.chisq_dof > chisq_dof_old) {
            jabs_message(MSG_WARNING, "Chisq increased, this probably shouldn't happen.\n");
        }
        if(chisq_change < chisq_tol) {
            return FIT_SUCCESS_CHISQ;
        }
    }
    return FIT_ERROR_MAXITER;
}

void fit_report_results(const fit_data *fit, const gsl_multifit_nlinear_workspace *w, const gsl_multifit_nlinear_fdf *fdf) {
    jabs_message(MSG_INFO, "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    jabs_message(MSG_INFO, "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    jabs_message(MSG_INFO, "function evaluations: %zu\n", fit->stats.n_evals);
#ifdef DEBUG
    jabs_message(MSG_INFO, "function evaluations (GSL): %zu\n", fdf->nevalf);
#endif
    jabs_message(MSG_INFO, "Jacobian evaluations: %zu\n", fdf->nevaldf);
    jabs_message(MSG_INFO, "number of spectra simulated: %zu\n", fit->stats.n_spectra);
    jabs_message(MSG_INFO, "reason for stopping: %s\n", fit_error_str(fit->stats.error));
    jabs_message(MSG_INFO, "initial |f(x)| = %f\n", sqrt(fit->stats.chisq0));
    jabs_message(MSG_INFO, "final   |f(x)| = %f\n", sqrt(fit->stats.chisq));
}


void fit_covar_print(const gsl_matrix *covar, jabs_msg_level msg_level) {
    jabs_message(msg_level, "\nCorrelation coefficients matrix:\n       | ");
    for(size_t i = 0; i < covar->size1; i++) {
        jabs_message(msg_level, " %4zu  ", i + 1);
    }
    jabs_message(msg_level, "\n");
    for(size_t i = 0; i < covar->size1; i++) {
        jabs_message(msg_level, "%6zu | ", i + 1);
        for(size_t j = 0; j <= i && j < covar->size2; j++) {
            jabs_message(msg_level, " %6.3f", gsl_matrix_get(covar, i, j) / sqrt(gsl_matrix_get(covar, i, i) * gsl_matrix_get(covar, j, j)));
        }
        jabs_message(msg_level, "\n");
    }
}

int fit_uncertainty_print(const fit_data *fit, const gsl_matrix *J, const gsl_matrix *covar, const gsl_vector *f, const gsl_vector *w, const char *filename) {
    FILE *f_errvec = fopen(filename, "w");
    if(!f_errvec) {
        return EXIT_FAILURE;
    }
    gsl_matrix *JC = gsl_matrix_alloc(J->size1, covar->size2); /* Product of J*C */
    gsl_matrix *JCJT = gsl_matrix_alloc(J->size1, J->size1); /* Product J*C*JT */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, J, covar, 0.0, JC);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, JC, J, 0.0, JCJT);
    gsl_vector_view err_vec = gsl_matrix_diagonal(JCJT);
    size_t i_vec = 0;
    for(size_t i_roi = 0; i_roi < fit->n_fit_ranges; i_roi++) {
        const roi *roi = &fit->fit_ranges[i_roi];
        const result_spectra *s = &(fit->spectra[roi->i_det]);
        for(size_t ch = roi->low; ch <= roi->high; ch++) {
            double sim_counts = jabs_histogram_get(result_spectra_simulated_histo(s), ch);
            double exp_counts = jabs_histogram_get(result_spectra_experimental_histo(s), ch);
            double residuals_weighted = gsl_vector_get(f, i_vec);  /* Residuals (with weights sqrt(W)) */
            double residuals_unweighted = exp_counts - sim_counts;
            double weight = gsl_vector_get(w, i_vec); /* Fit weight (1/variance, i.e. 1 / N_exp in bin) */
            double err = gsl_vector_get(&err_vec.vector, i_vec); /* Something(J x C x J^T diagonal), should be standard error of the fit, i.e. how much do the (statistical) errors in parameters propagated to this particular bin give us uncertainty. This is probably variance. */
            double error_fit_final = 2.0 * sqrt(err) * sqrt(sim_counts); /* TODO: sqrt(sim_counts)? */
            double error_prediction_final = 2.0 * sqrt(err + weight) * sqrt(sim_counts);
            double error_positive = sim_counts + error_prediction_final;
            double error_negative = sim_counts - error_prediction_final;
            if(error_negative < 0.0) {
                error_negative = 0.0;
            }
            fprintf(f_errvec, "%3zu %3zu %4zu %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
                    i_vec, i_roi, ch,
                    sim_counts, exp_counts,
                    residuals_unweighted, residuals_weighted,
                    weight,
                    err, error_fit_final,
                    error_prediction_final,
                    error_negative, error_positive);
            i_vec++;
        }
    }
    fclose(f_errvec);
    gsl_matrix_free(JCJT);
    gsl_matrix_free(JC);
    return EXIT_SUCCESS;
}

int fit(fit_data *fit) {
    struct fit_params *fit_params = fit->fit_params;
    if(!fit_params || fit_params->n_active == 0) {
        jabs_message(MSG_ERROR, "No parameters to fit.\n");
        return EXIT_FAILURE;
    }
    if(!fit->exp) {
        jabs_message(MSG_ERROR, "No experimental spectrum to fit.\n");
        return EXIT_FAILURE;
    }
    if(!fit->n_fit_ranges) {
        jabs_message(MSG_ERROR, "No fit range(s) given, can not fit.\n");
        return EXIT_FAILURE;
    }
    for(size_t i_range = 0; i_range < fit->n_fit_ranges; i_range++) {
        roi *range = &fit->fit_ranges[i_range];
        assert(range);
        detector *det = sim_det(fit->sim, range->i_det);
        jabs_histogram *exp = fit_data_exp(fit, range->i_det);
        if(!det) {
            jabs_message(MSG_ERROR, "Detector %zu (fit range %zu) does not exist.\n", range->i_det + 1, i_range + 1);
            return EXIT_FAILURE;
        }
        if(!exp) {
            jabs_message(MSG_ERROR, "Experimental spectrum for detector %zu (fit range %zu) does not exist.\n", range->i_det + 1, i_range + 1);
            return EXIT_FAILURE;
        }
    }
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    fdf_params.trs = gsl_multifit_nlinear_trs_lm;
    fdf_params.solver = gsl_multifit_nlinear_solver_qr;
    fdf_params.h_df = sqrt(GSL_DBL_EPSILON) * 10.0;
    gsl_multifit_nlinear_fdf *fdf = malloc(sizeof(gsl_multifit_nlinear_fdf));
    fdf->params = fit;
    fit->fdf = fdf;
    fdf->f = &fit_function;
    fdf->df = &fit_deriv_function; /*  &fit_deriv_function to use our own finite difference to determine Jacobian, NULL to use GSL. */
    if(fdf->df) {
        DEBUGSTR("Using our own function to calculate Jacobian.");
    }
    fdf->fvv = NULL; /* No geodesic acceleration */
    fdf->n = fit_data_ranges_calculate_number_of_channels(fit);
    fdf->p = fit_params->n_active;
    if(fdf->n < fdf->p) {
        jabs_message(MSG_ERROR, "Not enough data (%zu points) for given number of free parameters (%zu)\n", fdf->n, fdf->p);
        free(fdf);
        return -1;
    }
    fit->dof = fdf->n - fdf->p;
    fit->h_df = fdf_params.h_df;
    jabs_message(MSG_INFO, "%zu channels and %zu parameters in fit, %zu degrees of %s\n", fdf->n, fdf->p, fit->dof, fit->dof < 10000 ? "freedom.":"FREEDOOOOM!!!");
    if(fit_data_jspace_init(fit, fdf->n)) {
        jabs_message(MSG_ERROR, "Could not initialize Jacobian evaluation workspace for %zu parameters.\n", fit->fit_params->n_active);
    }
    gsl_vector *f;
    gsl_matrix *J;
    int status;
    double *weights = calloc(fdf->n, sizeof(double));
    if(!weights)
        return 1;
    size_t i_w = 0;
    for(size_t i_range = 0; i_range < fit->n_fit_ranges; i_range++) {
        roi *range = &fit->fit_ranges[i_range];
        assert(range);
        jabs_histogram *exp = fit_data_exp(fit, range->i_det);
        for(size_t i = range->low; i <= range->high && i < exp->n; i++) {
            if(exp->bin[i] > 1.0) {
                weights[i_w] = 1.0 / (exp->bin[i]);
            } else {
                weights[i_w] = 1.0; /* TODO: ?*/
            }
            i_w++;
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Set %zu weights.\n", i_w);
#endif
    if(fit_data_fdd_init(fit)) {
        return EXIT_FAILURE;
    }
    assert(i_w == fdf->n);
    fit->f_iter = gsl_vector_alloc(fdf->n);
    gsl_matrix *covar = gsl_matrix_alloc(fit_params->n_active, fit_params->n_active);
    gsl_vector *x = gsl_vector_alloc(fit_params->n_active);

    for(size_t i = 0; i < fit_params->n; i++) { /* Update all (including inactives) and set the (active) variable for each Jacobian parameter workspace */
        fit_variable *var = &(fit_params->vars[i]);
        var->err = 0.0;
        var->value_orig = *(var->value);
    }

    gsl_vector_view wts = gsl_vector_view_array(weights, i_w);

    w = gsl_multifit_nlinear_alloc(T, &fdf_params, fdf->n, fdf->p);
    sim_calc_params p_orig = *fit->sim->params; /* Store original values (will be used in final stage of fitting) */
    for(int phase = fit->phase_start; phase <= fit->phase_stop; phase++) { /* Phase 1 is "fast", phase 2 normal. */
        assert(phase >= FIT_PHASE_FAST && phase <= FIT_PHASE_SLOW);
        double xtol = fit->xtol;
        double chisq_tol = fit->chisq_tol;
        /* initialize solver with starting point and weights */
        fit->stats = fit_stats_init();
        fit->stats.phase = phase;
        if(fit->fit_iter_callback) { /* First call to callback quickly (before most initialization) */
            if(fit->fit_iter_callback(fit->stats)) {
                fit->stats.error = FIT_ERROR_ABORTED;
                break;
            }
        }
        if(phase == FIT_PHASE_FAST) {
            sim_calc_params_defaults_fast(fit->sim->params); /* Set current parameters to be faster in phase 0. */
            xtol *= FIT_FAST_XTOL_MULTIPLIER;
            chisq_tol = fit->chisq_fast_tol;
        } else {
            sim_calc_params_copy(&p_orig, fit->sim->params);
        }
        sim_calc_params_update(fit->sim->params);
        for(size_t i = 0; i < fit_params->n; i++) { /* Set active variables to vector */
            fit_variable *var = &(fit_params->vars[i]);
            if(var->active) {
                gsl_vector_set(x, var->i_v, *(var->value)/var->value_orig); /* We'll pass normalized values to GSL, so we are actually starting fit with always with vector full of 1.0. Next phase starts where previous ends. */
            }
        }
        jabs_message(MSG_INFO, "\nInitializing fit phase %i. Xtol = %e, chisq_tol %e\n", phase, xtol, chisq_tol);
        jabs_message(MSG_VERBOSE, "Simulation parameters for this phase:\n");
        sim_calc_params_print(fit->sim->params, MSG_VERBOSE);
        jabs_message(MSG_IMPORTANT, "Initializing fit...\n");
        status = gsl_multifit_nlinear_winit(x, &wts.vector, fdf, w);
        if(status != 0) {
            jabs_message(MSG_ERROR, "Fit aborted in initialization.\n");
            fit->stats.error = FIT_ERROR_INIT;
            break;
        }
        /* compute initial cost function */
        f = gsl_multifit_nlinear_residual(w);
        gsl_blas_ddot(f, f, &fit->stats.chisq0);
        jabs_message(MSG_IMPORTANT, "Done. Starting iteration...\n");
        status = jabs_gsl_multifit_nlinear_driver(fit->n_iters_max, xtol, chisq_tol, fit, w); /* Fit */
        fit->stats.error = status;
        if(status < 0) {
            jabs_message(MSG_ERROR, "Fit aborted in phase %i, reason: %s.\n", phase, fit_error_str(fit->stats.error));
            break;
        }
        jabs_message(MSG_IMPORTANT, "Phase %i finished. Time used for actual simulation so far: %.3lf s.\n", phase, fit->stats.cputime_cumul);
        fit_report_results(fit, w, fdf);
    }

    if(fit->stats.error < 0) { /* Revert changes on error */
        for(size_t i = 0; i < fit_params->n; i++) {
            fit_variable *var = &(fit_params->vars[i]);
            *(var->value) = var->value_orig;
        }
    } else { /* Do final calculations when fit was successful */
        /* compute covariance of best fit parameters */
        J = gsl_multifit_nlinear_jac(w);
        gsl_multifit_nlinear_covar(J, 0.0, covar);

        /* compute final cost */
        gsl_blas_ddot(f, f, &fit->stats.chisq);
        fit->stats.chisq_dof = fit->stats.chisq / fit->dof;

        fit_parameters_update(fit_params, w, covar, fit->stats.chisq_dof);
        if(sample_model_renormalize(fit->sm)) {
            jabs_message(MSG_WARNING, "Could not renormalize concentrations of sample model after the fit.\n");
        }
        fit_parameters_update_changed(fit_params); /* sample_model_renormalize() can and will change concentration values, this will recompute error (assuming relative error stays the same) */
        fit_params_print_final(fit_params);
        fit_covar_print(covar, MSG_VERBOSE);
        fit_data_print(fit, MSG_VERBOSE);
#ifdef DEBUG
        fit_uncertainty_print(fit, J, covar, f, &wts.vector, "errors.dat");
#endif
    }
    gsl_multifit_nlinear_free(w);
    gsl_vector_free(fit->f_iter);
    gsl_matrix_free(covar);
    gsl_vector_free(x);
    free(weights);
    free(fdf);
    fit_data_jspace_free(fit);
    return fit->stats.error;
}

int fit_set_roi_from_string(roi *r, const char *str) {
    const char *str_orig = str;
    if(!str)
        return EXIT_FAILURE;
    if(*str != '[') { /* Silent failure, intentionally, since this can signal end of valid roi range arguments. */
#ifdef DEBUG
        fprintf(stderr, "fit_set_roi_from_string() fails silently.\n");
#endif
        return EXIT_FAILURE;
    }
    str++; /* Skipping '[' */
    char *end;
    r->low = strtoull(str, &end, 10);
    str = end;
    if(*str != ':') {
        jabs_message(MSG_ERROR, "Can not parse range from \"%s\". Is ':' missing?\n", str_orig);
        return EXIT_FAILURE;
    }
    str++; /* Skipping ':' */
    r->high = strtoull(str, &end, 10);
    str = end;
    if(*str != ']') {
        jabs_message(MSG_ERROR, "Can not parse range from \"%s\". Is ']' missing near \"%s\"?\n", str_orig, str);
        return EXIT_FAILURE;
    }
    str++;
    if(*str != '\0') {
        jabs_message(MSG_ERROR, "Unexpected input when parsing a range, \"%s\" at end of \"%s\"\n", str, str_orig);
        return EXIT_FAILURE;
    }
    if(r->low > r->high) {
        jabs_message(MSG_ERROR, "Range from %zu to %zu is not valid, because %zu > %zu!\n", r->low, r->high, r->low, r->high);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

double fit_emin(struct fit_data *fit, size_t i_det) {
    double emin = fit->sim->beam_E;
    for(size_t i_range = 0; i_range < fit->n_fit_ranges; i_range++) {
        const roi *r = &(fit->fit_ranges[i_range]);
        if(r->i_det != i_det)
            continue;
        const detector *det = sim_det(fit->sim, i_det);
        for(int Z = JIBAL_ANY_Z; Z <= det->cal_Z_max; Z++) {
            double E = detector_calibrated(det, Z, r->low); /* TODO: assumes calibration is increasing monotonously. Is it guaranteed in all cases? */
            E -= 3.0 * det->calibration->resolution; /* TODO: bad approximation. */
            E *= 0.95;
            E -= 10.0 * C_KEV;
            if(E < emin) {
                emin = E;
            }
        }
    }
    return emin;
}

const char *fit_error_str(int error) {
    switch(error) {
        case FIT_SUCCESS_CHISQ:
            return "chi squared change below tolerance";
        case FIT_SUCCESS_DELTA:
            return "step size below tolerance";
        case FIT_ERROR_NONE:
            return "success";
        case FIT_ERROR_GENERIC:
            return "generic error";
        case FIT_ERROR_MAXITER:
            return "maximum number of iterations reached";
        case FIT_ERROR_NO_PROGRESS:
            return "iteration is not making progress";
        case FIT_ERROR_SANITY:
            return "simulation failed sanity check";
        case FIT_ERROR_WORKSPACE_INITIALIZATION:
            return "simulation workspace could not be initialized";
        case FIT_ERROR_IMPOSSIBLE:
            return "an impossible thing has happened";
        case FIT_ERROR_ABORTED:
            return "user requested abort";
        case FIT_ERROR_INIT:
            return "error during initialization";
        default:
            return "unknown";
    }
}

int fit_range_compare(const void *a, const void *b) {
    const roi *r_a = (const roi *) a;
    const roi *r_b = (const roi *) b;
    if(r_a->i_det < r_b->i_det) { /* i_det is size_t (unsigned) so we have to make this complicated */
        return -1;
    } else if(r_a->i_det > r_b->i_det){
        return 1;
    };
   /* Same detector, compare by low channel */
   if(r_a->low < r_b->low) {
       return -1;
   } else if(r_a->low > r_b->low) {
       return  1;
   }
   /* Same low channel, compare by high channel */
    if(r_a->high < r_b->high) {
        return -1;
    } else if(r_a->high > r_b->high) {
        return  1;
    }
   return 0; /* Ok, they are the same */
}
