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
#include <time.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "defaults.h"
#include "generic.h"
#include "jabs.h"
#include "spectrum.h"
#include "message.h"
#include "fit.h"

int fit_function(const gsl_vector *x, void *params, gsl_vector *f) {
    clock_t start, end;
    struct fit_data *fit_data = (struct fit_data *) params;
#ifdef DEBUG
    fprintf(stderr, "Fit iteration %zu, fit function evaluation %zu\n", fit_data->stats.iter, fit_data->stats.n_evals_iter);
#endif
    size_t i_v_active = fit_data->stats.iter_call - 1; /* Which variable is being varied by the fit algorithm in this particular call inside an iteration (compared to first call (== 0)) */
    fit_variable *var_active = NULL;
    for(size_t i_v = 0; i_v < fit_data->fit_params->n; i_v++) {
        fit_variable *var = &fit_data->fit_params->vars[i_v];
        if(var->i_v == i_v_active) {
            var_active = var;
            break;
        }
    }
    fit_data->stats.iter_call++;
#ifdef DEBUG
    fprintf(stderr, "Iter %zu (call %zu). var_active = %p (%s)\n", fit_data->stats.iter, fit_data->stats.iter_call, (void *)var_active, var_active?var_active->name:"none");
#endif

    if(fit_parameters_set_from_vector(fit_data, x)) {
        return GSL_FAILURE;
    }

    if(var_active && var_active->value == &fit_data->sim->fluence) {
        if(fit_scale_by_variable(fit_data, var_active)) {
            return GSL_FAILURE;
        }
        fit_data->stats.n_evals_iter++;
        fit_data->stats.error = fit_set_residuals(fit_data, f);
        if(fit_data->stats.error) {
            return GSL_FAILURE;
        } else {
            return GSL_SUCCESS;
        }
    }

    sample_free(fit_data->sim->sample);
    fit_data->sim->sample = NULL;
    fit_data_workspaces_free(fit_data);

    if(sample_model_sanity_check(fit_data->sm)) {
        fit_data->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }
    if(sim_sanity_check(fit_data->sim)) {
        fit_data->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }
    fit_data->sim->sample = sample_from_sample_model(fit_data->sm);
    sample_renormalize(fit_data->sim->sample);

    if(fit_data_workspaces_init(fit_data)) {
        fit_data_workspaces_free(fit_data); /* Some workspaces may have already been allocated */
        fit_data->stats.error = FIT_ERROR_SANITY;
        return GSL_FAILURE;
    }

#if 0 /* If and when this is enabled, make sure to calculate one complete spectrum with original emin at the end of fit */
    for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) { /* Sets the lowest energy in each simulation according to fit ranges */
        double emin = fit_emin(fit_data, i_det);
        sim_workspace *ws = fit_data_ws(fit_data, i_det);
        if(emin > ws->emin) { /* Only increase the emin, never reduce it. */
            ws->emin = emin;
        }
    }
#endif

    start = clock();
    for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) {
        simulate_with_ds(fit_data->ws[i_det]);
    }
    end = clock();
    fit_data->stats.cputime_iter += (((double) (end - start)) / CLOCKS_PER_SEC);
    fit_data->stats.n_evals_iter++;

    if(fit_data->stats.iter_call == 1) { /* First call of iter, store sum histograms to fit_data */
        fit_data_histo_sum_store(fit_data);
    }
    fit_data->stats.error = fit_set_residuals(fit_data, f);
    if(fit_data->stats.error) {
        return GSL_FAILURE;
    }
    return GSL_SUCCESS;
}

int fit_scale_by_variable(struct fit_data *fit, const fit_variable *var) {
    double scale = *(var->value) / var->value_iter;
#ifdef DEBUG
    fprintf(stderr, "Varying fluence (iter %zu, call %zu) by %12.10lf (variable %s).\n", fit->stats.iter, fit->stats.iter_call, scale, var->name);
#endif
    for(size_t i = 0; i < fit->n_ws; i++) {
        sim_workspace *ws = fit_data_ws(fit, i);
        sim_workspace_histograms_scale(ws, scale);
        sim_workspace_calculate_sum_spectra(ws);
    }
    return GSL_SUCCESS;
}

int fit_set_residuals(const struct fit_data *fit_data, gsl_vector *f) {
    size_t i_vec = 0;
    for(size_t i_range = 0; i_range < fit_data->n_fit_ranges; i_range++) {
        /* TODO: actually calculate which spectra we need to simulate when fitting. Calculate emin based on lowest energy in lowest range. */
        if(i_vec >= f->size) {
            jabs_message(MSG_ERROR, stderr, "Too many channels in fits for the residuals vector. This shouldn't happen.\n");
            return FIT_ERROR_IMPOSSIBLE;
        }

        roi *range = &fit_data->fit_ranges[i_range];
        if(range->i_det >= fit_data->sim->n_det) {
            jabs_message(MSG_ERROR, stderr, "Fit range %zu has detector %zu, but we're only supposed to have %zu detectors!\n", i_range + 1, range->i_det, fit_data->sim->n_det);
            return FIT_ERROR_IMPOSSIBLE;
        }
        sim_workspace *ws = fit_data->ws[range->i_det];
        assert(ws);
        gsl_histogram *exp = fit_data->exp[range->i_det];
        assert(exp);
        for(size_t i = range->low; i <= range->high; i++) {
            if(i >= ws->n_channels) { /* Outside range of simulated spectrum */
                gsl_vector_set(f, i_vec, exp->bin[i]);
            } else {
                gsl_vector_set(f, i_vec, exp->bin[i] - ws->histo_sum->bin[i]);
            }
            i_vec++;
        }
    }
    if(i_vec != f->size) {
        jabs_message(MSG_ERROR, stderr, "Not enough channels in fits for the residuals vector. This shouldn't happen.\n");
        return FIT_ERROR_IMPOSSIBLE;
    }
    return FIT_ERROR_NONE;
}

int fit_parameters_set_from_vector(struct fit_data *fit, const gsl_vector *x) {
    for(size_t i = 0; i < fit->fit_params->n; i++) {
        fit_variable *var = &fit->fit_params->vars[i];
        if(var->active) {
            if(!isfinite(*(var->value))) {
                fprintf(stderr, "Fit iteration %zu, call %zu, fit variable %zu (%s) is not finite.\n", fit->stats.iter, fit->stats.iter_call, var->i_v, var->name);
                fit->stats.error = FIT_ERROR_SANITY;
                return GSL_FAILURE;
            }
            *(var->value) = gsl_vector_get(x, var->i_v);
            if(fit->stats.iter_call == 1) { /* Store the value of fit parameters at the first function evaluation */
                var->value_iter = *(var->value);
            }
#ifdef DEBUG
            fprintf(stderr, "  %zu %12.10lf %12.10lf\n", var->i_v, *(var->value)/var->value_orig, *(var->value)/var->value_iter);
#endif
        }
    }
#ifdef DEBUG
    fprintf(stderr, "\n");
#endif
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
    fit_data->stats.cputime_cumul += fit_data->stats.cputime_iter;
}

void fit_iter_stats_print(const struct fit_stats *stats) {
    jabs_message(MSG_INFO, stderr, "%4zu | %12.6e | %14.8e | %12.7lf | %4zu | %13.3lf | %12.1lf |\n",
                 stats->iter, 1.0 / stats->rcond, stats->norm,
                 stats->chisq_dof, stats->n_evals, stats->cputime_cumul,
                 1000.0 * stats->cputime_iter / stats->n_evals_iter);
}

void fit_stats_print(FILE *f, const struct fit_stats *stats) {
    if(stats->chisq_dof > 0.0) {
        jabs_message(MSG_INFO, f, "Final chisq/dof = %.7lf\n", stats->chisq_dof);
    }
}

int fit_data_fit_range_add(struct fit_data *fit_data, const struct roi *range) { /* Makes a deep copy */
    if(range->low == 0 && range->high == 0) {
        return EXIT_FAILURE;
    }
    if(range->high < range->low) {
        jabs_message(MSG_ERROR, stderr, "Range from %zu to %zu is not valid!\n", range->low, range->high);
        return EXIT_FAILURE;
    }
    fit_data->n_fit_ranges++;
    fit_data->fit_ranges = realloc(fit_data->fit_ranges, fit_data->n_fit_ranges * sizeof(roi));
    if(!fit_data->fit_ranges) {
        fit_data->n_fit_ranges = 0;
        return EXIT_FAILURE;
    }
    fit_data->fit_ranges[fit_data->n_fit_ranges - 1] = *range;
    return EXIT_SUCCESS;
}

void fit_data_fit_ranges_free(struct fit_data *fit_data) {
    if(!fit_data)
        return;
    free(fit_data->fit_ranges);
    fit_data->fit_ranges = NULL;
    fit_data->n_fit_ranges = 0;
}

fit_data *fit_data_new(const jibal *jibal, simulation *sim) {
    struct fit_data *f = malloc(sizeof(struct fit_data));
    f->jibal = jibal;
    f->sim = sim;
    f->n_exp = 0;
    fit_data_exp_alloc(f);
    f->ref = NULL;
    f->sm = NULL; /* Can be set later. */
    f->ws = NULL; /* Initialized later */
    f->n_ws = 0; /* Number of allocated workspaces, initially same as number of detectors. */
    f->fit_params = NULL; /* Holds fit parameters AFTER a fit. */
    f->histo_sum_iter = NULL; /* Initialized later */
    f->fit_iter_callback = NULL; /* Optional */
    f->n_fit_ranges = 0;
    f->fit_ranges = NULL;
    fit_data_defaults(f);
    return f;
}

void fit_data_defaults(fit_data *f) {
    f->n_iters_max = FIT_ITERS_MAX;
    f->xtol = FIT_XTOL;
    f->gtol = FIT_GTOL;
    f->ftol = FIT_FTOL;
    f->chisq_tol = FIT_CHISQ_TOL;
    f->chisq_fast_tol = FIT_FAST_CHISQ_TOL;
    f->phase_start = FIT_PHASE_FAST;
    f->phase_stop = FIT_PHASE_SLOW;
}

void fit_data_free(fit_data *fit) {
    if(!fit)
        return;
    fit_params_free(fit->fit_params);
    fit_data_fit_ranges_free(fit);
    fit_data_histo_sum_free(fit);
    free(fit);
}

void fit_data_roi_print(FILE *f, const struct fit_data *fit_data, const struct roi *roi) {
    if(!fit_data) {
        return;
    }
    gsl_histogram *histo_sim = fit_data_histo_sum(fit_data, roi->i_det);
    gsl_histogram *histo_exp = fit_data_exp(fit_data, roi->i_det);
    gsl_histogram *histo_ref = fit_data->ref;
    size_t n_sim = spectrum_channels_in_range(histo_sim, roi->low, roi->high);
    size_t n_exp = spectrum_channels_in_range(histo_exp, roi->low, roi->high);
    size_t n_ref = spectrum_channels_in_range(histo_ref, roi->low, roi->high);
    double sim_cts = spectrum_roi(histo_sim, roi->low, roi->high);
    double exp_cts = spectrum_roi(histo_exp, roi->low, roi->high);
    double ref_cts = spectrum_roi(histo_ref, roi->low, roi->high);
    const detector *det = sim_det(fit_data->sim, roi->i_det);

    jabs_message(MSG_INFO, f, "          low = %12zu\n", roi->low);
    jabs_message(MSG_INFO, f, "         high = %12zu\n", roi->high);
    if(det) {
        jabs_message(MSG_INFO, f, "        E_low = %12.3lf keV (low energy edge of bin)\n", detector_calibrated(det, JIBAL_ANY_Z, roi->low) / C_KEV);
        jabs_message(MSG_INFO, f, "       E_high = %12.3lf keV (high energy edge of bin)\n", detector_calibrated(det, JIBAL_ANY_Z, roi->high + 1) / C_KEV);
    }
    if(histo_sim) {
        jabs_message(MSG_INFO, f, "        n_sim = %12zu\n", n_sim);
        jabs_message(MSG_INFO, f, "          sim  = %12.8g\n", sim_cts);
    }
    if(histo_exp) {
        jabs_message(MSG_INFO, f, "        n_exp = %12zu\n", n_exp);
        jabs_message(MSG_INFO, f, "          exp  = %12.8g\n", exp_cts);
        jabs_message(MSG_INFO, f, "    sqrt(exp)  = %12.5lf\n", sqrt(exp_cts));
    }
    if(histo_ref) {
        jabs_message(MSG_INFO, f, "        n_ref = %12zu\n", n_ref);
        jabs_message(MSG_INFO, f, "          ref  = %12.8g\n", ref_cts);
        jabs_message(MSG_INFO, f, "    sqrt(ref)  = %12.5lf\n", sqrt(ref_cts));
    }
    if(histo_sim && histo_exp) {
        jabs_message(MSG_INFO, f, "      exp-sim  = %12.8g\n", exp_cts - sim_cts);
        jabs_message(MSG_INFO, f, "      sim/exp  = %12.5lf\n", sim_cts / exp_cts);
        jabs_message(MSG_INFO, f, "      exp/sim  = %12.5lf\n", exp_cts / sim_cts);
        jabs_message(MSG_INFO, f, "  1/sqrt(exp)  = %12.5lf%%\n", 100.0 / sqrt(exp_cts));
        jabs_message(MSG_INFO, f, "(exp-sim)/exp  = %12.5lf%%\n", 100.0 * (exp_cts - sim_cts) / exp_cts);
    }
    if(histo_sim && histo_ref) {
        jabs_message(MSG_INFO, f, "      ref-sim  = %12g\n", ref_cts - sim_cts);
        jabs_message(MSG_INFO, f, "      sim/ref  = %12.5lf\n", sim_cts / ref_cts);
        jabs_message(MSG_INFO, f, "      ref/sim  = %12.5lf\n", ref_cts / sim_cts);
        jabs_message(MSG_INFO, f, "  1/sqrt(ref)  = %12.5lf%%\n", 100.0 / sqrt(ref_cts));
        jabs_message(MSG_INFO, f, "(ref-sim)/ref  = %12.5lf%%\n", 100.0 * (ref_cts - sim_cts) / ref_cts);
    }
}

gsl_histogram *fit_data_exp(const struct fit_data *fit_data, size_t i_det) {
    if(!fit_data || !fit_data->exp)
        return NULL;
    if(i_det >= fit_data->sim->n_det)
        return NULL;
    return fit_data->exp[i_det];
}

gsl_histogram *fit_data_sim(const struct fit_data *fit_data, size_t i_det) {
    sim_workspace *ws = fit_data_ws(fit_data, i_det);
    if(ws)
        return ws->histo_sum;
    else
        return NULL;
}


fit_params *fit_params_all(fit_data *fit) {
    simulation *sim = fit->sim;
    sample_model *sm = fit->sm;
    if(!sim)
        return NULL;
    size_t param_name_max_len = 256; /* Laziness. We use a fixed size temporary string. snprintf is used, so no overflows should occur, but very long names may be truncated. */
    char *param_name = malloc(sizeof(char) * param_name_max_len);
    fit_params *params = fit_params_new();
    fit_params_add_parameter(params, &sim->fluence, "fluence", "", 1.0); /* This must be the first parameter always, as there is a speedup in the fit routine */
    fit_params_add_parameter(params, &sim->channeling_offset, "channeling", "", 1.0);
    fit_params_add_parameter(params, &sim->channeling_slope, "channeling_slope", "1/keV", 1.0 / C_KEV);
    fit_params_add_parameter(params, &sim->sample_theta, "alpha", "deg", C_DEG);
    fit_params_add_parameter(params, &sim->beam_E, "energy", "keV", C_KEV);
    for(size_t i_det = 0; i_det < sim->n_det; i_det++) {
        detector *det = sim_det(sim, i_det);
        char *det_name = NULL;
        if(asprintf(&det_name, "det%zu_", i_det + 1) < 0) {
            return NULL;
        }
        snprintf(param_name, param_name_max_len, "%ssolid", det_name);
        fit_params_add_parameter(params, &det->solid, param_name, "msr", C_MSR);

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
                fit_params_add_parameter(params, calibration_get_param_ref(c, i), param_name, "keV", C_KEV);
            }
        }
        free(det_name);
    }
    if(sm) {
        for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) {
            sample_range *r = &(sm->ranges[i_range]);
            size_t range_index = i_range + 1; /* Human readable indexing */
            if(r->x > 0.0) {
                snprintf(param_name, param_name_max_len, "thick%zu", range_index);
                fit_params_add_parameter(params, &(r->x), param_name, "tfu", C_TFU);
            }

            snprintf(param_name, param_name_max_len, "yield%zu", range_index);
            fit_params_add_parameter(params, &(r->yield), param_name, "", 1.0);

            snprintf(param_name, param_name_max_len, "bragg%zu", range_index);
            fit_params_add_parameter(params, &(r->bragg), param_name, "", 1.0);

            snprintf(param_name, param_name_max_len, "stragg%zu", range_index);
            fit_params_add_parameter(params, &(r->stragg), param_name, "", 1.0);

            if(r->rough.model != ROUGHNESS_NONE && r->rough.x > 0.0) {
                snprintf(param_name, param_name_max_len, "rough%zu", range_index);
                fit_params_add_parameter(params, &(r->rough.x), param_name, "tfu", C_TFU);
            }
            for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
                if(*sample_model_conc_bin(sm, i_range, i_mat) < CONC_TOLERANCE) /* Don't add fit variables for negative, zero or very low concentrations. */
                    continue;
                snprintf(param_name, param_name_max_len, "conc%zu_%s", range_index, sm->materials[i_mat]->name);
                fit_params_add_parameter(params, sample_model_conc_bin(sm, i_range, i_mat), param_name, "%", C_PERCENT);
            }
        }
    }
    free(param_name);
    return params;
}

void fit_data_exp_alloc(fit_data *fit) {
    if(!fit) {
        return;
    }
    size_t n_alloc = fit->sim->n_det;
    if(fit->n_exp == 0) { /* This could be handled by realloc too, but this is an easy way to get null pointers as contents. */
        fit->exp = calloc(n_alloc, sizeof(gsl_histogram *));
    } else {
        fit->exp = realloc(fit->exp, sizeof(gsl_histogram *) * n_alloc);
        for(size_t i = fit->n_exp; i < n_alloc; i++) { /* Reset newly allocated space */
            fit->exp[i] = NULL;
        }
    }
}

void fit_data_exp_free(fit_data *fit) {
    if(!fit->exp)
        return;
    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        if(fit->exp[i_det]) {
            gsl_histogram_free(fit->exp[i_det]);
            fit->exp[i_det] = NULL;
        }
    }
    free(fit->exp);
    fit->exp = NULL;
    fit->n_exp = 0;
}

int fit_data_load_exp(struct fit_data *fit, size_t i_det, const char *filename) {
    gsl_histogram *h = spectrum_read_detector(filename, sim_det(fit->sim, i_det));
    if(!h) {
        jabs_message(MSG_ERROR, stderr, "Reading spectrum from file \"%s\" was not successful.\n", filename);
        return EXIT_FAILURE;
    }
    if(fit->exp[i_det]) {
        gsl_histogram_free(fit->exp[i_det]);
    }
    fit->exp[i_det] = h;
    return EXIT_SUCCESS;
}

void fit_data_histo_sum_free(struct fit_data *fit_data) {
    if(!fit_data->histo_sum_iter)
        return;
    for(size_t i = 0; i < fit_data->n_histo_sum; i++) {
        if(fit_data->histo_sum_iter[i]) {
            gsl_histogram_free(fit_data->histo_sum_iter[i]);
        }
    }
    free(fit_data->histo_sum_iter);
    fit_data->histo_sum_iter = NULL;
    fit_data->n_histo_sum = 0;
}

void fit_data_histo_sum_store(struct fit_data *fit_data) {
    fit_data_histo_sum_free(fit_data);
    fit_data->histo_sum_iter = calloc(fit_data->n_ws, sizeof(gsl_histogram *));
    for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) {
        sim_workspace *ws = fit_data_ws(fit_data, i_det);
        fit_data->histo_sum_iter[i_det] = gsl_histogram_clone(ws->histo_sum);
    }
    fit_data->n_histo_sum = fit_data->n_ws;
}

gsl_histogram *fit_data_histo_sum(const struct fit_data *fit_data, size_t i_det) {
    if(!fit_data || !fit_data->histo_sum_iter)
        return NULL;
    if(i_det >= fit_data->n_histo_sum)
        return NULL;
    return fit_data->histo_sum_iter[i_det];
}

int fit_data_add_det(struct fit_data *fit, detector *det) {
    if(!fit || !det)
        return EXIT_FAILURE;
    size_t n_det_old = fit->sim->n_det;
    if(sim_det_add(fit->sim, det)) {
        return EXIT_FAILURE;
    }
    fit_data_exp_alloc(fit); /* Number of detectors in sim changed, let these guys know it too */
    return EXIT_SUCCESS;
}

sim_workspace *fit_data_ws(const struct fit_data *fit_data, size_t i_det) {
    if(!fit_data || !fit_data->ws)
        return NULL;
    if(i_det >= fit_data->n_ws)
        return NULL;
    return fit_data->ws[i_det];
}


size_t fit_data_ranges_calculate_number_of_channels(const struct fit_data *fit_data) {
    size_t sum = 0;
    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *r = &fit_data->fit_ranges[i];
        detector *det = sim_det(fit_data->sim, r->i_det);
        if(!det) {
            continue;
        }
        if(r->high >= det->channels) { /* Limited by detector */
            sum += (det->channels - r->low);
        } else {
            sum += (r->high - r->low) + 1;
        }
    }
    return sum;
}

int fit_data_workspaces_init(struct fit_data *fit_data) {
    int status = EXIT_SUCCESS;
    fit_data_workspaces_free(fit_data);
    fit_data->n_ws = fit_data->sim->n_det;
    fit_data->ws = calloc(fit_data->n_ws, sizeof(sim_workspace *));
    if(!fit_data->ws)
        return EXIT_FAILURE;
    for(size_t i = 0; i < fit_data->n_ws; i++) {
        detector *det = sim_det(fit_data->sim, i);
        if(detector_sanity_check(det)) {
            jabs_message(MSG_ERROR, stderr, "Detector %zu failed sanity check!\n", i);
            status = EXIT_FAILURE;
            break;
        }
        detector_update(det);
        spectrum_set_calibration(fit_data_exp(fit_data, i), det, JIBAL_ANY_Z); /* Update the experimental spectra calibration (using default calibration) */
        sim_workspace *ws = sim_workspace_init(fit_data->jibal, fit_data->sim, det);
        if(!ws) {
            jabs_message(MSG_ERROR, stderr, "Workspace %zu failed to initialize!\n", i);
            status = EXIT_FAILURE;
        }
        fit_data->ws[i] = ws;
    }
    if(status == EXIT_FAILURE) {
        fit_data_workspaces_free(fit_data);
    }
    return status;
}

void fit_data_workspaces_free(struct fit_data *fit_data) {
    if(!fit_data->ws) {
        return;
    }
    for(size_t i_det = 0; i_det < fit_data->n_ws; i_det++) {
        sim_workspace_free(fit_data->ws[i_det]);
        fit_data->ws[i_det] = NULL;
    }
    free(fit_data->ws);
    fit_data->ws = NULL;
    fit_data->n_ws = 0;
}

struct fit_stats fit_stats_init() {
    struct fit_stats s;
    s.n_evals = 0;
    s.n_evals_iter = 0;
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

void fit_data_print(FILE *f, const struct fit_data *fit_data) {
    if(!fit_data) {
        return;
    }
    if(fit_data->n_fit_ranges == 0) {
        jabs_message(MSG_ERROR, f, "No fit ranges.\n");
        return;
    }
    jabs_message(MSG_INFO, f, "%zu fit ranges:\n", fit_data->n_fit_ranges);
    jabs_message(MSG_INFO, f, "  i |    low |   high |    exp cts |    sim cts | sim/exp |  sigmas |\n");
    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *range = &fit_data->fit_ranges[i];
        double exp_cts = spectrum_roi(fit_data_exp(fit_data, range->i_det), range->low, range->high);
        double sim_cts = spectrum_roi(fit_data_histo_sum(fit_data, range->i_det), range->low, range->high);
        if(exp_cts == 0.0) {
            jabs_message(MSG_INFO, f, "%3zu | %6lu | %6lu | %10.0lf | %10.1lf |         |         |\n", i + 1, range->low, range->high, exp_cts, sim_cts);
        } else {
            double ratio = sim_cts / exp_cts;
            double sigmas = (sim_cts - exp_cts) / sqrt(exp_cts);
            jabs_message(MSG_INFO, f, "%3zu | %6lu | %6lu | %10.0lf | %10.1lf | %7.5lf | %7.2lf |\n", i + 1, range->low, range->high, exp_cts, sim_cts, ratio, sigmas);
        }
    }
    jabs_message(MSG_INFO, f, "\nFit has %zu channels total.\n", fit_data_ranges_calculate_number_of_channels(fit_data));
}

int jabs_test_delta(const gsl_vector *dx, const gsl_vector *x, double epsabs, double epsrel) { /* test_delta() copied from GSL convergence.c and modified */
    int ok = TRUE;
    for(size_t i = 0; i < x->size; i++) {
        double xi = gsl_vector_get(x, i);
        double dxi = gsl_vector_get(dx, i);
        double tolerance = epsabs + epsrel * fabs(xi);
        double rel = fabs(dxi) / tolerance; /* "How many times over the acceptable tolerance are we */
#ifdef DEBUG
        fprintf(stderr, "Test delta: i %zu, xi %g, dxi %g, tolerance %g, rel %g\n", i, xi, dxi, tolerance, rel);
#endif
        if(rel >= 1.0) {
            ok = FALSE;
            break;
        }
    }
    if(ok)
        return GSL_SUCCESS;
    return GSL_CONTINUE;
}

int jabs_gsl_multifit_nlinear_driver(const size_t maxiter, const double xtol, const double chisq_tol, struct fit_data *fit_data, gsl_multifit_nlinear_workspace *w) {
    int status = 0;
    size_t iter;
    double chisq_dof_old;
    jabs_message(MSG_INFO, stderr, "iter |    cond(J)   |     |f(x)|     |   chisq/dof  | eval | cpu cumul (s) | cpu/eval (ms)|\n");
    for(iter = 0; iter <= maxiter; iter++) {
        fit_data->stats.iter_call = 0;
        fit_data->stats.iter = iter;
        if(iter) {
            chisq_dof_old = fit_data->stats.chisq_dof;
            fit_data->stats.cputime_iter = 0.0;
            fit_data->stats.n_evals_iter = 0;
            status = gsl_multifit_nlinear_iterate(w);
#ifdef DEBUG
            fprintf(stderr, "Iteration status %i (%s)\n", status, gsl_strerror(status));
#endif
        }
        if(fit_data->stats.error) {
            return fit_data->stats.error;
        }
        if(status == GSL_ENOPROG && iter == 1) {
            return FIT_ERROR_NO_PROGRESS;
        }
        fit_iter_stats_update(fit_data, w);
        fit_iter_stats_print(&fit_data->stats);
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
            jabs_message(MSG_WARNING, stderr, "Chisq increased, this probably shouldn't happen.\n");
        }
        if(chisq_change < chisq_tol) {
            return FIT_SUCCESS_CHISQ;
        }
    }
    return FIT_ERROR_MAXITER;
}

void fit_report_results(const fit_data *fit, const gsl_multifit_nlinear_workspace *w, const gsl_multifit_nlinear_fdf *fdf) {
    jabs_message(MSG_INFO, stderr, "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    jabs_message(MSG_INFO, stderr, "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    jabs_message(MSG_INFO, stderr, "function evaluations: %zu\n", fit->stats.n_evals);
#ifdef DEBUG
    jabs_message(MSG_INFO, stderr, "function evaluations (GSL): %zu\n", fdf->nevalf);
#endif
    jabs_message(MSG_INFO, stderr, "Jacobian evaluations: %zu\n", fdf->nevaldf);
    jabs_message(MSG_INFO, stderr, "reason for stopping: %s\n", fit_error_str(fit->stats.error));
    jabs_message(MSG_INFO, stderr, "initial |f(x)| = %f\n", sqrt(fit->stats.chisq0));
    jabs_message(MSG_INFO, stderr, "final   |f(x)| = %f\n", sqrt(fit->stats.chisq));
}


void fit_covar_print(const gsl_matrix *covar) {
    jabs_message(MSG_INFO, stderr, "\nCorrelation coefficients matrix:\n       | ");
    for(size_t i = 0; i < covar->size1; i++) {
        jabs_message(MSG_INFO, stderr, " %4zu  ", i + 1);
    }
    jabs_message(MSG_INFO, stderr, "\n");
    for(size_t i = 0; i < covar->size1; i++) {
        jabs_message(MSG_INFO, stderr, "%6zu | ", i + 1);
        for(size_t j = 0; j <= i && j < covar->size2; j++) {
            jabs_message(MSG_INFO, stderr, " %6.3f", gsl_matrix_get(covar, i, j) / sqrt(gsl_matrix_get(covar, i, i) * gsl_matrix_get(covar, j, j)));
        }
        jabs_message(MSG_INFO, stderr, "\n");
    }
}

int fit(struct fit_data *fit_data) {
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    fdf_params.trs = gsl_multifit_nlinear_trs_lm;
    fdf_params.solver = gsl_multifit_nlinear_solver_svd; /* Robust? */
    struct fit_params *fit_params = fit_data->fit_params;
    if(!fit_params || fit_params->n_active == 0) {
        jabs_message(MSG_ERROR, stderr, "No parameters to fit.\n");
        return EXIT_FAILURE;
    }
    if(!fit_data->exp) {
        jabs_message(MSG_ERROR, stderr, "No experimental spectrum to fit.\n");
        return EXIT_FAILURE;
    }
    gsl_multifit_nlinear_fdf fdf;
    fdf.params = fit_data;
    if(!fit_data->exp) {
        jabs_message(MSG_ERROR, stderr, "No experimental data, can not fit.\n");
        return EXIT_FAILURE;
    }
    if(!fit_data->n_fit_ranges) {
        jabs_message(MSG_ERROR, stderr, "No fit range(s) given, can not fit.\n");
        return EXIT_FAILURE;
    }

    for(size_t i = 0; i < fit_data->n_fit_ranges; i++) {
        roi *range = &fit_data->fit_ranges[i];
        jabs_message(MSG_INFO, stderr, "Fit range %zu [%zu:%zu]\n", i + 1, range->low, range->high);
    }

    fdf.f = &fit_function;
    fdf.df = NULL; /* Jacobian, with NULL using finite difference. */
    fdf.fvv = NULL; /* No geodesic acceleration */
    fdf.n = fit_data_ranges_calculate_number_of_channels(fit_data);
    fdf.p = fit_params->n_active;
    if(fdf.n < fdf.p) {
        jabs_message(MSG_ERROR, stderr, "Not enough data (%zu points) for given number of free parameters (%zu)\n", fdf.n, fdf.p);
        return -1;
    } else {
        jabs_message(MSG_INFO, stderr, "%zu channels and %zu parameters in fit, %zu degrees of freedom.\n", fdf.n, fdf.p, fdf.n - fdf.p);
    }
    gsl_vector *f;
    gsl_matrix *J;
    fit_data->dof = fdf.n - fdf.p;
    int status;

    double *weights = malloc(sizeof(double) * fdf.n);
    if(!weights)
        return 1;
    size_t i_w = 0;
    for(size_t i_range = 0; i_range < fit_data->n_fit_ranges; i_range++) {
        roi *range = &fit_data->fit_ranges[i_range];
        assert(range);
        detector *det = sim_det(fit_data->sim, range->i_det);
        gsl_histogram *exp = fit_data_exp(fit_data, range->i_det);
        if(!det) {
            jabs_message(MSG_ERROR, stderr, "Detector %zu (fit range %zu) does not exist.\n", range->i_det + 1, i_range + 1);
            free(weights);
            return 1;
        }
        if(!exp) {
            jabs_message(MSG_ERROR, stderr, "Experimental spectrum for detector %zu (fit range %zu) does not exist.\n", range->i_det + 1, i_range + 1);
            free(weights);
            return 1;
        }
        for(size_t i = range->low; i <= range->high && i < det->channels; i++) {
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
    assert(i_w == fdf.n);

    gsl_matrix *covar = gsl_matrix_alloc(fit_params->n_active, fit_params->n_active);
    gsl_vector *x = gsl_vector_alloc(fit_params->n_active);
    for(size_t i = 0; i < fit_params->n; i++) { /* Update all (including inactives) */
        fit_variable *var = &(fit_params->vars[i]);
        var->err = 0.0;
        var->value_orig = *(var->value);
    }

    gsl_vector_view wts = gsl_vector_view_array(weights, i_w);

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc(T, &fdf_params, fdf.n, fdf.p);

    sim_calc_params p_orig = *fit_data->sim->params; /* Store original values (will be used in final stage of fitting) */
    for(int phase = fit_data->phase_start; phase <= fit_data->phase_stop; phase++) { /* Phase 1 is "fast", phase 2 normal. */
        assert(phase >= FIT_PHASE_FAST && phase <= FIT_PHASE_SLOW);
        double xtol = fit_data->xtol;
        double chisq_tol = fit_data->chisq_tol;
        /* initialize solver with starting point and weights */
        fit_data->stats = fit_stats_init();
        fit_data->stats.phase = phase;
        if(fit_data->fit_iter_callback) { /* First call to callback quickly (before most initialization) */
            if(fit_data->fit_iter_callback(fit_data->stats)) {
                fit_data->stats.error = FIT_ERROR_ABORTED;
                break;
            }
        }
        if(phase == FIT_PHASE_FAST) {
            sim_calc_params_defaults_fast(fit_data->sim->params); /* Set current parameters to be faster in phase 0. */
            xtol *= FIT_FAST_XTOL_MULTIPLIER;
            chisq_tol = fit_data->chisq_fast_tol;
        } else if(phase == FIT_PHASE_SLOW) {
            sim_calc_params_copy(&p_orig, fit_data->sim->params);
        }
        sim_calc_params_update(fit_data->sim->params);
        for(size_t i = 0; i < fit_params->n; i++) { /* Set active variables to vector */
            fit_variable *var = &(fit_params->vars[i]);
            if(var->active) {
                gsl_vector_set(x, var->i_v, *(var->value)); /* Start phase from initial or fitted (previous phase) results */
            }
        }
        jabs_message(MSG_INFO, stderr, "\nInitializing fit phase %i. Xtol = %e, chisq_tol %e\n", phase, xtol, chisq_tol);
        jabs_message(MSG_INFO, stderr, "Simulation parameters for this phase:\n");
        sim_calc_params_print(fit_data->sim->params);
        gsl_multifit_nlinear_winit(x, &wts.vector, &fdf, w);

        /* compute initial cost function */
        f = gsl_multifit_nlinear_residual(w);
        gsl_blas_ddot(f, f, &fit_data->stats.chisq0);

        status = jabs_gsl_multifit_nlinear_driver(fit_data->n_iters_max, xtol, chisq_tol, fit_data, w); /* Fit */
        fit_data->stats.error = status;
        if(status < 0) {
            jabs_message(MSG_ERROR, stderr, "Fit aborted in phase %i, reason: %s.\n", phase, fit_error_str(fit_data->stats.error));
            break;
        }
        jabs_message(MSG_INFO, stderr, "Phase %i finished. CPU time used for actual simulation so far: %.3lf s.\n", phase, fit_data->stats.cputime_cumul);
        fit_report_results(fit_data, w, &fdf);
    }

    if(fit_data->stats.error < 0) { /* Revert changes on error */
        for(size_t i = 0; i < fit_params->n; i++) {
            fit_variable *var = &(fit_params->vars[i]);
            *(var->value) = var->value_orig;
        }
    } else { /* Do final calculations when fit was successful */
        /* compute covariance of best fit parameters */
        J = gsl_multifit_nlinear_jac(w);
        gsl_multifit_nlinear_covar(J, 0.0, covar);

        /* compute final cost */
        gsl_blas_ddot(f, f, &fit_data->stats.chisq);
        fit_data->stats.chisq_dof = fit_data->stats.chisq / fit_data->dof;

        fit_parameters_update(fit_params, w, covar, fit_data->stats.chisq_dof);
        sample_model_renormalize(fit_data->sm);
        fit_parameters_update_changed(fit_params); /* sample_model_renormalize() can and will change concentration values, this will recompute error (assuming relative error stays the same) */
        fit_params_print_final(fit_params);
        fit_covar_print(covar);
    }
    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);
    gsl_vector_free(x);
    free(weights);
    return fit_data->stats.error;
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
        jabs_message(MSG_ERROR, stderr, "Can not parse range from \"%s\". Is ':' missing?\n", str_orig);
        return EXIT_FAILURE;
    }
    str++; /* Skipping ':' */
    r->high = strtoull(str, &end, 10);
    str = end;
    if(*str != ']') {
        jabs_message(MSG_ERROR, stderr, "Can not parse range from \"%s\". Is ']' missing near \"%s\"?\n", str_orig, str);
        return EXIT_FAILURE;
    }
    str++;
    if(*str != '\0') {
        jabs_message(MSG_ERROR, stderr, "Unexpected input when parsing a range, \"%s\" at end of \"%s\"\n", str, str_orig);
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
        case FIT_ERROR_IMPOSSIBLE:
            return "an impossible thing has happened";
        case FIT_ERROR_ABORTED:
            return "user requested abort";
        default:
            return "unknown";
    }
}
