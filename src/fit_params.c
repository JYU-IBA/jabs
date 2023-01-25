#include <gsl/gsl_multifit_nlinear.h>
#include <assert.h>
#include "defaults.h"
#include "message.h"
#include "generic.h"
#include "fit_params.h"


fit_params *fit_params_new() {
    fit_params *p = malloc(sizeof(fit_params));
    p->n = 0;
    p->n_active = 0;
    p->vars = NULL;
    return p;
}

void fit_params_free(fit_params *p) {
    if(!p)
        return;
    for(size_t i = 0; i < p->n; i++) {
        free(p->vars[i].name); /* This should be fit_variable_free(), but then p->vars should probably be an array of pointers, too). */
    }
    free(p->vars);
    free(p);
}

void fit_params_update(fit_params *p) {
    if(!p)
        return;

    /* TODO: check if some ACTIVE parameters are duplicates. We don't care about inactive ones. */

    p->n_active = 0;
    for(size_t i = 0; i < p->n; i++) {
        if(p->vars[i].active) {
            p->vars[i].i_v = p->n_active;
            p->n_active++;
        } else {
            p->vars[i].i_v = p->n; /* Intentionally invalid index */
        }
    }
}

int fit_params_add_parameter(fit_params *p, double *value, const char *name, const char *unit, double unit_factor) {
    if(!value || !name) {
#ifdef DEBUG
        fprintf(stderr, "Didn't add a fit parameter since a NULL pointer was passed. Value = %p, name = %p (%s).\n", (void *)value, (void *)name, name);
#endif
        return EXIT_FAILURE;
    }
    for(size_t i = 0; i < p->n; i++) {
        if(p->vars[i].value == value) {
#ifdef DEBUG
            fprintf(stderr, "Didn't add fit parameter %s that points to value %p since it already exists.\n", name, (void *)value);
#endif
            return EXIT_SUCCESS; /* Parameter already exists, don't add. */ /* TODO: maybe this requirement could be relaxed? */
        }
    }
    p->n++;
    p->vars = realloc(p->vars, sizeof(fit_variable) * p->n);
    fit_variable *var = &(p->vars[p->n - 1]);
    var->value = value;
    var->name = strdup_non_null(name);
    var->unit = unit;
    var->unit_factor = unit_factor;
    var->active = FALSE;
#ifdef DEBUG
    fprintf(stderr, "Fit parameter %s added successfully (total %zu).\n", var->name, p->n);
#endif
    return EXIT_SUCCESS;
}

void fit_params_print(const fit_params *params, int active, const char *pattern) { /* Prints current values of all possible fit variables matching pattern. Pattern can be NULL too. */
    if(!params)
        return;
    if(params->n) {
        if(pattern) {
            jabs_message(MSG_INFO, stderr, "All possible fit variables matching pattern \"%s\":\n", pattern);
        } else {
            jabs_message(MSG_INFO, stderr, "All possible fit variables (use 'show fit variables <pattern>' to see variables matching pattern, wildcards are '*' and '?'):\n");
        }
    } else {
        jabs_message(MSG_INFO, stderr, "No fit variables.\n");
    }
    for(size_t i = 0; i < params->n; i++) {
        const fit_variable *var = &params->vars[i];
        if(active && !var->active)
            continue;
        if(pattern && !is_match(var->name, pattern))
            continue;
        jabs_message(MSG_INFO, stderr, "%s %24s = %g %s\n", var->active ? "X" : " ", var->name, *(var->value) / var->unit_factor, var->unit);
    }
}

void fit_params_print_final(const fit_params *params) { /* Prints final values of active fit variables. */
    if(!params)
        return;
    if(params->n_active) {
        jabs_message(MSG_INFO, stderr, "Final fit variables (%zu/%zu):\n", params->n_active, params->n);
    } else {
        jabs_message(MSG_INFO, stderr, "No fitted variables of total %zu.\n", params->n);
        return;
    }
    jabs_message(MSG_INFO, stderr, "  i |                 variable |  unit |        value |      error | rel %% |  orig. value | multipl. | sigmas |\n");
    for(size_t i = 0; i < params->n; i++) {
        const fit_variable *var = &params->vars[i];
        if(!var->active)
            continue;
        //jabs_message(MSG_INFO, stderr, "%24s(%3s) = %12g +- %12g (%.1lf%%)\n", var->name, var->unit, var->value_final/var->unit_factor, var->err/var->unit_factor, var->err/var->value_final*100.0);
        jabs_message(MSG_INFO, stderr, "%3zu | %24s | %5s | %12.6g | %10.4g | %5.1lf | %12.6g | %8.4lf | %6.1lf |\n", var->i_v + 1, var->name, var->unit,
                     var->value_final / var->unit_factor,
                     var->err / var->unit_factor,
                     100.0 * var->err_rel,
                     var->value_orig / var->unit_factor,
                     var->value_final / var->value_orig,
                     var->sigmas
        );
    }
}

size_t fit_params_enable(fit_params *params, const char *s, int enable) {
    size_t n_match = 0;
    for(size_t i = 0; i < params->n; i++) {
        fit_variable *var = &params->vars[i];
        if(is_match(var->name, s)) {
            var->active = enable;
            n_match++;
        }
    }
    return n_match;
}


void fit_parameters_update(const fit_params *fit_params, const gsl_multifit_nlinear_workspace *w, const gsl_matrix *covar, double chisq_dof) {
    double c = GSL_MAX_DBL(1, sqrt(chisq_dof));
    for(size_t i = 0; i < fit_params->n; i++) { /* Update final fitted values to the table (same as used for initial guess) */
        fit_variable *var = &(fit_params->vars[i]);
        if(!var->active)
            continue;
        assert(var->i_v < fit_params->n_active);
        var->value_final = gsl_vector_get(w->x, var->i_v);
        *(var->value) = var->value_final;
        var->err = c * sqrt(gsl_matrix_get(covar, var->i_v, var->i_v));
        var->err_rel = fabs(var->err / var->value_final);
        var->sigmas = fabs(var->value_final - var->value_orig) / var->value_orig / var->err_rel;
    }
}

void fit_parameters_update_changed(const fit_params *fit_params) {
    for(size_t i = 0; i < fit_params->n; i++) {
        fit_variable *var = &(fit_params->vars[i]);
        if(!var->active)
            continue;
        if(*(var->value) != var->value_final) { /* Values changed by something (renormalization) */
            double scale = *(var->value) / var->value_final;
            var->value_final = (*var->value);
            var->err *= scale;
            /* var->err_rel stays the same */
        }
    }
}