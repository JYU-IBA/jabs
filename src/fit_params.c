#include <gsl/gsl_multifit_nlinear.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include "jabs_debug.h"
#include "defaults.h"
#include "message.h"
#include "generic.h"
#include "fit_params.h"
#include "jibal_generic.h"


fit_params *fit_params_new(void) {
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

int fit_params_update(fit_params *p) {
    if(!p)
        return EXIT_FAILURE;
    p->n_active = 0;
    for(size_t i = 0; i < p->n; i++) { /* Count number of active variables and assign index numbers */
        fit_variable *var = &p->vars[i];
        if(var->active) {
            for(size_t j = 0; j < i; j++) { /* Check if this *active* variable is a duplicate (earlier active variable has the same value) */
                if(!p->vars[j].active)
                    continue;
                if(var->value == p->vars[j].value) {
                    jabs_message(MSG_ERROR, "Fit parameters %s and %s can not be used simultaneously, as they are actually the same variable.\n", var->name, p->vars[j].name);
                    return EXIT_FAILURE;
                }
            }
            if(*(var->value) == 0.0) {
                jabs_message(MSG_ERROR, "Fit parameter %s has a zero value, which is not allowed for fitting. Try changing it a little bit or don't fit it.\n", var->name);
                return EXIT_FAILURE;
            }
            var->i_v = p->n_active;
            p->n_active++;
        } else {
            var->i_v = p->n; /* Intentionally invalid index */
        }
    }
    return EXIT_SUCCESS;
}

int fit_params_add_parameter(fit_params *p, fit_variable_type type, double *value, const char *name, const char *unit, double unit_factor, size_t i_det) {
    if(!value || !name) {
        DEBUGMSG("Didn't add a fit parameter since a NULL pointer was passed. Value = %p, name = %p (%s).", (void *)value, (void *)name, name);
        return EXIT_FAILURE;
    }
    for(size_t i = 0; i < p->n; i++) {
        if(p->vars[i].value == value) {
            DEBUGMSG("Fit parameter %s that points to value %p already exists. This is intentional, but beware.", name, (void *)value);
        }
    }
    p->n++;
    p->vars = realloc(p->vars, sizeof(fit_variable) * p->n);
    fit_variable *var = &(p->vars[p->n - 1]);
    var->type = type;
    var->value = value;
    var->value_orig = 0.0;
    var->value_final = 0.0;
    var->err = 0.0;
    var->err_rel = 0.0;
    var->sigmas = 0.0;
    var->name = strdup_non_null(name);
    var->unit = unit;
    var->unit_factor = unit_factor;
    var->active = FALSE;
    var->i_v = SIZE_MAX;
    if(type == FIT_VARIABLE_DETECTOR) {
        var->i_det = i_det;
    } else {
        var->i_det = SIZE_MAX;
    }
    DEBUGMSG("Fit parameter %s added successfully (total %zu).", var->name, p->n);
    return EXIT_SUCCESS;
}

void fit_params_print(const fit_params *params, int active, const char *pattern, jabs_msg_level msg_level) { /* Prints current values of all possible fit variables matching pattern. Pattern can be NULL too. */
    if(!params)
        return;
    const char *whatkind = active ? "active" : "possible";
    if(params->n) {
        if(pattern) {
            jabs_message(msg_level, "All %s fit variables matching pattern \"%s\":\n", whatkind, pattern);
        } else {
            jabs_message(msg_level, "All %s fit variables (use 'show fit variables <pattern>' to see variables matching pattern, wildcards are '*' and '?'):\n", whatkind);
        }
    } else {
        jabs_message(msg_level, "No %s fit variables.\n", whatkind);
    }
    for(size_t i = 0; i < params->n; i++) {
        const fit_variable *var = &params->vars[i];
        if(active && !var->active)
            continue;
        if(pattern && !is_match(var->name, pattern))
            continue;
        jabs_message(msg_level, "%s %24s = %g %s\n", var->active && !active ? "X" : "  ", var->name, *(var->value) / var->unit_factor, var->unit);
        DEBUGMSG("This one is %p, raw value %.12g, i_v = %zu", (void *)var->value, *var->value, var->i_v);
    }
}

void fit_params_print_final(const fit_params *params) { /* Prints final values of active fit variables. */
    if(!params)
        return;
    if(params->n_active) {
        jabs_message(MSG_INFO, "Final fit variables (%zu/%zu):\n", params->n_active, params->n);
    } else {
        jabs_message(MSG_INFO, "No fitted variables of total %zu.\n", params->n);
        return;
    }
    size_t maxlen_name = 8;
    for(size_t i = 0; i < params->n; i++) {
        const fit_variable *var = &params->vars[i];
        if(!var->active)
            continue;
        maxlen_name = JABS_MAX(maxlen_name, strlen(var->name));
    }

    jabs_message(MSG_INFO, "  i | %*s |  unit |        value |     error |  rel %% |  orig. value | multipl. | sigmas |\n", maxlen_name, "variable");
    for(size_t i = 0; i < params->n; i++) {
        const fit_variable *var = &params->vars[i];
        if(!var->active)
            continue;
        //jabs_message(MSG_INFO, stderr, "%24s(%3s) = %12g +- %12g (%.1lf%%)\n", var->name, var->unit, var->value_final/var->unit_factor, var->err/var->unit_factor, var->err/var->value_final*100.0);
        jabs_message(MSG_INFO, "%3zu | %*s | %5s | %12.6g | %9.3g | %6.1lf | %12.6g | %8.4lf | %6.1lf |\n", var->i_v + 1, maxlen_name, var->name, var->unit,
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
        double val_norm = gsl_vector_get(w->x, var->i_v); /* This is the actual fitted value, multiplier to original */
        var->value_final = val_norm * var->value_orig; /* So we get final fitted value by multiplication */
        *(var->value) = var->value_final;
        double err_norm = c * sqrt(gsl_matrix_get(covar, var->i_v, var->i_v)); /* This is the absolute error of the actual fitted value */
        var->err_rel =  err_norm / val_norm; /* Relative error is simple */
        var->err = err_norm * var->value_orig; /* Absolute error of the fitted value is also multiplied */
        var->sigmas = (val_norm - 1.0) / var->err_rel; /* How many standard deviations was the change in value */
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


int fit_params_enable_using_string(fit_params *params, const char *fit_vars) {
    DEBUGMSG("fitvars = %s", fit_vars);
    int status = EXIT_SUCCESS;
    if(!fit_vars)
        return EXIT_FAILURE;
    char *token, *s, *s_orig;
    s_orig = s = strdup(fit_vars);
    while((token = jibal_strsep_with_quotes(&s, ",")) != NULL) { /* parse comma separated list of parameters to fit */
        if(fit_params_enable(params, token, TRUE) == 0) {
            jabs_message(MSG_ERROR, "No matches for %s. See 'show fit variables' for a list of possible fit variables.\n", token);
            status = EXIT_FAILURE;
        }
        if(status == EXIT_FAILURE)
            break;
    }
    free(s_orig);
    return status;
}

fit_variable *fit_params_find_active(const fit_params *params, size_t i_v) {
    for(size_t i = 0; i < params->n; i++) {
        if(params->vars[i].i_v == i_v) {
            return &params->vars[i];
        }
    }
    DEBUGMSG("Could not find active fit parameter number i_v = %zu. This is a bug.", i_v);
    return NULL;
}

const char *fit_variable_type_str(const fit_variable *var) {
    if(!var) {
        return "null";
    }
    switch(var->type) {
        case FIT_VARIABLE_NONE:
            return "none";
        case FIT_VARIABLE_GENERIC:
            return "generic";
        case FIT_VARIABLE_BEAM:
            return "beam";
        case FIT_VARIABLE_GEOMETRY:
            return "geometry";
        case FIT_VARIABLE_DETECTOR:
            return "detector";
        case FIT_VARIABLE_SAMPLE:
            return "sample";
        default:
            return "unknown";
    }
}
