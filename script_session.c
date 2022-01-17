/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <stdlib.h>
#include <string.h>
#include "fit.h"
#include "message.h"
#include "script_session.h"

script_session *script_session_init(jibal *jibal, simulation *sim) {
    if(!jibal)
        return NULL;
    struct script_session *s = malloc(sizeof(struct script_session));
    s->jibal = jibal;
    if(!sim) { /* Sim shouldn't be NULL. If it is, we make a new one. */
        sim = sim_init(jibal);
    }
    s->fit = fit_data_new(jibal, sim); /* Not just fit, but this conveniently holds everything we need. */
    s->cf = jibal_config_file_init(jibal->units);
    if(!s->fit || !s->cf) {
        jabs_message(MSG_ERROR, stderr,"Script session initialization failed.\n");
        free(s);
        return NULL;
    }
    script_session_reset_vars(s);
    s->output_filename = NULL;
    s->bricks_out_filename = NULL;
    s->sample_out_filename = NULL;
    s->detector_out_filename = NULL;
    s->file_depth = 0;
    s->files[0] = NULL;
    return s;
}

int script_session_reset_vars(script_session *s) {
    free(s->cf->vars);
    s->cf->vars = NULL;
    return jibal_config_file_set_vars(s->cf, script_make_vars(s)); /* Loading and resetting things can reset some pointers (like fit->det, so we need to update those to the vars */
}

jibal_config_var *script_make_vars(script_session *s) {
    struct fit_data *fit = s->fit;
    if(!fit)
        return NULL;
    simulation *sim = fit->sim;
    if(!sim)
        return NULL;
    jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT,   "fluence",                &sim->fluence,                       NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy",                 &sim->beam_E,                        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy_broad",           &sim->beam_E_broad,                  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "emin",                   &sim->emin,                          NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "alpha",                  &sim->sample_theta,                  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "sample_azi",             &sim->sample_phi,                    NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling",             &sim->channeling_offset,             NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling_slope",       &sim->channeling_slope,              NULL},
            {JIBAL_CONFIG_VAR_STRING, "output",                 &s->output_filename,                 NULL},
            {JIBAL_CONFIG_VAR_STRING, "bricks_out",             &s->bricks_out_filename,             NULL},
            {JIBAL_CONFIG_VAR_STRING, "sample_out",             &s->sample_out_filename,             NULL},
            {JIBAL_CONFIG_VAR_STRING, "det_out",                &s->detector_out_filename,           NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "erd",                    &sim->erd,                           NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rbs",                    &sim->rbs,                           NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "fit_maxiter",            &fit->n_iters_max,                   NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "fit_xtol",               &fit->xtol,                          NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "fit_gtol",               &fit->gtol,                          NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "fit_ftol",               &fit->ftol,                          NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "ds",                     &sim->params.ds,                     NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rk4",                    &sim->params.rk4,                    NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_incident",      &sim->params.stop_step_incident,     NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_exiting",       &sim->params.stop_step_exiting,      NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "stop_step_fudge",        &sim->params.stop_step_fudge_factor, NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "nucl_stop_accurate",     &sim->params.nucl_stop_accurate,     NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "mean_conc_and_energy",   &sim->params.mean_conc_and_energy,   NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "geostragg",              &sim->params.geostragg,              NULL},
            {JIBAL_CONFIG_VAR_OPTION, "beam_aperture",          &sim->beam_aperture.type, aperture_option},
            {JIBAL_CONFIG_VAR_UNIT,   "beam_aperture_diameter", &sim->beam_aperture.diameter,        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "beam_aperture_width",    &sim->beam_aperture.width,           NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "beam_aperture_height",   &sim->beam_aperture.height,          NULL},
            {JIBAL_CONFIG_VAR_NONE, NULL, NULL,                                                      NULL}
    };
    int n_vars;
    for(n_vars = 0; vars[n_vars].type != 0; n_vars++);
    size_t var_size = sizeof(jibal_config_var)*(n_vars + 1); /* +1 because the null termination didn't count */
    jibal_config_var *vars_out = malloc(var_size);
    if(vars_out) {
        memcpy(vars_out, vars, var_size);
    }
    return vars_out;
}

int script_session_load_script(script_session *s, const char *filename) {
    if(s->file_depth >= SCRIPT_NESTED_MAX) {
        jabs_message(MSG_ERROR, stderr, "Script files nested too deep.\n");
        return EXIT_FAILURE;
    }
    script_file *sfile = script_file_open(filename);
    if(!sfile) {
        jabs_message(MSG_ERROR, stderr, "Can not open file \"%s\".\n", filename);
        return EXIT_FAILURE;
    }
    s->files[s->file_depth] = sfile;

    s->file_depth++;
#ifdef DEBUG
    fprintf(stderr, "Successfully opened file %s, depth now %zu.\n", sfile->filename, s->file_depth);
#endif
    return EXIT_SUCCESS;
}

void script_session_free(script_session *s) {
    if(!s)
        return;
    free(s->output_filename);
    free(s->bricks_out_filename);
    free(s->sample_out_filename);
    free(s->detector_out_filename);
    jibal_config_file_free(s->cf);
    fit_data_workspaces_free(s->fit);
    fit_data_exp_free(s->fit);
    sim_free(s->fit->sim);
    sample_model_free(s->fit->sm);
    fit_data_free(s->fit);
    free(s);
}

int script_get_detector_number(const simulation *sim, int allow_empty, int * const argc, char * const ** const argv, size_t *i_det) {
    char *end;
    if(!argc || !argv || !i_det) {
#ifdef DEBUG
        fprintf(stderr, "Null pointer passed to script_get_detector_number()\n");
#endif
        return EXIT_FAILURE;
    }
    if(*argc < 1) {
        return EXIT_SUCCESS;
    }
    char *s = (*argv)[0];
    if(*s == '\0') {
        return EXIT_FAILURE;
    }
    size_t number = strtoul(s, &end, 10);
    if(end == s) { /* No digits at all! */
        if(allow_empty) {
            return EXIT_SUCCESS; /* First argument was not a number, don't change i_det! */
        } else {
            return EXIT_FAILURE;
        }
    }
    if(*end == '\0') { /* Entire string was valid */
        *i_det = number - 1;
        if(*i_det > sim->n_det) {
            jabs_message(MSG_ERROR, stderr, "Detector number %zu is not valid (n_det = %zu).\n", number, sim->n_det);
            return EXIT_FAILURE;
        }
        *argc -= 1;
        (*argv)++;
        return EXIT_SUCCESS;
    }
#ifdef DEBUG
    fprintf(stderr, "Unknown failure! End points to %p, (== '%c')\n", (void *)end, *end);
#endif
    return EXIT_FAILURE;
}
