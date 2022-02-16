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
#include "script_command.h"
#include "script_file.h"
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
    if(!s->fit) {
        jabs_message(MSG_ERROR, stderr,"Script session initialization failed.\n");
        free(s);
        return NULL;
    }
    s->output_filename = NULL;
    s->file_depth = 0;
    s->files[0] = NULL;
    s->commands = script_commands_create(s);
    return s;
}

int script_session_load_script(script_session *s, const char *filename) {
    if(s->file_depth >= SCRIPT_FILES_NESTED_MAX) {
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
    fit_data_workspaces_free(s->fit);
    fit_data_exp_free(s->fit);
    sim_free(s->fit->sim);
    sample_model_free(s->fit->sm);
    fit_data_free(s->fit);
    script_commands_free(s->commands);
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
#ifdef DEBUG
    fprintf(stderr, "Trying to determine if %s is a detector number.\n", s);
#endif
    if(strcmp(s, "first") == 0) {
        if(sim->n_det) {
            *i_det = 0;
        } else {
            jabs_message(MSG_ERROR, stderr, "No detectors.\n");
            return EXIT_FAILURE;
        }
    } else if(strcmp(s, "last") == 0) {
        if(sim->n_det) {
            *i_det = sim->n_det - 1;
        } else {
            jabs_message(MSG_ERROR, stderr, "No detectors.\n");
            return EXIT_FAILURE;
        }
    } else {
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
            if(*i_det >= sim->n_det) {
                jabs_message(MSG_ERROR, stderr, "Detector number %zu is not valid (n_det = %zu).\n", number, sim->n_det);
                return EXIT_FAILURE;
            }
        } else {
#ifdef DEBUG
            fprintf(stderr, "Unknown failure! End points to %p, (== '%c')\n", (void *)end, *end);
#endif
            return EXIT_FAILURE;
        }
    }
    /* Branches that consumed an argument successfully lead here */
    (*argc)--;
    (*argv)++;
    return EXIT_SUCCESS;
}
