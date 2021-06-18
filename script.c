#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <jibal_units.h>
#include <jibal_generic.h>
#include "defaults.h"
#include "fit.h"
#include "generic.h"
#include "spectrum.h"
#include "script.h"

static const struct script_command commands[] = {
        {"help",    &script_help,           "Print help."},
        {"show",    &script_show,           "Show information on things."},
        {"set",     &script_set,            "Set variables."},
        {"add",     &script_add,            "Add things."},
        {"simulate",    &script_simulate,   "Run a simulation."},
        {"load",    &script_load,           "Load something."},
        {"reset",   &script_reset,           "Reset something."},
        {"fit",     &script_fit,            "Do a fit."},
        {"roi",     &script_roi,            "Show information from a region of interest."},
        {"save",    &script_save,           "Save something."},
        {"exit", NULL, "Exit."},
        {"quit", NULL, NULL},
        {NULL, NULL, NULL}
}; /* TODO: more commands... */

void script_print_commands(FILE *f) {
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        if(!c->help_text)
            continue;
        fprintf(f, "%22s    %s\n", c->name, c->help_text);
    }
}

int script_load(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc == 0) {
        fprintf(stderr, "Usage: load [script|sample|detector|exp|reaction] [file]\n");
        return -1;
    }
    if(strcmp(argv[0], "script") == 0) {
        if(argc == 2) {
            return script_process(s, argv[1]);
        } else {
            fprintf(stderr, "Usage: load script [file]\n");
        }
        return 0;
    } else if(strcmp(argv[0], "sample") == 0) {
        if(argc == 2) {
            sample_model *sm = sample_model_from_file(fit->jibal, argv[1]);
            if(!sm) {
                fprintf(stderr, "Sample load from \"%s\" failed.\n", argv[1]);
                return -1;
            }
            sample_model_free(fit->sm);
            fit->sm = sm;
            return 0;
        } else {
            fprintf(stderr, "Usage: load sample [file]\n");
        }
        return 0;
    } else if(strcmp(argv[0], "det") == 0) {
        int status;
        size_t i_det = 0;
        if(argc == 3) {
            i_det = strtoul(argv[1], NULL, 10);
            argc--;
            argv++;
        }
        if(argc != 2) {
            fprintf(stderr, "Usage: load det [number] file\n");
            return EXIT_SUCCESS;
        }
        detector *det = detector_from_file(fit->jibal, argv[1]);
        if(!det)
            return EXIT_FAILURE;
        if(fit->sim->n_det == 0 && i_det == 0) {
            status = fit_data_add_det(fit, det); /* Adds a new detector (and space for experimental spectrum) */
        } else {
            status = sim_det_set(fit->sim, det, i_det);
        }
        return status;
    } else if(strcmp(argv[0], "exp") == 0) {
        size_t i_det = 0;
        if(argc == 3) {
            i_det = strtoul(argv[1], NULL, 10);
            argc--;
            argv++;
        }
        if(argc != 2) {
            fprintf(stderr, "Usage: load exp [number] file\n");
            return EXIT_SUCCESS;
        }
        fit->exp[i_det] = spectrum_read(argv[1], sim_det(fit->sim, i_det));
        if(!fit->exp[i_det]) {
            fprintf(stderr, "Reading spectrum from file \"%s\" was not successful.\n", argv[1]);
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    } else if(strcmp(argv[0], "reaction") == 0) {
        fprintf(stderr, "Loading reactions from files not implemented yet, sorry.\n");
        return EXIT_SUCCESS;
    }
    fprintf(stderr, "I don't know what to load (%s?)\n", argv[0]);
    return EXIT_FAILURE;
}

int script_reset(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(!fit_data) {
        return -1;
    }
    fit_data_fit_ranges_free(fit_data);
    fit_params_free(fit_data->fit_params);
    fit_data->fit_params = NULL;
    fit_data_workspaces_free(fit_data);
    fit_data->ws = NULL;
    sim_free(fit_data->sim);
    fit_data->sim = sim_init(NULL);
    return 0;
}

int script_show(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc == 0) {
        fprintf(stderr, "Usage show [sim|fit|sample|detector].\n");
        return 0;
    }
    if(strcmp(argv[0], "sim") == 0) {
         simulation_print(stderr, fit->sim);
         return 0;
     }
    if(strcmp(argv[0], "fit") == 0) {
        fit_data_print(stderr, fit);
        return 0;
    }
    if(strcmp(argv[0], "sample") == 0) {
        return sample_model_print(NULL, fit->sm);
    }
    if(strcmp(argv[0], "det") == 0) {

        for(size_t i_det = 0 ; i_det < fit->sim->n_det; i_det++) {  /* TODO: prettier output, maybe a table */
            fprintf(stderr, "DETECTOR %zu\n", i_det);
            detector_print(NULL, fit->sim->det[i_det]);
            fprintf(stderr, "\n");
        }
        return EXIT_SUCCESS;
    }
    fprintf(stderr, "Don't know what \"%s\" is.\n", argv[0]);
    return -1;
}

int script_set(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        fprintf(stderr, "Nothing to set. See \"help set\" for more information.\n");
        return 0;
    }
    if(strcmp(argv[0], "ion") == 0) {
        if(argc != 2) {
            fprintf(stderr, "Usage: set ion [ion]\nExample: set ion 4He\n");
            return -1;
        }
        fit->sim->beam_isotope = jibal_isotope_find(fit->jibal->isotopes, argv[1], 0, 0);
        if(!fit->sim->beam_isotope) {
            fprintf(stderr, "No such isotope: %s\n", argv[1]);
            return -1;
        }
        return 0;
    } else if(strcmp(argv[0], "sample") == 0) {
        if(argc < 2) {
            fprintf(stderr, "Usage: set sample [sample]\nExample: set sample TiO2 1000tfu Si 10000tfu\n");
            return -1;
        }
        sample_model *sm_new = sample_model_from_argv(fit->jibal, argc-1, argv+1);
        if(sm_new) {
            sample_model_free(fit->sm);
            fit->sm = sm_new;
        } else {
            fprintf(stderr, "Sample is not valid.\n");
            return -1;
        }
        return 0;
    } else if(strcmp(argv[0], "det") == 0) {
        size_t i_det = 0;
        if(argc == 4) {
            i_det = strtoul(argv[1], NULL, 10);
            argc--;
            argv++;
        }
        if(argc != 3) {
            fprintf(stderr, "Usage: set det [number] variable value\n");
            return EXIT_FAILURE;
        }
        if(detector_set_var(s->jibal, sim_det(fit->sim, i_det), argv[1], argv[2])) {
            fprintf(stderr, "Can't set \"%s\" to be \"%s\"!\n", argv[1], argv[2]);
        }
        return EXIT_SUCCESS;
    }
    for(jibal_config_var *var = s->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
        if(strcmp(argv[0], var->name) == 0) {
            if(argc != 2) {
                fprintf(stderr, "Usage: set %s [value]\n", var->name);
                return -1;
            }
            jibal_config_var_set(fit->jibal->units, var, argv[1], NULL);
#ifdef DEBUG
            fprintf(stderr, "%s = %g\n", var->name, *((double *)var->variable));
#endif
            return 0;
        }
    }
    fprintf(stderr, "Don't know what \"%s\" is.\n", argv[0]);
    return -1;
}

int script_add(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc < 1) {
        fprintf(stderr, "Nothing to add. See \"help add\" for more information.\n");
        return 0;
    }
    if(strcmp(argv[0], "fit_range") == 0) {
        roi range = {.i_det = 0};
        if(argc == 4) {
            range.i_det = strtoul(argv[1], NULL, 10);
            argc--;
            argv++;
        }
        if(argc == 3) {
            range.low = strtoul(argv[1], NULL, 10);
            range.high = strtoul(argv[2], NULL, 10);
        } else {
            fprintf(stderr, "Usage: add fit_range [detector] low high\n");
            return -1;
        }
        fit_data_fit_range_add(fit_data, &range);
        return 0;
    } else if(strcmp(argv[0], "det") == 0) {
        if(argc != 2) {
            fprintf(stderr, "Usage: add det filename\n");
            return EXIT_FAILURE;
        }
        detector *det = detector_from_file(s->jibal, argv[1]);
        if(!det) {
            return EXIT_FAILURE;
        }
        return fit_data_add_det(fit_data, det);
    }
    fprintf(stderr, "Don't know what \"%s\" is.\n", argv[0]);
    return -1;
}


int script_help(script_session *s, int argc, char * const *argv) {
    (void) fit; /* Unused */
    static const struct help_topic topics[] = {
            {"help", "This is help on help. How meta. Help is available on following topics:\n"},
            {"commands", "I recognize the following commands:\n"},
            {"version", "JaBS version: "},
            {"set", "The following variables can be set (unit optional, SI units assumed otherwise):\n"},
            {"show", "Show things (print to screen).\n"},
            {NULL, NULL}
    };
    if(argc == 0) {
        fprintf(stderr, "Type help [topic] for information on a particular topic or \"help help\" for help on help.\n");
        return 0;
    }
    for(const struct help_topic *t = topics; t->name != NULL; t++) {
        if(strcmp(t->name, argv[0]) == 0) {
            fputs(t->help_text, stderr);
            if(strcmp(t->name, "help") == 0) {
                size_t i = 0;
                for(const struct help_topic *t2 = topics; t2->name != NULL; t2++) {
                    i++;
                    fprintf(stderr, "%18s", t2->name);
                    if(i % 4 == 0) {
                        fputc('\n', stderr);
                    }
                }
                fprintf(stderr, "\n");
            } else if(strcmp(t->name, "commands") == 0) {
                script_print_commands(stderr);
            } else if(strcmp(t->name, "version") == 0) {
                fprintf(stderr, "%s\n", jabs_version());
            } else if(strcmp(t->name, "set") == 0) {
                size_t i = 0;
                for(jibal_config_var *var = s->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    if(var->type != JIBAL_CONFIG_VAR_UNIT)
                        continue;
                    i++;
                    fprintf(stderr, "%18s", var->name);
                    if(i % 4 == 0) {
                        fputc('\n', stderr);
                    }
                }
                i = 0;
                fprintf(stderr, "\nThe following variables are not in SI units:\n");
                for(jibal_config_var *var = s->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    if(var->type == JIBAL_CONFIG_VAR_UNIT)
                        continue;
                    i++;
                    fprintf(stderr, "%18s", var->name);
                    if(i % 4 == 0) {
                        fputc('\n', stderr);
                    }
                }
                fprintf(stderr, "\nAlso the following things can be set: ion, sample, foil. Special syntax applies for each.\n");
            }
            return 0;
        }
    }
    fprintf(stderr, "Sorry, no help for '%s'.\n", argv[0]);
    return -1;
}

void script_make_vars(script_session *s) {
    free(s->vars);
    struct fit_data *fit = s->fit;
    if(!fit)
        return;
    simulation *sim = fit->sim;
    jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT, "fluence", &sim->fluence, NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy",         &sim->beam_E,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy_broad",   &sim->beam_E_broad,NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "alpha",          &sim->sample_theta,NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "sample_azi",     &sim->sample_phi,  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling",     &sim->channeling_offset, NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling_slope",&sim->channeling_slope, NULL},
            {JIBAL_CONFIG_VAR_STRING, "output",         &s->output_filename,      NULL},
            {JIBAL_CONFIG_VAR_STRING, "bricks_out",     &s->bricks_out_filename,   NULL},
            {JIBAL_CONFIG_VAR_STRING, "sample_out",     &s->sample_out_filename,   NULL},
            {JIBAL_CONFIG_VAR_STRING, "det_out",        &s->detector_out_filename, NULL},
            {JIBAL_CONFIG_VAR_NONE,NULL,NULL,                              NULL}
    };
    int n_vars;
    for(n_vars = 0; vars[n_vars].type != 0; n_vars++);
    size_t var_size = sizeof(jibal_config_var)*(n_vars + 1); /* +1 because the null termination didn't count */
    jibal_config_var *vars_out = malloc(var_size);
    if(vars_out) {
        memcpy(vars_out, vars, var_size);
    }
    s->vars = vars_out;
}

int script_simulate(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(script_prepare_sim_or_fit(s)) {
        return EXIT_FAILURE;
    }
    if(fit_data_workspaces_init(fit)) {
        fprintf(stderr, "Could not initialize simulation workspace(s).\n");
        return EXIT_FAILURE;
    }
    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        simulate_with_ds(fit->ws[i_det]);
    }
    script_finish_sim_or_fit(s);
    return 0;
}

int script_fit(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc != 1) {
        fprintf(stderr, "Usage: fit [fitvar1,fitvar2,...]\n");
        return -1;
    }
    fit_params_free(fit_data->fit_params);
    fit_data->fit_params = fit_params_new();
    fit_params_add(fit_data->sim, fit_data->sm, fit_data->fit_params, argv[0]);
    if(fit_data->fit_params->n == 0) {
        fprintf(stderr, "No parameters for fit.\n");
        return -1;
    }
    if(!fit_data->exp) { /* TODO: not enough to check this */
        fprintf(stderr, "No experimental spectrum set.\n");
        return -1;
    }
    if(fit_data->n_fit_ranges == 0) {
        fprintf(stderr, "No fit range(s) given.\n");
        return -1;
    }
    if(script_prepare_sim_or_fit(s)) {
        return -1;
    }
    if(fit(fit_data)) {
        fprintf(stderr, "Fit failed!\n");
        return -1;
    }
    script_finish_sim_or_fit(s);
    fprintf(stderr, "\nFinal parameters:\n");
    simulation_print(stderr, fit_data->sim);
    fprintf(stderr, "\nFinal profile:\n");
    sample_print(NULL, fit_data->sim->sample, FALSE);
    sample_areal_densities_print(stderr, fit_data->sim->sample, FALSE);
    fprintf(stderr, "\nFinal sample model:\n");
    sample_model_print(NULL, fit_data->sm);
    fprintf(stderr, "\n");
    fit_stats_print(stderr, &fit_data->stats);
    return 0;
}

int script_save(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc < 1) {
        fprintf(stderr, "Nothing to save. See \"help save\" for more information.\n");
        return -1;
    }
    if(strcmp(argv[0], "spectra") == 0) {
        size_t i_det = 0;
        if(argc == 3) {
            i_det = strtoul(argv[1], NULL, 10);
            argc--;
            argv++;
        }
        if(argc != 2) {
            fprintf(stderr, "Usage: save spectra [detector] file\n");
            return EXIT_FAILURE;
        }
        if(print_spectra(argv[1], fit_data_ws(fit_data, i_det), fit_data_exp(fit_data, i_det))) {
            fprintf(stderr, "Could not save spectra of detector %zu to file \"%s\"!\n", i_det, argv[1]);
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    } else if(strcmp(argv[0], "sample") == 0) {
        if(argc != 2) {
            fprintf(stderr, "Usage: save sample [file]\n");
        }
        if(!fit_data->sm) {
            fprintf(stderr, "No sample set.\n");
            return -1;
        }
        if(sample_model_print(argv[1], fit_data->sm)) {
            fprintf(stderr, "Could not write sample to file \"%s\".\n", argv[1]);
            return -1;
        }
        return 0;
    } else if(strcmp(argv[0], "det") == 0) {
        size_t i_det = 0;
        if(argc == 3) {
            i_det = strtoul(argv[1], NULL, 10);
            argc--;
            argv++;
        }
        if(argc != 2) {
            fprintf(stderr, "Usage: save det [detector] file\n");
            return EXIT_FAILURE;
        }
        if(detector_print(argv[1], sim_det(fit_data->sim, i_det))) {
            fprintf(stderr, "Could not write detector %zu to file \"%s\".\n", i_det, argv[1]);
            return -1;
        }
        return 0;
    }
    fprintf(stderr, "I don't know what to save (%s?)\n", argv[0]);
    return -1;
}

int script_roi(script_session *s, int argc, char * const *argv) {
    struct roi r;
    if(argc == 3) {
        r.i_det = strtoul(argv[0], NULL, 10);
        argv++;
        argc--;
    }
    if (argc == 2) {
        r.low = strtoul(argv[0], NULL, 10);
        r.high = strtoul(argv[1], NULL, 10);
    } else {
        fprintf(stderr, "Usage: roi [det] low high");
    }
    fit_data_roi_print(stderr, s->fit, &r);
    return EXIT_SUCCESS;
}

script_session *script_session_init(jibal *jibal, simulation *sim) {
    if(!jibal)
        return NULL;
    struct script_session *s = malloc(sizeof(struct script_session));
    s->jibal = jibal;
    if(!sim) { /* Sim shouldn't be NULL */
        sim = sim_init(NULL);
    }
    s->fit = fit_data_new(jibal, sim); /* Not just fit, but this conveniently holds everything we need. */
    s->vars = NULL;
    s->output_filename = NULL;
    s->bricks_out_filename = NULL;
    s->sample_out_filename = NULL;
    s->detector_out_filename = NULL;
    script_make_vars(s);
    return s;
}
void script_session_free(script_session *s) {
    free(s->output_filename);
    free(s->bricks_out_filename);
    free(s->sample_out_filename);
    free(s->detector_out_filename);
    free(s->vars);
    fit_data_workspaces_free(s->fit);
    fit_data_exp_free(s->fit);
    sim_free(s->fit->sim);
    sample_model_free(s->fit->sm);
    fit_data_free(s->fit);
    free(s);
}


int script_process(script_session *s, const char *filename) {
    char *line=NULL;
    size_t line_size=0;
    size_t lineno=0;
    FILE *f = fopen_file_or_stream(filename, "r");
    if(!f) {
        return EXIT_FAILURE;
    }
    int interactive = (f == stdin && isatty(fileno(stdin)));
    const char *prompt = "jabs> ";
    if(interactive) {
        fputs(prompt, stderr);
    } else if(filename) {
        fprintf(stderr, "\nRunning script \"%s\"\n\n", filename);
    }
    int exit = FALSE;
    int status = 0;
    while(getline(&line, &line_size, f) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strip newlines */
        if(*line == '#') /* Comment */
            continue;
        char **argv = string_to_argv(line);
        if(!argv) {
            fprintf(stderr, "Something went wrong in parsing arguments.\n");
            continue;
        }
        char **a = argv;
        int argc = 0;
        while(*a != NULL) {
            a++;
            argc++;
        }
#ifdef DEBUG
        for(int i = 0; i < argc; i++) {
            fprintf(stderr, "args: %i: \"%s\"\n", i, argv[i]);
        }
#endif
        if(argc) {
            int found = FALSE;
            for(const struct script_command *c = commands; c->name != NULL; c++) {
                if(strcmp(c->name, argv[0]) == 0) {
                    found = TRUE;
                    if(c->f == NULL) {
                        exit = TRUE;
                        break;
                    }
                    status = c->f(s, argc - 1, argv + 1);
                    if(c->f == script_load || c->f == script_reset) {
                        script_make_vars(s); /* Loading and resetting things can reset some pointers (like fit->det, so we need to update those to the vars */
                    }
                    break;
                }
            }
            if(!found) {
                if(interactive) {
                    fprintf(stderr, "Command \"%s\" not recognized. See 'help commands'.\n", line);
                } else {
                    fprintf(stderr, "Command \"%s\" not recognized on line %zu\n", line, lineno);
                    return EXIT_FAILURE; /* TODO: leaks memory */
                }
            }
        }
        free(argv[0]);
        free(argv);
        if(exit)
            break;
        if(interactive) {
            fputs(prompt, stderr);
        } else if(status) {
            fprintf(stderr, "Error (%i) on line %zu. Aborting.\n", status, lineno);
            break;
        }
    }
    free(line);
    fclose_file_or_stream(f);
    if(interactive) {
        fprintf(stderr, "Bye.\n");
    } else if(filename) {
        if(status) {
            fprintf(stderr, "Error running script \"%s\"\n", filename);
        } else {
            fprintf(stderr, "Finished running script \"%s\"\n", filename);
        }
    }
    return status;
}

int script_prepare_sim_or_fit(script_session *s) {
    /* TODO: option to disable RBS or ERD */
    /* TODO: reactions from files? Should "sim_reactions_add" be somewhere else? */
    fit_data *fit = s->fit;
    if(!fit->sm) {
        fprintf(stderr, "No sample has been defined!\n");
        return -1;
    }
    if(!fit->sim->beam_isotope) {
        fprintf(stderr, "No ion has been defined!\n");
        return -1;
    }
    if(!fit->sim->det) { /* Shouldn't be possible */
        fprintf(stderr, "No detector has been defined!\n");
        return -1;
    }
    fit_data_workspaces_free(s->fit);
    sample_free(fit->sim->sample);
    fit->sim->sample = sample_from_sample_model(fit->sm);
    if(!fit->sim->sample) {
        fprintf(stderr, "Could not make a sample based on model description. This should never happen.\n");
        return -1;
    }
    fprintf(stderr, "Simplified sample model for simulation:\n");
    sample_print(NULL, fit->sim->sample, TRUE);

    for(size_t i = 0; i < fit->sim->n_reactions; i++) {
        reaction_free(&fit->sim->reactions[i]);
    }
    free(fit->sim->reactions);
    fit->sim->reactions = NULL;
    fit->sim->n_reactions = 0;
    sim_reactions_add(fit->sim, REACTION_RBS, fit->jibal->config->cs_rbs, 0.0); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
    sim_reactions_add(fit->sim, REACTION_ERD, fit->jibal->config->cs_erd, 0.0);

    if(fit->sim->n_reactions == 0) {
        fprintf(stderr, "No reactions, nothing to do.\n");
        return -1;
    }
    fprintf(stderr, "\nReactions:\n");
    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);

    jibal_gsto_assign_clear_all(fit->jibal->gsto); /* Is it necessary? No. Here? No. Does it clear old stuff? Yes. */
    if(assign_stopping(fit->jibal->gsto, fit->sim)) {
        fprintf(stderr, "Could not assign stopping.\n");
        return -1;
    }
    jibal_gsto_print_assignments(fit->jibal->gsto);
    jibal_gsto_print_files(fit->jibal->gsto, TRUE);
    fprintf(stderr, "Loading stopping data.\n");
    jibal_gsto_load_all(fit->jibal->gsto);
    simulation_print(stderr, fit->sim);
    s->start = clock();
    return 0;
}

int script_finish_sim_or_fit(script_session *s) {
    s->end = clock();
    double cputime_total = (((double) (s->end - s->start)) / CLOCKS_PER_SEC);
    fprintf(stderr, "...finished!\n\n");
    fprintf(stderr, "Total CPU time: %.3lf s.\n", cputime_total);
    struct fit_data *fit = s->fit;
#ifdef MULTIDET_FIXED
    if(s->output_filename) {
        print_spectra(s->output_filename, fit->ws, fit->exp);
    }
    if(s->bricks_out_filename) {
        print_bricks(s->bricks_out_filename, fit->ws);
    }
#endif
    if(s->sample_out_filename) {
        sample_model_print(s->sample_out_filename, fit->sm);
    }
#ifdef MULTIDET_FIXED
    if(s->detector_out_filename) {
        detector_print(s->detector_out_filename, fit->ws->det);
    }
#endif
    return 0;
}
