#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <jibal_units.h>
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
        fprintf(stderr, "Usage: load [sample|detector|exp|reaction] [file]\n");
        return -1;
    }
    if(strcmp(argv[0], "sample") == 0) {
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
    } else if(strcmp(argv[0], "detector") == 0) {
        if(argc == 2) {
            detector *det = detector_from_file(fit->jibal, argv[1]);
            if(!det) {
                return -1;
            }
            detector_free(fit->sim->det);
            fit->sim->det = det;
            return 0;
        } else {
            fprintf(stderr, "Usage: load detector [file]\n");
        }
    } else if(strcmp(argv[0], "exp") == 0) {
        if(argc == 2) {
            if(!fit->sim->det) {
                fprintf(stderr, "No detector has been set, experimental spectrum can not be read.\n");
                return -1;
            }
            fit->exp = read_experimental_spectrum(argv[1], fit->sim->det);
            if(!fit->exp) {
                return -1;
            }
            return 0;
        } else {
            fprintf(stderr, "Usage: load exp [file]\n");
        }
    } else if(strcmp(argv[0], "reaction") == 0) {
        fprintf(stderr, "Loading reactions from files not implemented yet, sorry.\n");
        return 0;
    }
    fprintf(stderr, "I don't know what to load (%s?)\n", argv[0]);
    return -1;
}

int script_reset(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(!fit) {
        return -1;
    }
    free(fit->fit_ranges);
    fit->fit_ranges = NULL;
    fit->n_fit_ranges = 0;
    fit_params_free(fit->fit_params);
    fit->fit_params = NULL;
    sim_workspace_free(fit->ws);
    fit->ws = NULL;
    sim_free(fit->sim);
    fit->sim = sim_init();
    return 0;
}

int script_show(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc == 0) {
        fprintf(stderr, "Usage show [sim|fit|sample|spectra|detector].\n");
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
        sample_model_print(stderr, fit->sm);
        return 0;
    }
    if(strcmp(argv[0], "spectra") == 0) {
        return print_spectra(NULL, fit->ws, fit->exp);
    }
    if(strcmp(argv[0], "detector") == 0) {
        detector_print(stderr, fit->sim->det);
        return 0;
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
    } else if(strcmp(argv[0], "foil") == 0) {
        if(argc < 2) {
            fprintf(stderr, "Usage: set foil [foil]\nExample: set foil Si 500tfu\n");
            return -1;
        }
        if(!fit->sim->det) {
            fprintf(stderr, "No detector has been set.\n");
            return -1;
        }
        free(fit->sim->det->foil_description);
        fit->sim->det->foil_description = argv_to_string(argc-1, argv+1);
        sample_model *sm = sample_model_from_argv(fit->jibal, argc-1, argv+1);
        fit->sim->det->foil = sample_from_sample_model(sm);
        sample_model_free(sm);
        if(!fit->sim->det->foil) {
            fprintf(stderr, "Could not set foil.\n");
            return -1;
        }
        return 0;
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
        if(argc != 3) {
            fprintf(stderr, "Usage: add fit_range [low] [high]\n");
            return -1;
        }
        fit_range range = {.low = strtoul(argv[1], NULL, 10),
                           .high = strtoul(argv[2], NULL, 10)
        };
        fit_range_add(fit_data, &range);
        return 0;
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
    size_t i = 0;
    for(const struct help_topic *t = topics; t->name != NULL; t++) {
        if(strcmp(t->name, argv[0]) == 0) {
            fputs(t->help_text, stderr);
            if(strcmp(t->name, "help") == 0) {
                i = 0;
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
                i = 0;
                for(jibal_config_var *var = s->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    i++;
                    fprintf(stderr, "%18s", var->name);
                    if(i % 4 == 0) {
                        fputc('\n', stderr);
                    }
                }
                fputc('\n', stderr);
                fprintf(stderr, "\nAlso the following things can be set: ion, sample, foil. Special syntax applies for each.\n");
            }
            return 0;
        }
    }
    fprintf(stderr, "Sorry, no help for '%s'.\n", argv[0]);
    return -1;
}

jibal_config_var *script_make_vars(struct fit_data *fit) {
    if(!fit)
        return NULL;
    simulation *sim = fit->sim;
    detector *det = fit->sim->det;
    jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT,   "fluence",        &sim->p_sr,        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy",         &sim->beam_E,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy_broad",   &sim->beam_E_broad,NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "alpha",          &sim->sample_theta,NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "sample_azi",     &sim->sample_phi,  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling",     &sim->channeling_offset, NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling_slope",&sim->channeling_slope, NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "slope",          &det->slope,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "offset",         &det->offset,     NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "resolution",     &det->resolution, NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "theta",          &det->theta,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "phi",            &det->phi,        NULL},
            {JIBAL_CONFIG_VAR_INT,    "number",         &det->number,     NULL},
            {JIBAL_CONFIG_VAR_INT,    "channels",       &det->channels,   NULL},
            {JIBAL_CONFIG_VAR_INT,    "compress",       &det->compress,   NULL},
            {JIBAL_CONFIG_VAR_NONE,NULL,NULL,NULL}
    };
    int n_vars;
    for(n_vars = 0; vars[n_vars].type != 0; n_vars++);
    size_t s = sizeof(jibal_config_var)*(n_vars + 1); /* +1 because the null termination didn't count */
    jibal_config_var *vars_out = malloc(s);
    if(vars_out) {
        memcpy(vars_out, vars, s);
    }
    return vars_out;
}

int script_simulate(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(script_prepare_sim_or_fit(fit)) {
        return -1;
    }
    fit->ws = sim_workspace_init(fit->sim, fit->jibal);
    simulate_with_ds(fit->ws);
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
    if(script_prepare_sim_or_fit(fit_data)) {
        return -1;
    }
    if(fit(fit_data)) {
        fprintf(stderr, "Fit failed!\n");
        return -1;
    }
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
        if(argc != 2) {
            fprintf(stderr, "Usage: save spectra [file]\n");
            return -1;
        }
        if(!fit_data->ws) {
            fprintf(stderr, "No simulation or fit has been made, cannot save spectra.\n");
            return -1;
        }
        print_spectra(argv[1], fit_data->ws, fit_data->exp);
        return 0;
    } else if(strcmp(argv[0], "sample") == 0) {
        if(argc != 2) {
            fprintf(stderr, "Usage: save sample [file]\n");
        }
        if(!fit_data->sm) {
            fprintf(stderr, "No sample set.\n");
            return -1;
        }
        FILE *f_sout;
        if((f_sout = fopen(argv[1], "w"))) {
            sample_model_print(f_sout, fit_data->sm);
        } else {
            fprintf(stderr, "Could not write sample to file \"%s\".\n", argv[1]);
            return -1;
        }
        return 0;
    } else if(strcmp(argv[0], "det") == 0) {
        if(argc != 2) {
            fprintf(stderr, "Usage: save det [file]\n");
            return -1;
        }
        if(!fit_data->sim->det) {
            fprintf(stderr, "No detector set.\n");
            return -1;
        }
        FILE *f_det;
        if((f_det = fopen(argv[1], "w"))) {
            detector_print(f_det, fit_data->sim->det); /* TODO: is this the fitted detector? */
        } else {
            fprintf(stderr, "Could not write detector to file \"%s\".\n", argv[1]);
            return -1;
        }
        return 0;
    }
    fprintf(stderr, "I don't know what to save (%s?)\n", argv[0]);
    return -1;
}

script_session *script_session_init(jibal *jibal) {
    struct script_session *s = malloc(sizeof(struct script_session));
    s->jibal = jibal;
    s->fit = fit_data_new(jibal, sim_init(), NULL, NULL); /* Not just fit, but this conveniently holds everything we need. */
    s->vars = script_make_vars(s->fit);
    return s;
}
void script_session_free(script_session *s) {
    free(s->vars);
    sim_workspace_free(s->fit->ws);
    sim_free(s->fit->sim);
    fit_data_free(s->fit);
    free(s);
}


int script_process(script_session *s, FILE *f) {
    clock_t start, end;
    struct fit_data *fit = fit_data_new(s->jibal, sim_init(), NULL, NULL); /* Not just fit, but this conveniently holds everything we need. */
    jibal_config_var *vars = script_make_vars(fit);
    char *line=NULL;
    size_t line_size=0;
    size_t lineno=0;
    int interactive = (f == stdin && isatty(fileno(stdin)));
    const char *prompt = "jabs> ";
    if(interactive) {
        fputs(COPYRIGHT_STRING, stderr);
        fprintf(stderr, "Welcome to interactive mode.\nType \"help\" for help or run \"jabs -h\" for command line help.\n\n");
        fputs(prompt, stderr);
    }
    int exit = FALSE;
    while(getline(&line, &line_size, f) > 0) {
        int status = 0;
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
                    start = clock();
                    status = c->f(s, argc - 1, argv + 1);
                    end  = clock();
                    if(c->f == script_load || c->f == script_reset) {
                        free(vars);
                        vars = script_make_vars(fit); /* Loading and resetting things can reset some pointers (like fit->det, so we need to update those to the vars */
                    }
                    if((c->f == script_simulate || c->f == script_fit) && status == 0) {
                        double cputime_total =(((double) (end - start)) / CLOCKS_PER_SEC);
                        fprintf(stderr, "...finished!\n\n");
                        fprintf(stderr, "Total CPU time: %.3lf s.\n", cputime_total);
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
        } else if(status < 0) {
            fprintf(stderr, "Error (%i) on line %zu. Aborting.\n", status, lineno);
            break;
        }
    }
    free(line);

    if(interactive) {
        fprintf(stderr, "Bye.\n");
    }
    return EXIT_SUCCESS;
}

int script_prepare_sim_or_fit(struct fit_data *fit) {
    /* TODO: option to disable RBS or ERD */
    /* TODO: reactions from files? Should "make_reactions" be somewhere else? */
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
    sim_workspace_free(fit->ws);
    fit->ws = NULL;
    sample_free(fit->sim->sample);
    fit->sim->sample = sample_from_sample_model(fit->sm);
    if(!fit->sim->sample) {
        fprintf(stderr, "Could not make a sample based on model description. This should never happen.\n");
        return -1;
    }
    fprintf(stderr, "Simplified sample model for simulation:\n");
    sample_print(stderr, fit->sim->sample, TRUE);

    for(size_t i = 0; i < fit->sim->n_reactions; i++) {
        reaction_free(&fit->sim->reactions[i]);
    }
    free(fit->sim->reactions);
    fit->sim->reactions = NULL;
    fit->sim->n_reactions = 0;
    sim_reactions_add(fit->sim, REACTION_RBS, fit->jibal->config->cs_rbs);
    sim_reactions_add(fit->sim, REACTION_ERD, fit->jibal->config->cs_erd);

    if(fit->sim->n_reactions == 0) {
        fprintf(stderr, "No reactions, nothing to do.\n");
        return -1;
    }
    fprintf(stderr, "Reactions:\n");
    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);

    jibal_gsto_assign_clear_all(fit->jibal->gsto); /* Is it necessary? No. Here? No. Does it clear old stuff? Yes. */
    if(assign_stopping(fit->jibal->gsto, fit->sim)) {
        fprintf(stderr, "Could not assign stopping.\n");
        return -1;
    }
    jibal_gsto_print_assignments(fit->jibal->gsto);
    jibal_gsto_print_files(fit->jibal->gsto, TRUE);
    jibal_gsto_load_all(fit->jibal->gsto);
    simulation_print(stderr, fit->sim);
    return 0;
}
