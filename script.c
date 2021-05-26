#include <string.h>
#include <jibal_units.h>
#include "defaults.h"
#include "script.h"
#include "fit.h"
#include "generic.h"

static const struct script_command commands[] = {
        {"help",    &script_help,           "Print help."},
        {"show",    &script_show,           "Show information on things."},
        {"set",     &script_set,            "Set variables."},
        {"simulate",    &script_simulate,   "Run a simulation."},
        {"reset",   &script_show,           "Reset something."},
        {"exit", NULL, "Exit."},
        {"quit", NULL, NULL},
        {NULL, NULL, NULL}
}; /* TODO: commands: "load", "save" etc... */

void script_print_commands(FILE *f) {
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        if(!c->help_text)
            continue;
        fprintf(f, "%22s    %s\n", c->name, c->help_text);
    }
}

int script_reset(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv) {
    if(!fit) {
        return 0; /* TODO: or error? */
    }
    fit_params_free(fit->fit_params);
    fit->fit_params = NULL;
    sim_workspace_free(fit->ws);
    fit->ws = NULL;
    sample_free(fit->sample);
    fit->sample = NULL;
    return 0;
}

int script_show(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv) {
    if(argc == 0) {
        fprintf(stderr, "Nothing to show.\n");
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
    fprintf(stderr, "Don't know what \"%s\" is.\n", argv[0]);
    return -1;
}

int script_set(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv) {
    if(argc == 0) {
        fprintf(stderr, "Nothing to set.\n");
        return 0;
    }
    if(strcmp(argv[0], "ion") == 0) {
        if(argc != 2) {
            fprintf(stderr, "Usage: set ion [ion]\n");
            return -1;
        }
        fit->sim->beam_isotope = jibal_isotope_find(fit->jibal->isotopes, argv[1], 0, 0);
        if(!fit->sim->beam_isotope) {
            fprintf(stderr, "No such isotope: %s\n", argv[1]);
            return -1;
        }
        return 0;
    }
    if(strcmp(argv[0], "sample") == 0) {
        if(argc < 2) {
            fprintf(stderr, "Usage: set sample [sample]\n");
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
    }
    if(argc == 2) {
        for(jibal_config_var *var = vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
            if(strcmp(argv[0], var->name) == 0) {
                jibal_config_var_set(fit->jibal->units, var, argv[1], NULL);
#ifdef DEBUG
                fprintf(stderr, "%s = %g\n", var->name, *((double *)var->variable));
#endif
                return 0;
            }
        }
    }
    fprintf(stderr, "Don't know what \"%s\" is.\n", argv[0]);
    return -1;
}


int script_help(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv) {
    (void) fit; /* Unused */
    static const struct help_topic topics[] = {
            {"help", "This is help on help. How meta. Help is available on following topics:\n"},
            {"commands", "I recognize the following commands:\n"},
            {"version", "JaBS version: "},
            {"set", "The following variables can be set (unit optional, SI units assumed otherwise):\n"},
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
                    fprintf(stderr, "%16s", t2->name);
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
                for(jibal_config_var *var = vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    i++;
                    fprintf(stderr, "%16s", var->name);
                    if(i % 4 == 0) {
                        fputc('\n', stderr);
                    }
                }
                fputc('\n', stderr);
            }
            return 0;
        }
    }
    fprintf(stderr, "Sorry, no help for '%s'.\n", argv[0]);
    return -1;
}


jibal_config_var *script_make_vars(struct fit_data *fit) {
    simulation *sim = fit->sim;
    detector *det = fit->sim->det;
    jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT,   "fluence",        &sim->p_sr,        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy",         &sim->beam_E,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy_broad",   &sim->beam_E_broad,NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "alpha",          &sim->sample_theta,NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "sample_azi",     &sim->sample_phi,  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling",     &sim->channeling,  NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "fit_low",        &fit->low_ch,      NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "fit_high",        &fit->high_ch,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "slope",          &det->slope,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "offset",         &det->offset,     NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "resolution",     &det->resolution, NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "theta",          &det->theta,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "phi",            &det->phi,        NULL},
            {JIBAL_CONFIG_VAR_INT,    "number",         &det->number,     NULL},
            {JIBAL_CONFIG_VAR_INT,    "channels",       &det->channels,   NULL},
            {JIBAL_CONFIG_VAR_INT,    "compress",       &det->compress,   NULL},
            {JIBAL_CONFIG_VAR_STRING, "foil",           &det->foil_description,       NULL}, /* TODO: det->foil must also be set somewhere */
            {JIBAL_CONFIG_VAR_NONE,NULL,NULL,NULL}
    };
    int n_vars;
    for(n_vars=0; vars[n_vars].type != 0; n_vars++);
    size_t s=sizeof(jibal_config_var)*(n_vars+1); /* +1 because the null termination didn't count */
    jibal_config_var *vars_out=malloc(s);
    memcpy(vars_out, vars, s);
    return vars_out;
}

int script_simulate(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv) {
    /* TODO: free existing data from fit */
    /* TODO: option to disable RBS or ERD */
    /* TODO: reactions from files? Should "make_reactions" be somewhere else? */
    if(!fit->sm) {
        fprintf(stderr, "Cannot simulate, no sample has been defined.\n");
        return -1;
    }
    if(!fit->sim->beam_isotope) {
        fprintf(stderr, "Cannot simulate, no ion has been defined.\n");
        return -1;
    }
    sim_workspace_free(fit->ws);
    fit->ws = NULL;
    sample_free(fit->sample);
    fit->sample = sample_from_sample_model(fit->sm);
    if(!fit->sample) {
        fprintf(stderr, "Could not make a sample based on model description. This should never happen.\n");
        return -1;
    }
    fprintf(stderr, "Simplified sample model for simulation:\n");
    sample_print(stderr, fit->sample, TRUE);

    fit->reactions = make_reactions(fit->sample, fit->sim, fit->jibal->config->cs_rbs, fit->jibal->config->cs_erd);
    if(!fit->reactions || fit->reactions[0] == NULL ) {
        fprintf(stderr, "No reactions, nothing to do.\n");
        return -1;
    }
    fprintf(stderr, "Reactions:\n");
    reactions_print(stderr, fit->reactions);

    jibal_gsto_assign_clear_all(fit->jibal->gsto); /* Is it necessary? No. Here? No. Does it clear old stuff? Yes. */
    if(assign_stopping(fit->jibal->gsto, fit->sim, fit->sample, fit->reactions)) {
        fprintf(stderr, "Could not assign stopping.\n");
        return -1;
    }
    jibal_gsto_print_assignments(fit->jibal->gsto);
    jibal_gsto_print_files(fit->jibal->gsto, TRUE);
    jibal_gsto_load_all(fit->jibal->gsto);
    simulation_print(stderr, fit->sim);
    fit->ws = sim_workspace_init(fit->sim, fit->reactions, fit->sample, fit->jibal);
    simulate_with_ds(fit->ws);
    return 0;
}

int script_process(jibal *jibal, FILE *f) { /* TODO: pass initial fit_data (includes settings in sim!) */
    struct fit_data *fit = fit_data_new(jibal, sim_init(), NULL, NULL, NULL, NULL, 0, 0, 0); /* Not just fit, but this conveniently holds everything we need. */
    jibal_config_var *vars = script_make_vars(fit);
    char *line=NULL;
    size_t line_size=0;
    size_t lineno=0;
    int interactive = (f == stdin);
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
                    status = c->f(fit, vars, argc - 1, argv + 1);
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
        }
    }
    free(line);
    free(vars);
    sim_workspace_free(fit->ws);
    sim_free(fit->sim);
    fit_data_free(fit);
    if(interactive) {
        fprintf(stderr, "Bye.\n");
    }
    return EXIT_SUCCESS;
}
