#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#include <jibal_units.h>
#include "fit.h"
#include "generic.h"
#include "spectrum.h"
#include "script.h"
#include "jabs.h"
#include "message.h"
#include "geostragg.h"


void script_print_commands(FILE *f, const struct script_command *commands) {
    if(!commands)
        return;
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        if(!c->help_text)
            continue;
        jabs_message(MSG_INFO, f, " %16s    %s\n", c->name, c->help_text);
    }
}

script_command_status script_load_script(script_session *s, int argc, char * const *argv) {
    if(argc < 1) {
        jabs_message(MSG_INFO, stderr, "Usage: load script [file]\n");
    }
     return script_process(s, argv[1]);
}

script_command_status script_load_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load sample [file]\n");
        return EXIT_FAILURE;
    }
    sample_model *sm = sample_model_from_file(fit->jibal, argv[1]);
    if(!sm) {
        jabs_message(MSG_INFO, stderr, "Sample load from \"%s\" failed.\n", argv[1]);
        return -1;
    }
    sample_model_free(fit->sm);
    fit->sm = sm;
    return 0;
}

script_command_status script_load_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR,  stderr,"Usage: load detector [detector] filename\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return sim_det_set(fit->sim, detector_from_file(fit->jibal, argv[0]), i_det);
}

script_command_status script_load_experimental(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR,  stderr,"Usage: load experimental [detector] filename\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    gsl_histogram *h = spectrum_read(argv[0], sim_det(fit->sim, i_det));
    if(!h) {
        jabs_message(MSG_ERROR,  stderr,"Reading spectrum from file \"%s\" was not successful.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit->exp[i_det]) {
        gsl_histogram_free(fit->exp[i_det]);
    }
    fit->exp[i_det] = h;
    return EXIT_SUCCESS;
}

script_command_status script_load_reaction(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load reaction [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return sim_reactions_add_r33(fit->sim, fit->jibal->isotopes, argv[0]);
}


script_command_status script_load(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc == 0) {
        jabs_message(MSG_INFO, stderr, "Usage: load [script|sample|det|exp|reaction] [file]\n");
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    if(strcmp(argv[0], "script") == 0) {

    } else if(strcmp(argv[0], "sample") == 0) {

        return 0;
    } else if(strcmp(argv[0], "det") == 0) {

    } else if(strcmp(argv[0], "detectors") == 0) {
        if(argc != 2) {
            jabs_message(MSG_INFO, stderr, "Usage: load detectors file\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        size_t n_d = 0;
        detector **d = detectors_from_file(fit->jibal, argv[1], &n_d);
        if(d && n_d) {
            /* TODO: set detectors, number of detectors etc... Clear spectra? */
            fit->sim->det = d;
        }
    } else if(strcmp(argv[0], "exp") == 0) {

    } else if(strcmp(argv[0], "reaction") == 0) {

    }
    return SCRIPT_COMMAND_NOT_FOUND;
}

script_command_status script_reset(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(!fit_data) {
        return -1;
    }
    if(argc >= 1 && strcmp(argv[0], "reactions") == 0) {
        sim_reactions_free(fit_data->sim);
        return SCRIPT_COMMAND_SUCCESS;
    } else if(argc >= 1 && strcmp(argv[0], "detectors") == 0) {
        for(size_t i_det = 0; i_det < s->fit->sim->n_det; i_det++) {
            detector_free(sim_det(s->fit->sim, i_det));
        }
        s->fit->sim->n_det = 0;
        return SCRIPT_COMMAND_SUCCESS;
    } else if(argc >= 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: reset [reactions|detectors]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_data_fit_ranges_free(fit_data);
    fit_params_free(fit_data->fit_params);
    fit_data->fit_params = NULL;
    fit_data_exp_free(s->fit);
    fit_data_workspaces_free(fit_data);
    sample_model_free(fit_data->sm);
    fit_data->sm = NULL;
    sim_free(fit_data->sim);
    fit_data->sim = sim_init(s->jibal);
    fit_data->exp = calloc(fit_data->sim->n_det, sizeof(gsl_histogram *));
    jibal_gsto_assign_clear_all(fit_data->jibal->gsto);
    script_session_reset_vars(s);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_WARNING, stderr, "No sample has been set.\n");
        return SCRIPT_COMMAND_SUCCESS;
    } else {
        return sample_model_print(NULL, fit->sm);
    }
}

script_command_status script_show_simulation(script_session *s, int argc, char * const *argv) {
    simulation_print(stderr, s->fit->sim);
    return SCRIPT_COMMAND_SUCCESS;
}
script_command_status script_show_fit(script_session *s, int argc, char * const *argv) {
    fit_data_print(stderr, s->fit);
    return SCRIPT_COMMAND_SUCCESS;
}
script_command_status script_show_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    detector_print(NULL, fit->sim->det[i_det]);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_reactions(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;

    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_INFO, stderr, "No reactions.\n");
        return SCRIPT_COMMAND_SUCCESS;
    }
    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_variables(script_session *s, int argc, char * const *argv) {
    const jibal_config_var *var;
    for(var = s->cf->vars; var->type != 0; var++) {
        if(var->variable == NULL)
            continue;
        switch(var->type) {
            case JIBAL_CONFIG_VAR_NONE:
                break;
            case JIBAL_CONFIG_VAR_PATH:
            case JIBAL_CONFIG_VAR_STRING:
                if(*((void **) var->variable) == NULL)
                    continue;
                jabs_message(MSG_INFO, stderr, "%s = %s\n", var->name, *((char **) var->variable));
                break;
            case JIBAL_CONFIG_VAR_BOOL:
                jabs_message(MSG_INFO, stderr, "%s = %s\n", var->name, *((int *) var->variable) ? "true" : "false");
                break;
            case JIBAL_CONFIG_VAR_INT:
                jabs_message(MSG_INFO, stderr, "%s = %i\n", var->name, *((int *) var->variable));
                break;
            case JIBAL_CONFIG_VAR_DOUBLE:
                jabs_message(MSG_INFO, stderr, "%s = %g\n", var->name, *((double *) var->variable));
                break;
            case JIBAL_CONFIG_VAR_UNIT:
                jabs_message(MSG_INFO, stderr, "%s = %g\n", var->name, *((double *) var->variable));
                break;
            case JIBAL_CONFIG_VAR_OPTION:
                jabs_message(MSG_INFO, stderr, "%s = %s\n", var->name, jibal_option_get_string(var->option_list, *((int *) var->variable)));
                break;
            case JIBAL_CONFIG_VAR_SIZE:
                jabs_message(MSG_INFO, stderr, "%s = %zu\n", var->name, *((size_t *) var->variable));
                break;
        }
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Nothing to set. See 'help set' for more information.\n");
        return 0;
    }
    if(strcmp(argv[0], "ion") == 0) {
        if(argc != 2) {
            jabs_message(MSG_ERROR, stderr, "Usage: set ion [ion]\nExample: set ion 4He\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        fit->sim->beam_isotope = jibal_isotope_find(fit->jibal->isotopes, argv[1], 0, 0);
        if(!fit->sim->beam_isotope) {
            jabs_message(MSG_ERROR, stderr,"No such isotope: %s\n", argv[1]);
            return SCRIPT_COMMAND_FAILURE;
        }
        return 0;
    } else if(strcmp(argv[0], "sample") == 0) {
        if(argc < 2) {
            jabs_message(MSG_ERROR, stderr, "Usage: set sample [sample]\nExample: set sample TiO2 1000tfu Si 10000tfu\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        sample_model *sm_new = sample_model_from_argv(fit->jibal, argc-1, argv+1);
        if(sm_new) {
            sample_model_free(fit->sm);
            fit->sm = sm_new;
        } else {
            jabs_message(MSG_ERROR, stderr, "Sample is not valid.\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        return SCRIPT_COMMAND_SUCCESS;
    } else if(strcmp(argv[0], "det") == 0) {
        size_t i_det = 1;
        if(argc == 4) {
            i_det = strtoul(argv[1], NULL, 10);
            if(i_det == 0) {
                jabs_message(MSG_ERROR, stderr, "Usage: set det [number] variable value\n");
            }
            argc--;
            argv++;
        }
        if(i_det == 1 && fit->sim->n_det == 0) { /* This happens on first "set det" */
            jabs_message(MSG_VERBOSE, stderr, "No detectors were defined. Detector was added.\n");
            fit_data_add_det(fit, detector_default(NULL));
        }
        if(i_det == 0 || i_det > fit->sim->n_det) {
            jabs_message(MSG_ERROR, stderr, "Detector number not valid.\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        i_det--;
        if(argc != 3) {
            jabs_message(MSG_ERROR, stderr, "Usage: set det [number] variable value\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        if(detector_set_var(s->jibal, sim_det(fit->sim, i_det), argv[1], argv[2])) {
            jabs_message(MSG_ERROR, stderr, "Can't set \"%s\" to be \"%s\"!\n", argv[1], argv[2]);
        }
        return EXIT_SUCCESS;
    }
    if(argc != 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: set variable [value]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(jibal_config_file_var_set(s->cf, argv[0], argv[1])) {
        jabs_message(MSG_ERROR, stderr,"Error in setting \"%s\" to \"%s\". Does the variable exist? Use show vars.\n", argv[0], argv[1]);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

script_command_status script_add_reaction(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: add reaction TYPE isotope\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    reaction *r = reaction_make_from_argv(fit->jibal, fit->sim->beam_isotope, argc - 1, argv + 1);
    if(!r) {
        jabs_message(MSG_ERROR, stderr, "Could not make a reaction based on given description.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(r->cs == JIBAL_CS_NONE) {
        jibal_cross_section_type cs = sim_cs(fit->sim, r->type);
        r->cs = cs;
        jabs_message(MSG_VERBOSE, stderr, "Reaction cross section not given or not valid, assuming default for %s: %s.\n",
                     reaction_name(r), jibal_cs_types[cs].s);
    }
    return sim_reactions_add_reaction(fit->sim, r);
}
script_command_status script_add_reactions(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR, stderr, "Cannot add reactions before sample has been set (I need to know which reactions to add!).\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc > 0) {
        if(strcmp(argv[0], "RBS") == 0) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS));
        } else if(strcmp(argv[0], "ERD") == 0) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
        } else {
            jabs_message(MSG_ERROR, stderr, "What is %s anyways? Did you mean \"RBS\" or \"ERD\"?\n", argv[0]);
            return SCRIPT_COMMAND_FAILURE;
        }
        return SCRIPT_COMMAND_SUCCESS;
    }
    if(fit->sim->rbs) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS)); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
    }
    if(fit->sim->erd) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_add_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR,  stderr,"Usage: add detector filename\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return fit_data_add_det(fit, detector_from_file(fit->jibal, argv[0])); /* Adds a new detector (and space for experimental spectrum) */
}

script_command_status script_add_fit_range(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
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
        jabs_message(MSG_ERROR, stderr, "Usage: add fit_range [detector] low high\n");
        return -1;
    }
    fit_data_fit_range_add(fit, &range);
    return 0;
}

script_command_status script_help(script_session *s, int argc, char * const *argv) {
    (void) fit; /* Unused */
    static const struct help_topic topics[] = {
            {"help", "This is help on help. How meta.\nHelp is available on following topics:\n"},
            {"commands", "I recognize the following commands (try 'help' followed by command name):\n"},
            {"version", "JaBS version: "},
            {"set", "The following variables can be set (unit optional, SI units assumed otherwise):\n"},
            {"show", "Show things. Possible things: sim, fit, sample, det, vars.\n"},
            {"fit", "Make a fit. Provide list of variables to fit.\n"},
            {NULL, NULL}
    };
    if(argc == 0) {
        jabs_message(MSG_INFO, stderr, "Type help [topic] for information on a particular topic or \"help help\" for help on help.\n");
        return 0;
    }

    int found = 0;
    for(const struct help_topic *t = topics; t->name != NULL; t++) {
        if(strcmp(t->name, argv[0]) == 0) {
            found++;
            jabs_message(MSG_INFO, stderr, "%s", t->help_text);
            if(strcmp(t->name, "help") == 0) {
                size_t i = 0;
                for(const struct help_topic *t2 = topics; t2->name != NULL; t2++) {
                    i++;
                    jabs_message(MSG_INFO, stderr, "%18s", t2->name);
                    if(i % 4 == 0) {
                        jabs_message(MSG_INFO, stderr, "\n");
                    }
                }
                jabs_message(MSG_INFO, stderr, "\n");
            } else if(strcmp(t->name, "commands") == 0) {
                script_print_commands(stderr, script_commands);
            } else if(strcmp(t->name, "version") == 0) {
                jabs_message(MSG_INFO, stderr, "%s\n", jabs_version());
            } else if(strcmp(t->name, "set") == 0) {
                if(!s->cf || !s->cf->vars)
                    break;
                size_t i = 0;
                for(jibal_config_var *var = s->cf->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    if(var->type != JIBAL_CONFIG_VAR_UNIT)
                        continue;
                    i++;
                    jabs_message(MSG_INFO, stderr," %25s", var->name);
                    if(i % 3 == 0) {
                        jabs_message(MSG_INFO, stderr,"\n");
                    }
                }
                fprintf(stderr, "\n\nThe following variables are not in SI units:\n");
                for(jibal_config_var *var = s->cf->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    if(var->type == JIBAL_CONFIG_VAR_UNIT)
                        continue;
                    jabs_message(MSG_INFO, stderr, " %25s: %s\n", var->name, jibal_config_var_type_name(var->type));
                }
                jabs_message(MSG_INFO, stderr,"\n\nAlso the following things can be set: ion, sample, det. Special syntax applies for each.\n");
            }
            return 0;
        }
    }

    for(const struct script_command *c = script_commands; c->name != NULL; c++) {
        if(strcmp(c->name, argv[0]) == 0) {
            if(!found) { /* There wasn't a help topic  */
                jabs_message(MSG_INFO, stderr, "\"%s\" is a valid command, but no additional help is available!\n\n", c->name);
            }
            found++;
            if(c->subcommands) {
                jabs_message(MSG_INFO, stderr, "At least the following sub-commands are recognized:\n", c->name);
                script_print_commands(stderr, c->subcommands);
            }
            break;
        }
    }

    if(!found) {
        jabs_message(MSG_ERROR, stderr,"Sorry, no help for '%s' available.\n", argv[0]);
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    return SCRIPT_COMMAND_SUCCESS;
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

script_command_status script_simulate(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(script_prepare_sim_or_fit(s)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data_workspaces_init(fit)) {
        jabs_message(MSG_ERROR, stderr, "Could not initialize simulation workspace(s).\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        if(simulate_with_ds(fit->ws[i_det])) {
            jabs_message(MSG_ERROR, stderr, "Simulation failed.\n");
            return SCRIPT_COMMAND_FAILURE;
        }
    }
    script_finish_sim_or_fit(s);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_fit(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc != 1) {
        fprintf(stderr, "Usage: fit [fitvar1,fitvar2,...]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_params_free(fit_data->fit_params);
    fit_data->fit_params = fit_params_new();
    if(fit_params_add(fit_data->sim, fit_data->sm, fit_data->fit_params, argv[0])) {
        jabs_message(MSG_ERROR, stderr, "Could not add some fit parameters.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data->fit_params->n == 0) {
        jabs_message(MSG_ERROR, stderr, "No parameters for fit.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!fit_data->exp) { /* TODO: not enough to check this */
        jabs_message(MSG_ERROR, stderr, "No experimental spectrum set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data->n_fit_ranges == 0) {
        jabs_message(MSG_ERROR, stderr, "No fit range(s) given.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(script_prepare_sim_or_fit(s)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit(fit_data)) {
        jabs_message(MSG_ERROR, stderr, "Fit failed!\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    script_finish_sim_or_fit(s);
    jabs_message(MSG_INFO, stderr, "\nFinal parameters:\n");
    simulation_print(stderr, fit_data->sim);
    jabs_message(MSG_INFO, stderr, "\nFinal profile:\n");
    sample_print(NULL, fit_data->sim->sample, FALSE);
    sample_areal_densities_print(stderr, fit_data->sim->sample, FALSE);
    jabs_message(MSG_INFO, stderr, "\nFinal sample model:\n");
    sample_model_print(NULL, fit_data->sm);
    jabs_message(MSG_INFO, stderr, "\n");
    fit_stats_print(stderr, &fit_data->stats);
    return 0;
}

script_command_status script_save_spectra(script_session *s, int argc, char * const *argv) {
    size_t i_det = 0;
    struct fit_data *fit_data = s->fit;
    if(script_get_detector_number(fit_data->sim, &argc, &argv, &i_det) || argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save spectra [detector] file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(print_spectra(argv[0], fit_data_ws(fit_data, i_det), fit_data_exp(fit_data, i_det))) {
        jabs_message(MSG_ERROR, stderr, "Could not save spectra of detector %zu to file \"%s\"! There should be %zu detector(s).\n", i_det + 1, argv[0], fit_data->sim->n_det);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_save_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save sample [file]\n");
    }
    if(!fit_data->sm) {
        jabs_message(MSG_ERROR, stderr, "No sample set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sample_model_print(argv[0], fit_data->sm)) {
        jabs_message(MSG_ERROR, stderr, "Could not write sample to file \"%s\".\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_save_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit_data->sim, &argc, &argv, &i_det) || argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save detector [detector] file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(detector_print(argv[0], sim_det(fit_data->sim, i_det))) {
        jabs_message(MSG_ERROR, stderr,"Could not write detector %zu to file \"%s\".\n", i_det, argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_remove(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: remove [reaction] ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(strcmp(argv[0], "reaction") == 0) {
        if(argc != 3) {
            jabs_message(MSG_ERROR, stderr, "Usage: remove reaction [TYPE] [target_isotope]   OR   remove reaction number [number]\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        if(strcmp(argv[1], "number") == 0) {
            size_t i = strtoull(argv[2], NULL, 10);
            return sim_reactions_remove_reaction(fit_data->sim, i - 1);
        }
        reaction_type type = reaction_type_from_string(argv[1]);
        const jibal_isotope *target = jibal_isotope_find(fit_data->jibal->isotopes, argv[2], 0, 0);
        if(type == REACTION_NONE) {
            jabs_message(MSG_ERROR, stderr, "This is not a valid reaction type: \"%s\".\n", argv[1]);
            return SCRIPT_COMMAND_FAILURE;
        }
        if(!target) {
            jabs_message(MSG_ERROR, stderr, "This is not a valid isotope: \"%s\".\n", argv[2]);
            return SCRIPT_COMMAND_FAILURE;
        }
        for(size_t i = 0; i < fit_data->sim->n_reactions; i++) {
            if(fit_data->sim->reactions[i]->type == type && fit_data->sim->reactions[i]->target == target) {
                return sim_reactions_remove_reaction(fit_data->sim, i);
            }
        }
        jabs_message(MSG_ERROR, stderr, "No matching reaction found!\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_NOT_FOUND;
}

script_command_status script_roi(script_session *s, int argc, char * const *argv) {
    struct roi r;
    if(argc == 3) {
        size_t i = strtoul(argv[0], NULL, 10);
        if(i == 0) {
            jabs_message(MSG_ERROR, stderr, "Detector number must be > 0\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        if(i > s->fit->sim->n_det) {
            jabs_message(MSG_WARNING, stderr, "Warning: Detector number %zu > %zu.\n", i, s->fit->sim->n_det);
        }
        r.i_det = i - 1;
        argv++;
        argc--;
    }
    if (argc == 2) {
        r.i_det = 0;
        r.low = strtoul(argv[0], NULL, 10);
        r.high = strtoul(argv[1], NULL, 10);
    } else {
        jabs_message(MSG_ERROR, stderr, "Usage: roi [det] low high\n");
    }
    fit_data_roi_print(stderr, s->fit, &r);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_exit(script_session *s, int argc, char * const *argv) {
    return SCRIPT_COMMAND_EXIT;
}

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
    return s;
}
int script_session_reset_vars(script_session *s) {
    free(s->cf->vars);
    s->cf->vars = NULL;
    return jibal_config_file_set_vars(s->cf, script_make_vars(s)); /* Loading and resetting things can reset some pointers (like fit->det, so we need to update those to the vars */
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

const script_command *script_command_find(const script_command *commands, const char *cmd_string) {
    int found = 0;
    const script_command *c_found = NULL;
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
            found++;
            c_found = c;
            if(strlen(cmd_string) == strlen(c->name)) { /* Exact match, can not be ambiguous */
                found = 1;
                break;
            }
        }
    }
    if(found == 1) {
        return c_found;
    }
    if(found > 1) {
        jabs_message(MSG_ERROR, stderr, "\"%s\" is ambiguous (%i matches):", cmd_string, found);
        for(const struct script_command *c = commands; c->name != NULL; c++) {
            if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
                jabs_message(MSG_ERROR, stderr, " %s", c->name);
            }
        }
        jabs_message(MSG_ERROR, stderr, "\n");
        return NULL;
    }
    return NULL;
}

script_command_status script_process_command(script_session *s, const script_command *commands, int argc, char **argv) {
    if(!s || !commands || !argv)
        return SCRIPT_COMMAND_FAILURE;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "No arguments given.\n");
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    const script_command *cmds = commands;
    const script_command *c_parent = NULL;
    while(1) {
        if(argc && cmds) { /* Arguments and subcommands remain. Try to find the right one, if possible. */
            const script_command *c = script_command_find(cmds, argv[0]);
            if(c) { /* Subcommand found */
#ifdef DEBUG
                fprintf(stderr, "Debug: Found command %s. Going deeper.\n", c->name);
#endif
                cmds = c->subcommands;
                argc--;
                argv++;
                c_parent = c;
                continue;
            } else {
#ifdef DEBUG
                fprintf(stderr, "Debug: Didn't find subcommand %s.\n", argv[0]);
#endif
            }

        }
        if(c_parent) { /* This is a subcommand */
            if(c_parent->f) { /* Fallback function exists and should handle this case */
#ifdef DEBUG
                fprintf(stderr, "There is a function %p in command %s. Calling it with %i arguments.\n", (void *) c_parent->f, c_parent->name, argc);
#endif
                return c_parent->f(s, argc, argv);
            } else {
#ifdef DEBUG
                fprintf(stderr, "There is no function to call. Parent was %s.\n", c_parent->name);
#endif
                script_command_not_found(argv[0], c_parent);
            }
        } else {
            jabs_message(MSG_ERROR, stderr, "Command not recognized! Try 'help commands'.\n");
            return SCRIPT_COMMAND_NOT_FOUND;
        }
        return SCRIPT_COMMAND_NOT_FOUND;
    }
}

void script_command_not_found(const char *cmd, const script_command *parent) {
    if(parent && parent->subcommands) {
        if(cmd) {
            jabs_message(MSG_ERROR, stderr, "Sub-command \"%s\" not found!\n\n", cmd);
        } else {
            jabs_message(MSG_ERROR, stderr, "Not enough arguments!\n\n");
        }
        jabs_message(MSG_ERROR, stderr, "Following subcommands are recognized:\n");
        script_print_commands(stderr, parent->subcommands);
    } else {
        if(cmd) {
            jabs_message(MSG_ERROR, stderr, "\"%s\" not understood\n", cmd);
        }
    }
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
        jabs_message(MSG_INFO, stderr,"\nRunning script \"%s\"\n\n", filename);
    }
    script_command_status status = SCRIPT_COMMAND_SUCCESS;
    while(getline(&line, &line_size, f) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strip newlines */
        if(*line == '#') {/* Comment */
            continue;
        }
        if(!interactive) {
            jabs_message(MSG_INFO, stderr, "%s%s\n", prompt, line);
        }
        int argc = 0;
        char **argv = string_to_argv(line, &argc);
        if(!argv) {
            jabs_message(MSG_ERROR, stderr, "Something went wrong in parsing arguments.\n");
            continue;
        }
#ifdef DEBUG
        for(int i = 0; i < argc; i++) {
            fprintf(stderr, "args: %i: \"%s\"\n", i, argv[i]);
        }
#endif
        if(argc) {
            status = script_process_command(s, script_commands, argc, argv);
        }
        argv_free(argv, argc);
        if(status == SCRIPT_COMMAND_EXIT) {
            break;
        }
        if(interactive) {
            fputs(prompt, stderr);
        } else if(status != SCRIPT_COMMAND_SUCCESS) {
            jabs_message(MSG_ERROR, stderr, "Error (%i) on line %zu. Aborting.\n", status, lineno);
            break;
        }
    }
    free(line);
    fclose_file_or_stream(f);
    if(interactive) {
        jabs_message(MSG_INFO, stderr,"Bye.\n");
    } else if(filename) {
        if(status != SCRIPT_COMMAND_SUCCESS) {
            jabs_message(MSG_ERROR, stderr,"Error running script \"%s\"\n", filename);
        } else {
            jabs_message(MSG_ERROR, stderr, "Finished running script \"%s\"\n", filename);
        }
    }
    fflush(stderr);
    return status;
}

int script_prepare_sim_or_fit(script_session *s) {
    /* TODO: option to disable RBS or ERD */
    /* TODO: reactions from files? Should "sim_reactions_add" be somewhere else? */
    fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR, stderr,"No sample has been defined!\n");
        return -1;
    }
    if(!fit->sim->beam_isotope) {
        jabs_message(MSG_ERROR, stderr,"No ion has been defined!\n");
        return -1;
    }
    if(!fit->sim->det || fit->sim->n_det == 0) {
        jabs_message(MSG_ERROR, stderr,"No detector has been defined!\n");
        return -1;
    }
    fit_data_workspaces_free(s->fit);
    sample_free(fit->sim->sample);
#ifdef DEBUG
    fprintf(stderr, "Original sample model:\n");
    sample_model_print(NULL, fit->sm);
#endif
    fit->sim->sample = sample_from_sample_model(fit->sm);
    if(!fit->sim->sample) {
        jabs_message(MSG_ERROR, stderr, "Could not make a sample based on model description. This should never happen.\n");
        return -1;
    }
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_WARNING, stderr, "No reactions, adding some automatically. Please be aware there are commands called \"reset reactions\" and \"add reactions\".\n");
        if(fit->sim->rbs) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS)); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
        }
        if(fit->sim->erd) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
        }
    }
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_ERROR, stderr, "No reactions. Nothing to do.\n");
        return EXIT_FAILURE;
    }
    sim_sort_reactions(fit->sim);
    jabs_message(MSG_INFO, stderr, "Simplified sample model for simulation:\n");
    sample_print(NULL, fit->sim->sample, TRUE);

    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);

    jibal_gsto_assign_clear_all(fit->jibal->gsto); /* Is it necessary? No. Here? No. Does it clear old stuff? Yes. */
    if(assign_stopping(fit->jibal->gsto, fit->sim)) {
        jabs_message(MSG_ERROR, stderr,"Could not assign stopping or straggling. Failure. Provide more data, check that JIBAL Z2_max is sufficiently large (currently %i) or disable unwanted reactions (e.g. ERD).\n", s->jibal->config->Z_max);
        return -1;
    }
    jibal_gsto_print_assignments(fit->jibal->gsto);
    jibal_gsto_print_files(fit->jibal->gsto, TRUE);
    jabs_message(MSG_VERBOSE, stderr, "Loading stopping data.\n");
    jibal_gsto_load_all(fit->jibal->gsto);
    simulation_print(stderr, fit->sim);
    s->start = clock();
    return 0;
}

int script_finish_sim_or_fit(script_session *s) {
    s->end = clock();
    double cputime_total = (((double) (s->end - s->start)) / CLOCKS_PER_SEC);
    jabs_message(MSG_INFO, stderr, "...finished! Total CPU time: %.3lf s.\n", cputime_total);

    struct fit_data *fit = s->fit;

    if(fit->sim->n_det == 1) { /* TODO: multidetector automatic spectra saving! */
        size_t i_det = 0;
        sim_workspace *ws = fit_data_ws(fit, i_det);
        if(ws) {
            if(s->output_filename) {
                if(print_spectra(s->output_filename, ws, fit_data_exp(fit, i_det))) {
                    jabs_message(MSG_ERROR, stderr, "Could not save spectra of detector %zu to file \"%s\"\n", i_det, s->output_filename);
                    return EXIT_FAILURE;
                }
            }
            if(s->bricks_out_filename) {
                print_bricks(s->bricks_out_filename, ws);
            }
            if(s->detector_out_filename) {
                detector_print(s->detector_out_filename, ws->det);
            }
        }
    }
    if(s->sample_out_filename) {
        sample_model_print(s->sample_out_filename, fit->sm);
    }
    return 0;
}

int script_get_detector_number(const simulation *sim, int * const argc, char * const ** const argv, size_t *i_det) {
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
        return EXIT_SUCCESS; /* First argument was not a number, don't change i_det! */
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
