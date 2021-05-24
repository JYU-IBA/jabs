#include <string.h>
#include "defaults.h"
#include "script.h"
#include "fit.h"
#include "generic.h"

static const struct script_command commands[] = {
        {"help", &script_help, "Print help."},
        {"show", &script_show, "Show information on things."},
        {"reset", &script_show, "Reset something."},
        {"exit", NULL, "Exit."},
        {"quit", NULL, NULL},
        {NULL, NULL, NULL}
};

void script_print_commands(FILE *f) {
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        if(!c->help_text)
            continue;
        fprintf(f, "%22s    %s\n", c->name, c->help_text);
    }
}

int script_reset(struct fit_data *fit, int argc, char * const *argv) {
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

int script_show(struct fit_data *fit, int argc, char * const *argv) {
    if(argc == 0) {
        fprintf(stderr, "Nothing to show.\n");
        return 0;
    }
     if(strcmp(argv[0], "sim") == 0) {
         simulation_print(stderr, fit->sim);
         return 0;
     } else {
         fprintf(stderr, "Don't know what \"%s\" is.\n", argv[0]);
         return -1;
     }
    return 0;
}


int script_help(struct fit_data *fit, int argc, char * const *argv) {
    (void) fit; /* Unused */
    static const struct help_topic topics[] = {
            {"help", "This is help on help. How meta. Help is available on following topics:\n"},
            {"commands", "I recognize the following commands:\n"},
            {"version", "JaBS version: "},
            {NULL, NULL}
    };
    if(argc == 0) {
        fprintf(stderr, "Type help [topic] for information on a particular topic.\n");
        return 0;
    }
    for(const struct help_topic *t = topics; t->name != NULL; t++) {
        if(strcmp(t->name, argv[0]) == 0) {
            fputs(t->help_text, stderr);
            if(strcmp(t->name, "help") == 0) {
                size_t i = 0;
                for(const struct help_topic *t2 = topics; t2->name != NULL; t2++) {
                    fprintf(stderr, "%16s", t2->name);
                    if(i % 4 != 0) {
                        fputc('\n', stderr);
                    }
                }
                fprintf(stderr, "\n");
            } else if(strcmp(t->name, "commands") == 0) {
                script_print_commands(stderr);
            } else if(strcmp(t->name, "version") == 0) {
                fprintf(stderr, "%s\n", jabs_version());
            }
            return 0;
        }
    }
    fprintf(stderr, "Sorry, no help for '%s'.\n", argv[0]);
    return -1;
}


int script_process(jibal *jibal, FILE *f) {
    struct fit_data *fit = fit_data_new(jibal, sim_init(), NULL, NULL, NULL, NULL, 0, 0, 0); /* Not just fit, but this conveniently holds everything we need. */
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
                    status = c->f(fit, argc - 1, argv + 1);
                    break;
                }
            }
            if(!found) {
                if(interactive) {
                    fprintf(stderr, "Command \"%s\" not recognized. See 'help commands'.\n", line);
                } else {
                    fprintf(stderr, "Command \"%s\" not recognized on line %zu\n", line, lineno);
                    return EXIT_FAILURE;
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
    if(interactive) {
        fprintf(stderr, "Bye.\n");
    }
    return EXIT_SUCCESS;
}
