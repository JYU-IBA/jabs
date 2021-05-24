#include <string.h>
#include "defaults.h"
#include "script.h"
#include "fit.h"
#include "generic.h"

struct script_command {
    const char *name;
    int (*f)(struct fit_data *, int, char * const *);
    const char *help_text;
};

void script_print_commands(FILE *f, const struct script_command *commands) {
    fprintf(f, "I recognize the following commands: \n");
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        if(!c->help_text)
            continue;
        fprintf(f, "%22s    %s\n", c->name, c->help_text);
    }
}

int script_help(struct fit_data *fit, int argc, char * const *argv) {
    (void) fit; /* Unused */
    fprintf(stderr, "HEEEEEeeeeeeeellllppp!\n");
    return 0;
}

int script_process(jibal *jibal, FILE *f) {
    static const struct script_command commands[] = {
            {"help", &script_help, "Print help."},
            {"exit", NULL, "Exit."},
            {"quit", NULL, NULL},
            {NULL, NULL, NULL}
    };
    struct fit_data *fit = fit_data_new(jibal, sim_init(), NULL, NULL, NULL, NULL, 0, 0, 0); /* Not just fit, but this conveniently holds everything we need. */
    char *line=NULL, *arguments;
    size_t line_size=0;
    size_t lineno=0;
    int interactive = (f == stdin);
    const char *prompt = "jabs> ";
    if(interactive) {
        fputs(COPYRIGHT_STRING, stderr);
        fprintf(stderr, "Welcome to interactive mode.\nType \"help\" for help or run \"jabs -h\" for command line help.\n");
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
                    fprintf(stderr, "Command \"%s\" not recognized. Type \"exit\" to exit or \"help\" for help.\n", line);
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
    return EXIT_SUCCESS;
}
