#include <string.h>
#include "defaults.h"
#include "script.h"


int script_process(jibal *jibal, FILE *f) {
    char *line=NULL;
    size_t line_size=0;
    size_t lineno=0;
    int interactive = (f == stdin);
    const char *prompt = "\njabs> ";
    if(interactive) {
        fputs(COPYRIGHT_STRING, stderr);
        fprintf(stderr, "Welcome to interactive mode.\nType \"help\" for help or run \"jabs -h\" for command line help.\n");
        fputs(prompt, stderr);
    }
    while(getline(&line, &line_size, f) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strip newlines */
        if(*line == '#') /* Comment */
            continue;
        if(strcmp(line, "exit") == 0 || strcmp(line, "quit") == 0) {
            fprintf(stderr, "Exiting due to user request after %zu lines processed.\n", lineno);
            break;
        }
        if(strcmp(line, "help") == 0) {
            fprintf(stderr, "No help available yet.\n");
            continue;
        }
        fprintf(stderr, "Command \"%s\" not recognized. Type \"exit\" to exit or \"help\" for help.\n", line);
        if(interactive) {
            fputs(prompt, stderr);
        }
    }
    free(line);
    return EXIT_SUCCESS;
}
