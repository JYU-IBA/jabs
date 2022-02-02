#ifndef JABS_SCRIPT_GENERIC_H
#define JABS_SCRIPT_GENERIC_H

#include <time.h>
#include <jibal.h>
#include "defaults.h"

struct script_command;

typedef struct script_file {
    FILE *f;
    char *filename;
    char *line;
    size_t line_size;
    size_t lineno;
    int interactive;
} script_file;

typedef struct script_session {
    jibal *jibal;
    struct fit_data *fit;
    char *output_filename; /* File name for automatic spectra saving */ /* TODO: multidetector! */
    clock_t start, end;
    script_file *files[SCRIPT_NESTED_MAX];
    size_t file_depth;
    struct script_command *commands;
    size_t i_det_active;
} script_session;

#define SCRIPT_COMMAND_SUCCESS (0) /* Anything above and including zero is success */
#define SCRIPT_COMMAND_FAILURE (-1)
#define SCRIPT_COMMAND_NOT_FOUND (-2)
#define SCRIPT_COMMAND_EXIT (-3)
#define SCRIPT_COMMAND_EOF (-4)

typedef int script_command_status; /* Script commands should return negative on error (see defines above) and number of arguments (zero or positive) consumed successfully */

typedef struct script_command {
    char *name;
    script_command_status (*f)(struct script_session *, int, char * const *); /* Function to process argument vectors. Called if it is non-NULL, even if subcommands exist! Return value is important.*/
    script_command_status (*f_var)(struct script_session *, jibal_config_var *var, int, char * const *); /* Function to process variables (from subcommands). */
    script_command_status (*f_val)(struct script_session *, int, int, char * const *); /* Function to process vals (from subcommands). */
    char *help_text; /* Short help text */
    struct script_command *subcommands;
    jibal_config_var *var;
    int val;
    struct script_command *next;
} script_command;

struct help_topic {
    const char *name;
    const char *help_text;
};

#endif // JABS_SCRIPT_GENERIC_H
