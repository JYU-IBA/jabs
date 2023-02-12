/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include <jibal_defaults.h>

#include "version.h"
#include "defaults.h"
#include "options.h"
#include "git.h"

#define USAGE_STRING "Usage: jabs [OPTION [<argument>]] [OPTION2 ...] ... [<FILE> [<FILE2>] ...] | sample <material1> <thickness1> [<material2> <thickness2> ...]]\n\nRun without arguments or \"-i\" to use JaBS interactively, give a sample or script file(s) otherwise.\n\nExample: jabs -E 2MeV --alpha=10deg --theta=170deg --out=spectrum.csv sample Au 500tfu SiO2 1000tfu Si 10000tfu\nExample: jabs example.jbs\n"

const char *jabs_version() {
    if(git_populated()) {
        return git_describe();
    }
    return jabs_version_simple();
}

const char *jabs_version_simple() {
    return jabs_VERSION;
}

void usage() {
    fprintf(stderr, USAGE_STRING);
}

void greeting(int interactive) {
    fprintf(stderr, "JaBS version %s. Copyright (C) 2021 - 2023 Jaakko Julin.\n", jabs_version()); /* These are printed when running non-interactively with just command line parameters */
    fprintf(stderr, "Compiled using JIBAL %s, current library version %s.\n", JIBAL_VERSION, jibal_version());
    fputs(COPYRIGHT_STRING, stderr);
    if(interactive) {
        fprintf(stderr, "Welcome to interactive mode.\nType \"help\" for help or run \"jabs -h\" for command line help.\n\n");
    }
}

void read_options(cmdline_options *cmd_opt, int *argc, char *const **argv) {
    static struct option long_options[] = {
            {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'V'},
            {"interactive", no_argument,       NULL, 'i'},
            {"verbose",     optional_argument, NULL, 'v'},
            {NULL, 0,                          NULL, 0}
    };
    static const char *help_texts[] = {
            "Print this message.",
            "Print version number.",
            "Interactive mode. If script file(s) are given, they will be run first.",
            "Increase or give verbosity level.",
            NULL
    }; /* It is important to have the elements of this array correspond to the elements of the long_options[] array to avoid confusion. */
    while (1) {
        int option_index = 0;
        int c = getopt_long(*argc, *argv, "ihvV", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case 'h':
                usage();
                fprintf(stderr, "\nThe following options (prefix with --) are supported: \n");
                struct option *o = long_options;
                const char **h = help_texts;
                while (o->name != NULL) {
                    fprintf(stderr, "  %14s  ", o->name);
                    if (*h != NULL) {
                        fprintf(stderr, "%s\n", *h);
                        h++;
                    } else {
                        fprintf(stderr, "\n");
                    }
                    o++;
                }
                exit(EXIT_SUCCESS);
            case 'V':
                printf("%s\n", jabs_version());
                exit(EXIT_SUCCESS);
            case 'v':
                if (optarg)
                    cmd_opt->verbose = atoi(optarg);
                else
                    cmd_opt->verbose++;
                break;
            case 'i':
                cmd_opt->interactive = TRUE;
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }
    *argc -= optind;
    *argv += optind;
}

cmdline_options *cmdline_options_init() {
    cmdline_options *cmd_opt = malloc(sizeof(cmdline_options));
    memset(cmd_opt, 0, sizeof(cmdline_options)); /* Everything not listed below are zero or NULL by default */
    return cmd_opt;
}

void cmdline_options_free(cmdline_options *cmd_opt) {
    free(cmd_opt);
}
