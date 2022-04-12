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
    fprintf(stderr, "JaBS version %s. Copyright (C) 2021 - 2022 Jaakko Julin.\n", jabs_version()); /* These are printed when running non-interactively with just command line parameters */
    fprintf(stderr, "Compiled using JIBAL %s, current library version %s.\n", JIBAL_VERSION, jibal_version());
    fputs(COPYRIGHT_STRING, stderr);
    if(interactive) {
        fprintf(stderr, "Welcome to interactive mode.\nType \"help\" for help or run \"jabs -h\" for command line help.\n\n");
    }
}

void read_options(const jibal *jibal, simulation *sim, cmdline_options *cmd_opt,  int *argc, char * const **argv) {
    static struct option long_options[] = {
            {"help",          no_argument,       NULL, 'h'},
            {"version",       no_argument,       NULL, 'V'},
            {"interactive",   no_argument,       NULL, 'i'},
            {"verbose",       optional_argument, NULL, 'v'},
            {"out",           required_argument, NULL, 'o'},
            {"ion",           required_argument, NULL, 'I'},
            {"energy",        required_argument, NULL, 'E'},
            {"broad",         required_argument, NULL, 'B'},
            {"alpha",         required_argument, NULL, 'a'},
            {"theta",         required_argument, NULL, 't'},
            {"phi",           required_argument, NULL, 'p'},
            {"fluence",       required_argument, NULL, '3'},
            {"reaction",      required_argument, NULL, 'r'},
            {"resolution",    required_argument, NULL, 'R'},
            {"step_incident", required_argument, NULL, 2  },
            {"step_exiting",  required_argument, NULL, 1  },
            {"sample",        required_argument, NULL, 's'},
            {"sample_out",    required_argument, NULL, 'S'},
            {"slope",         required_argument, NULL, '1'},
            {"offset",        required_argument, NULL, '2'},
            {"compress",      required_argument, NULL, 'c'},
            {"fast",          optional_argument, NULL, 'f'},
            {"ds",            no_argument,       NULL,  0},
            {"exp",           required_argument, NULL, 'e'},
            {"fit",           no_argument,       NULL, 'F'},
            {"fit_low",       required_argument, NULL, '4'},
            {"fit_high",      required_argument, NULL, '5'},
            {"fit_vars",      required_argument, NULL, '6'},
            {"bricks_out",    required_argument, NULL, '7'},
            {"depthsteps",    required_argument, NULL, '8'},
            {"norbs",         no_argument,       NULL, '#'},
            {"noerd",         no_argument,       NULL, '9'},
            {"channeling",    required_argument, NULL, 'C'},
            {"channeling_slope",    required_argument, NULL, 3},
            {NULL, 0,                NULL,   0}
    };
    static const char *help_texts[] = {
            "Print this message.",
            "Print version number.",
            "Interactive mode. If script file(s) are given, they will be run first.",
            "Increase or give verbosity level.",
            "Output to file instead of standard output. For CSV use .csv suffix.",
            "Incident ion (without charge state), e.g. 4He",
            "Incident beam energy. Please give units as well, e.g. 2MeV or \"2 MeV\" or 2000keV",
            "Incident beam broadening (FWHM)",
            "Incident angle (from sample normal).",
            "Detector or scattering angle (from beam).",
            "Detector azimuthal angle (0 deg or 180 deg = IBM, +/-90 deg or 270 deg = Cornell).",
            "Fluence (or actually particles * sr).",
            "Reaction file to load (in R33 format).",
            "Resolution of detector (FHWM of Gaussian)",
            "Incident ion step size. Zero is automatic.",
            "Exiting particle step size.",
            "Sample file.",
            "Sample file (output).",
            "Slope of energy calibration.",
            "Offset of energy calibration.",
            "Compress channels in spectra by an integer factor.",
            "Make things faster, but worse.",
            "Dual scattering (broken!)",
            "Load experimental spectrum from file.",
            "Fit",
            "Fit range, low",
            "Fit range, high",
            "Comma separated list of parameters to fit, e.g. \"calib,fluence,thickness1\"",
            "Save intermediate raw data to file",
            "Maximum number of depth steps",
            "Don't make an RBS spectrum",
            "Don't make an ERD spectrum (ERD automatically turns on forward angles)",
            "Ad-hoc substrate channeling yield correction (constant)",
            "Ad-hoc substrate channeling yield correction (energy slope 1/keV)",
            NULL
    }; /* It is important to have the elements of this array correspond to the elements of the long_options[] array to avoid confusion. */
    while (1) {
        int option_index = 0;
        char c = getopt_long(*argc, *argv, "ihvVE:o:a:t:prob:I:r:R:s:fe:Fd:B:c:C:", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case 0:
                if(strcmp(long_options[option_index].name, "ds") == 0) {
                    cmd_opt->ds = TRUE;
                }
                break;
            case 'f':
                if (optarg)
                    cmd_opt->fast = atoi(optarg);
                else
                    cmd_opt->fast++;
                break;
            case 1:
                cmd_opt->stop_step_exiting = jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case '1':
                calibration_set_param(sim->det[0]->calibration, CALIBRATION_PARAM_SLOPE, jibal_get_val(jibal->units, UNIT_TYPE_ANY, optarg));
                break;
            case '2':
                calibration_set_param(sim->det[0]->calibration, CALIBRATION_PARAM_OFFSET, jibal_get_val(jibal->units, UNIT_TYPE_ANY, optarg));
                break;
            case 'c':
                sim->det[0]->compress = atoi(optarg);
                break;
            case 'C':
                sim->channeling_offset = strtod(optarg, NULL);
                break;
            case 3:
                sim->channeling_slope = strtod(optarg, NULL)/C_KEV;
                break;
            case '3':
                sim->fluence = strtod(optarg, NULL);
                break;
            case '4':
                cmd_opt->fit_low = atoi(optarg);
                break;
            case '5':
                cmd_opt->fit_high = atoi(optarg);
                break;
            case '6':
                cmd_opt->fit_vars = strdup(optarg);
                break;
            case '8':
                cmd_opt->depthsteps_max = atoi(optarg);
                break;
            case '9':
                cmd_opt->erd = 0;
                break;
            case '#':
                cmd_opt->rbs = 0;
                break;
            case 'F':
                cmd_opt->fit = 1;
                break;
            case 2:
                cmd_opt->stop_step_incident = jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 's':
                cmd_opt->sample_filename = strdup(optarg);
                break;
            case 'a':
                sim->sample_theta = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 't':
                sim->det[0]->theta = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 'p':
                sim->det[0]->phi = jibal_get_val(jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 'h':
                usage();
                fprintf(stderr, "\nThe following options (prefix with --, e.g. --energy=2MeV) are supported: \n");
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
                fprintf(stderr, "\nPlease note that many variables require units or they are assumed to be in relevant SI units.\n");
                exit(EXIT_SUCCESS);
                break;
            case 'V':
                printf("%s\n", jabs_version());
                exit(EXIT_SUCCESS);
                break; /* Unnecessary */
            case 'v':
                if (optarg)
                    cmd_opt->verbose = atoi(optarg);
                else
                    cmd_opt->verbose++;
                break;
            case 'i':
                cmd_opt->interactive = TRUE;
                break;
            case 'o':
                cmd_opt->output_filename = strdup(optarg);
                break;
            case 'e':
                cmd_opt->exp_filename = strdup(optarg);
                break;
            case 'E':
                sim->beam_E = jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'B':
                sim->beam_E_broad = pow2(jibal_get_val(jibal->units, UNIT_TYPE_ENERGY, optarg) / C_FWHM);
                break;
            case 'I':
                sim->beam_isotope = jibal_isotope_find(jibal->isotopes, optarg, 0, 0);
                if(!sim->beam_isotope) {
                    fprintf(stderr, "%s is not a valid isotope.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'r':
                cmd_opt->n_reaction_filenames++;
                cmd_opt->reaction_filenames=realloc(cmd_opt->reaction_filenames, cmd_opt->n_reaction_filenames * sizeof(char *));
                cmd_opt->reaction_filenames[cmd_opt->n_reaction_filenames - 1] = strdup(optarg);
                break;
            case 'R':
                calibration_set_param(sim->det[0]->calibration, CALIBRATION_PARAM_RESOLUTION, jibal_get_val(jibal->units, UNIT_TYPE_ANY, optarg));
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
                break;
        }
    }
    *argc -= optind;
    *argv += optind;
}

cmdline_options *cmdline_options_init() {
    cmdline_options *cmd_opt = malloc(sizeof(cmdline_options));
    memset(cmd_opt, 0, sizeof(cmdline_options)); /* Everything not listed below are zero or NULL by default */
    cmd_opt->rbs = TRUE;
    cmd_opt->erd = TRUE;
    return cmd_opt;
}

void cmdline_options_free(cmdline_options *cmd_opt) {
    free(cmd_opt->output_filename);
    free(cmd_opt->exp_filename);
    free(cmd_opt->fit_vars);
    free(cmd_opt->sample_filename);
    if(cmd_opt->reaction_filenames) {
        for(size_t i = 0; i < cmd_opt->n_reaction_filenames; i++) {
            free(cmd_opt->reaction_filenames[i]);
        }
    }
    free(cmd_opt->reaction_filenames);
}
