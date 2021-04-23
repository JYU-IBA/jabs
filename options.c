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

#include "version.h"
#include "options.h"



#define USAGE_STRING "Usage: jabs [-E <energy>] <material1> <thickness1> [<material2> <thickness2> ...]\n\nExample: jabs -E 2MeV --alpha=10deg --theta=170deg --out=spectrum.csv Au 500tfu SiO2 1000tfu Si 10000tfu\n"
#define COPYRIGHT_STRING "    Jaakko's Backscattering Simulator (JaBS)\n    Copyright (C) 2021 Jaakko Julin\n\n    This program is free software; you can redistribute it and/or modify \n    it under the terms of the GNU General Public License as published by\n    the Free Software Foundation; either version 2 of the License, or\n    (at your option) any later version.\n\n   See LICENSE.txt for the full license.\n\n"

const char *jabs_version() {
    return jabs_VERSION;
}

void usage() {
    fprintf(stderr, USAGE_STRING);
}

void read_options(global_options *global, simulation *sim, int *argc, char ***argv) {
    static struct option long_options[] = {
            {"help",          no_argument,       NULL, 'h'},
            {"version",       no_argument,       NULL, 'V'},
            {"verbose",       optional_argument, NULL, 'v'},
            {"out",           required_argument, NULL, 'o'},
            {"ion",           required_argument, NULL, 'I'},
            {"energy",        required_argument, NULL, 'E'},
            {"alpha",         required_argument, NULL, 'a'},
            {"theta",         required_argument, NULL, 't'},
            {"phi",           required_argument, NULL, 'p'},
            {"fluence",       required_argument, NULL, '3'},
            {"resolution",    required_argument, NULL, 'R'},
            {"step_incident", required_argument, NULL, 'S'},
            {"step_exiting",  required_argument, NULL, '0'},
            {"sample",        required_argument, NULL, 's'},
            {"detector",      required_argument, NULL, 'd'},
            {"detector_out",  required_argument, NULL, 'D'},
            {"slope",         required_argument, NULL, '1'},
            {"offset",        required_argument, NULL, '2'},
            {"compress",      required_argument, NULL, 'c'},
            {"fast",          optional_argument, NULL, 'f'},
            {"exp",           required_argument, NULL, 'e'},
            {"fit",           no_argument,       NULL, 'F'},
            {"fit_low",       required_argument, NULL, '4'},
            {"fit_high",      required_argument, NULL, '5'},
            {"fit_vars",      required_argument, NULL, '6'},
            {"bricks_out",    required_argument, NULL, '7'},
            {"depthsteps",    required_argument, NULL, '8'},
            {"norbs",         no_argument,       NULL, '#'},
            {"noerd",         no_argument,       NULL, '9'},
            {"isotopes",      no_argument,       NULL, 0},
            {"channeling",    required_argument, NULL, 'C'},
            {NULL, 0,                NULL,   0}
    };
    static const char *help_texts[] = {
            "Print this message.",
            "Print version number.",
            "Increase or give verbosity level.",
            "Output to file instead of standard output. For CSV use .csv suffix.",
            "Incident ion (without charge state), e.g. 4He",
            "Incident beam energy. Please give units as well, e.g. 2MeV or \"2 MeV\" or 2000keV",
            "Incident angle (from sample normal).",
            "Detector or scattering angle (from beam).",
            "Detector azimuthal angle (0 deg or 180 deg = IBM, +/-90 deg or 270 deg = Cornell).",
            "Fluence (or actually particles * sr).",
            "Resolution of detector (FHWM of Gaussian)",
            "Incident ion step size. Zero is automatic.",
            "Exiting particle step size.",
            "Sample file.",
            "Detector file.",
            "Detector file (output).",
            "Slope of energy calibration.",
            "Offset of energy calibration.",
            "Compress channels in spectra by an integer factor.",
            "Make things faster, but worse.",
            "Load experimental spectrum from file.",
            "Fit",
            "Fit range, low",
            "Fit range, high",
            "Comma separated list of parameters to fit, e.g. \"calib,fluence,thickness1\"",
            "Save intermediate raw data to file",
            "Maximum number of depth steps",
            "Don't make an RBS spectrum",
            "Don't make an ERD spectrum (ERD automatically turns on forward angles)",
            "Print isotopes (in concentration table)",
            "Ad-hoc substrate channeling yield correction",
            NULL
    }; /* It is important to have the elements of this array correspond to the elements of the long_options[] array to avoid confusion. */
    while (1) {
        int option_index = 0;
        char c = getopt_long(*argc, *argv, "hvVE:o:a:t:prob:I:R:S:s:fe:Fd:D:c:C:", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case 0:
                if(strcmp(long_options[option_index].name, "isotopes") == 0) {
                    global->print_isotopes = TRUE;
                }
            case 'f':
                if (optarg)
                    sim->fast = atoi(optarg);
                else
                    sim->fast++;
                break;
            case '0':
                sim->stop_step_exiting = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case '1':
                sim->det.slope = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case '2':
                sim->det.offset = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'c':
                sim->det.compress = atoi(optarg);
                break;
            case 'C':
                sim->channeling = strtod(optarg, NULL);
                break;
            case '3':
                sim->p_sr = strtod(optarg, NULL);
                break;
            case '4':
                global->fit_low = atoi(optarg);
                break;
            case '5':
                global->fit_high = atoi(optarg);
                break;
            case '6':
                global->fit_vars = optarg;
                break;
            case '7':
                global->bricks_filename = optarg;
                break;
            case '8':
                sim->depthsteps_max = atoi(optarg);
                break;
            case '9':
                global->erd = 0;
                break;
            case '#':
                global->rbs = 0;
                break;
            case 'F':
                global->fit = 1;
                break;
            case 'S':
                sim->stop_step_incident = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 's':
                global->sample_filename = optarg;
                break;
            case 'a':
                sim->sample_theta = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 't':
                sim->det.theta = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 'p':
                sim->det.phi = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 'h':
                fputs(COPYRIGHT_STRING, stderr);
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
                    global->verbose = atoi(optarg);
                else
                    global->verbose++;
                break;
            case 'o':
                global->out_filename = optarg;
                break;
            case 'e':
                global->exp_filename = optarg;
                break;
            case 'E':
                sim->beam_E = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'I':
                sim->beam_isotope = jibal_isotope_find(global->jibal->isotopes, optarg, 0, 0);
                if(!sim->beam_isotope) {
                    fprintf(stderr, "%s is not a valid isotope.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'R':
                sim->det.resolution = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg)/C_FWHM;
                sim->det.resolution *= sim->det.resolution; /* square */
                break;
            case 'd':
                sim->det = detector_from_file(global->jibal->units, optarg);
                break;
            case 'D':
                global->detector_out_filename = optarg;
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
