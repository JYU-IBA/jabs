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

#include "version.h"
#include "options.h"



#define USAGE_STRING "Usage: jabs [-E <energy>] <material1> <thickness1> [<material2> <thickness2> ...]\n\nExample: jabs -E 2MeV --alpha=10deg --beta=0deg --theta=170deg --out=spectrum.csv Au 500tfu SiO2 1000tfu Si 10000tfu\n"
#define COPYRIGHT_STRING "    Jaakko's Backscattering Simulator (JaBS)\n    Copyright (C) 2021 Jaakko Julin\n\n    This program is free software; you can redistribute it and/or modify \n    it under the terms of the GNU General Public License as published by\n    the Free Software Foundation; either version 2 of the License, or\n    (at your option) any later version.\n\n   See LICENSE.txt for the full license.\n\n"

const char *jabs_version() {
    return jabs_VERSION;
}

void usage() {
    fprintf(stderr, USAGE_STRING);
}

void read_options(global_options *global, simulation *sim, int *argc, char ***argv) {
    static struct option long_options[] = {
            {"help",      no_argument,       NULL, 'h'},
            {"version",   no_argument,       NULL, 'V'},
            {"verbose",   optional_argument, NULL, 'v'},
            {"out",       required_argument, NULL, 'o'},
            {"ion",       required_argument, NULL, 'I'},
            {"energy",    required_argument, NULL, 'E'},
            {"alpha",     required_argument, NULL, 'a'},
            {"beta",      required_argument, NULL, 'b'},
            {"theta",     required_argument, NULL, 't'},
            {"fluence",   required_argument, NULL, '3'},
            {"resolution",required_argument, NULL, 'R'},
            {"step_incident",required_argument, NULL, 'S'},
            {"step_exiting",required_argument, NULL, '0'},
            {"slope",     required_argument, NULL, '1'},
            {"offset",    required_argument, NULL, '2'},
            {"fast",      optional_argument, NULL, 'f'},
            {"exp",       required_argument, NULL, 'e'},
            {"fit",        no_argument,      NULL, 'F'},
            {"fit_low",   required_argument, NULL, '4'},
            {"fit_high",   required_argument, NULL, '5'},
            {"fit_vars",   required_argument, NULL, '6'},
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
            "Exit angle (from sample norma).",
            "Scattering angle (from beam).",
            "Fluence (or actually particles * sr).",
            "Resolution of detector (FHWM of Gaussian)",
            "Incident ion step size.",
            "Exiting particle step size.",
            "Slope of energy calibration.",
            "Offset of energy calibration.",
            "Make things faster, but worse.",
            "Load experimental spectrum from file.",
            "Fit",
            "Fit range, low",
            "Fit range, high",
            "Comma separated list of parameters to fit, e.g. \"calib,fluence,thickness1\"",
            NULL
    }; /* It is important to have the elements of this array correspond to the elements of the long_options[] array to avoid confusion. */
    while (1) {
        int option_index = 0;
        char c = getopt_long(*argc, *argv, "hvVE:o:a:b:t:I:R:S:fe:F", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
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
                sim->energy_slope = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case '2':
                sim->energy_offset = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
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
            case 'F':
                global->fit = 1;
                break;
            case 'S':
                sim->stop_step_incident = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'a':
                sim->alpha = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 'b':
                sim->beta = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
                break;
            case 't':
                sim->theta = jibal_get_val(global->jibal->units, UNIT_TYPE_ANGLE, optarg);
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
                sim->ion.E = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg);
                break;
            case 'I':
                ion_set_isotope(&sim->ion, jibal_isotope_find(global->jibal->isotopes, optarg, 0, 0));
                if(!sim->ion.isotope) {
                    fprintf(stderr, "%s is not a valid isotope.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'R':
                sim->energy_resolution = jibal_get_val(global->jibal->units, UNIT_TYPE_ENERGY, optarg)/C_FWHM;
                sim->energy_resolution *= sim->energy_resolution; /* square */
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