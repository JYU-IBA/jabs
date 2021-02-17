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



#define USAGE_STRING "Usage: jabs [-E <energy>] <material1> <thickness1> [<material2> <thickness2> ...]\n\nExample: jabs -E 2MeV --alpha 10deg --beta 0deg -theta 170deg Au 500tfu SiO2 1000tfu Si 10000tfu\n"
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
            {"fluence",   required_argument, NULL, 'F'},
            {"resolution",required_argument, NULL, 'R'},
            {"step_incident",required_argument, NULL, 'S'},
            {"step_exiting",required_argument, NULL, '0'},
            {"fast", optional_argument, NULL, 'f'},
            {NULL, 0,                NULL,   0}
    };
    while (1) {
        int option_index = 0;
        char c = getopt_long(*argc, *argv, "hvVE:o:a:b:t:I:F:R:S:f:", long_options, &option_index);
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
            case 'F':
                sim->p_sr = strtod(optarg, NULL);
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
