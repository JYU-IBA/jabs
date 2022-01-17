/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <jibal.h>
#include "simulation.h"
#include "script_file.h"
#include "defaults.h"

typedef struct script_session {
    jibal *jibal;
    struct fit_data *fit;
    jibal_config_file *cf; /* a pseudo-configuration "file" for our session */
    char *output_filename; /* File name for automatic spectra saving */ /* TODO: multidetector! */
    char *bricks_out_filename; /* File name for automatic bricks saving */
    char *sample_out_filename; /* File name for automatic sample saving */
    char *detector_out_filename; /* File name for automatic detector saving */
    clock_t start, end;
    script_file *files[SCRIPT_NESTED_MAX];
    size_t file_depth;
} script_session;

script_session *script_session_init(jibal *jibal, simulation *sim); /* sim can be NULL or a previously initialized sim can be given. Note that it will be free'd by script_session_free()! */
jibal_config_var *script_make_vars(script_session *s);
int script_session_reset_vars(script_session *s);
void script_session_free(script_session *s);
int script_session_load_script(script_session *s, const char *filename);
int script_get_detector_number(const simulation *sim, int *argc, char *const **argv, size_t *i_det);