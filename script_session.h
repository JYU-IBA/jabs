/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SCRIPT_SESSION_H
#define JABS_SCRIPT_SESSION_H


#include <jibal.h>
#include "simulation.h"
#include "defaults.h"
#include "script_generic.h"

script_session *script_session_init(jibal *jibal, simulation *sim); /* sim can be NULL or a previously initialized sim can be given. Note that it will be free'd by script_session_free()! */
void script_session_free(script_session *s);
int script_session_load_script(script_session *s, const char *filename);
int script_get_detector_number(const simulation *sim, int allow_empty, int *argc, char * const **argv, size_t *i_det);

#endif //JABS_SCRIPT_SESSION_H
