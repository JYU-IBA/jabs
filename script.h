/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_SCRIPT_H
#define JABS_SCRIPT_H

#include "script_generic.h"
int script_process(script_session *s, const char *filename);
#endif // JABS_SCRIPT_H
