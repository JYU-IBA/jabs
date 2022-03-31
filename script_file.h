/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <stdio.h>
#include "win_compat.h"
#include "script_generic.h"

script_file *script_file_open(const char *filename);
void script_file_close(script_file *sfile);

ssize_t script_file_getline(script_file *sfile);
