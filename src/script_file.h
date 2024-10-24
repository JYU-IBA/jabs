/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_SCRIPT_FILE_H
#define JABS_SCRIPT_FILE_H
#include <stdio.h>
#include "win_compat.h"
#include "script_generic.h"
#ifdef __cplusplus
extern "C" {
#endif

script_file *script_file_open(const char *filename);
void script_file_close(script_file *sfile);

ssize_t script_file_getline_no_rl(script_file *sfile); /* Implementation that does not use readline. */
ssize_t script_file_getline(script_file *sfile);
#ifdef __cplusplus
}
#endif
#endif
