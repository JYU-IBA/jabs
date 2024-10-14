/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_IDF2JBS_H
#define JABS_IDF2JBS_H
#include "idfparse.h"
#ifdef __cplusplus
extern "C" {
#endif
idf_error idf_parse(const char *filename, char **filename_out);
/* Reads an IDF file, parses it and saves output as an JaBS script file (automatic file name stored in newly allocated filename_out if it is not NULL). */
#ifdef __cplusplus
}
#endif
#endif
