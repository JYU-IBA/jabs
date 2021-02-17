/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */


#ifndef JABS_LAYERS_H
#define JABS_LAYERS_H
#include <jibal.h>

#define LAYERS_INITIAL_ALLOC 8

jibal_layer **read_layers(jibal *jibal, int argc, char **argv, size_t *n_layers);
void layers_free(jibal_layer **layers, int n_layers);
#endif //JABS_LAYERS_H
