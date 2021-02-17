/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include "layers.h"
#include <stdlib.h>

jibal_layer **read_layers(jibal *jibal, int argc, char **argv, size_t *n_layers) {
    size_t n = 0; /* How many layers have been allocated */
    jibal_layer **layers = NULL;
    *n_layers=0;
    while (argc >= 2) {
        if(*n_layers == n) {
            if(n == 0) {
                n = LAYERS_INITIAL_ALLOC;
            } else {
                n *= 2;
            }
#ifdef DEBUG
            fprintf(stderr, "(Re)allocating space for up to %lu layers.\n", n);
#endif
            layers = realloc(layers, n*sizeof(jibal_layer *));
            if(!layers)
                return NULL;
        }
        jibal_layer *layer = jibal_layer_new(jibal_material_create(jibal->elements, argv[0]),
                                             jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, argv[1]));
        if (!layer) {
            fprintf(stderr, "Not a valid layer: %s!\n", argv[0]);
            free(layers);
            return NULL;
        }

        layers[*n_layers] = layer;
        argc -= 2;
        argv += 2;
        (*n_layers)++;
    }
    return layers;
}

void layers_free(jibal_layer **layers, int n_layers) {
    int i;
    for(i = 0; i < n_layers; i++) {
        jibal_layer_free(layers[i]);
    }
    free(layers);
}
