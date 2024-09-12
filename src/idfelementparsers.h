/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

    */

#ifndef IDFELEMENTPARSERS_H
#define IDFELEMENTPARSERS_H
#include "idfparse.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CALIB_PARAMS_MAX (10)

idf_error idf_parse_note(idf_parser *idf, xmlNode *note);
idf_error idf_parse_sample(idf_parser *idf, xmlNode *sample);
idf_error idf_parse_layerelement(idf_parser *idf, xmlNode *element);
idf_error idf_parse_layerelements(idf_parser *idf, xmlNode *elements);
idf_error idf_parse_layer(idf_parser *idf, xmlNode *layer);
idf_error idf_parse_layers(idf_parser *idf, xmlNode *layers);
idf_error idf_parse_simple_data(xmlNode *simple_data, const char *filename);
idf_error idf_parse_spectrum(idf_parser *idf, xmlNode *spectrum);
idf_error idf_parse_calibrations(idf_parser *idf, xmlNode *calibrations); /* inside spectrum */
idf_error idf_parse_beam(idf_parser *idf, xmlNode *beam); /* inside spectrum */
idf_error idf_parse_data(idf_parser *idf, xmlNode *data); /* inside spectrum */
idf_error idf_parse_geometry(idf_parser *idf, xmlNode *geometry); /* inside spectrum */
idf_error idf_parse_spectra(idf_parser *idf, xmlNode *spectra);
idf_error idf_parse_detector(idf_parser *idf, xmlNode *detector);
idf_error idf_parse_detectorshape(idf_parser *idf, xmlNode *detectorshape);
idf_error idf_parse_detector_foil(idf_parser *idf, xmlNode *stoppingfoil);
idf_error idf_parse_spot(idf_parser *idf, xmlNode *spot);
idf_error idf_detectorshape_or_spot(idf_parser *idf, xmlNode *node, const char *prefix); /* used by idf_parse_detectorshape and idf_parse_spot, since parsing is so similar */
idf_error idf_parse_process(idf_parser *idf, xmlNode *process);
idf_error idf_parse_physicsdefaults(idf_parser *idf, xmlNode *physicsdefaults);
idf_error idf_parse_simulations(idf_parser *idf, xmlNode *simulations);
idf_error idf_parse_simulation(idf_parser *idf, xmlNode *simulation);
#ifdef __cplusplus
}
#endif
#endif //IDFELEMENTPARSERS_H
