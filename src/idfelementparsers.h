#ifndef IDFELEMENTPARSERS_H
#define IDFELEMENTPARSERS_H
#include "idfparse.h"

#define CALIB_PARAMS_MAX (10)

idf_error idf_parse_note(idf_parser *idf, xmlNode *note);
idf_error idf_parse_sample(idf_parser *idf, xmlNode *sample);
idf_error idf_parse_layerelement(idf_parser *idf, xmlNode *element);
idf_error idf_parse_layerelements(idf_parser *idf, xmlNode *elements);
idf_error idf_parse_layer(idf_parser *idf, xmlNode *layer);
idf_error idf_parse_layers(idf_parser *idf, xmlNode *layers);
idf_error idf_parse_spectrum(idf_parser *idf, xmlNode *spectrum);
idf_error idf_parse_spectra(idf_parser *idf, xmlNode *spectra);
idf_error idf_parse_detector(idf_parser *idf, xmlNode *detector);
idf_error idf_parse_detectorshape(idf_parser *idf, xmlNode *detectorshape);
idf_error idf_parse_spot(idf_parser *idf, xmlNode *spot);
idf_error idf_detectorshape_or_spot(idf_parser *idf, xmlNode *node, const char *prefix); /* used by idf_parse_detectorshape and idf_parse_spot, since parsing is so similar */
idf_error idf_parse_process(idf_parser *idf, xmlNode *process);
idf_error idf_parse_physicsdefaults(idf_parser *idf, xmlNode *physicsdefaults);
#endif //IDFELEMENTPARSERS_H
