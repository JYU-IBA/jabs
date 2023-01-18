#ifndef IDFELEMENTPARSERS_H
#define IDFELEMENTPARSERS_H
#include "idfparse.h"

#define CALIB_PARAMS_MAX (10)

idf_error idf_parse_sample(idf_parser *idf, xmlNode *sample);
idf_error idf_parse_layerelement(idf_parser *idf, xmlNode *elements);
idf_error idf_parse_layerelements(idf_parser *idf, xmlNode *elements);
idf_error idf_parse_layer(idf_parser *idf, xmlNode *layer);
idf_error idf_parse_layers(idf_parser *idf, xmlNode *layers);
idf_error idf_parse_spectrum(idf_parser *idf, xmlNode *spectrum);
idf_error idf_parse_spectra(idf_parser *idf, xmlNode *spectra);
idf_error idf_parse_detector(idf_parser *idf, xmlNode *spectra);
#endif //IDFELEMENTPARSERS_H
