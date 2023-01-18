#ifndef IDFELEMENTPARSERS_H
#define IDFELEMENTPARSERS_H
#include "idfparse.h"

#define CALIB_PARAMS_MAX (10)

int idf_parse_sample(idfparser *idf, xmlNode *sample);
void idf_parse_layerelements(idfparser *idf, xmlNode *elements);
int idf_parse_layer(idfparser *idf, xmlNode *layer);
int idf_parse_layers(idfparser *idf, xmlNode *layers);
int idf_parse_spectrum(idfparser *idf, xmlNode *spectrum);
int idf_parse_spectra(idfparser *idf, xmlNode *spectra);
#endif //IDFELEMENTPARSERS_H
