/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_SIMULATION2IDF_H
#define JABS_SIMULATION2IDF_H
#include "simulation.h"
#include "fit.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct idf_writer {
        xmlDocPtr doc;
        xmlNodePtr root;
        xmlNsPtr ns; /* IDF namespace (default) */
        xmlNsPtr simnra; /* SIMNRA namespace. Set to NULL if SIMNRA extensions are not to be added. */
} idf_writer;

int simulation2idf(const fit_data *fit, const char *filename, int simnra_extensions);
xmlNodePtr simulation2idf_notes(void);
xmlNodePtr simulation2idf_attributes(const char *filename);
xmlNodePtr simulation2idf_elementsandmolecules(const sample_model *sm);
xmlNodePtr simulation2idf_structure(const idf_writer *idf, const sample_model *sm);
xmlNodePtr simulation2idf_layers(const idf_writer *idf, const sample_model *sm);
xmlNodePtr simulation2idf_spectra(const idf_writer *idf, const fit_data *fit);
xmlNodePtr simulation2idf_beam(const simulation *sim);
xmlNodePtr simulation2idf_geometry(const simulation *sim, const detector *det);
xmlNodePtr simulation2idf_detection(const idf_writer *idf, const detector *det);
xmlNodePtr simulation2idf_calibrations(const jibal_element *elements, const simulation *sim, const detector *det);
xmlNodePtr simulation2idf_detectorresolution(const detector *det, const calibration *cal, const char  *ion_name);
xmlNodePtr simulation2idf_energycalibration(const detector *det, const calibration *cal, const char *ion_name);
xmlNodePtr simulation2idf_aperture(const char *name, const aperture *aperture);
xmlNodePtr simulation2idf_reactions(const simulation *sim, const detector *det);
xmlNodePtr simulation2idf_process(const simulation *sim, const result_spectra *spectra);
xmlNodePtr simulation2idf_data(const jabs_histogram *h);
xmlNodePtr simulation2idf_simpledata(const jabs_histogram *h);

#ifdef __cplusplus
}
#endif
#endif
