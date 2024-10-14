/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <string.h>
#include <libxml/parser.h>
#include "defaults.h"
#include "options.h"
#include "geostragg.h"
#include "idfparse.h"
#include "jabs_debug.h"
#include "simulation2idf.h"


int simulation2idf(const struct fit_data *fit, const char *filename) {
    if(!fit->sim) {
        return EXIT_FAILURE;
    }
    FILE *f = fopen(filename, "w");
    if(!f) {
        return EXIT_FAILURE;
    }
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
    xmlNodePtr n = xmlNewNode(NULL, BAD_CAST "idf");
    xmlSetProp(n, BAD_CAST "xmlns:xsi", BAD_CAST "http://www.w3.org/2001/XMLSchema-instance");
    xmlSetProp(n, BAD_CAST "xmlns", BAD_CAST "http://idf.schemas.itn.pt");
    xmlDocSetRootElement(doc, n);

    xmlNodePtr notes = simulation2idf_notes();
    if(notes) {
        xmlAddChild(n, notes);
    }

    xmlNodePtr attributes = simulation2idf_attributes(filename);
    if(attributes) {
        xmlAddChild(n, attributes);
    }

    xmlNodePtr sample = xmlNewChild(n, NULL, BAD_CAST "sample", NULL);

    sample_model *sm2 = sample_model_split_elements(fit->sm); /* TODO: this splits molecules. We might not necessarily want it always, but it is a reasonable default. */
    xmlNodePtr elementsandmolecules = simulation2idf_elementsandmolecules(sm2);
    xmlNodePtr structure = simulation2idf_structure(sm2);
    if(elementsandmolecules && structure) {
        xmlAddChild(sample, elementsandmolecules);
        xmlAddChild(sample, structure);
    }
    sample_model_free(sm2);

    xmlNodePtr spectra = simulation2idf_spectra(fit);
    if(spectra) {
        xmlAddChild(sample, spectra);
    }

    xmlChar *xmlbuff;
    int buffersize;
    xmlDocDumpFormatMemory(doc, &xmlbuff, &buffersize, 1);
    fprintf(f, "%s", (char *) xmlbuff);
    fclose(f);
    xmlFreeDoc(doc);
    return EXIT_SUCCESS;
}

xmlNodePtr simulation2idf_notes(void) {
    xmlNodePtr notes = xmlNewNode(NULL, BAD_CAST "notes");
    xmlNodePtr note = xmlNewChild(notes, NULL, BAD_CAST "note", NULL);

    char *jabs_note = NULL;
    if(asprintf(&jabs_note, "File created by JaBS %s", jabs_version()) > 0) {
        xmlNodeSetContent(note, BAD_CAST jabs_note);
        free(jabs_note);
    }
    return notes;
}

xmlNodePtr simulation2idf_attributes(const char *filename) {
    xmlNodePtr attributes = xmlNewNode(NULL, BAD_CAST "attributes");
    xmlNewChild(attributes, NULL, BAD_CAST "idfversion", BAD_CAST "1.02");
    /* xmlNewChild(attributes, NULL, BAD_CAST "filename", BAD_CAST filename); */
    return attributes;
}

xmlNodePtr simulation2idf_elementsandmolecules(const sample_model *sm) {
    if(!sm) {
        return NULL;
    }
    xmlNodePtr elementsandmolecules = xmlNewNode(NULL, BAD_CAST "elementsandmolecules");
    xmlNodePtr elements = xmlNewChild(elementsandmolecules, NULL, BAD_CAST "elements", NULL);
    xmlAddChild(elements, idf_new_node_printf(BAD_CAST "nelements", "%zu", sm->n_materials));
    for(size_t i = 0; i < sm->n_materials; i++) {
        xmlNodePtr element = xmlNewChild(elementsandmolecules, NULL, BAD_CAST "element", NULL);
        xmlNewChild(element, NULL, BAD_CAST "name", BAD_CAST sm->materials[i]->name);
    }
    return elementsandmolecules;
}

xmlNodePtr simulation2idf_structure(const sample_model *sm) {
    if(!sm) {
        return NULL;
    }
    xmlNodePtr structure = xmlNewNode(NULL, BAD_CAST "structure");
    if(sm->type == SAMPLE_MODEL_LAYERED) {
        xmlNodePtr layeredstructure = xmlNewChild(structure, NULL, BAD_CAST "layeredstructure", NULL);
        xmlAddChild(layeredstructure, idf_new_node_printf(BAD_CAST "nlayers", "%zu", sm->n_ranges));
        xmlNodePtr layers = simulation2idf_layers(sm);
        if(layers) {
            xmlAddChild(layeredstructure, layers);
        }
    }
    /* TODO: support point-by-point profiles using IDF "pointbypointstructure" */
    return structure;
}

xmlNodePtr simulation2idf_layers(const sample_model *sm) {
    xmlNodePtr layers = xmlNewNode(NULL, BAD_CAST "layers");
    for(size_t i = 0; i < sm->n_ranges; i++) {
        const sample_range *r = &sm->ranges[i];
        xmlNodePtr layer = xmlNewChild(layers, NULL, BAD_CAST "layer", NULL);
        xmlNodePtr layerthickness = idf_new_node_units(BAD_CAST "layerthickness", BAD_CAST IDF_UNIT_TFU, NULL, r->x);
        xmlAddChild(layer, layerthickness);
        if(r->rough.model == ROUGHNESS_GAMMA && r->rough.x > 0.0) {
            xmlNodePtr layeruniformity = idf_new_node_units(BAD_CAST "layeruniformity", BAD_CAST IDF_UNIT_TFU, BAD_CAST "FWHM", r->rough.x);
            xmlAddChild(layer, layeruniformity);
        }
        xmlNodePtr layerelements = xmlNewChild(layer, NULL, BAD_CAST "layerelements", NULL);
        for(size_t j = 0; j < sm->n_materials; j++) {
            double conc = *sample_model_conc_bin(sm, i, j);
            if(conc == 0.0) {
                continue;
            }
            xmlNodePtr layerelement = xmlNewNode(NULL, BAD_CAST "layerelement");
            xmlNewChild(layerelement, NULL, BAD_CAST "name", BAD_CAST sm->materials[j]->name);
            xmlNodePtr concentration = idf_new_node_units(BAD_CAST "concentration", BAD_CAST IDF_UNIT_FRACTION, NULL, conc);
            xmlAddChild(layerelement, concentration);
            xmlAddChild(layerelements, layerelement);
        }
    }
    return layers;
}

xmlNodePtr simulation2idf_spectra(const struct fit_data *fit) {
    if(!fit || !fit->sim || !fit->sim->beam_isotope) {
        return NULL;
    }
    xmlNodePtr spectra = xmlNewNode(NULL, BAD_CAST "spectra");
    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        const detector *det = sim_det(fit->sim, i_det);
        xmlNodePtr spectrum = xmlNewNode(NULL, BAD_CAST "spectrum");
        xmlAddChild(spectrum, simulation2idf_beam(fit->sim));
        xmlAddChild(spectrum, simulation2idf_geometry(fit->sim, det));
        xmlAddChild(spectrum, simulation2idf_detection(det));
        xmlAddChild(spectrum, simulation2idf_calibrations(fit->jibal->elements, fit->sim, det));
        xmlAddChild(spectrum, simulation2idf_reactions(fit->sim, det));
        if(i_det < fit->n_det_spectra) {
            xmlNodePtr data = simulation2idf_data(result_spectra_experimental_histo(&(fit->spectra[i_det])));
            if(data) {
                xmlAddChild(spectrum, data);
            }
            xmlNodePtr process = simulation2idf_process(fit->sim, &fit->spectra[i_det]);
            if(process) {
                xmlAddChild(spectrum, process);
            }
        }
        xmlAddChild(spectra, spectrum);
    }
    return spectra;
}

xmlNodePtr simulation2idf_beam(const simulation *sim) {
    xmlNodePtr beam = xmlNewNode(NULL, BAD_CAST "beam");
    xmlAddChild(beam, idf_new_node_printf(BAD_CAST "beamparticle", "%s", sim->beam_isotope->name));
    xmlAddChild(beam, idf_new_node_printf(BAD_CAST "beamZ", "%i", sim->beam_isotope->Z));
    xmlAddChild(beam, idf_new_node_units(BAD_CAST "beammass", BAD_CAST IDF_UNIT_AMU, NULL, sim->beam_isotope->mass));
    xmlAddChild(beam, idf_new_node_units(BAD_CAST "beamenergy", BAD_CAST IDF_UNIT_KEV, NULL, sim->beam_E));
    xmlAddChild(beam, idf_new_node_units(BAD_CAST "beamenergyspread", BAD_CAST IDF_UNIT_KEV, BAD_CAST IDF_MODE_FWHM, sim->beam_E_broad));
    xmlAddChild(beam, idf_new_node_units(BAD_CAST "beamfluence", BAD_CAST IDF_UNIT_PARTICLES, NULL, sim->fluence));
    if(sim->beam_aperture) {
        xmlAddChild(beam, simulation2idf_aperture("beamshape", sim->beam_aperture));
    }
    return beam;
}

xmlNodePtr simulation2idf_geometry(const simulation *sim, const detector *det) {
    xmlNodePtr geometry = xmlNewNode(NULL, BAD_CAST "geometry");
    double incidenceangle = sim_alpha_angle(sim);
    double scatteringangle = det->theta;
    double exitangle = sim_exit_angle(sim, det);

    double ibm1 = exit_angle(incidenceangle, 0.0, scatteringangle, 0.0); /* Exit angle in IBM geometry */
    double ibm2 = exit_angle(-1.0 * incidenceangle, 0.0, scatteringangle, 0.0); /* Exit angle in IBM geometry (another solution) */
    double cornell = exit_angle(incidenceangle, 0.0, scatteringangle, 90.0 * C_DEG); /* Exit angle in Cornell geometry */
    char *geometrytype;
    if(fabs(exitangle - ibm1) < 0.1 * C_DEG || fabs(exitangle - ibm2) < 0.1 * C_DEG) {
        geometrytype = "IBM";
    } else if(fabs(exitangle - cornell) < 0.1 * C_DEG) {
        geometrytype = "Cornell";
    } else {
        geometrytype = "general";
    }
    xmlNewChild(geometry, NULL, BAD_CAST "geometrytype", BAD_CAST geometrytype);
    xmlAddChild(geometry, idf_new_node_units(BAD_CAST "incidenceangle", BAD_CAST IDF_UNIT_DEGREE, NULL, incidenceangle));
    xmlAddChild(geometry, idf_new_node_units(BAD_CAST "scatteringangle", BAD_CAST IDF_UNIT_DEGREE, NULL, scatteringangle));
    xmlAddChild(geometry, idf_new_node_units(BAD_CAST "exitangle", BAD_CAST IDF_UNIT_DEGREE, NULL, exitangle));
    return geometry;
}

xmlNodePtr simulation2idf_detection(const detector *det) {
    xmlNodePtr detection = xmlNewNode(NULL, BAD_CAST "detection");
    xmlNodePtr detector = xmlNewChild(detection, NULL, BAD_CAST "detector", NULL);
    char *detectortype = NULL;
    switch(det->type) {
        case DETECTOR_ENERGY:
            detectortype = "SSB";
            break;
        case DETECTOR_TOF:
            detectortype = "ToF";
            break;
        case DETECTOR_ELECTROSTATIC:
            detectortype = "ElSpec";
            break;
        default:
            break;
    }
    if(detectortype) {
        xmlAddChild(detector, idf_new_node_printf(BAD_CAST "detectortype", "%s", detectortype));
    }
    xmlAddChild(detector, idf_new_node_units(BAD_CAST "solidangle", BAD_CAST IDF_UNIT_MSR, NULL, det->solid));
    if(det->aperture) {
        xmlAddChild(detector, simulation2idf_aperture("detectorshape", det->aperture));
    }
    if(det->type == DETECTOR_TOF) {
        xmlNodePtr tof = xmlNewChild(detector, NULL, BAD_CAST "tof", NULL);
        xmlAddChild(tof, idf_new_node_units(BAD_CAST "toflength", BAD_CAST IDF_UNIT_MM, NULL, det->length));
        xmlAddChild(tof, idf_new_node_units(BAD_CAST "toftimeresolution", BAD_CAST IDF_UNIT_NS, BAD_CAST IDF_MODE_FWHM, sqrt(det->calibration->resolution_variance)));
    }
    if(det->distance > 0.0) {
        xmlAddChild(detector, idf_new_node_units(BAD_CAST "distancedetectortosample", BAD_CAST IDF_UNIT_MM, NULL, det->distance));
    }
    if(det->foil_sm) {
        xmlNodePtr stoppingfoil = xmlNewChild(detection, NULL, BAD_CAST "stoppingfoil", NULL);
        xmlNodePtr foillayers = xmlNewChild(stoppingfoil, NULL, BAD_CAST "foillayers", NULL);
        xmlAddChild(foillayers, simulation2idf_layers(det->foil_sm));
        xmlAddChild(detection, stoppingfoil);
    }
    return detection;
}

xmlNodePtr simulation2idf_calibrations(const jibal_element *elements, const simulation *sim, const detector *det) {
    xmlNodePtr calibrations = xmlNewNode(NULL, BAD_CAST "calibrations");
    xmlNodePtr detectorresolutions = xmlNewChild(calibrations, NULL, BAD_CAST "detectorresolutions", NULL);
    xmlNodePtr energycalibrations = xmlNewChild(calibrations, NULL, BAD_CAST "energycalibrations", NULL);
    xmlAddChild(energycalibrations, simulation2idf_energycalibration(det, det->calibration, NULL));
    for(int Z = 0; Z <= det->cal_Z_max; Z++) {
        const calibration *cal = detector_get_calibration(det, Z);
        if(cal == det->calibration) {
            continue;
        }
        const char *ion_name = jibal_element_name(elements, Z);
        xmlAddChild(detectorresolutions, simulation2idf_detectorresolution(det, cal, ion_name));
        xmlAddChild(energycalibrations, simulation2idf_energycalibration(det, cal, ion_name));
    }
    return calibrations;
}

xmlNodePtr simulation2idf_detectorresolution(const detector *det, const calibration *cal, const char  *ion_name) {
    xmlNodePtr detectorresolution = xmlNewNode(NULL, BAD_CAST "detectorresolution");
    if(ion_name) {
        xmlNewChild(detectorresolution, NULL, BAD_CAST "resolutionion", BAD_CAST ion_name); /* This can be either ion or isotope name, (e.g. "He" or "4He"). If NULL, we omit this field (default resolution). This is not exactly correct way to do it with NDF if multiple species are detected. */
    }
    xmlNodePtr resolutionparameters = xmlNewChild(detectorresolution, NULL, BAD_CAST "resolutionparameters", NULL);
    double resolution;
    if(det->type ==  DETECTOR_ENERGY) {
        resolution = cal->resolution;
    } else {
        resolution = DETECTOR_RESOLUTION_DEFAULT; /* TODO: SIMNRA saves energy resolution even for ToF detector. We can't. Save default. */
    }
    resolution /=  C_FWHM;
    xmlAddChild(resolutionparameters, idf_new_node_units(BAD_CAST "resolutionparameter", BAD_CAST IDF_UNIT_KEV, BAD_CAST IDF_MODE_FWHM, resolution));
    return detectorresolution;
}

xmlNodePtr simulation2idf_energycalibration(const detector *det, const calibration *cal, const char *ion_name) {
    xmlNodePtr energycalibration = xmlNewNode(NULL, BAD_CAST "energycalibration");
    if(ion_name) {
        xmlNewChild(energycalibration, NULL, BAD_CAST "calibrationion", BAD_CAST ion_name); /* This can be either ion or isotope name, (e.g. "He" or "4He"). If NULL, we omit this field (default calibration). This is not exactly correct way to do it with NDF if multiple species are detected. */
    }
    xmlNewChild(energycalibration, NULL, BAD_CAST "calibrationmode", BAD_CAST "energy");
    xmlNodePtr calibrationparameters = xmlNewChild(energycalibration, NULL, BAD_CAST "calibrationparameters", NULL);
    if(cal->type == CALIBRATION_LINEAR || cal->type == CALIBRATION_POLY) {
        xmlAddChild(calibrationparameters, idf_new_node_units(BAD_CAST "calibrationparameter", BAD_CAST IDF_UNIT_KEV, NULL, calibration_get_param(cal, CALIBRATION_PARAM_OFFSET)));
        xmlAddChild(calibrationparameters, idf_new_node_units(BAD_CAST "calibrationparameter", BAD_CAST IDF_UNIT_KEVCH, NULL, calibration_get_param(cal, CALIBRATION_PARAM_SLOPE)));
        if(calibration_get_number_of_params(cal) > 2) {
            xmlAddChild(calibrationparameters, idf_new_node_units(BAD_CAST "calibrationparameter", BAD_CAST IDF_UNIT_KEVCH2, NULL, calibration_get_param(cal, CALIBRATION_PARAM_QUAD)));
        }
        /* TODO: higher degree polynomials could be supported, but then some solution for the units "channel^n" needs to be figured out */
    }
    return energycalibration;
}

xmlNodePtr simulation2idf_aperture(const char *name, const aperture *aperture) {
    if(!aperture) {
        return NULL;
    }
    xmlNodePtr n = xmlNewNode(NULL, BAD_CAST name);
    if(aperture->type == APERTURE_CIRCLE) {
        xmlNewChild(n, NULL, BAD_CAST "shape", BAD_CAST "circular");
        xmlAddChild(n, idf_new_node_units(BAD_CAST "l1", BAD_CAST IDF_UNIT_MM, NULL, aperture->width));
    } else if(aperture->type == APERTURE_ELLIPSE) {
        xmlNewChild(n, NULL, BAD_CAST "shape", BAD_CAST "ellipse");
        xmlAddChild(n, idf_new_node_units(BAD_CAST "l1", BAD_CAST IDF_UNIT_MM, NULL, aperture->width));
        xmlAddChild(n, idf_new_node_units(BAD_CAST "l2", BAD_CAST IDF_UNIT_MM, NULL, aperture->height));
    } else if(aperture->type == APERTURE_RECTANGLE) {
        xmlNewChild(n, NULL, BAD_CAST "shape", BAD_CAST "rectangular");
        xmlAddChild(n, idf_new_node_units(BAD_CAST "l1", BAD_CAST IDF_UNIT_MM, NULL, aperture->width));
        xmlAddChild(n, idf_new_node_units(BAD_CAST "l2", BAD_CAST IDF_UNIT_MM, NULL, aperture->height));
    } else if(aperture->type == APERTURE_SQUARE) {
        xmlNewChild(n, NULL, BAD_CAST "shape", BAD_CAST "square");
        xmlAddChild(n, idf_new_node_units(BAD_CAST "l1", BAD_CAST IDF_UNIT_MM, NULL, aperture->width));
    }
    return n;
}

xmlNodePtr simulation2idf_reactions(const simulation *sim, const detector *det) {
    xmlNodePtr reactions = xmlNewNode(NULL, BAD_CAST "reactions");
    /* TODO: technique? We could see if any RBS reactions are used and mark (primary) technique as RBS? */

    reaction_type default_reaction_type = REACTION_NONE;
    for(size_t i = 0; i < sim->n_reactions; i++) {
        const reaction *r = sim->reactions[i];
        if(r->type == REACTION_ERD) {
            default_reaction_type = REACTION_ERD;
            break; /* Even one ERDA reaction => this is definitively ERDA (primary technique) */
        }
        if(r->type == REACTION_RBS) {
            default_reaction_type = REACTION_RBS; /* One RBS reaction => this could be RBS (note no break after this assignment) */
        }
    }
    char *technique;
    switch(default_reaction_type) {
        case REACTION_RBS:
            technique = "RBS";
            break;
        case REACTION_ERD:
            technique = "ERDA";
            break;
        default:
            technique = NULL;
    }
    if(technique) {
        xmlNewChild(reactions, NULL, BAD_CAST "technique", BAD_CAST technique);
    }

    xmlNodePtr reactionlist = xmlNewNode(NULL, BAD_CAST "reactionlist");
    xmlAddChild(reactions, reactionlist);
    for(size_t i = 0; i < sim->n_reactions; i++) {
        const reaction *r = sim->reactions[i];
        if(!reaction_is_possible(r, sim->params, det->theta)) { /* This is what is used to figure out which reaction is possible with this particular detector. We should probably use the result of an earlier calculation that actually was used to determine the viability. */
            continue;
        }
        if(default_reaction_type == r->type) { /* Don't list reactions that are part of the "technique", usually RBS or ERDA */
            continue;
        }
        xmlNodePtr reaction = xmlNewNode(NULL, BAD_CAST "reaction");
        xmlNewChild(reaction, NULL, BAD_CAST "initialtargetparticle", BAD_CAST r->target->name);
        xmlNewChild(reaction, NULL, BAD_CAST "incidentparticle", BAD_CAST r->incident->name);
        xmlNewChild(reaction, NULL, BAD_CAST "exitparticle", BAD_CAST r->product->name);
        xmlNewChild(reaction, NULL, BAD_CAST "finaltargetparticle", BAD_CAST r->residual->name);
        xmlAddChild(reaction, idf_new_node_units(BAD_CAST "reactionQ", BAD_CAST IDF_UNIT_MEV, NULL, r->Q));
        xmlAddChild(reactionlist, reaction);
    }
    return reactions;
}

xmlNodePtr simulation2idf_data(const jabs_histogram *h) {
    if(!h) {
        return NULL;
    }
    xmlNodePtr data = xmlNewNode(NULL, BAD_CAST "data");
    xmlNewChild(data, NULL, BAD_CAST "datamode", BAD_CAST "simple");
    xmlNewChild(data, NULL, BAD_CAST "channelmode", BAD_CAST "left");
    xmlNodePtr simpledata = simulation2idf_simpledata(h);
    if(simpledata) { /* Note that if there is no experimental spectrum, exp is NULL, simpledata becomes NULL and no <simpledata> is added at all! (We could make an all-zero spectrum if we wanted to) */
        xmlAddChild(data, simpledata);
    }
    return data;
}

xmlNodePtr simulation2idf_process(const simulation *sim, const result_spectra *spectrum) {
    xmlNodePtr process = xmlNewNode(NULL, BAD_CAST "process");
    xmlNodePtr simulations = xmlNewChild(process, NULL, BAD_CAST "simulations", NULL);
    xmlNodePtr simulation = xmlNewChild(simulations, NULL, BAD_CAST "simulation", NULL);
    xmlNewChild(simulation, NULL, BAD_CAST "simulationtype", BAD_CAST "total");
    xmlNewChild(simulation, NULL, BAD_CAST "datamode", BAD_CAST "simple");
    xmlNewChild(simulation, NULL, BAD_CAST "channelmode", BAD_CAST "left");
    xmlNodePtr simpledata = simulation2idf_simpledata(result_spectra_simulated_histo(spectrum));
    if(simpledata) {
        xmlAddChild(simulation, simpledata);
    }
    return process;
}


xmlNodePtr simulation2idf_simpledata(const jabs_histogram *h) {
    if(!h) {
        return NULL;
    }
    xmlNodePtr simpledata = xmlNewNode(NULL, BAD_CAST "simpledata");
    xmlNodePtr xaxis = xmlNewChild(simpledata, NULL, BAD_CAST "xaxis", NULL);
    xmlNewChild(xaxis, NULL, BAD_CAST "axisname", BAD_CAST "channel");
    xmlNewChild(xaxis, NULL, BAD_CAST "axisunit", BAD_CAST "#");
    xmlNodePtr yaxis = xmlNewChild(simpledata, NULL, BAD_CAST "yaxis", NULL);
    xmlNewChild(yaxis, NULL, BAD_CAST "axisname", BAD_CAST "yield");
    xmlNewChild(yaxis, NULL, BAD_CAST "axisunit", BAD_CAST "counts");

    int prec = 10; /* Precision we ask */
    size_t l_max = prec + 10; /* In addition to "prec" digits we needs space for leading sign, decimal separator, 'e', sign of exponent, and for the exponent (max 4 digits), and a space. This is nine characters. I probably forgot something, so ten it is.  We use the same math for channel numbers, because we are lazy. */
    size_t x_alloc = l_max * h->n;
    size_t y_alloc = l_max * h->n;
    char *x_data_orig = calloc(x_alloc, sizeof(char));
    char *y_data_orig = calloc(y_alloc, sizeof(char));
    if(!x_data_orig || !y_data_orig) {
        free(x_data_orig);
        free(y_data_orig);
        return NULL;
    }
    DEBUGMSG("IDF simpledata output, space for %zu and %zu characters for x and y, respectively based on %zu bins.", x_alloc, y_alloc, h->n);
    char *x_data = x_data_orig;
    char *y_data = y_data_orig;
    char *x_data_last = x_data_orig + x_alloc - 1;
    char *y_data_last = y_data_orig + y_alloc - 1;
    for(size_t i = 0; i <= h->n; i++) {
        int l = snprintf(x_data, x_data_last - x_data, "%zu ", i);
        if(l < 0) {
            /* TODO: error handling */
            break;
        }
        x_data += l;
        l = snprintf(y_data, y_data_last - y_data , "%.*g ", prec, i < h->n ? h->bin[i] : 0.0);
        if(l < 0) {
            /* TODO: error handling */
            break;
        }
        y_data += l;
    }
    *x_data = '\0';
    *y_data = '\0';
    xmlNewChild(simpledata, NULL, BAD_CAST "x", BAD_CAST x_data_orig);
    xmlNewChild(simpledata, NULL, BAD_CAST "y", BAD_CAST y_data_orig);
    free(x_data_orig);
    free(y_data_orig);
    return simpledata;
}
