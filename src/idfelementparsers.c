/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

    Some parts of this source file under different license, see below!

 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <string.h>
#include "idfelementparsers.h"
#include "generic.h"
#ifdef WIN32
#include "win_compat.h"
#endif
#include "geostragg.h"

idf_error idf_parse_note(idf_parser *idf, xmlNode *note) {
    char *s = idf_node_content_to_str(note);
    if(*s == '\0') { /* Empty */
        return IDF2JBS_SUCCESS;
    }
    idf_output_printf(idf, "#Note in file: %s\n", jabs_strip_newline(s)); /* TODO: Support multiline notes. This truncates the string. */
    free(s);
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_sample(idf_parser *idf, xmlNode *sample) {
    idf->i_sample++; /* Numbering from 1 .. idf->n_samples inclusive */
    free(idf->sample_basename);
    if(idf->n_samples == 1) {
        idf->sample_basename = strdup(idf->basename);
    } else {
        int len = asprintf(&idf->sample_basename, "%s-sample%zu", idf->basename, idf->i_sample);
        if(len < 0) {
            return IDF2JBS_FAILURE;
        }
    }
    idf_output_printf(idf, "#Sample number %zu\n", idf->i_sample);
    char *description = idf_node_content_to_str(idf_findnode(sample, "description"));
    if(strlen(description) > 0) {
        idf_output_printf(idf, "#Description: %s\n", jabs_strip_newline(description));
    }
    free(description);
    idf_foreach(idf, idf_findnode(sample, "notes"), "note", idf_parse_note);
    idf_parse_layers(idf, idf_findnode(sample, "structure/layeredstructure/layers"));
    idf_parse_spectra(idf, idf_findnode(sample, "spectra"));
    idf_output_printf(idf, "#End of sample %zu\n", idf->i_sample);
    if(idf->i_sample != idf->n_samples) {
        idf_output_printf(idf, "reset\n\n", idf->i_sample);
    }
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_layerelement(idf_parser *idf, xmlNode *element) {
    if(!element) {
        return IDF2JBS_FAILURE;
    }
    char *name = idf_node_content_to_str(idf_findnode(element, "name"));
    if(!name) {
        return IDF2JBS_FAILURE;
    }
    double conc = idf_node_content_to_double(idf_findnode(element, "concentration"));
    if(conc == 1.0) {
        idf_output_printf(idf, "%s", name);
    } else {
        idf_output_printf(idf, "%s%g", name, conc);
    }
    free(name);
    return IDF2JBS_SUCCESS;
}

int idf_parse_layerelements(idf_parser *idf, xmlNode *elements) {
    idf_output_printf(idf, " ");
    return idf_foreach(idf, elements, "layerelement", idf_parse_layerelement);
}

idf_error idf_parse_layer(idf_parser *idf, xmlNode *layer) {
    if(idf_parse_layerelements(idf, idf_findnode(layer, "layerelements")) == IDF2JBS_FAILURE) {
        return IDF2JBS_FAILURE;
    }
    idf_output_printf(idf, " %gtfu", idf_node_content_to_double(idf_findnode(layer, "layerthickness")) / C_TFU);
    double roughness = idf_node_content_to_double(idf_findnode(layer, "layeruniformity"));
    if(roughness > 0.0) {
        idf_output_printf(idf, " rough %gtfu", roughness / C_TFU); /* Note that output is *not* FWHM */
    }
    return IDF2JBS_SUCCESS;
}


idf_error idf_parse_energycalibrations(idf_parser *idf, xmlNode *energycalibrations) {
    xmlNode *energycalibration = idf_findnode(energycalibrations, "energycalibration");
    if(!energycalibration) {
        fprintf(stderr, "No energy calibration in energycalibrations.\n");
        return IDF2JBS_FAILURE;
    }
#ifdef IDF_DEBUG
    fprintf(stderr, "There is energy calibration.\n");
#endif
    xmlNode *calibrationparameters = idf_findnode(energycalibration, "calibrationparameters");
#ifdef IDF_DEBUG
    fprintf(stderr, "There are calibration parameters.\n");
#endif
    xmlNode *cur;
    int n_params = 0;
    double params[CALIB_PARAMS_MAX]; /* First CALIB_PARAMS_MAX parameters are stored here */
    for (cur = calibrationparameters->children; cur; cur = cur->next) {
        if(cur->type == XML_ELEMENT_NODE) {
            if(idf_stringeq(cur->name, "calibrationparameter")) {
                double param = idf_node_content_to_double(cur);
                if(n_params < 3) {
                    params[n_params] = param;
                }
                n_params++;
                if(n_params == CALIB_PARAMS_MAX) {
                    break;
                }
            }
        }
    }
#ifdef DEBUG
    fprintf(stderr, "There are %i params (or more than %i).\n", n_params, CALIB_PARAMS_MAX);
#endif
    if(n_params < 2) {
        return IDF2JBS_FAILURE;
    }
    if(n_params == 3 && params[2] == 0.0) { /* Quad term is zero, this is linear */
        n_params--;
    }
    if(n_params == 2) {
        idf_output_printf(idf, "set det calib linear slope %gkeV offset %gkeV\n", params[1]/C_KEV, params[0]/C_KEV);
    } else {
        idf_output_printf(idf, "set det calib poly %i", n_params - 1);
        for(int i = 0; i < n_params; i++) {
            idf_output_printf(idf, " %gkeV", params[i] / C_KEV);
        }
        idf_output_printf(idf, "\n");
    }
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_simple_data(xmlNode *simple_data, const char *filename) {
    if(!simple_data) {
        return IDF2JBS_FAILURE;
    }
    char *x_raw = idf_node_content_to_str(idf_findnode(simple_data, "x"));
    char *y_raw = idf_node_content_to_str(idf_findnode(simple_data, "y"));
    return idf_write_simple_data_to_file(filename, x_raw, y_raw);
}


idf_error idf_parse_spectrum(idf_parser *idf, xmlNode *spectrum) {
    idf->i_spectrum++;
#ifdef DEBUG
    fprintf(stderr, "parse spectrum called, i_spectrum = %zu.\n", idf->i_spectrum);
#endif
    idf_parse_beam(idf, idf_findnode(spectrum, "beam"));
    idf_parse_process(idf, idf_findnode(spectrum, "process"));
    idf_parse_geometry(idf, idf_findnode(spectrum, "geometry"));
    idf_parse_calibrations(idf, idf_findnode(spectrum, "calibrations"));
    idf_parse_detector(idf, idf_findnode(spectrum, "detection/detector"));
    idf_parse_detector_foil(idf, idf_findnode(spectrum, "detection/stoppingfoil"));
    /* TODO: parse reactions */
    idf_parse_data(idf, idf_findnode(spectrum, "data"));

    idf_output_printf(idf, "simulate\n");

    char *savespectrafilename = idf_spectrum_out_name(idf, idf->i_spectrum);
    if(savespectrafilename) {
        idf_output_printf(idf, "save spectra \"%s\"\n", savespectrafilename);
        free(savespectrafilename);
    }
    if(idf->i_spectrum != idf->n_spectra) {
        idf_output_printf(idf, "reset detectors\nadd detector default\n", idf->i_sample);
    }
    return IDF2JBS_SUCCESS;
}


idf_error idf_parse_calibrations(idf_parser *idf, xmlNode *calibrations) {
    if(!calibrations) {
        return IDF2JBS_FAILURE;
    }
#ifdef DEBUG
    fprintf(stderr, "There are calibrations.\n");
#endif
    idf_parse_energycalibrations(idf, idf_findnode(calibrations, "energycalibrations"));
    double resolution = idf_node_content_to_double(idf_findnode(calibrations, "detectorresolutions/detectorresolution/resolutionparameters/resolutionparameter"));
    if(resolution > 0.0) {
        idf_output_printf(idf, "set det resolution %gkeV\n", resolution / C_KEV * C_FWHM);
    }
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_beam(idf_parser *idf, xmlNode *beam) {
    if(!beam) {
        return IDF2JBS_FAILURE;
    }
#ifdef DEBUG
    fprintf(stderr, "There is beam.\n");
#endif
    char *particle = idf_node_content_to_str(idf_findnode(beam, "beamparticle"));
    idf_output_printf(idf, "set ion %s\n", particle);
    free(particle);
    double energy = idf_node_content_to_double(idf_findnode(beam, "beamenergy"));
    idf_output_printf(idf, "set energy %gkeV\n", energy/C_KEV);
    double energyspread = idf_node_content_to_double(idf_findnode(beam, "beamenergyspread"));
    if(energyspread > 0.0) {
        idf_output_printf(idf, "set energy_broad %gkeV\n", energyspread / C_KEV * C_FWHM);
    }
    double fluence = idf_node_content_to_double(idf_findnode(beam, "beamfluence"));
    idf_output_printf(idf, "set fluence %e\n", fluence);
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_data(idf_parser *idf, xmlNode *data) {
    if(!data) {
        return IDF2JBS_FAILURE;
    }
    char *datamode = idf_node_content_to_str(idf_findnode(data, "datamode"));
    if(idf_stringeq(datamode, "simple")) {
        char *simpledata_filename = idf_exp_name(idf, idf->i_spectrum);
        if(idf_parse_simple_data(idf_findnode(data, "simpledata"), simpledata_filename) == IDF2JBS_SUCCESS) {
            idf_output_printf(idf, "load exp \"%s\"\n", simpledata_filename);
        }
        free(simpledata_filename);
    }
    free(datamode);
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_geometry(idf_parser *idf, xmlNode *geometry) {
    if(!geometry) {
        return IDF2JBS_FAILURE;
    }
#ifdef DEBUG
    fprintf(stderr, "There is geometry.\n");
#endif
    char *geotype = idf_node_content_to_str(idf_findnode(geometry, "geometrytype"));
    const char *phi_str = "";
    double incidenceangle = idf_node_content_to_double(idf_findnode(geometry, "incidenceangle"));
    double scatteringangle = idf_node_content_to_double(idf_findnode(geometry, "scatteringangle"));
    double exitangle = idf_node_content_to_double(idf_findnode(geometry, "exitangle"));
    double ibm1 = exit_angle(incidenceangle, 0.0, scatteringangle, 0.0);
    double ibm2 = exit_angle(-1.0 * incidenceangle, 0.0, scatteringangle, 0.0);
    double cornell = exit_angle(incidenceangle, 0.0, scatteringangle, 90.0 * C_DEG);
#ifdef DEBUG
    fprintf(stderr, "Possible angles: IBM1: %g deg, IBM2: %g deg,  Cornell: %g deg\n", ibm1 / C_DEG, ibm2 / C_DEG, cornell / C_DEG);
#endif
    if(idf_stringeq(geotype, "IBM")) {
        idf_output_printf(idf, "#IBM geometry set in file.\n");
        if(fabs(exitangle - ibm1) < 0.1 * C_DEG) {
            idf_output_printf(idf, "#Geometry in file matches with IBM exit angle %g deg.\n", ibm1 / C_DEG);
        } else if(fabs(exitangle - ibm2) < 0.1 * C_DEG) {
            idf_output_printf(idf, "#Geometry in file matches with IBM exit angle %g deg.\n", ibm2 / C_DEG);
            incidenceangle *= -1.0;
        } else {
            idf_output_printf(idf, "#Geometry in file (exit angle %g deg) does not match with IBM geometry. Check it!\n", exitangle / C_DEG);
        }
    } else if(idf_stringeq(geotype, "Cornell")) {
        idf_output_printf(idf, "#Cornell geometry set in file.\n");
        if(fabs(exitangle - cornell) < 0.1 * C_DEG) {
            idf_output_printf(idf, "#Geometry in file matches with Cornell exit angle %g deg.\n", cornell / C_DEG);
            phi_str = " phi 90.0deg";
        } else {
            idf_output_printf(idf, "#Geometry in file (exit angle %g deg) does not match with Cornell geometry! Check it.\n", exitangle / C_DEG);
        }
    } else {
        idf_output_printf(idf, "#Not IBM or Cornell geometry (exit angle %g deg)? Use set det phi XXXdeg to tell JaBS where the detector is.\n", exitangle / C_DEG);
        phi_str = " phi 0.0deg";
    }
    free(geotype);
    idf_output_printf(idf, "set alpha %gdeg\n", incidenceangle / C_DEG);
    idf_output_printf(idf, "set det theta %gdeg%s\n", scatteringangle / C_DEG, phi_str);
    idf_parse_spot(idf, idf_findnode(geometry, "spot"));
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_spectra(idf_parser *idf, xmlNode *spectra) {
    idf->i_spectrum = 0;
    idf->n_spectra = idf_foreach(idf, spectra, "spectrum", NULL);
    return idf_foreach(idf, spectra, "spectrum", idf_parse_spectrum);
}

idf_error idf_parse_layers(idf_parser *idf, xmlNode *layers) {
    idf_output_printf(idf, "set sample");
    idf_foreach(idf, layers, "layer", idf_parse_layer);
    idf_output_printf(idf, "\n");
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_detector(idf_parser *idf, xmlNode *detector) {
    if(!detector) {
        return IDF2JBS_FAILURE;
    }
    char *type = idf_node_content_to_str(idf_findnode(detector, "detectortype"));
    if(idf_stringeq(type, "SSB")) {
        idf_output_printf(idf, "set detector type energy\n");
    } else if(idf_stringeq(type, "ToF")) {
        double length = idf_node_content_to_double(idf_findnode(detector, "tof/toflength"));
        double timing_resolution = idf_node_content_to_double(idf_findnode(detector, "tof/toftimeresolution"));;
        idf_output_printf(idf, "set det type tof length %gmm resolution %gps\n", length / C_MM, timing_resolution / C_PS * C_FWHM);
    }
    free(type);
    double solid = idf_node_content_to_double(idf_findnode(detector, "solidangle"));
    if(solid > 0.0) {
        idf_output_printf(idf, "set det solid %gmsr\n", solid / C_MSR);
    }
    double distance = idf_node_content_to_double(idf_findnode(detector, "distancedetectortosample"));
    if(distance > 0.0) {
        idf_output_printf(idf, "set det distance %gmm\n", distance / C_MM);
    }
    idf_parse_detectorshape(idf, idf_findnode(detector, "detectorshape"));



    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_detectorshape(idf_parser *idf, xmlNode *detectorshape) {
    return idf_detectorshape_or_spot(idf, detectorshape, "set det aperture");
}

idf_error idf_parse_detector_foil(idf_parser *idf, xmlNode *stoppingfoil) {
    if(!stoppingfoil) {
        return IDF2JBS_FAILURE;
    }
    idf_output_printf(idf, "set det foil ");
    xmlNode *foillayers = idf_findnode(stoppingfoil, "foillayers");
    idf_foreach(idf, foillayers, "layer", idf_parse_layer);
    idf_output_printf(idf, "\n");
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_spot(idf_parser *idf, xmlNode *spot) {
    return idf_detectorshape_or_spot(idf, spot, "set aperture");
}

idf_error idf_detectorshape_or_spot(idf_parser *idf, xmlNode *node, const char *prefix) {
    if(!node) {
        return IDF2JBS_FAILURE;
    }
    double l1 = idf_node_content_to_double(idf_findnode(node, "l1"));
    double l2 = idf_node_content_to_double(idf_findnode(node, "l2"));
    char *shape = idf_node_content_to_str(idf_findnode(node, "shape"));
    if(idf_stringneq(shape, "circ", 4)) {
        if(l1 > 0.0) {
            idf_output_printf(idf, "%s circle diameter %gmm\n", prefix, l1 / C_MM);
        }
    } else if(idf_stringneq(shape, "rect", 4)) {
        if(l1 > 0.0 && l2 > 0.0) {
            idf_output_printf(idf, "%s rectangle width %gmm height %gmm\n", prefix, l1 / C_MM, l2 / C_MM);
        }
    } else if(idf_stringneq(shape, "ellip", 4)) {
        if(l1 > 0.0 && l2 > 0.0) {
            idf_output_printf(idf, "%s ellipse width %gmm height %gmm\n", prefix, l1 / C_MM, l2 / C_MM);
        }
    } else if(idf_stringneq(shape, "squ", 4)) {
        if(l1 > 0.0) {
            idf_output_printf(idf, "%s square width %gmm\n", prefix, l1 / C_MM);
        }
    }
    free(shape);
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_process(idf_parser *idf, xmlNode *process) {
    if(!process) {
        return IDF2JBS_FAILURE;
    }
    idf_parse_physicsdefaults(idf, idf_findnode(process, "physicsdefaults"));
    idf_parse_simulations(idf, idf_findnode(process, "simulations"));
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_physicsdefaults(idf_parser *idf, xmlNode *physicsdefaults) {
    if(!physicsdefaults) {
        return IDF2JBS_FAILURE;
    }
    xmlNode *crosssectiondefault = idf_findnode(physicsdefaults, "crosssectiondefault");
    if(crosssectiondefault) {
        char *screening = idf_node_content_to_str(idf_findnode(crosssectiondefault, "screening"));
        const char *cs;
        if(idf_stringeq(screening, "Andersen") || idf_stringeq(screening, "Universal")) {
            cs = screening;
        } else if(idf_stringeq(screening, "Ecuyer")) {
            cs = "LEcuyer";
        } else if(idf_stringeq(screening, "none")) {
            cs = "Rutherford";
        } else {
            cs = NULL;
        }
        if(cs) {
            idf_output_printf(idf, "set cs_rbs %s\n", cs);
            idf_output_printf(idf, "set cs_erd %s\n", cs);
        }
    }


    xmlNode *energyspreaddefault = idf_findnode(physicsdefaults, "energyspreaddefault");
    if(energyspreaddefault) {
        int beamsizespread = idf_node_content_to_boolean(idf_findnode(energyspreaddefault, "geometricspread/beamsize"));
        int detectoraperturespread = idf_node_content_to_boolean(idf_findnode(energyspreaddefault, "geometricspread/beamsize"));
        if(beamsizespread == TRUE && detectoraperturespread == TRUE) {
            idf_output_printf(idf, "set geostragg true\n");
        } else if(beamsizespread == TRUE) {
            idf_output_printf(idf, "#In IDF file only beam size spread is set to true. The following line may also enable detector size spread.\n");
            idf_output_printf(idf, "set geostragg true\n");
        } else if(detectoraperturespread == TRUE) {
            idf_output_printf(idf, "#In IDF file only detector size spread is set to true. The following line may also enable beam size spread.\n");
            idf_output_printf(idf, "set geostragg true\n");
        }
    }
    return IDF2JBS_SUCCESS;
}


idf_error idf_parse_simulations(idf_parser *idf, xmlNode *simulations) {
    idf_foreach(idf, simulations, "simulation", idf_parse_simulation);
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_simulation(idf_parser *idf, xmlNode *simulation) {
    char *type = idf_node_content_to_str(idf_findnode(simulation, "simulationtype"));
    if(!idf_stringeq(type, "total")) { /* Ignore other spectra, just check the total one. */
        free(type);
        return IDF2JBS_SUCCESS;
    }
    char *datamode = idf_node_content_to_str(idf_findnode(simulation, "datamode"));
    if(idf_stringeq(datamode, "simple")) {
        char *simpledata_filename = idf_sim_name(idf, idf->i_spectrum);
        if(idf_parse_simple_data(idf_findnode(simulation, "simpledata"), simpledata_filename) == IDF2JBS_SUCCESS) {
            idf_output_printf(idf, "load ref \"%s\"\n", simpledata_filename);
        }
        free(simpledata_filename);
    }
    free(datamode);
    return IDF2JBS_SUCCESS;
}

