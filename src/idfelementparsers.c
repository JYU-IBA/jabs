#include "idfelementparsers.h"

int idf_parse_sample(idf_parser *idf, xmlNode *sample) {
    xmlNode *spectra_node = findnode(sample, "spectra");
    if(spectra_node) {
        idf_parse_spectra(idf, spectra_node);
    }
    xmlNode *layers_node = findnode(sample, "structure/layeredstructure/layers");
    if(layers_node) {
        idf_parse_layers(idf, layers_node);
    }
    return IDF2JBS_SUCCESS;
}

int idf_parse_layerelement(idf_parser *idf, xmlNode *element) {
    if(!element) {
        return IDF2JBS_FAILURE;
    }
    xmlNode *name = findnode(element, "name");
    xmlNode *concentration = findnode(element, "concentration");
    if(name && concentration) {
        xmlChar *name_str = xmlNodeGetContent(name);
        double conc = idf_node_content_to_double(concentration);
        if(conc == 1.0) {
            idf_output_printf(idf, "%s", name_str);
        } else {
            idf_output_printf(idf, "%s%g", name_str, conc);
        }
        free(name_str);
    } else {
        fprintf(stderr, "Layer ignored.\n");
    }
    return IDF2JBS_SUCCESS;
}

int idf_parse_layerelements(idf_parser *idf, xmlNode *elements) {
    idf_output_printf(idf, " ");
    return idf_foreach(idf, elements, "layerelement", idf_parse_layerelement);
}

int idf_parse_layer(idf_parser *idf, xmlNode *layer) {
    idf_parse_layerelements(idf, findnode(layer, "layerelements"));
    idf_output_printf(idf, " %gtfu", idf_node_content_to_double(findnode(layer, "layerthickness"))/C_TFU);
    return IDF2JBS_SUCCESS;
}


int idf_parse_energycalibrations(idf_parser *idf, xmlNode *energycalibrations) {
    xmlNode *energycalibration = findnode(energycalibrations, "energycalibration");
    if(!energycalibration) {
        fprintf(stderr, "No energy calibration in energycalibrations.\n");
        return IDF2JBS_FAILURE;
    }
#ifdef DEBUG
    fprintf(stderr, "There is energy calibration.\n");
#endif
    xmlNode *calibrationparameters = findnode(energycalibration, "calibrationparameters");
#ifdef DEBUG
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
    }
    return IDF2JBS_SUCCESS;
}

int idf_parse_simple_data(idf_parser *idf, xmlNode *simple_data) {
    if(!simple_data) {
        return IDF2JBS_FAILURE;
    }
    char *x_raw = idf_node_content_to_str(findnode(simple_data, "x"));
    char *y_raw = idf_node_content_to_str(findnode(simple_data, "y"));
    char *filename = idf_file_name_with_suffix(idf, SPECTRUM_FILE_SUFFIX);
    if(idf_write_simple_data_to_file(filename, x_raw, y_raw) == IDF2JBS_SUCCESS) {
        idf_output_printf(idf, "load exp \"%s\"\n", filename);
    }
    free(filename);
    return IDF2JBS_SUCCESS;
}


int idf_parse_spectrum(idf_parser *idf, xmlNode *spectrum) {
#ifdef DEBUG
    fprintf(stderr, "parse spectrum called.\n");
#endif
    xmlNode *beam = findnode(spectrum, "beam");
    if(beam) {
#ifdef DEBUG
        fprintf(stderr, "There is beam.\n");
#endif
        char *particle = idf_node_content_to_str(findnode(beam, "beamparticle"));
        idf_output_printf(idf, "set ion %s\n", particle);
        free(particle);
        double energy = idf_node_content_to_double(findnode(beam, "beamenergy"));
        idf_output_printf(idf, "set energy %gkeV\n", energy/C_KEV);
        double energyspread = idf_node_content_to_double(findnode(beam, "beamenergyspread"));
        if(energyspread > 0.0) {
            idf_output_printf(idf, "set energy_broad %gkeV\n", energyspread / C_KEV);
        }
        double fluence = idf_node_content_to_double(findnode(beam, "beamfluence"));
        idf_output_printf(idf, "set fluence %e\n", fluence);
    }
    xmlNode *geometry = findnode(spectrum, "geometry");
    if(geometry) {
#ifdef DEBUG
        fprintf(stderr, "There is geometry.\n");
#endif
        char *geotype = idf_node_content_to_str(findnode(geometry, "geometrytype"));
        if(idf_stringeq(geotype, "IBM")) {
            /* do nothing */
        } else if(idf_stringeq(geotype, "Cornell")) {
            /* TODO: do something */
        }
        free(geotype);
        double incidenceangle = idf_node_content_to_double(findnode(geometry, "incidenceangle"));
        idf_output_printf(idf, "set alpha %gdeg\n", incidenceangle/C_DEG);
        double scatteringangle = idf_node_content_to_double(findnode(geometry, "scatteringangle"));
        idf_output_printf(idf, "set det theta %gdeg\n", scatteringangle/C_DEG);
    }
    xmlNode *calibrations = findnode(spectrum, "calibrations");
    if(calibrations) {
#ifdef DEBUG
        fprintf(stderr, "There are calibrations.\n");
#endif
        idf_parse_energycalibrations(idf, findnode(calibrations, "energycalibrations"));
        double resolution = idf_node_content_to_double(findnode(calibrations, "detectorresolutions/detectorresolution/resolutionparameters/resolutionparameter"));
        if(resolution > 0.0) {
            idf_output_printf(idf, "set det resolution %gkeV\n", resolution / C_KEV);
        }
    }

    idf_parse_detector(idf, findnode(spectrum, "detection/detector"));

    /* TODO: parse reactions */


    xmlNode *data_node = findnode(spectrum, "data");
    if(data_node) {
        xmlNode *simple_data_node = findnode(data_node, "simpledata");
        idf_parse_simple_data(idf, simple_data_node);
    }
    return IDF2JBS_SUCCESS;
}

int idf_parse_spectra(idf_parser *idf, xmlNode *spectra) {
    return idf_foreach(idf, spectra, "spectrum", idf_parse_spectrum);
}

int idf_parse_layers(idf_parser *idf, xmlNode *layers) {
    idf_output_printf(idf, "set sample");
    idf_foreach(idf, layers, "layer", idf_parse_layer);
    idf_output_printf(idf, "\n");
    return IDF2JBS_SUCCESS;
}

int idf_parse_detector(idf_parser *idf, xmlNode *detector) {
    if(!detector) {
        return IDF2JBS_FAILURE;
    }
    char *type = idf_node_content_to_str(findnode(detector, "detectortype"));
    if(idf_stringeq(type, "SSB")) {
#if 0
        idf_output_printf(idf, "set detector type energy\n");
#endif
    }
    double solid = idf_node_content_to_double(findnode(detector, "solidangle"));
    if(solid > 0.0) {
        idf_output_printf(idf, "set det solid %gmsr\n", solid / C_MSR);
    }
    return IDF2JBS_SUCCESS;
}
