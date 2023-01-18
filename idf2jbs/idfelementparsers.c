#include "idfelementparsers.h"

int idf_parse_sample(idfparser *idf, xmlNode *sample) {
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

int idf_parse_layerelement(idfparser *idf, xmlNode *element) {
    if(!element) {
        return IDF2JBS_FAILURE;
    }
    xmlNode *name = findnode(element, "name");
    xmlNode *concentration = findnode(element, "concentration");
    if(name && concentration) {
        xmlChar *name_str = xmlNodeGetContent(name);
        double conc = idf_node_content_to_double(concentration);
        if(conc == 1.0) {
            fprintf(stdout, "%s", name_str);
        } else {
            fprintf(stdout, "%s%g", name_str, conc);
        }
        free(name_str);
    } else {
        fprintf(stderr, "Layer ignored.\n");
    }
    return IDF2JBS_SUCCESS;
}

int idf_parse_layerelements(idfparser *idf, xmlNode *elements) {
    fprintf(stdout, " ");
    return idf_foreach(idf, elements, "layerelement", idf_parse_layerelement);
}

int idf_parse_layer(idfparser *idf, xmlNode *layer) {
    idf_parse_layerelements(idf, findnode(layer, "layerelements"));
    fprintf(stdout, " thick %gtfu", idf_node_content_to_double(findnode(layer, "layerthickness"))/C_TFU);
    return IDF2JBS_SUCCESS;
}


int idf_parse_energycalibrations(idfparser *idf, xmlNode *energycalibrations) {
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
        fprintf(stdout, "set det calib linear slope %gkeV offset %gkeV\n", params[1]/C_KEV, params[0]/C_KEV);
    } else {
        fprintf(stdout, "set det calib poly %i", n_params - 1);
        for(int i = 0; i < n_params; i++) {
            fprintf(stdout, " %gkeV", params[i] / C_KEV);
        }
    }
    return IDF2JBS_SUCCESS;
}

int idf_parse_simple_data(idfparser *idf, xmlNode *simple_data) {
    if(!simple_data) {
        return IDF2JBS_FAILURE;
    }
    char *x_raw = idf_node_content_to_str(findnode(simple_data, "x"));
    char *y_raw = idf_node_content_to_str(findnode(simple_data, "y"));
    idf_write_simple_data_to_file("out.dat", x_raw, y_raw); /* TODO: figure out the spectrum output filename */
    return IDF2JBS_SUCCESS;
}


int idf_parse_spectrum(idfparser *idf, xmlNode *spectrum) {
#ifdef DEBUG
    fprintf(stderr, "parse spectrum called.\n");
#endif
    xmlNode *beam = findnode(spectrum, "beam");
    if(beam) {
#ifdef DEBUG
        fprintf(stderr, "There is beam.\n");
#endif
        char *particle = idf_node_content_to_str(findnode(beam, "beamparticle"));
        fprintf(stdout, "set ion %s\n", particle);
        free(particle);
        double energy = idf_node_content_to_double(findnode(beam, "beamenergy"));
        fprintf(stdout, "set energy %gkeV\n", energy/C_KEV);
        double energyspread = idf_node_content_to_double(findnode(beam, "beamenergyspread"));
        if(energyspread > 0.0) {
            fprintf(stdout, "set energy_broad %gkeV\n", energyspread / C_KEV);
        }
        double fluence = idf_node_content_to_double(findnode(beam, "beamfluence"));
        fprintf(stdout, "set fluence %e\n", fluence);
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
        fprintf(stdout, "set alpha %gdeg\n", incidenceangle/C_DEG);
        double scatteringangle = idf_node_content_to_double(findnode(geometry, "scatteringangle"));
        fprintf(stdout, "set det theta %gdeg\n", scatteringangle/C_DEG);
    }
    xmlNode *calibrations = findnode(spectrum, "calibrations");
    if(calibrations) {
#ifdef DEBUG
        fprintf(stderr, "There are calibrations.\n");
#endif
        idf_parse_energycalibrations(idf, findnode(calibrations, "energycalibrations"));
        double resolution = idf_node_content_to_double(findnode(calibrations, "detectorresolutions/detectorresolution/resolutionparameters/resolutionparameter"));
        if(resolution > 0.0) {
            fprintf(stdout, "set det resolution %gkeV\n", resolution / C_KEV);
        }
    }
    /* TODO: parse reactions */


    xmlNode *data_node = findnode(spectrum, "data");
    if(data_node) {
        xmlNode *simple_data_node = findnode(data_node, "simpledata");
        idf_parse_simple_data(idf, simple_data_node);
    }
    return IDF2JBS_SUCCESS;
}

int idf_parse_spectra(idfparser *idf, xmlNode *spectra) {
    return idf_foreach(idf, spectra, "spectrum", idf_parse_spectrum);
}

int idf_parse_layers(idfparser *idf, xmlNode *layers) {
    fprintf(stdout, "set sample");
    idf_foreach(idf, layers, "layer", idf_parse_layer);
    fprintf(stdout, "\n");
    return IDF2JBS_SUCCESS;
}
