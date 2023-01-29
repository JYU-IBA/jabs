#include <string.h>
#include "idfelementparsers.h"
#include "generic.h"
#ifdef WIN32
#include "win_compat.h"
#endif

idf_error idf_parse_note(idf_parser *idf, xmlNode *note) {
    char *s = idf_node_content_to_str(note);
    if(*s == '\0') { /* Empty */
        return IDF2JBS_SUCCESS;
    }
    idf_output_printf(idf, "#Note in IDF file: %s\n", jabs_strip_newline(s)); /* TODO: Support multiline notes. This truncates the string. */
    free(s);
    return IDF2JBS_SUCCESS;
}
int idf_parse_sample(idf_parser *idf, xmlNode *sample) {
    idf->i_sample++; /* Numbering from 1 .. idf->n_samples inclusive */
    free(idf->sample_basename);
    if(idf->n_samples == 1) {
        idf->sample_basename = strdup(idf->basename);
    } else {
        asprintf(&idf->sample_basename, "%s-sample%zu", idf->basename, idf->i_sample);
    }
    idf_output_printf(idf, "#Sample number %zu in file %s\n", idf->i_sample, idf->filename);
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

int idf_parse_layerelement(idf_parser *idf, xmlNode *element) {
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

int idf_parse_layer(idf_parser *idf, xmlNode *layer) {
    if(idf_parse_layerelements(idf, idf_findnode(layer, "layerelements")) == IDF2JBS_FAILURE) {
        return IDF2JBS_FAILURE;
    }
    idf_output_printf(idf, " %gtfu", idf_node_content_to_double(idf_findnode(layer, "layerthickness")) / C_TFU);
    return IDF2JBS_SUCCESS;
}


int idf_parse_energycalibrations(idf_parser *idf, xmlNode *energycalibrations) {
    xmlNode *energycalibration = idf_findnode(energycalibrations, "energycalibration");
    if(!energycalibration) {
        fprintf(stderr, "No energy calibration in energycalibrations.\n");
        return IDF2JBS_FAILURE;
    }
#ifdef DEBUG
    fprintf(stderr, "There is energy calibration.\n");
#endif
    xmlNode *calibrationparameters = idf_findnode(energycalibration, "calibrationparameters");
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
    char *x_raw = idf_node_content_to_str(idf_findnode(simple_data, "x"));
    char *y_raw = idf_node_content_to_str(idf_findnode(simple_data, "y"));
    char *filename = idf_exp_name(idf, idf->i_spectrum);
    if(idf_write_simple_data_to_file(filename, x_raw, y_raw) == IDF2JBS_SUCCESS) {
        idf_output_printf(idf, "load exp \"%s\"\n", filename);
    }
    free(filename);
    return IDF2JBS_SUCCESS;
}


int idf_parse_spectrum(idf_parser *idf, xmlNode *spectrum) {
    idf->i_spectrum++;
#ifdef DEBUG
    fprintf(stderr, "parse spectrum called, i_spectrum = %zu.\n", idf->i_spectrum);
#endif
    xmlNode *beam = idf_findnode(spectrum, "beam");
    if(beam) {
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
    }
    idf_parse_process(idf, idf_findnode(spectrum, "process"));

    xmlNode *geometry = idf_findnode(spectrum, "geometry");
    if(geometry) {
#ifdef DEBUG
        fprintf(stderr, "There is geometry.\n");
#endif
        char *geotype = idf_node_content_to_str(idf_findnode(geometry, "geometrytype"));
        const char *phi_str = "";
        if(idf_stringeq(geotype, "IBM")) {
            idf_output_printf(idf, "#IBM geometry set in file.\n");
        } else if(idf_stringeq(geotype, "Cornell")) {
            phi_str = " phi 90.0deg";
            idf_output_printf(idf, "#Cornell geometry set in file.\n");
        } else {
            idf_output_printf(idf, "#Not IBM or Cornell geometry? Use set det phi XXXdeg to tell JaBS where the detector is.\n");
        }
        free(geotype);
        double incidenceangle = idf_node_content_to_double(idf_findnode(geometry, "incidenceangle"));
        idf_output_printf(idf, "set alpha %gdeg\n", incidenceangle/C_DEG);
        double scatteringangle = idf_node_content_to_double(idf_findnode(geometry, "scatteringangle"));
        idf_output_printf(idf, "set det theta %gdeg%s\n", scatteringangle/C_DEG, phi_str);
        double exitangle = idf_node_content_to_double(idf_findnode(geometry, "exitangle"));
        idf_output_printf(idf, "#Check output to confirm JaBS uses correct exit angle (%g deg?)\n", exitangle / C_DEG);
        idf_parse_spot(idf, idf_findnode(geometry, "spot"));
    }
    xmlNode *calibrations = idf_findnode(spectrum, "calibrations");
    if(calibrations) {
#ifdef DEBUG
        fprintf(stderr, "There are calibrations.\n");
#endif
        idf_parse_energycalibrations(idf, idf_findnode(calibrations, "energycalibrations"));
        double resolution = idf_node_content_to_double(idf_findnode(calibrations, "detectorresolutions/detectorresolution/resolutionparameters/resolutionparameter"));
        if(resolution > 0.0) {
            idf_output_printf(idf, "set det resolution %gkeV\n", resolution / C_KEV * C_FWHM);
        }
    }

    idf_parse_detector(idf, idf_findnode(spectrum, "detection/detector"));

    /* TODO: parse reactions */


    xmlNode *data_node = idf_findnode(spectrum, "data");
    if(data_node) {
        xmlNode *simple_data_node = idf_findnode(data_node, "simpledata");
        idf_parse_simple_data(idf, simple_data_node);
    }
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

int idf_parse_spectra(idf_parser *idf, xmlNode *spectra) {
    idf->i_spectrum = 0;
    idf->n_spectra = idf_foreach(idf, spectra, "spectrum", NULL);
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
    }
    free(shape);
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_process(idf_parser *idf, xmlNode *process) {
    if(!process) {
        return IDF2JBS_FAILURE;
    }
    idf_parse_physicsdefaults(idf, idf_findnode(process, "physicsdefaults"));
    return IDF2JBS_SUCCESS;
}

idf_error idf_parse_physicsdefaults(idf_parser *idf, xmlNode *physicsdefaults) {
    if(!physicsdefaults) {
        return IDF2JBS_FAILURE;
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
