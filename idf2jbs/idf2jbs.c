#include <stdio.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "idf2jbs.h"

int stringeq(const void *a, const void *b) {
    if(!a || !b)
        return -1;
    return (strcmp(a,b) == 0);
}

int nodename_equals(const xmlNode *node, const char *s) {
    if(!node || !s)
        return -1;
    return (strcmp((char *)node->name, s) == 0);
}

double unit_string_to_SI(xmlChar *unit) {
    for(idfunit *u = idfunits; u->unit; u++) {
        if(xmlStrEqual(unit, xmlstr(u->unit))) {
            return u->factor;
        }
    }
    return 0.0;
}

double node_content_to_double(const xmlNode *node) {
    xmlChar *content = xmlNodeGetContent(node);
    double out = strtod((char *)content, NULL); /* TODO: error checking */
    xmlChar *unit = xmlGetProp(node, xmlstr("units"));
    if(unit) {
        out *= unit_string_to_SI(unit);
        free(unit);
    }
    free(content);
    return out;
}

char *node_content_to_str(const xmlNode *node) {
    if(!node) {
        return strdup("");
    }
    xmlChar *content = xmlNodeGetContent(node);
    if(!content) {
        return strdup("");
    }
    return (char *) content;
}

const xmlChar *xmlstr(const char *s) {
    return (const xmlChar *)s;
}

xmlNode *findnode(xmlNode *root, const char *path) {
#ifdef DEBUG
    fprintf(stderr, "findnode called with path = \"%s\".\n", path);
#endif
    size_t offset = strcspn(path, "/"); /* We use strcspn to tokenize because we don't have to mutilate path */
    if(offset == 0) {
        fprintf(stderr, "Weirdness. Offset is zero.\n");
        return NULL;
    }
    int last = (path[offset] == '\0');
    const char *next_path = path + offset + 1;
    xmlNode *cur_node = NULL;
#ifdef DEBUG
    fprintf(stderr, "Should find %.*s in %s\n", (int)offset, path, root->name);
#endif
    for (cur_node = root->children; cur_node; cur_node = cur_node->next) {
        if(cur_node->type == XML_ELEMENT_NODE) {
            size_t len = xmlStrlen(cur_node->name);
            if(len != offset) /* Prevents partial matches in strncmp */
                continue;
            if(strncmp((char *)cur_node->name, path, offset) == 0) {
#ifdef DEBUG
                fprintf(stderr, "Found something.\n");
#endif
                if(last) {
#ifdef DEBUG
                    fprintf(stderr, "Found the last one.\n");
#endif
                    return cur_node;
                } else {
                    return findnode(cur_node, next_path);
                }
            }
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Returning NULL\n");
#endif
    return NULL;
}

void parse_idf(xmlNode *idfnode) {
#ifdef DEBUG
    fprintf(stderr, "Parsing root node.\n");
#endif
    xmlNode *cur_node = NULL;
    for (cur_node = idfnode->children; cur_node; cur_node = cur_node->next) {
        if(cur_node->type == XML_ELEMENT_NODE) {

        }
    }
}

void parse_layerelements(xmlNode *elements) {
    if(!elements)
        return;
    xmlNode *cur = NULL;
    fprintf(stdout, " ");
    for (cur = elements->children; cur; cur = cur->next) {
        if(cur->type == XML_ELEMENT_NODE) {
            if(nodename_equals(cur, "layerelement")) {
                double conc = 0.0;
                xmlNode *name = findnode(cur, "name");
                xmlNode *concentration = findnode(cur, "concentration");
                if(name && concentration) {
                    xmlChar *name_str = xmlNodeGetContent(name);
                    conc = node_content_to_double(concentration);
                    if(conc == 1.0) {
                        fprintf(stdout, "%s", name_str);
                    } else {
                        fprintf(stdout, "%s%g", name_str, conc);
                    }
                    free(name_str);
                } else {
                    fprintf(stderr, "Layer ignored.\n");
                }
            }
        }
    }
}

int write_simple_data_to_file(const char *filename, const char *x, const char *y) {
#if 0
    FILE *f = fopen(filename, "w");
    if(!f) {
        return IDF2JBS_FAILURE;
    }
#endif
    while(1) {
        char *x_end, *y_end;
        double x_dbl = strtod(x, &x_end);
        double y_dbl = strtod(y, &y_end);
        if(y == y_end || x == x_end)  { /* No conversion */
            break;
        }
        fprintf(stderr, "%g %g\n", floor(x_dbl), y_dbl);
        y = y_end;
        x = x_end;
    }
#if 0
    fclose(f);
#endif
}

void parse_simple_data(xmlNode *simple_data) {
    if(!simple_data) {
        return;
    }
    char *x_raw = node_content_to_str(findnode(simple_data, "x"));
    char *y_raw = node_content_to_str(findnode(simple_data, "y"));
    write_simple_data_to_file("out.dat", x_raw, y_raw); /* TODO: figure out the spectrum output filename */
}

void parse_layer(xmlNode *layer) {
    xmlNode *cur = NULL;
    double thickness = 0.0;
    parse_layerelements(findnode(layer, "layerelements"));

    for (cur = layer->children; cur; cur = cur->next) {
        if(cur->type == XML_ELEMENT_NODE) {
            if(nodename_equals(cur, "layerthickness")) {
                thickness = node_content_to_double(cur);
            }
        }
    }
    fprintf(stdout, " thick %gtfu", thickness/C_TFU);
}

void parse_energycalibrations(xmlNode *energycalibrations) {
    xmlNode *energycalibration = findnode(energycalibrations, "energycalibration");
    if(!energycalibration) {
        fprintf(stderr, "No energy calibration in energycalibrations.\n");
        return;
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
            if(stringeq(cur->name, "calibrationparameter")) {
                double param = node_content_to_double(cur);
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
        return;
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
}

void parse_spectrum(xmlNode *spectrum) {
#ifdef DEBUG
    fprintf(stderr, "parse spectrum called.\n");
#endif
    xmlNode *beam = findnode(spectrum, "beam");
    if(beam) {
#ifdef DEBUG
        fprintf(stderr, "There is beam.\n");
#endif
        char *particle = node_content_to_str(findnode(beam, "beamparticle"));
        fprintf(stdout, "set ion %s\n", particle);
        free(particle);
        double energy = node_content_to_double(findnode(beam, "beamenergy"));
        fprintf(stdout, "set energy %gkeV\n", energy/C_KEV);
        double energyspread = node_content_to_double(findnode(beam, "beamenergyspread"));
        if(energyspread > 0.0) {
            fprintf(stdout, "set energy_broad %gkeV\n", energyspread / C_KEV);
        }
        double fluence = node_content_to_double(findnode(beam, "beamfluence"));
        fprintf(stdout, "set fluence %e\n", fluence);
    }
    xmlNode *geometry = findnode(spectrum, "geometry");
    if(geometry) {
#ifdef DEBUG
        fprintf(stderr, "There is geometry.\n");
#endif
        char *geotype = node_content_to_str(findnode(geometry, "geometrytype"));
        if(stringeq(geotype, "IBM")) {
            /* do nothing */
        } else if(stringeq(geotype, "Cornell")) {
            /* TODO: do something */
        }
        free(geotype);
        double incidenceangle = node_content_to_double(findnode(geometry, "incidenceangle"));
        fprintf(stdout, "set alpha %gdeg\n", incidenceangle/C_DEG);
        double scatteringangle = node_content_to_double(findnode(geometry, "scatteringangle"));
        fprintf(stdout, "set det theta %gdeg\n", scatteringangle/C_DEG);
    }
    xmlNode *calibrations = findnode(spectrum, "calibrations");
    if(calibrations) {
#ifdef DEBUG
        fprintf(stderr, "There are calibrations.\n");
#endif
        parse_energycalibrations(findnode(calibrations, "energycalibrations"));
        double resolution = node_content_to_double(findnode(calibrations, "detectorresolutions/detectorresolution/resolutionparameters/resolutionparameter"));
        if(resolution > 0.0) {
            fprintf(stdout, "set det resolution %gkeV\n", resolution / C_KEV);
        }
    }
    /* TODO: parse reactions */


    xmlNode *data_node = findnode(spectrum, "data");
    if(data_node) {
        xmlNode *simple_data_node = findnode(data_node, "simpledata");
        parse_simple_data(simple_data_node);
    }

}

void parse_spectra(xmlNode *spectra) {
    xmlNode *cur_node = NULL;
    for(cur_node = spectra->children; cur_node; cur_node = cur_node->next) {
        if(cur_node->type == XML_ELEMENT_NODE) {
            if(stringeq(cur_node->name, "spectrum")) {
                parse_spectrum(cur_node);
            }
        }
    }
}

void parse_layers(xmlNode *layers) {
    xmlNode *cur_node = NULL;
    fprintf(stdout, "set sample");
    for (cur_node = layers->children; cur_node; cur_node = cur_node->next) {
        if(cur_node->type == XML_ELEMENT_NODE) {
            if(stringeq(cur_node->name, "layer")) {
                parse_layer(cur_node);
            }
        }
    }
    fprintf(stdout, "\n");
}

void parse_sample(xmlNode *sample) {
    xmlNode *spectra_node = findnode(sample, "spectra");
    if(spectra_node) {
        parse_spectra(spectra_node);
    }
    xmlNode *layers_node = findnode(sample, "structure/layeredstructure/layers");
    if(layers_node) {
        parse_layers(layers_node);
    }
}

int parse_xml(const char *filename) {
    xmlDoc *doc = NULL;
    xmlNode *root_element = NULL;
    xmlNode *cur = NULL;
    doc = xmlReadFile(filename, NULL, 0);
    if(!doc) {
        return IDF2JBS_FAILURE;
    }
    root_element = xmlDocGetRootElement(doc);
    if(stringeq(root_element->name, "idf")) {
        parse_idf(root_element);
    }
    for (cur = root_element->children; cur; cur = cur->next) {
        if(cur->type == XML_ELEMENT_NODE) {
            if(stringeq(cur->name, "sample")) {
                parse_sample(cur); /* TODO: multiple samples? */
            }
        }
    }
    xmlFreeDoc(doc);
    xmlCleanupParser();
    return IDF2JBS_SUCCESS;
}
