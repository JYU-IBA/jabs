/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

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
#include "simulation2idf.h"
#include "options.h"


int simulation2idf(struct fit_data *fit, const char *filename) {
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

    xmlNodePtr sample = xmlNewChild(n, NULL, BAD_CAST "sample", NULL);

    sample_model *sm2 = sample_model_split_elements(fit->sm); /* TODO: this splits molecules. We might not necessarily want it always, but it is a reasonable default. */
    xmlNodePtr elementsandmolecules = simulation2idf_elementsandmolecules(sm2);
    if(elementsandmolecules) {
        xmlAddChild(sample, elementsandmolecules);
    }
    sample_model_free(sm2);

    xmlChar *xmlbuff;
    int buffersize;
    xmlDocDumpFormatMemory(doc, &xmlbuff, &buffersize, 1);
    fprintf(f, "%s", (char *) xmlbuff);
    fclose(f);
    xmlFreeDoc(doc);
    jabs_message(MSG_WARNING, "Saving simulations to IDF files support is incomplete, output file was created but it might not work as intended.");
    return EXIT_SUCCESS;
}

xmlNodePtr simulation2idf_notes() {
    xmlNodePtr notes = xmlNewNode(NULL, BAD_CAST "notes");
    xmlNodePtr note = xmlNewChild(notes, NULL, BAD_CAST "note", NULL);

    char *jabs_note = NULL;
    if(asprintf(&jabs_note, "File created by JaBS %s", jabs_version()) > 0) {
        xmlNodeSetContent(note, BAD_CAST jabs_note);
        free(jabs_note);
    }
    return notes;
}

xmlNodePtr simulation2idf_elementsandmolecules(const sample_model *sm) {
    if(!sm) {
        return NULL;
    }
    xmlNodePtr elementsandmolecules = xmlNewNode(NULL, BAD_CAST "elementsandmolecules");
    xmlNodePtr elements = xmlNewChild(elementsandmolecules, NULL, BAD_CAST "elements", NULL);
    char *nelements;
    asprintf(&nelements, "%zu", sm->n_materials);
    xmlNewChild(elements, NULL, BAD_CAST "nelements", BAD_CAST nelements);
    free(nelements);
    for(size_t i = 0; i < sm->n_materials; i++) {
        xmlNodePtr element = xmlNewChild(elementsandmolecules, NULL, BAD_CAST "element", NULL);
        xmlNewChild(element, NULL, BAD_CAST "name", BAD_CAST sm->materials[i]->name);
    }
    return elementsandmolecules;
}
