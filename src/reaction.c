/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "jabs_debug.h"
#include "generic.h"
#include "defaults.h"
#include "reaction.h"
#include "message.h"


void reactions_print(reaction * const * reactions, size_t n_reactions) {
    if(!reactions)
        return;
    jabs_message(MSG_INFO, stderr, "Reactions (%zu):\n", n_reactions);
    for(size_t i = 0; i < n_reactions; i++) {
        const reaction *r = reactions[i];
        if(!r) {
#ifdef DEBUG
            fprintf(stderr, "Reaction i=%zu (%zu) is NULL\n", i, i + 1);
#endif
            continue;
        }
        jabs_message(MSG_INFO, stderr, "%3zu: %s", i + 1, reaction_name(r));
#ifdef DEBUG
        jabs_message(MSG_INFO, stderr, " Reaction product stopping emin %g keV ", r->ion_gsto->emin / C_KEV); /* This is only valid after sim has been prepared */
#endif
        if(r->E_min > E_MIN || r->E_max < E_MAX) {
            jabs_message(MSG_INFO, stderr, ", E = [%.6g MeV, %.6g MeV]", r->E_min / C_MEV, r->E_max / C_MEV);
        }
        if(r->Q != 0.0) {
            jabs_message(MSG_INFO, stderr, ", Q = %g MeV ", r->Q / C_KEV);
        }
        if(r->type == REACTION_FILE) {
            jabs_message(MSG_INFO, stderr, ", Theta = %g deg,  Data from file \"%s\".\n", r->theta/C_DEG, r->filename);
        }
#ifdef JABS_PLUGINS
        else if(r->type == REACTION_PLUGIN) {
            jabs_message(MSG_INFO, stderr, ", Plugin \"%s\" filename \"%s\".\n", r->plugin->name, r->filename);
        }
#endif
        else if(r->type == REACTION_RBS || r->type == REACTION_RBS_ALT || r->type == REACTION_ERD){
            assert(r->Q == 0.0);
            jabs_message(MSG_INFO, stderr, ", %s cross sections (built-in).", jabs_reaction_cs_to_string(r->cs));
            jabs_message(MSG_INFO, stderr, "\n");
        }
    }
}

const char *reaction_name(const reaction *r) {
    if(!r)
        return "NULL";
    return r->name;
}

const char *reaction_type_to_string(reaction_type type) {
    switch (type) {
        case REACTION_NONE:
            return "NONE";
        case REACTION_RBS:
            return "RBS";
        case REACTION_RBS_ALT:
            return "RBS-";
        case REACTION_ERD:
            return "ERD";
        case REACTION_FILE:
            return "FILE";
        case REACTION_PLUGIN:
            return "PLUGIN";
        default:
            return "???";
    }
}

reaction_type reaction_type_from_string(const char *s) {
    if(!s)
        return REACTION_NONE;
    if(strcmp(s, "RBS") == 0)
        return REACTION_RBS;
    if(strcmp(s, "RBS-") == 0)
        return REACTION_RBS_ALT;
    if(strcmp(s, "ERD") == 0)
        return REACTION_ERD;
    if(strcmp(s, "FILE") == 0)
        return REACTION_FILE;
    if(strcmp(s, "PLUGIN") == 0)
        return REACTION_PLUGIN;
#ifdef DEBUG
    fprintf(stderr, "Reaction type string %s is not valid.\n", s);
#endif
    return REACTION_NONE;
}

reaction *reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, jabs_reaction_cs cs) {
    if(!target || !incident) {
        return NULL;
    }
    reaction *r = malloc(sizeof(reaction));
    r->type = type;
    r->cs = cs;
    r->incident = incident;
    r->target = target;
    r->filename = NULL;
    r->cs_table = NULL;
    r->n_cs_table = 0;
#ifdef JABS_PLUGINS
    r->plugin = NULL;
    r->plugin_r = NULL;
#endif
    r->nucl_stop = NULL; /* Will be handled when reaction is prepared */
    r->E_min = E_MIN;
    r->E_max = E_MAX;
    r->Q = 0.0;
    switch(type) {
        case REACTION_ERD:
            r->product = target;
            r->residual = incident;
            break;
        case REACTION_RBS: /* Falls through */
        case REACTION_RBS_ALT:
            r->product = incident;
            r->residual = target;
            break;
        case REACTION_FILE:
        case REACTION_PLUGIN:
        default:
            r->product = NULL;
            r->residual = NULL;
            break;
    }
    reaction_generate_name(r);
    return r;
}

reaction *reaction_make_from_argv(const jibal *jibal, const jibal_isotope *incident, int *argc, char * const **argv) {
    if((*argc) < 2) {
        jabs_message(MSG_ERROR, stderr, "Not enough arguments\n");
        return NULL;
    }
    reaction_type type = reaction_type_from_string((*argv)[0]);
    const jibal_isotope *target = jibal_isotope_find(jibal->isotopes, (*argv)[1], 0, 0);
    if(type == REACTION_NONE) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid reaction type: \"%s\".\n", (*argv)[0]);
        return NULL;
    }
    if(!target) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid isotope: \"%s\".\n", (*argv)[1]);
        return NULL;
    }
    reaction *r = reaction_make(incident, target, type, JABS_CS_NONE); /* Warning: JABS_CS_NONE used here, something sane must be supplied after this somewhere! */
    (*argc) -= 2;
    (*argv) += 2;
    while ((*argc) >= 2) {
        if(strcmp((*argv)[0], "max") == 0) {
            if(jabs_unit_convert(jibal->units, JIBAL_UNIT_TYPE_ENERGY, (*argv)[1], &r->E_max) < 0) {
                reaction_free(r);
                return NULL;
            }
        } else if(strcmp((*argv)[0], "min") == 0) {
            if(jabs_unit_convert(jibal->units, JIBAL_UNIT_TYPE_ENERGY, (*argv)[1], &r->E_min) < 0) {
                reaction_free(r);
                return NULL;
            }
        } else if(strcmp((*argv)[0], "cs") == 0) {
            r->cs = jibal_option_get_value(jibal_cs_types, (*argv)[1]);
            if(r->cs == 0) {
                jabs_message(MSG_ERROR, stderr, "Cross section type \"%s\" not recognized.\n", (*argv)[1]);
            }
        }
        (*argc) -= 2;
        (*argv) += 2;
    }
    return r;
}

void reaction_free(reaction *r) {
    if(!r)
        return;
#ifdef JABS_PLUGINS
    jabs_plugin_reaction_free(r->plugin, r->plugin_r);
    jabs_plugin_close(r->plugin);
#endif
    nuclear_stopping_free(r->nucl_stop);
    ion_gsto_free(r->ion_gsto);
    free(r->cs_table);
    free(r->filename);
    free(r->name);
    free(r);
}


int reaction_is_possible(const reaction *r, const sim_calc_params *p, double theta) {
    reaction_type type = r->type;
    const jibal_isotope *incident = r->incident;
    const jibal_isotope *target = r->target;
    if(type == REACTION_RBS || type == REACTION_RBS_ALT) {
        if(incident->mass >= target->mass && theta > asin(target->mass / incident->mass)) {
            DEBUGMSG("RBS with %s is not possible (theta %g deg > %g deg)", target->name, theta/C_DEG, asin(target->mass / incident->mass)/C_DEG);
            return FALSE;
        }
        if(type == REACTION_RBS_ALT) {
            if(incident->mass <= target->mass) {
                DEBUGMSG("RBS(-) with %s is not possible (target must be lighter than incident)", target->name);
                return FALSE;
            }
        }
    } else if(type == REACTION_ERD) {
        if(theta >= C_PI / 2.0) {
            return FALSE;
        }
    } else if(type == REACTION_FILE) {
            return reaction_theta_within_tolerance(r, p, theta);
    }
    return TRUE;
}

int reaction_theta_within_tolerance(const reaction *r, const sim_calc_params *p, double theta) {
    return (fabs(theta - r->theta) <= p->reaction_file_angle_tolerance);
}

reaction *r33_file_to_reaction(const jibal_isotope *isotopes, const r33_file *rfile) {
    const jibal_isotope *nuclei[R33_N_NUCLEI];
    if(rfile->n_data < 2) {
        jabs_message(MSG_ERROR, stderr, "Not enough data in file (n = %zu).\n", rfile->n_data);
        return NULL;
    }
    for(size_t i = 0; i < R33_N_NUCLEI; i++) {
#ifdef R33_IGNORE_REACTION_STRING
        nuclei[i] = jibal_isotope_find(isotopes, NULL, r33_double_to_int(rfile->zeds[i]), r33_double_to_int(rfile->masses[i]));
#else
        nuclei[i] = jibal_isotope_find(isotopes, rfile->reaction_nuclei[i], 0, 0);
#endif
        if(!nuclei[i]) {
#ifdef R33_IGNORE_REACTION_STRING
            fprintf(stderr, "Could not parse an isotope from Z=%g, mass=%g.", rfile->zeds[i], rfile->masses[i]);
#else
            jabs_message(MSG_ERROR, stderr, "Could not parse an isotope from \"%s\".\n", rfile->reaction_nuclei[i]);
#endif
            return NULL;
        }
    }
    if(rfile->composition) {
        jabs_message(MSG_ERROR, stderr, "This program does not currently support \"Composition\" in R33 files.\n");
        return NULL;
    }


    reaction *r = calloc(1, sizeof(reaction));
    r->cs = JABS_CS_NONE;
    r->type = REACTION_FILE;
#ifdef R33_IGNORE_REACTION_STRING
    r->incident = nuclei[0];
    r->target = nuclei[1];
#else
    r->incident = nuclei[1];
    r->target = nuclei[0];
#endif
    r->product = nuclei[2];
    r->residual = nuclei[3];
    if(rfile->unit == R33_UNIT_RR && r->incident != r->product) {
        jabs_message(MSG_ERROR, stderr, "R33 file is in units of ratio to Rutherford, but reaction product (%s) is not the same as target (%s). I don't know what to do.\n", r->product->name, r->target->name);
        return NULL;
    }
    r->theta = rfile->theta * C_DEG;
    r->Q = rfile->Qvalues[0] * C_KEV; /* TODO: other values? */
    if(rfile->filename) {
        r->filename = strdup(rfile->filename);
    } else {
        r->filename = NULL;
    }
    /* TODO: copy and convert data (check units etc) */
    r->n_cs_table = rfile->n_data;
    r->cs_table = calloc(r->n_cs_table, sizeof(struct reaction_point));
    for(size_t i = 0; i < rfile->n_data; i++) {
        struct reaction_point *rp = &r->cs_table[i];
        rp->E = rfile->data[i][0] * rfile->enfactors[0] * C_KEV; /* TODO: other factors? */
        if(rfile->unit == R33_UNIT_RR) {
            rp->sigma = rfile->data[i][2] * rfile->sigfactors[0] * jibal_cross_section_rbs(r->incident, r->target, r->theta, rp->E, JIBAL_CS_RUTHERFORD);
        } else if (rfile->unit == R33_UNIT_MB){
            rp->sigma = rfile->data[i][2] * rfile->sigfactors[0] * C_MB_SR;
        } else {
            rp->sigma = 0.0;
        }
    }
    r->E_min = r->cs_table[0].E;
    r->E_max = r->cs_table[r->n_cs_table - 1].E;
    r->nucl_stop = NULL; /* Will be handled when reaction is added */
    reaction_generate_name(r);
    return r;
}
int reaction_compare(const void *a, const void *b) {
    const reaction *r_a = *((const reaction **)a);
    const reaction *r_b = *((const reaction **)b);
    if(r_a == NULL)
        return 1; /* TODO: or 1? */
    if(r_b == NULL)
        return -1; /* TODO: or -1? */
    const jibal_isotope *isotope_a = r_a->target;
    const jibal_isotope *isotope_b = r_b->target;
    if(r_a->type < r_b->type)
        return -1;
    if(r_a->type > r_b->type)
        return 1;
    if(isotope_a->Z == isotope_b->Z) { /* Same Z, compare by A */
        return isotope_a->A - isotope_b->A;
    } else {
        return isotope_a->Z - isotope_b->Z;
    }
}

double reaction_product_energy(const reaction *r, double theta, double E) { /* Hint: call with E == 1.0 to get kinematic factor. Note that this function does not check if reaction is possible and may return "nan". */
    const double Q = r->Q;
    if(Q == 0.0) {
        if(r->type == REACTION_RBS || r->type == REACTION_RBS_ALT) {
            if(r->type == REACTION_RBS_ALT) {
                return jibal_kin_rbs(r->incident->mass, r->target->mass, theta, '-') * E;
            } else {
                return jibal_kin_rbs(r->incident->mass, r->target->mass, theta, '+') * E;
            }
        }  else if(r->type == REACTION_ERD) {
            return jibal_kin_erd(r->incident->mass, r->target->mass, theta)*E;
        }
    }
    const double m1 = r->incident->mass;
    const double m2 = r->target->mass;
    const double m3 = r->product->mass;
    const double m4 = r->residual->mass;
    double E_total = E + Q;
    double a13 = ((m1 * m3)/((m1 + m2)*(m3 + m4))) * (E / E_total);
    double a24 = ((m2 * m4)/((m1 + m2)*(m3 + m4))) * (1 + (m1/m2) * Q / E_total);
    if(a13 > a24) {
        double theta_max = sqrt(asin(a24/a13));
        if(theta > theta_max) {
            return 0.0;
        }
    }
    double E_out = E_total * a13 * pow2(cos(theta) + sqrt(a24/a13 - pow2(sin(theta)))); /* Solution with '-' is ignored. */
    return E_out;
}

const char *jabs_reaction_cs_to_string(jabs_reaction_cs cs) {
    switch(cs) {
        case JABS_CS_NONE:
            return "none";
        case JABS_CS_RUTHERFORD:
            return "Rutherford";
        case JABS_CS_ANDERSEN:
            return "Andersen";
        default:
            return "???";
    }
}

int reaction_generate_name(reaction *r) {
    int len = asprintf(&(r->name), "%-4s %s(%s,%s)%s", reaction_type_to_string(r->type), r->target->name, r->incident->name, r->product->name, r->residual->name);
    if(len < 0) {
        DEBUGSTR("Could not generate reaction name.");
        r->name = NULL;
    }
    return len;
}

jabs_reaction_cs jabs_reaction_cs_from_jibal_cs(jibal_cross_section_type jcs) {
    switch(jcs) {
        case JIBAL_CS_NONE:
            return JABS_CS_NONE;
        case JIBAL_CS_RUTHERFORD:
            return JABS_CS_RUTHERFORD;
        case JIBAL_CS_ANDERSEN:
            return JABS_CS_ANDERSEN;
        default:
#ifdef DEBUG
            fprintf(stderr, "Unknown JIBAL cs type: %i\n", jcs);
#endif
            return JABS_CS_NONE;
    }
}
