#include <string.h>
#include <assert.h>
#include <math.h>
#include <jibal_units.h>
#include <jibal_kin.h>


#include "defaults.h"
#include "rotate.h"
#include "roughness.h"
#include "jabs.h"

double stop_sample(sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, const depth depth, double E) {
    double em=E/incident->mass;
    double S1 = 0.0;
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        double c;
        if(sample->no_conc_gradients) {
            c = *sample_conc_bin(sample, depth.i, i_isotope);
        } else {
            c = get_conc(sample, depth, i_isotope);
        }
        if(c < ABUNDANCE_THRESHOLD)
            continue;
        if (type == GSTO_STO_TOT) {
            S1 += c * (
                    jibal_gsto_get_em(ws->gsto, GSTO_STO_ELE, incident->Z, sample->isotopes[i_isotope]->Z, em)
                    #ifdef NUCLEAR_STOPPING_FROM_JIBAL
                    +jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, sample->isotopes[i_isotope]->Z, sample->isotopes[i_isotope]->mass)
                    #else
                    + ion_nuclear_stop(incident, sample->isotopes[i_isotope], ws->isotopes, ws->nucl_stop_accurate)
                    #endif
                    );
        } else {
            S1 += c * (
                    jibal_gsto_get_em(ws->gsto, type, incident->Z, sample->isotopes[i_isotope]->Z, em)
            );
        }
    }
    //assert(S1 > 0.0);
    return S1;
}

double next_crossing(const ion *incident, const sample  *sample, depth *depth) { /* This will also update the depth to the right bin! */
    if(incident->inverse_cosine_theta > 0) { /* Going deeper */
        while(depth->i < sample->n_ranges - 1 && depth->x >= sample->ranges[depth->i + 1].x) {
            depth->i++;
        }
        assert(depth->i < sample->n_ranges - 1); /* There is a bug elsewhere in the code if you try to go this deep (deeper than last depth bin). */
        return sample->ranges[depth->i + 1].x - depth->x;
    } else if(incident->inverse_cosine_theta < 0.0) { /* Going towards the surface */
        while(depth->i > 0 && depth->x <= sample->ranges[depth->i].x) {
            depth->i--;
        }
        return sample->ranges[depth->i].x - depth->x;
    } else {
        fprintf(stderr, "WARNING: Inverse cosine is exactly zero. This is an issue!\n");
        return 0.0;
    }
}

depth stop_step(sim_workspace *ws, ion *incident, const sample *sample, depth depth, double step) {
    double k1, k2, k3, k4, stop, dE, E;
    double h_max_perp = next_crossing(incident, sample, &depth);
#ifdef DEBUG_STOP_STEP
    fprintf(stderr, "stop_step depth %g tfu (i=%zu) distance to next crossing %g tfu.\n", depth.x/C_TFU, depth.i, h_max_perp/C_TFU);
#endif
    /* k1...k4 are slopes of energy loss (stopping) at various x (depth) and E. Note convention: positive values, i.e. -dE/dx! */
    E = incident->E;
    k1 = stop_sample(ws, incident, sample, ws->stopping_type, depth, E);
    if(k1 < 0.001*C_EV_TFU) { /* Fail on positive values, zeroes (e.g. due to zero concentrations) and too small negative values */
#ifdef DEBUG_STOP_STEP
        fprintf(stderr, "stop_step returns no progress, because k1 = %g eV/tfu (x = %.3lf tfu, E = %.3lg keV)\n", k1/C_EV_TFU, depth.x/C_TFU, E/C_KEV);
#endif
        return depth;
    }
    double h_max = h_max_perp * incident->inverse_cosine_theta; /*  we can take bigger steps since we are going sideways. Note that inverse_cosine_theta can be negative and in this case h_max should also be negative so h_max is always positive! */
    assert(h_max >= 0.0);
    double h =  (step / k1); /* (energy) step should always be positive, as well as k1, so depth step h (not perpendicular, but "real" depth) is always positive  */
    assert(h > 0.0);
    double h_perp; /* has a sign (same as h_max_perp ) */
    if(h >= h_max) {
        h = h_max;
        h_perp = h_max_perp;
    } else {
        h_perp = h*incident->cosine_theta; /* x + h_perp is the actual perpendicular depth */
    }
    const struct depth halfdepth = {.x = depth.x + h_perp/2, .i = depth.i};
    const struct depth fulldepth = {.x = depth.x + h_perp, .i = depth.i};
    if(ws->rk4) {
        k2 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k1);
        k3 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k2);
        k4 = stop_sample(ws, incident, sample, ws->stopping_type, fulldepth, E - h * k3);
        stop = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    } else {
        stop = k1;
    }
    assert(stop > 0.0);
    dE =  -1.0* h * stop; /* Energy change in thickness "h". It is always negative! */
#if 0
    fprintf(stderr, "%s stop = %.3lf eV/tfu ( x = %.3lf, h = %.3lf, h_max = %.3lf), E = %.3lf keV, h = %6.3lf  dE = %.5lf keV\n", incident->isotope->name, stop/C_EV_TFU, x/C_TFU, h/C_TFU, h_max/C_TFU, E/C_KEV, h/C_TFU, dE/C_KEV);
#endif
#ifdef DEBUG_STOP_STEP
    if(fabs(stop) < 0.1*C_EV_TFU) {
        fprintf(stderr, "Not good!\n");
        return depth;
    }
#endif
#ifndef STATISTICAL_STRAGGLING
    double s_ratio = stop_sample(ws, incident, sample, ws->stopping_type, depth, E + dE) / k1; /* Ratio of stopping for non-statistical broadening. TODO: at x? */
#ifdef DEBUG
    //if((s_ratio)*(s_ratio) < 0.9 || (s_ratio)*(s_ratio) > 1.1) { /* Non-statistical broadening. */
    //   fprintf(stderr, "YIKES, s_ratio = %g, sq= %g\n", s_ratio, (s_ratio)*(s_ratio));
    //}
#endif
    incident->S *= (s_ratio)*(s_ratio);
#endif
    incident->S += h* stop_sample(ws, incident, sample, GSTO_STO_STRAGG, halfdepth, (E + dE / 2)); /* Straggling, calculate at mid-energy */

    assert(isnormal(incident->S));
    incident->E += dE;
    return fulldepth; /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}

double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, size_t i_isotope, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after) {
   if(ws->mean_conc_and_energy) {
        const depth d_halfdepth = {.x = (d_before->x + d_after->x)/2.0, .i = d_after->i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double c = get_conc(sample, d_halfdepth, i_isotope);
        if(c < ABUNDANCE_THRESHOLD)
            return 0.0;
        const double E_mean = (E_front + E_back) / 2.0;
        double sigma = sim_reaction_cross_section_rbs(sim_r, E_mean);
        return sigma*c;
    } else {

   }
    return 0.0;
}

void simulate(const ion *incident, const double x_0, sim_workspace *ws, const sample *sample) { /* Ion is expected to be in the sample system at depth x_0 */
    assert(sample->n_ranges);
    double thickness = sample->ranges[sample->n_ranges-1].x;
    size_t i_depth;
    ion ion1 = *incident; /* Shallow copy of the incident ion */
    double theta, phi; /* Generic polar and azimuth angles */
    double scatter_theta, scatter_phi;
    rotate(ws->sim.det.theta, ws->sim.det.phi, ws->sim.sample_theta, ws->sim.sample_phi, &theta, &phi); /* Detector in sample coordinate system */
    rotate(ws->sim.det.theta, ws->sim.det.phi, ion1.theta, ion1.phi, &scatter_theta, &scatter_phi); /* Detector in ion system */
    rotate(scatter_theta, scatter_phi, -ws->sim.sample_theta, -ws->sim.sample_phi, &scatter_theta, &scatter_phi); /* Counter sample rotation. Detector in lab (usually). If ion was somehow "deflected" then this is the real scattering angle. Compare to sim->theta.  */
    double K_min = 1.0;
    depth d_before = depth_seek(sample, x_0);
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = &ws->reactions[i];
        r->theta = scatter_theta;
        sim_reaction_recalculate_internal_variables(r);
        ion *p = &r->p;
        p->E = ion1.E * r->K;
        p->S = 0.0;
        r->max_depth = sample_isotope_max_depth(sample, r->i_isotope);
        ion_set_angle(p, theta, phi); /* Reaction products travel towards the detector (in the sample system), calculated above */
        brick *b = &r->bricks[0];
        b->E = ion1.E * r->K;
        b->S = 0.0;
        b->d = d_before;
        b->E_0 = ion1.E;
        if(r->K < K_min)
            K_min = r->K;
        if(r->i_isotope >= sample->n_isotopes) { /* No target isotope for reaction. */
            r->stop = TRUE;
            b->Q = -1.0;
        } else {
            r->stop = FALSE;
            b->Q = 0.0;
        }
#ifdef DEBUG
        fprintf(stderr, "Simulation reaction %zu: %s. Max depth %g tfu. i_isotope=%zu, stop = %i.\n", i, r->r->target->name, r->max_depth / C_TFU, r->i_isotope, r->stop);
#endif
    }
    assert(K_min > 0.0);
    i_depth=1;

#ifdef DEBUG
    ion_print(stderr, incident);
    fprintf(stderr, "Reaction product angles in sample system: (%g deg, %g deg)\n", theta/C_DEG, phi/C_DEG);
    fprintf(stderr, "Detector angles in ion system: (%g deg, %g deg). Sim theta is %g deg\n", scatter_theta/C_DEG, scatter_phi/C_DEG, ws->sim.theta/C_DEG);
#endif
    while(d_before.x < thickness) {
        if (ion1.E < ws->sim.emin) {
#ifdef DEBUG
            fprintf(stderr, "Break due to low energy (%.3lf keV < %.3lf keV), x = %.3lf, i_range = %lu.\n", ion1.E/C_KEV, ws->sim.emin/C_KEV, d_before.x/C_TFU, d_before.i);
#endif
            break;
        }
        const double E_front = ion1.E;
        const double S_front = ion1.S;
        depth d_after = stop_step(ws, &ion1, sample, d_before, ws->sim.stop_step_incident == 0.0?sqrt(ws->sim.det.resolution+K_min*(ion1.S)):ws->sim.stop_step_incident);
#ifdef DEBUG
        fprintf(stderr, "After:  %g tfu in range %zu\n", d_after.x/C_TFU, d_after.i);
#endif
        const double d_diff = depth_diff(d_before, d_after);
        if(d_diff < 0.001*C_TFU) {
            fprintf(stderr, "Error: no progress was made (E = %g keV, depth = %g tfu), check stopping.\n", ion1.E/C_KEV, d_before.x/C_TFU);
            break;
        }
        /* DEPTH BIN [x, x+h) */
        const double E_back = ion1.E;
        const double S_back = ion1.S;
#ifdef DEBUG_VERBOSE
        double E_diff = E_front-E_back;
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3lf keV to  %8.3lf keV (diff %6.4lf keV)\n", x/C_TFU, (x+h)/C_TFU, E_front/C_KEV, ws->ion.E/C_KEV, E_diff/C_KEV);
#endif
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(ion1.S) / C_KEV);
#endif
        for (size_t i = 0; i < ws->n_reactions; i++) {
            sim_reaction *r = &ws->reactions[i];
            if (r->stop)
                continue;
            if (i_depth >= r->n_bricks) {
                fprintf(stderr, "Too many bricks. Data partial.\n");
                r->stop = TRUE;
                continue;
            }
            brick *b = &r->bricks[i_depth];
            r->p.E = ion1.E * r->K;
            r->p.S = ion1.S * r->K;
            b->d = d_after;
            b->E_0 = ion1.E; /* Sort of energy just before the reaction. */
            assert(r->p.E > 0.0);
            if (d_before.x >= r->max_depth) {
#ifdef DEBUG
                fprintf(stderr, "Reaction %lu with %s stops, because maximum depth is reached at x = %.3lf tfu.\n",
                        i, r->r->target->name, d_before.x / C_TFU); /* TODO: give reactions a name */
#endif
                b->Q = -1.0;
                r->stop = TRUE;
                continue;
            }
            depth d = d_after;
            depth d_exit;
#ifdef DEBUG_REACTION
            fprintf(stderr, "Reaction %s (%zu): %s\n", reaction_name(r->r), i, r->r->target->name);
#endif
            while(1) {
#ifdef DEBUG_REACTION
                fprintf(stderr, "  Exiting... depth = %g tfu (i = %zu)\n", d.x, d.i);
#endif
                d_exit = stop_step(ws, &r->p, sample, d, ws->sim.stop_step_exiting == 0.0?r->p.E*0.1+sqrt(r->p.S):ws->sim.stop_step_exiting); /* TODO: 10% of energy plus straggling is a weird rule. Automatic stop size should be based more on required accuracy in stopping. */
                if(r->p.E < ws->sim.emin) {
#ifdef DEBUG_REACTION
                    fprintf(stderr,
                            "  Reaction %lu with %s: Energy below EMIN when surfacing from %.3lf tfu, break break.\n",
                            i, r->r->target->name, d_after.x / C_TFU);
#endif
                    break;
                }
                assert(d_exit.x < d.x /*|| (d_exit.x == d.x && d_exit.i != d.i)*/); /* Going towards the surface */
                if(d_exit.x <= 0.0) {
                    break;
                }
                d = d_exit;
            }
            b->E = r->p.E; /* Now exited from sample */
            b->S = r->p.S;
            if (r->p.E > ws->sim.emin) {
                double sigma_conc = cross_section_concentration_product(ws, sample, r->i_isotope, r, E_front, E_back,
                                                                        &d_before, &d_after); /* Product of concentration and sigma for isotope i_isotope target and this reaction. */
                if(sigma_conc > 0.0) {
                    if(d_after.i == sample->n_ranges - 2) {
                        sigma_conc *= ws->sim.channeling;
                    }
                    b->Q = fabs(incident->inverse_cosine_theta) * sigma_conc * d_diff; /* Note that ion is not the same as incident anymore. Incident has the original angles. */
#ifdef DEBUG
                    fprintf(stderr, "    %s: type=%i, E_front = %.3lf, E_after = %.3lf, E_out = %.3lf (sigma*conc = %g mb/sr, Q = %g (thickness = %.4lf tfu)\n",
                                 r->r->target->name, r->r->type, E_front/C_KEV, (ion1.E * r->K)/C_KEV, r->p.E/C_KEV, sigma_conc/C_MB_SR, b->Q, d_diff/C_TFU);
#endif
                 } else {
                    b->Q = 0.0;
                }
            } else {
                r->stop = TRUE;
                r->bricks[i_depth].Q = -1.0;
#ifdef DEBUG
                fprintf(stderr, "This was last depth for this reaction (and it didn't count anymore)\n");
#endif
            }
        }
        d_before = d_after;
        ion1.S = S_back;
        ion1.E = E_back;
        i_depth++;
    }
#ifdef DEBUG
    fprintf(stderr, "Last depth bin (brick) i_depth = %lu\n", i_depth);
#endif
    for (size_t i = 0; i < ws->n_reactions; i++) {
        if (ws->reactions[i].stop)
            continue;
        if (i_depth < ws->reactions[i].n_bricks)
            ws->reactions[i].bricks[i_depth].Q = -1.0; /* Set the last counts to negative to indicate end of calculation */
    }
    convolute_bricks(ws);
}

reaction *make_reactions(const sample *sample, const simulation *sim, jibal_cross_section_type cs_rbs, jibal_cross_section_type cs_erd) { /* Note that sim->ion needs to be set! */
    int rbs = (cs_rbs != JIBAL_CS_NONE);
    int erd = (cs_erd != JIBAL_CS_NONE);
    if(sim->theta > C_PI/2.0) {
        erd = FALSE;
        cs_erd = JIBAL_CS_NONE; /* Default when ERD is not possible :) */
    }
    size_t n_reactions = (sample->n_isotopes*rbs + sample->n_isotopes*erd + 1); /* TODO: we can predict this more accurately */
    reaction *reactions = malloc(n_reactions*sizeof(reaction));
    reaction *r = reactions;
    if(rbs) {
        for (size_t i = 0; i < sample->n_isotopes; i++) {
            *r = reaction_make(sim->beam_isotope, sample->isotopes[i], REACTION_RBS, cs_rbs, sim->theta);
            if (r->type == REACTION_NONE) {
                fprintf(stderr, "Failed to make an RBS reaction %lu (with %s)\n", i, sample->isotopes[i]->name);
            } else {
                r++;
            }
        };
    }
    if(erd) {
        for (size_t i = 0; i < sample->n_isotopes; i++) {
            *r = reaction_make(sim->beam_isotope, sample->isotopes[i], REACTION_ERD, cs_erd, sim->theta);
            if (r->type == REACTION_NONE) {
                fprintf(stderr, "Failed to make an ERD reaction %lu (with %s)\n", i, sample->isotopes[i]->name);
            } else {
                r++;
            }
        };
    }
    r->type = REACTION_NONE; /* Last reaction is a dummy one */
    return reactions;
}

int assign_stopping(jibal_gsto *gsto, const simulation *sim, const sample *sample, const reaction *reactions) {
    for(size_t i = 0; i < sample->n_isotopes; i++) {
        int Z2 = sample->isotopes[i]->Z;
        if (!jibal_gsto_auto_assign(gsto, sim->beam_isotope->Z, Z2)) { /* This should handle RBS */
            fprintf(stderr, "Can not assign stopping.\n");
            return 1;
        }
        for(const reaction *r = reactions; r->type != REACTION_NONE; r++) {
            if(r->type == REACTION_ERD) {
                if (!jibal_gsto_auto_assign(gsto, r->target->Z, Z2)) {
                    fprintf(stderr, "Can not assign stopping.\n");
                    return 1;
                }
            }
        }
    }
    return 0;
}

int print_spectra(const char *filename, const sim_workspace *ws, const gsl_histogram *exp) {
    char sep = ' ';
    FILE *f;
    if(filename) {
        f = fopen(filename, "w");
        if(!f) {
            fprintf(stderr, "Can't open file \"%s\" for output.\n", filename);
            return EXIT_FAILURE;
        }
        size_t l = strlen(filename);
        if(l > 4 && strncmp(filename+l-4, ".csv", 4) == 0) { /* For CSV: print header line */
            sep = ','; /* and set the separator! */
            fprintf(f, "\"Channel\",\"Simulated\"");
            if(exp) {
                fprintf(f, ",\"Experimental\",\"Energy (keV)\"");
            }
            for(size_t j = 0; j < ws->n_reactions; j++) {
                const reaction *r = ws->reactions[j].r;
                fprintf(f, ",\"%s (%s)\"", r->target->name, reaction_name(r));
            }
            fprintf(f, "\n");
        }
    } else {
        f = stdout;
    }
    for(size_t i = 0; i < ws->n_channels; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < ws->n_reactions; j++) { /* Sum comes always first, which means we have to compute it first. */
            if(i < ws->reactions[j].histo->n)
                sum += ws->reactions[j].histo->bin[i];
        }
        if(sum == 0.0) {
            fprintf(f, "%lu%c0", i, sep); /* Tidier output with a clean zero */
        } else {
            fprintf(f, "%lu%c%e", i, sep, sum);
        }
        if(exp) {
            if(i < exp->n) {
                fprintf(f, "%c%g", sep, exp->bin[i]);
            } else {
                fprintf(f, "%c0", sep);
            }
            fprintf(f,"%c%g", sep, exp->range[i]/C_KEV);
        }
        for (size_t j = 0; j < ws->n_reactions; j++) {
            if(i >= ws->reactions[j].histo->n || ws->reactions[j].histo->bin[i] == 0.0) {
                fprintf(f,"%c0", sep);
            } else {
                fprintf(f, "%c%e", sep, ws->reactions[j].histo->bin[i]);
            }
        }
        fprintf(f, "\n");
    }
    if(f != stdout) {
        fclose(f);
    } else {
        fprintf(f, "\n\n");
    }
    return EXIT_SUCCESS;
}

void add_fit_params(global_options *global, simulation *sim, const sample_model *sm, fit_params *params) {
#ifdef DEBUG
    fprintf(stderr, "fitvars = %s\n", global->fit_vars);
#endif
    if(!global->fit_vars)
        return;
    char *token, *s, *s_orig;
    s_orig = s = strdup(global->fit_vars);
    assert(s != NULL);
    while ((token = strsep(&s, ",")) != NULL) { /* parse comma separated list of parameters to fit */
#ifdef DEBUG
        fprintf(stderr, "Thing to fit: \"%s\"\n", token);
#endif
        if(strncmp(token, "calib", 5) == 0) {
            fit_params_add_parameter(params, &sim->det.slope); /* TODO: prevent adding already added things */
            fit_params_add_parameter(params, &sim->det.offset);
            fit_params_add_parameter(params, &sim->det.resolution);
        }
        if(strcmp(token, "slope") == 0) {
            fit_params_add_parameter(params, &sim->det.slope);
        }
        if(strcmp(token, "offset") == 0) {
            fit_params_add_parameter(params, &sim->det.offset);
        }
        if(strncmp(token, "reso", 4) == 0) {
            fit_params_add_parameter(params, &sim->det.resolution);
        }
        if(strcmp(token, "fluence") == 0) {
            fit_params_add_parameter(params, &sim->p_sr);
        }
        if(strcmp(token, "channeling") == 0) {
            fit_params_add_parameter(params, &sim->channeling);
        }
        if(strncmp(token, "rough", 5) == 0 && strlen(token) > 5) {
            size_t i_layer = strtoul(token+5, NULL, 10);
            if(i_layer >= 1 && i_layer <= sm->n_ranges) {
                fit_params_add_parameter(params, &sm->ranges[i_layer-1].rough.x);
            } else {
                fprintf(stderr, "No layer %zu (parsed from \"%s\")\n", i_layer, token);
            }
        }
        if(strncmp(token, "thickness", 9) == 0 && strlen(token) > 9) {
            size_t i_layer = strtoul(token+9, NULL, 10);
            if(i_layer >= 1 && i_layer <= sm->n_ranges) {
                fit_params_add_parameter(params, &sm->ranges[i_layer-1].x);
            } else {
                fprintf(stderr, "No layer %zu (parsed from \"%s\")\n", i_layer, token);
            }
        }
        size_t i,j;
        if(sscanf(token, "conc%lu_%lu", &i, &j) == 2) {
            if (i >= 1 && i <= sm->n_ranges && j >= 1  && j <= sm->n_materials) {
                fit_params_add_parameter(params, sample_model_conc_bin(sm, i-1, j-1));
            } else {
                fprintf(stderr, "No element %lu in layer %lu\n", j, i);
            }
        }
    }
    free(s_orig);
}

void output_bricks(const char *filename, const sim_workspace *ws) {
    FILE *f;
    if(!filename)
        return;
    if(strcmp(filename, "-") == 0)
        f=stdout;
    else {
        f = fopen(filename, "w");
    }
    if(!f)
        return;
    for(size_t i = 0; i < ws->n_reactions; i++) {
        const sim_reaction *r = &ws->reactions[i];
        fprintf(f, "#%s %s\n", reaction_name(r->r), r->r->target->name);
        for(size_t j = 0; j < r->n_bricks; j++) {
            brick *b = &r->bricks[j];
            if(b->Q < 0.0)
                break;
            fprintf(f, "%2lu %2lu %8.3lf %8.3lf %8.3lf %8.3lf %12.3lf\n",
                    i, j, b->d.x/C_TFU, b->E_0/C_KEV, b->E/C_KEV, sqrt(b->S)/C_KEV, b->Q * ws->sim.p_sr);
        }
        fprintf(f, "\n\n");
    }
    if(f != stdout)
        fclose(f);
}



void no_ds(sim_workspace *ws, const sample *sample) {
    double p_sr = ws->sim.p_sr;
    size_t n_rl = 0; /* Number of rough layers */
    for(size_t i = 0; i < sample->n_ranges; i++) {
        if(sample->ranges[i].rough.model != ROUGHNESS_NONE)
            n_rl++;
    }
#ifdef DEBUG
    fprintf(stderr, "%zu rough layers\n", n_rl);
#endif
    if(!n_rl) {
        simulate(&ws->ion, 0.0, ws, sample);
        return;
    }
    struct sample *sample_rough = sample_copy(sample);
    size_t *index = malloc(sizeof(size_t) * n_rl);
    size_t *modulos = malloc(sizeof(size_t) * n_rl);
    size_t j = 0;
    thick_prob_dist **tpd = malloc(sizeof(thick_prob_dist *) * n_rl);
    for(size_t i = 0; i < sample->n_ranges; i++) {
        if(sample->ranges[i].rough.model != ROUGHNESS_NONE) {
#ifdef DEBUG
            fprintf(stderr, "Range %zu is rough, model %i, amount %g tfu, n = %zu spectra\n", i, sample->ranges[i].rough.model, sample->ranges[i].rough.x/C_TFU, sample->ranges[i].rough.n);
#endif
            assert(sample->ranges[i].rough.n > 0);
            tpd[j] = thickness_probability_table_gen(sample->ranges[i].x, sample->ranges[i].rough.x, sample->ranges[i].rough.n);
            index[j] = i;
            if(j)
                modulos[j] = modulos[j-1] *  tpd[j-1]->n;
            else
                modulos[j] = 1;
            j++;
        }
    }
    size_t iter_total = modulos[n_rl-1] * tpd[n_rl-1]->n;
    for(size_t i_iter = 0; i_iter < iter_total; i_iter++) {
#ifdef DEBUG
        fprintf(stderr, "Gamma roughness step %zu/%zu\n", i_iter+1, iter_total);
#endif
        double p = 1.0;
        for(size_t i_range = 0; i_range < sample->n_ranges; i_range++) { /* Reset ranges for every iter */
            sample_rough->ranges[i_range].x = sample->ranges[i_range].x;
        }
        for(size_t i = 0; i < n_rl; i++) {
            j = (i_iter / modulos[i]) % tpd[i]->n; /* "j"th roughness element */
            //fprintf(stderr, " %zu", j);

            size_t i_range = index[i];
            p *= tpd[i]->p[j].prob; /* Probability is multiplied by the "i"th roughness, element "j" */
            double x_diff = tpd[i]->p[j].x - sample->ranges[i_range].x; /* Amount to change thickness of this and and all subsequent layers */
            for(; i_range < sample->n_ranges; i_range++) {
                sample_rough->ranges[i_range].x += x_diff;
            }
#ifdef DEBUG
            fprintf(stderr, "Gamma roughness diff %g tfu (from %g tfu, index i_range=%zu), probability %.3lf%%)\n", x_diff/C_TFU, sample->ranges[i_range].x/C_TFU, i_range, tpd[i]->p[j].prob*100.0);
            fprintf(stderr, "Gamma roughness, ranges");
            for(i_range = 0; i_range < sample->n_ranges; i_range++) {
                fprintf(stderr, ", %zu: %g tfu ", i_range, sample_rough->ranges[i_range].x/C_TFU);
            }
            fprintf(stderr, "\n");
#endif
        }
        //fprintf(stderr, "\n");
        ws->sim.p_sr = p * p_sr;
        ion_set_angle(&ws->ion, 0.0, 0.0);
        ion_rotate(&ws->ion, ws->sim.sample_theta, ws->sim.sample_phi);
        simulate(&ws->ion, 0.0, ws, sample_rough);
    }
    for(size_t i = 0; i < n_rl; i++) {
        thickness_probability_table_free(tpd[i]);
    }
    free(modulos);
    free(index);
    free(tpd);
}

void ds(sim_workspace *ws, const sample *sample) { /* TODO: the DS routine is more pseudocode at this stage... */
    ion_set_angle(&ws->ion, 0.0, 0.0);
    ion_rotate(&ws->ion, ws->sim.sample_theta, ws->sim.sample_phi);
    ion ion1 = ws->ion;
    depth d = depth_seek(sample, 0.0);
    while(1) {
        depth d_after = stop_step(ws, &ws->ion, sample, d, ws->sim.stop_step_incident);
        ion ion2 = ion1;

        int i_ds;
        for(i_ds = 0; i_ds < ws->sim.n_ds; i_ds++) {
            int i_polar = i_ds / ws->sim.ds_steps_azi;
            double cosine = (ws->sim.ds_steps_polar - 2 * i_polar) / (1.0 * ws->sim.ds_steps_polar);
            double ds_polar = acos(cosine);
            int i_azi = i_ds % ws->sim.ds_steps_azi;
            double ds_azi = 2 * M_PI * (1.0 * i_azi) / (ws->sim.ds_steps_azi * 1.0);

            ion_rotate(&ion2, ds_polar, ds_azi); /* Dual scattering: first scattering to some angle (scattering angle: ds_polar). Note that this does not follow SimNRA conventions. */
            /*TODO: is this correct ?*/

            /* TODO: reset/init ws? */
            simulate(&ion2, d.x, ws, sample); /* d.x or d_after.x?
 * TODO: simulate() cannot handle something that goes towards the surface. Improve stop_step() to handle both cases and remove i_range magic from elsewhere.
 * TODO: ion_rotate(ion, -ds_polar, -ds_azi); aka Undo scattering rotation. Is simulate() able to figure out where the detector is? */
            convolute_bricks(ws); /* TODO: does this work */
        }
        if(depth_diff(d,d_after) == 0.0)
            break;
        d = d_after;
    }
}
