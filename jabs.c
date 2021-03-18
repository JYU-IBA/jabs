#include <string.h>
#include <assert.h>
#include <math.h>
#include <jibal_units.h>
#include <jibal_kin.h>


#include "defaults.h"
#include "rotate.h"
#include "jabs.h"

double stop_sample(sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, double x, double E, size_t *range_hint) {
    double em=E/incident->mass;
    double S1 = 0.0;
    get_concs(sample, x, ws->c, range_hint);
    for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
        if(ws->c[i_isotope] < ABUNDANCE_THRESHOLD)
            continue;
        if (type == GSTO_STO_TOT) {
            S1 += ws->c[i_isotope] * (
                    jibal_gsto_get_em(ws->gsto, GSTO_STO_ELE, incident->Z, sample->isotopes[i_isotope]->Z, em)
                    #ifdef NUCLEAR_STOPPING_FROM_JIBAL
                    +jibal_gsto_stop_nuclear_universal(E, incident->Z, incident->mass, sample->isotopes[i_isotope]->Z, sample->isotopes[i_isotope]->mass)
                    #else
                    + ion_nuclear_stop(incident, sample->isotopes[i_isotope], ws->isotopes, ws->nucl_stop_accurate)
                    #endif
                    );
        } else {
            S1 += ws->c[i_isotope] * (
                    jibal_gsto_get_em(ws->gsto, type, incident->Z, sample->isotopes[i_isotope]->Z, em)
            );
        }
    }
    //assert(S1 > 0.0);
    return S1;
}

double stop_step(sim_workspace *ws, ion *incident, const sample *sample, double x, double step, const size_t range_hint) {
    double k1, k2, k3, k4, stop, dE, E;
    double h_max_perp;
    assert(range_hint < sample->n_ranges-1);
    if(incident->inverse_cosine_theta < 0.0) {
        h_max_perp = sample->cranges[range_hint] - x;
    } else {
        h_max_perp = sample->cranges[range_hint + 1] - x;
    }
#if 0
    if(h_max_perp != h_max_orig) {
        fprintf(stderr, "h_max_orig = %g, h_max_perp = %g\n", h_max_orig, h_max_perp);
    }
#endif
    size_t range = range_hint;
    /* k1...k4 are slopes of energy loss (stopping) at various x (depth) and E. Note convention: positive values, i.e. -dE/dx! */
    E = incident->E;
    k1 = stop_sample(ws, incident, sample, ws->stopping_type, x, E, &range);
    if(k1 < 0.001*C_EV_TFU) { /* Fail on positive values, zeroes (e.g. due to zero concentrations) and too small negative values */
#ifdef DEBUG
        fprintf(stderr, "stop_step returns 0.0, because k1 = %g eV/tfu (x = %.3lf tfu, E = %.3lg keV)\n", k1/C_EV_TFU, x/C_TFU, E/C_KEV);
#endif
        return 0.0;
    }
    double h_max = h_max_perp * incident->inverse_cosine_theta; /*  we can take bigger steps since we are going sideways. Note that inverse_cosine_theta can be negative and in this case h_max should also be negative so h_max is always positive! */
    assert(h_max > 0.0);
    double h =  (step / k1); /* (energy) step should always be positive, as well as k1, so depth step h (not perpendicular, but "real" depth) is always positive  */
    assert(h > 0.0);
    double h_perp; /* has a sign (same as h_max_perp ) */
    if(h >= h_max) {
        h = h_max;
        h_perp = h_max_perp;
    } else {
        h_perp = h*incident->cosine_theta; /* x + h_perp is the actual perpendicular depth */
    }
    if(ws->rk4) {
        k2 = stop_sample(ws, incident, sample, ws->stopping_type, x + (h_perp / 2.0), E - (h / 2.0) * k1, &range);
        k3 = stop_sample(ws, incident, sample, ws->stopping_type, x + (h_perp / 2.0), E - (h / 2.0) * k2, &range);
        k4 = stop_sample(ws, incident, sample, ws->stopping_type, x + h_perp, E - h * k3, &range);
        stop = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    } else {
        stop = k1;
    }
    assert(stop > 0.0);
    dE =  -1.0* h * stop; /* Energy change in thickness "h". It is always negative! */
#if 0
    fprintf(stderr, "%s stop = %.3lf eV/tfu ( x = %.3lf, h = %.3lf, h_max = %.3lf), E = %.3lf keV, h = %6.3lf  dE = %.5lf keV\n", incident->isotope->name, stop/C_EV_TFU, x/C_TFU, h/C_TFU, h_max/C_TFU, E/C_KEV, h/C_TFU, dE/C_KEV);
#endif
#ifdef DEBUG
    if(fabs(stop) < 0.1*C_EV_TFU) {
        fprintf(stderr, "Not good!\n");
        return 0.0;
    }
#endif
#ifndef STATISTICAL_STRAGGLING
    double s_ratio = stop_sample(ws, incident, sample, ws->stopping_type, x, E + dE, &range) / k1; /* Ratio of stopping for non-statistical broadening. TODO: at x? */
#ifdef DEBUG
    //if((s_ratio)*(s_ratio) < 0.9 || (s_ratio)*(s_ratio) > 1.1) { /* Non-statistical broadening. */
    //   fprintf(stderr, "YIKES, s_ratio = %g, sq= %g\n", s_ratio, (s_ratio)*(s_ratio));
    //}
#endif
    incident->S *= (s_ratio)*(s_ratio);
#endif
    incident->S += h* stop_sample(ws, incident, sample, GSTO_STO_STRAGG, x + (h_perp / 2.0), (E + dE / 2), &range); /* Straggling, calculate at mid-energy */

    assert(isnormal(incident->S));
    incident->E += dE;
    return h_perp; /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}

void simulate(const ion *incident, const double x_0, sim_workspace *ws, const sample *sample) {
    double x;
    assert(sample->n_ranges);
    double thickness = sample->cranges[sample->n_ranges-1];
    double next_crossing = sample->cranges[1];
    double h_max;
    size_t i_range = 0;
    size_t i_depth;
    ion ion1 = *incident; /* Shallow copy of the incident ion */
    double theta, phi;
    rotate(ws->sim.det.theta, ws->sim.det.theta, ws->sim.sample_theta, ws->sim.sample_phi, &theta, &phi); /* Detector in sample coordinate system */
    double K_min = 1.0;
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = &ws->reactions[i];
        ion *p = &r->p;
        p->E = ion1.E * r->r->K;
        p->S = 0.0;
        ion_set_angle(p, theta, phi);
        r->stop = 0;
        brick *b = &r->bricks[0];
        b->E = ion1.E * r->r->K;
        b->S = 0.0;
        b->d = 0.0;
        b->Q = 0.0;
        b->E_0 = ion1.E;
        if(r->r->K < K_min)
            K_min = r->r->K;
    }
    assert(K_min > 0.0);
    i_depth=1;

#ifdef DEBUG
    ion_print(stderr, incident);
#endif
    for (x = x_0; x < thickness;) {
        while (i_range < sample->n_ranges - 1 && x >= sample->cranges[i_range + 1]) {
            i_range++;
#ifdef DEBUG
            fprintf(stderr, "Crossing to range %lu = [%g, %g)\n", i_range, sample->cranges[i_range] / C_TFU,
                    sample->cranges[i_range + 1] / C_TFU);
#endif
            next_crossing = sample->cranges[i_range + 1];
        }
        if (ion1.E < ws->sim.emin) {
#ifdef DEBUG
            fprintf(stderr, "Break due to low energy (%.3lf keV < %.3lf keV), x = %.3lf, i_range = %lu.\n", ion1.E/C_KEV, ws->sim.emin/C_KEV, x/C_TFU, i_range);
#endif
            break;
        }

        h_max = next_crossing - x;
        assert(h_max > 0.001 * C_TFU);
        double E_front = ion1.E;
        double h = stop_step(ws, &ion1, sample, x, ws->sim.stop_step_incident == 0.0?sqrt(ws->sim.det.resolution+K_min*(ion1.S)):ws->sim.stop_step_incident, i_range);
        assert(h > 0.0);
        /* DEPTH BIN [x, x+h) */
        double E_back = ion1.E;
#ifdef DEBUG_VERBOSE
        double E_diff = E_front-E_back;
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3lf keV to  %8.3lf keV (diff %6.4lf keV)\n", x/C_TFU, (x+h)/C_TFU, E_front/C_KEV, ws->ion.E/C_KEV, E_diff/C_KEV);
#endif
        double S_back = ion1.S;
        double E_mean = (E_front + E_back) / 2.0;
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(ion1.S) / C_KEV);
#endif

        //ion ion2 = *incident; /* Make a shallow copy of ion */
        //ion_rotate(&ion2, -ws->sim.sample_theta, -ws->sim.sample_phi); /* Rotate to lab TODO: all this rotating is unnecessary, angles stay fixed */
        //ion_rotate(&ion2, ws->sim.detector_theta, ws->sim.detector_phi); /* Rotate to detector (we can determine the scattering angle from this!) */


        for (size_t i = 0; i < ws->n_reactions; i++) {
            sim_reaction *r = &ws->reactions[i];
            if (r->stop)
                continue;
            if (i_depth >= r->n_bricks) {
                r->stop = 1;
                continue;
            }
            brick *b = &r->bricks[i_depth];
            //ion_set_angle(&r->p, theta_after_second, phi_after_second);

            r->p.E = ion1.E * r->r->K;
            r->p.S = ion1.S * r->r->K;
            b->d = x + h;
            b->E_0 = ion1.E; /* Sort of energy just before the reaction. */

            assert(r->p.E > 0.0);

            if (x >= r->r->max_depth) {
#ifdef DEBUG
                fprintf(stderr, "Reaction %lu with %s stops, because maximum depth is reached at x = %.3lf tfu.\n",
                        i, r->r->isotope->name, x / C_TFU); /* TODO: give reactions a name */
#endif
                b->Q = -1.0;
                r->stop = 1;
                continue;
            }

            double x_out;
            size_t i_range_out = i_range;
            for (x_out = x + h; x_out > 0.0;) { /* Calculate energy and straggling of backside of slab */
                double h_out = stop_step(ws, &r->p, sample, x_out, ws->sim.stop_step_exiting == 0.0?r->p.E*0.1+sqrt(r->p.S):ws->sim.stop_step_exiting, i_range_out); /* TODO: 10% of energy plus straggling is a weird rule. Automatic stop size should be based more on required accuracy in stopping. */
                x_out += h_out;
                if (r->p.E < ws->sim.emin) {
#ifdef DEBUG
                    fprintf(stderr,
                            "Reaction %lu with %s: Energy below EMIN when surfacing from %.3lf tfu, break break.\n",
                            i, r->r->isotope->name, (x + h) / C_TFU);
#endif
                    break;
                }
                assert(h_out < 0.0);
                while (i_range_out > 0 && x_out <= sample->cranges[i_range_out]) {
                    i_range_out--;
#ifdef DEBUG_VERBOSE
                    fprintf(stderr, "Outgoing from reaction %i crossing to range %i = [%g, %g) when x_out = %g tfu.\n", i_reaction, i_range_out, sample->cranges[i_range_out]/C_TFU, sample->cranges[i_range_out+1]/C_TFU, x_out/C_TFU);
#endif
                }
            }
            double c = get_conc(sample, x + (h / 2.0), r->r->i_isotope, &i_range); /* TODO: x+h/2.0 is actually exact for linearly varying concentration profiles. State this clearly somewhere. */
            b->E = r->p.E; /* Now exited from sample */
            b->S = r->p.S;
            if (c > ABUNDANCE_THRESHOLD && r->p.E > ws->sim.emin) {/* TODO: concentration cutoff? TODO: it->E should be updated when we start calculating it again?*/
                double sigma;
                switch (r->r->type) {
                    case REACTION_RBS:
                        sigma = jibal_cross_section_rbs(ion1.isotope, r->r->isotope, ws->sim.theta, E_mean, ws->jibal_config->cs_rbs);
                        break;
                    case REACTION_ERD:
                        sigma = jibal_cross_section_erd(ion1.isotope, r->r->isotope, ws->sim.theta, E_mean, ws->jibal_config->cs_erd);
                        break;
                    default:
                        sigma = 0.0;
                }
                double Q = c * fabs(incident->inverse_cosine_theta) * sigma * h; /* TODO: worst possible approximation... */ /* Note that ion is not the same as incident anymore. Incident has the original angles. */
#ifdef DEBUG
                fprintf(stderr, "    %s: type=%i, E_mean = %.3lf, E_after = %.3lf, E_out = %.3lf (sigma = %g mb/sr, Q = %g (c = %.4lf%%)\n",
                                 r->r->isotope->name, r->r->type, E_mean/C_KEV, (ion1.E * r->r->K)/C_KEV, r->p.E/C_KEV, sigma/C_MB_SR, Q, c*100.0);
#endif
                assert(sigma > 0.0);
                assert(r->p.S > 0.0);
                assert(Q < 1.0e7 || Q > 0.0);
                assert(i_depth < r->n_bricks);
                b->Q = Q;
            } else {
                if (r->p.E < ws->sim.emin) {
                    r->stop = 1;
                    r->bricks[i_depth].Q = -1.0;
                } else {
                    b->Q = 0.0;
                }
            }
        }
        x += h;
        ion1.S = S_back;
        ion1.E = E_back;
        i_depth++;
    }
#ifdef DEBUG
    fprintf(stderr, "Last depth bin %lu\n", i_depth);
#endif
    for (size_t i = 0; i < ws->n_reactions; i++) {
        if (ws->reactions[i].stop)
            continue;
        if (i_depth < ws->reactions[i].n_bricks)
            ws->reactions[i].bricks[i_depth].Q = -1.0; /* Set the last counts to negative to indicate end of calculation */
    }
    convolute_bricks(ws);
}

reaction make_reaction(const sample *sample, const simulation *sim, const size_t i_isotope, reaction_type type) {
#ifdef DEBUG
    fprintf(stderr, "Attempting to make reaction with sample isotope %lu (%s)\n", i_isotope, sample->isotopes[i_isotope]->name);
#endif
    reaction r;
    r.type = type;
    r.isotope = sample->isotopes[i_isotope];
    r.i_isotope = i_isotope;
    if(!r.isotope) {
        r.type = REACTION_NONE;
        return r;
    }
    r.i_isotope = i_isotope;
    if(type == REACTION_RBS) {
        double theta_max=asin(r.isotope->mass/sim->beam_isotope->mass);
        if(sim->beam_isotope->mass >= r.isotope->mass && sim->theta > theta_max) {
            fprintf(stderr, "RBS with %s is not possible (theta max %g deg, sim theta %g deg)\n", r.isotope->name, theta_max/C_DEG, sim->theta/C_DEG);
            r.type = REACTION_NONE;
            return r;
        }
        r.K = jibal_kin_rbs(sim->beam_isotope->mass, r.isotope->mass, sim->theta, '+');
    } else if (type == REACTION_ERD) {
        if(sim->theta > C_PI/2.0) {
            fprintf(stderr, "ERD with %s is not possible (theta %g deg > 90.0 deg)", r.isotope->name, sim->theta);
            r.type = REACTION_NONE;
            return r;
        }
        r.K = jibal_kin_erd(sim->beam_isotope->mass, r.isotope->mass, sim->theta);
    } else {
        r.type = REACTION_NONE;
        return r;
    }
    r.max_depth = sample_isotope_max_depth(sample, r.i_isotope);
    return r;
}

reaction *make_reactions(const sample *sample, const simulation *sim, int rbs, int erd) { /* Note that sim->ion needs to be set! */
    if(sim->theta > C_PI/2.0)
        erd = 0; /* Default when ERD is not possible :) */
    size_t n_reactions = (sample->n_isotopes*rbs + sample->n_isotopes*erd + 1); /* TODO: we can predict this more accurately */
    reaction *reactions = malloc(n_reactions*sizeof(reaction));
    reaction *r = reactions;
    if(rbs) {
        for (size_t i = 0; i < sample->n_isotopes; i++) {
            *r = make_reaction(sample, sim, i, REACTION_RBS);
            if (r->type == REACTION_NONE) {
                fprintf(stderr, "Failed to make an RBS reaction %lu (with %s)\n", i, sample->isotopes[i]->name);
            } else {
                r++;
            }
        };
    }
    if(erd) {
        for (size_t i = 0; i < sample->n_isotopes; i++) {
            *r = make_reaction(sample, sim, i, REACTION_ERD);
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
                if (!jibal_gsto_auto_assign(gsto, r->isotope->Z, Z2)) {
                    fprintf(stderr, "Can not assign stopping.\n");
                    return 1;
                }
            }
        }
    }
    return 0;
}

void print_spectra(FILE *f, const global_options *global, const sim_workspace *ws, const sample *sample, const gsl_histogram *exp) {
    char sep = ' ';
    if(global->out_filename) {
        size_t l = strlen(global->out_filename);
        if(l > 4 && strncmp(global->out_filename+l-4, ".csv", 4) == 0) { /* For CSV: print header line */
            sep = ','; /* and set the separator! */
            fprintf(f, "\"Channel\",\"Simulated\"");
            if(exp) {
                fprintf(f, ",\"Experimental\"");
            }
            for(size_t j = 0; j < ws->n_reactions; j++) {
                const reaction *r = ws->reactions[j].r;
                fprintf(f, ",\"%s (%s)\"", sample->isotopes[r->i_isotope]->name, reaction_name(r));
            }
            fprintf(f, "\n");
        }
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
}

void add_fit_params(global_options *global, simulation *sim, jibal_layer **layers, const int n_layers, fit_params *params) {
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
        if(strncmp(token, "thickness", 9) == 0 && strlen(token) > 9) {
            int i_layer = atoi(token+9);
            if(i_layer >= 1 && i_layer <= n_layers) {
                fit_params_add_parameter(params, &layers[i_layer-1]->thickness);
            } else {
                fprintf(stderr, "No layer %i (parsed from \"%s\")\n", i_layer, token);
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
        fprintf(f, "#%s %s\n", reaction_name(r->r), r->r->isotope->name);
        for(size_t j = 0; j < r->n_bricks; j++) {
            brick *b = &r->bricks[j];
            if(b->Q < 0.0)
                break;
            fprintf(f, "%2lu %2lu %8.3lf %8.3lf %8.3lf %8.3lf %12.3lf\n",
                    i, j, b->d/C_TFU, b->E_0/C_KEV, b->E/C_KEV, sqrt(b->S)/C_KEV, b->Q * ws->sim.p_sr);
        }
        fprintf(f, "\n\n");
    }
    if(f != stdout)
        fclose(f);
}



void no_ds(sim_workspace *ws, const sample *sample) {
    ion_set_angle(&ws->ion, 0.0, 0.0);
    ion_rotate(&ws->ion, ws->sim.sample_theta, ws->sim.sample_phi);
    simulate(&ws->ion, 0.0, ws, sample);
}

void ds(sim_workspace *ws, const sample *sample) { /* TODO: the DS routine is more pseudocode at this stage... */
    double thickness = sample->cranges[sample->n_ranges-1];
    /* Go deeper */
    double x;
    ion_set_angle(&ws->ion, 0.0, 0.0);
    ion_rotate(&ws->ion, ws->sim.sample_theta, ws->sim.sample_phi);
    ion ion1 = ws->ion;
    size_t i_range = 0;
    for(x = 0.0; x < thickness;) {
        /* Go deeper and at every step start making new spectra. */
        while (i_range < sample->n_ranges - 1 && x >= sample->cranges[i_range + 1]) {
            i_range++;
        }
        double h = stop_step(ws, &ws->ion, sample, x, ws->sim.stop_step_incident, i_range);
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
            simulate(&ion2, x, ws, sample); /*
 * TODO: simulate() cannot handle something that goes towards the surface. Improve stop_step() to handle both cases and remove i_range magic from elsewhere.
 * TODO: ion_rotate(ion, -ds_polar, -ds_azi); aka Undo scattering rotation. Is simulate() able to figure out where the detector is? */
            convolute_bricks(ws); /* TODO: does this work */
        }
        x += h;
    }
}
