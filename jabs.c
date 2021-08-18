#include <string.h>
#include <assert.h>
#include <math.h>
#include <jibal_units.h>
#include <jibal_kin.h>
#include <jibal_r33.h>

#include "generic.h"
#include "rotate.h"
#include "roughness.h"
#include "jabs.h"
#include "defaults.h"

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
                    + ion_nuclear_stop(incident, sample->isotopes[i_isotope], ws->isotopes, ws->params.nucl_stop_accurate)
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

depth next_crossing(const ion *incident, const sample *sample, const depth *d_from) {
    depth d = *d_from;
    if(incident->inverse_cosine_theta > 0) { /* Going deeper */
        while(d.i < sample->n_ranges - 1 && d.x >= sample->ranges[d.i + 1].x) {
            d.i++;
        }
        assert(d.i < sample->n_ranges - 1); /* There is a bug elsewhere in the code if you try to go this deep (deeper than last depth bin). */
        d.x = sample->ranges[d.i + 1].x;
    } else if(incident->inverse_cosine_theta < 0.0) { /* Going towards the surface */
        while(d.i > 0 && d.x <= sample->ranges[d.i].x) {
            d.i--;
        }
        d.x = sample->ranges[d.i].x;
    } else {
        fprintf(stderr, "WARNING: Inverse cosine is exactly zero. This is an issue!\n");
    }
    return d;
}

depth stop_step(sim_workspace *ws, ion *incident, const sample *sample, depth depth, double step) {
    double k1, k2, k3, k4, stop, dE, E;
    struct depth depth_next = next_crossing(incident, sample, &depth);
    double h_max_perp = depth_next.x - depth.x;
#ifdef DEBUG_STOP_STEP
    if(depth_next.i != depth.i) {
        fprintf(stderr, "stop_step crossing depth range from %zu to %zu at depth %lf tfu. E = %.3lf keV, Inverse cosine %lf\n", depth.i, depth_next.i, depth.x/C_TFU, incident->E/C_KEV, incident->inverse_cosine_theta);
    } else {
        fprintf(stderr, "stop_step depth %g tfu (i=%zu) distance to next crossing %g tfu.\n", depth.x/C_TFU, depth.i, h_max_perp/C_TFU);
    }
    assert(fabs(h_max_perp) > 0.01 * C_TFU);
#endif
    depth.i = depth_next.i; /* depth may have the right depth (.x), but the index can be old. We need to cross the depth range somewhere, and it happens to be here. All calculations after this take place inside the same depth range (index depth.i). */
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
    double h = (step / k1); /* (energy) step should always be positive, as well as k1, so depth step h (not perpendicular, but "real" depth) is always positive  */
    assert(h > 0.0);
    struct depth halfdepth;
    struct depth fulldepth;
    halfdepth.i = depth.i;
    if(h >= h_max) { /* Depth step would take us beyond a depth range. We stop exactly on the boundary */
        h = h_max;
        halfdepth.x = depth.x + h_max_perp/2.0;
        fulldepth = depth_next;
    } else {
        double h_perp = h*incident->cosine_theta; /* x + h_perp is the actual perpendicular depth */
        halfdepth.x = depth.x + h_perp/2.0;
        fulldepth.i = depth.i;
        fulldepth.x = depth.x + h_perp;
    }



    if(ws->params.rk4) {
        k2 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k1);
        k3 = stop_sample(ws, incident, sample, ws->stopping_type, halfdepth, E - (h / 2.0) * k2);
        k4 = stop_sample(ws, incident, sample, ws->stopping_type, fulldepth, E - h * k3);
        stop = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    } else {
        stop = k1;
    }
    assert(stop > 0.0);
    dE =  -1.0* h * stop; /* Energy change in thickness "h". It is always negative! */
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
    incident->E += dE;
    return fulldepth; /*  Stopping is calculated in material the usual way, but we only report progress perpendicular to the sample. If incident->angle is 45 deg, cosine is 0.7-ish. */
}


double normal_pdf(double x, double mean, double sigma) {
    static const double inv_sqrt_2pi = 0.398942280401432703;
    double a = (x - mean) / sigma;
    return (inv_sqrt_2pi / sigma) * exp(-0.5 * a * a);
}

double cross_section_straggling(const sim_reaction *sim_r, int n_steps, double E, double S) {
    const double sigmas = 2.0; /* TODO: what is enough or too much? */
    const double std_dev = sqrt(S);
    const int half_n = n_steps/2;
    const double w = sigmas/(half_n);
    double cs_sum = 0.0;
    double prob_sum = 0.0;
    //static const double inv_sqrt_2pi = 0.398942280401432703;
    for(int i = 0; i < n_steps; i++) {
        double x = w*(i-half_n);
        double E_stragg = E + x * std_dev;
        double prob = normal_pdf(x, 0.0, 1.0)*w; /* TODO: if this is always a normal distribution and n_steps doesn't change, this function call could be replaced by a lookup table */
        prob_sum += prob;
        cs_sum += prob * sim_r->cross_section(sim_r, E_stragg);
    }
    cs_sum /= prob_sum;
#ifdef DEBUG_CS_WEIGHT
    double unweighted = sim_r->cross_section(sim_r, E);
    double diff = (cs_sum-unweighted)/unweighted;
    fprintf(stderr, "Got cs %.7lf mb/sr (unweighted by straggling %.7lf mb/sr) diff %.5lf%%. Sum of probs %lf%% (compensated for).\n", cs_sum/C_MB_SR, unweighted/C_MB_SR, 100.0*diff, prob_sum*100.0);
#endif
    return cs_sum;
}


double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back) {
   if(ws->params.mean_conc_and_energy) { /* This if-branch is slightly faster (maybe) and also serves as a testing branch, since it is a lot easier to understand... */
        const depth d_halfdepth = {.x = (d_before->x + d_after->x)/2.0, .i = d_after->i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double c = get_conc(sample, d_halfdepth, sim_r->i_isotope);
        if(c < ABUNDANCE_THRESHOLD)
            return 0.0;
        const double E_mean = (E_front + E_back) / 2.0;
        assert(sim_r->cross_section);
        double sigma = sim_r->cross_section(sim_r, E_mean);
        return sigma*c;
    } else {
        depth d;
        d.i = d_after->i;
        const double x_step = (d_after->x - d_before->x) * ws->params.cs_frac;
        const double E_step = (E_back - E_front) * ws->params.cs_frac;
        const double S_step = (S_back - S_front) * ws->params.cs_frac;
        double sum = 0.0;
        for(int i = 1; i <= ws->params.cs_n_steps; i++) { /* Compute cross section and concentration product in several "sub-steps" */
            d.x = d_before->x + x_step * i;
            double E = E_front + E_step * i;
#ifdef DEBUG
            fprintf(stderr, "i=%i, E = %g keV, (E_front = %g keV, E_back = %g keV)\n", i, E/C_KEV, E_front/C_KEV, E_back/C_KEV);
#endif
            double c = get_conc(sample, d, sim_r->i_isotope);
            double sigma;
            if(ws->params.cs_n_stragg_steps > 1) { /* Further weighted with straggling */
                double S = S_front + S_step * i;
                sigma = cross_section_straggling(sim_r, ws->params.cs_n_stragg_steps, E, S);
            } else {
                sigma = sim_r->cross_section(sim_r, E);
            }
            sum += sigma * c;
        }
       return sum/(ws->params.cs_n_steps*1.0);
   }
    return 0.0;
}

void post_scatter_exit(ion *p, const depth depth_start, sim_workspace *ws, const sample *sample) {
    depth d = depth_start;
    while(1) { /* Exit from sample (hopefully) */
#ifdef DEBUG_REACTION
        fprintf(stderr, "  Exiting... depth = %g tfu (i = %zu)\n", d.x, d.i);
#endif
        if(d.x <= DEPTH_TOLERANCE) {
            break;
        }
        depth d_after = stop_step(ws, p, sample, d, ws->params.stop_step_exiting == 0.0?(p->E*0.10+sqrt(p->S)+2.0*C_KEV):ws->params.stop_step_exiting); /* TODO: 10% of energy plus straggling plus 2 keV is a weird rule. Automatic stop size should be based more on required accuracy in stopping. */
        if(p->E < ws->params.emin) {
#ifdef DEBUG_REACTION
            fprintf(stderr,
                            "  Reaction %lu with %s: Energy below EMIN when surfacing from %.3lf tfu, break break.\n",
                            i, r->r->target->name, d_after.x / C_TFU);
#endif
            return;
        }
        assert(d_after.x <= d.x /*|| (d_exit.x == d.x && d_exit.i != d.i)*/); /* Going towards the surface */
        d = d_after;
    }
    const struct sample *foil = ws->det->foil;
    if(!foil)
        return;
    d.i = 0;
    d.x = 0.0;
    ion ion_foil = *p;
    ion_set_angle(&ion_foil, 0.0, 0.0); /* Foils are not tilted. We use a temporary copy of "p" to do this step. */
    while(1) {
        if(d.x >= foil->ranges[foil->n_ranges-1].x) {
            break;
        }
        depth d_after  = stop_step(ws, &ion_foil, foil, d, ws->params.stop_step_exiting == 0.0?p->E*0.1+sqrt(p->S):ws->params.stop_step_exiting);
        if(p->E < ws->params.emin) {
            break;
        }
        d = d_after;
    }
    p->E = ion_foil.E;
    p->S = ion_foil.S;
}

double scattering_angle(const ion *incident, sim_workspace *ws) { /* Calculate scattering angle necessary for ion (in sample coordinate system) to hit detector */

    double scatter_theta, scatter_phi;

    rotate(ws->det->theta, ws->det->phi, incident->theta, incident->phi, &scatter_theta, &scatter_phi);/* Detector in ion system */
#ifdef DEBUG
    fprintf(stderr, "Detector in ion system angles %g deg and %g deg.\n", scatter_theta/C_DEG, scatter_phi/C_DEG);
#endif
    rotate(scatter_theta, scatter_phi, -1.0*ws->sim->sample_theta, ws->sim->sample_phi, &scatter_theta, &scatter_phi); /* Counter sample rotation. Detector in lab (usually). If ion was somehow "deflected" then this is the real scattering angle. Compare to sim->theta.  */
#ifdef DEBUG
    fprintf(stderr, "Incident angles %g deg and %g deg (in sample)\n", incident->theta/C_DEG, incident->phi/C_DEG);
    fprintf(stderr, "Sample lab angles %g deg and %g deg, detector lab angles %g deg and %g deg\n", ws->sim->sample_theta/C_DEG, ws->sim->sample_phi/C_DEG, ws->det->theta/C_DEG, ws->det->phi/C_DEG);
    fprintf(stderr, "Final lab scattering angles %g deg and %g deg\n", scatter_theta/C_DEG, scatter_phi/C_DEG);
    if(!ws->params.ds) {
        assert(fabs(ws->det->theta - scatter_theta) < 0.01 * C_DEG); /* with DS this assert will fail */
    }
#endif
    return scatter_theta;
}

void set_ion_exit_angles(const simulation *sim, const detector *det, ion *p) {
    double theta, phi; /* Generic polar and azimuth angles */
    rotate(det->theta, det->phi, sim->sample_theta, sim->sample_phi, &theta, &phi); /* Detector in sample coordinate system */
    ion_set_angle(p, theta, phi);
#ifdef DEBUG
    fprintf(stderr, "Reaction angles (in sample) %g deg and %g deg\n", theta/C_DEG, phi/C_DEG);
#endif
}

void simulate(const ion *incident, const depth depth_start, sim_workspace *ws, const sample *sample) { /* Ion is expected to be in the sample coordinate system at starting depth. Also note that sample may be slightly different (e.g. due to roughness) to ws->sim->sample */
    assert(sample->n_ranges);
    int warnings = 0;
    double thickness = sample->ranges[sample->n_ranges-1].x;
    size_t i_depth;
    ion ion1 = *incident; /* Shallow copy of the incident ion */
    double scatter_theta = scattering_angle(incident, ws);
#ifdef DEBUG
    ion_print(stderr, incident);
    fprintf(stderr, "Simulate from depth %g tfu (index %zu), detector theta = %g deg, calculated theta = %g deg. %zu reactions.\n", depth_start.x/C_TFU, depth_start.i, ws->det->theta/C_DEG, scatter_theta/C_DEG, ws->n_reactions);
#endif
    depth d_before = depth_start;
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = &ws->reactions[i];
        r->last_brick = 0;
        r->stop = FALSE;
        r->theta = scatter_theta;
        sim_reaction_recalculate_internal_variables(r);
        if(r->stop) {
            r->max_depth = 0.0;
            continue;
        }
        ion *p = &r->p;
        p->E = ion1.E * r->K;
        p->S = ion1.S * r->K;
        r->max_depth = sample_isotope_max_depth(sample, r->i_isotope);
        set_ion_exit_angles(ws->sim, ws->det, p); /* Calculating this for every reaction is not necessary (the angles are the same), but what do we know... */
        sim_reaction_reset_bricks(r);
        brick *b = &r->bricks[0];
        b->d = d_before;
        b->E_0 = ion1.E;
        b->Q = 0.0;
        if(r->i_isotope >= sample->n_isotopes) { /* No target isotope for reaction. */
            r->stop = TRUE;
        }
        post_scatter_exit(&r->p, b->d, ws, sample); /* Calculates the exit energy if calculation doesn't start from the surface */
        b->E = p->E;
        b->S = p->S;
#ifdef DEBUG
        fprintf(stderr, "Simulation reaction %zu: %s %s. Max depth %g tfu. i_isotope=%zu, stop = %i.\nCross section at %g keV is %g mb/sr, exit %g keV\n",
                i, reaction_name(r->r), r->r->target->name, r->max_depth / C_TFU, r->i_isotope, r->stop, ion1.E/C_KEV, r->cross_section(r, ion1.E)/C_MB_SR, p->E/C_KEV);
#endif
    }
    i_depth=1;
    if(fabs(ion1.cosine_theta) < 1e-6) {
#ifdef DEBUG
        fprintf(stderr, "Ion was going sideways in the sample, we nudged it a little.\n");
#endif
        ion1.theta = 90.01 * C_DEG;
        ion1.cosine_theta = cos(ion1.theta);
        return;
    }
    while(1) {
        if(warnings > SIMULATE_WARNING_LIMIT) {
            fprintf(stderr, "Warning limit reached. Won't calculate anything.\n");
            break;
        }
        if(d_before.x >= thickness) /* We're in too deep. */
            break;
        if(ion1.inverse_cosine_theta < 0.0 && d_before.x < 0.001*C_TFU) /* We're coming out of the sample and about to exit */
            break;
        if (ion1.E < ws->params.emin) {
#ifdef DEBUG
            fprintf(stderr, "Break due to low energy (%.3lf keV < %.3lf keV), x = %.3lf, i_range = %lu.\n", ion1.E/C_KEV, ws->params.emin/C_KEV, d_before.x/C_TFU, d_before.i);
#endif
            break;
        }
        const double E_front = ion1.E;
        const double S_front = ion1.S;
        double E_step = ws->params.stop_step_incident == 0.0?STOP_STEP_AUTO_FUDGE_FACTOR*sqrt(ws->det->resolution+ion1.S):ws->params.stop_step_incident;
        depth d_after = stop_step(ws, &ion1, sample, d_before, E_step);
#ifdef DEBUG
        fprintf(stderr, "After:  %g tfu in range %zu\n", d_after.x/C_TFU, d_after.i);
#endif
        const double d_diff = depth_diff(d_before, d_after);
        /* DEPTH BIN [x, x+d_diff) */
        const double E_back = ion1.E;
        const double S_back = ion1.S;
        if(fabs(d_diff) < 0.001*C_TFU && E_front-E_back < 0.001*C_KEV) {
            fprintf(stderr, "Warning: no or very little progress was made (E step (goal) %g keV, E from %g keV to E = %g keV, depth = %g tfu, d_diff = %g tfu), check stopping or step size.\n", E_step/C_KEV, E_front/C_KEV , E_back/C_KEV, d_before.x/C_TFU, d_diff/C_TFU);
            //sample_print(stderr, sample, FALSE);
            d_before.x += incident->inverse_cosine_theta*0.0001*C_TFU;
            warnings++;
            continue;
        }

#ifdef DEBUG_VERBOSE
        double E_diff = E_front-E_back;
        fprintf(stderr, "x = %8.3lf, x+h = %6g, E = %8.3lf keV to  %8.3lf keV (diff %6.4lf keV)\n", x/C_TFU, (x+h)/C_TFU, E_front/C_KEV, ws->ion.E/C_KEV, E_diff/C_KEV);
#endif
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "For incident beam: E_front = %g MeV, E_back = %g MeV,  E_mean = %g MeV, sqrt(S) = %g keV\n",
                        E_front / C_MEV, E_back / C_MEV, E_mean / C_MEV, sqrt(ion1.S) / C_KEV);
#endif
        int all_stop = 1; /* Do we have any reactions left (or are we too deep in the sample) */
        for (size_t i = 0; i < ws->n_reactions; i++) {
            sim_reaction *r = &ws->reactions[i];
            if(r->stop)
                continue;
            all_stop = 0;
            if(i_depth >= r->n_bricks) {
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
                b->Q = 0.0;
                r->stop = TRUE;
                continue;
            }
#ifdef DEBUG_REACTION
            fprintf(stderr, "Reaction %s (%zu): %s\n", reaction_name(r->r), i, r->r->target->name);
#endif
            post_scatter_exit(&r->p, d_after, ws, sample);
            b->E = r->p.E; /* Now exited from sample */
            b->S = r->p.S;
            if (r->p.E < ws->params.emin) {
                r->stop = TRUE;
                b->Q = 0.0;
                continue;
            }
            double sigma_conc = cross_section_concentration_product(ws, sample, r, E_front, E_back, &d_before, &d_after, S_front, S_back); /* Product of concentration and sigma for isotope i_isotope target and this reaction. */
            if(sigma_conc > 0.0) {
                if(d_after.i == sample->n_ranges - 2) {
                    sigma_conc *= ws->sim->channeling_offset + ws->sim->channeling_slope * (E_front + E_back)/2.0;
                }
                b->Q = ion1.inverse_cosine_theta * sigma_conc * d_diff;
                assert(b->Q >= 0.0);
#ifdef DEBUG
                fprintf(stderr, "    %s: type=%i, E_front = %.3lf, E_after = %.3lf, E_out = %.3lf (sigma*conc = %g mb/sr, Q = %g (thickness = %.4lf tfu)\n",
                                 r->r->target->name, r->r->type, E_front/C_KEV, (ion1.E * r->K)/C_KEV, r->p.E/C_KEV, sigma_conc/C_MB_SR, b->Q, d_diff/C_TFU);
#endif
            } else {
                b->Q = 0.0;
            }
            r->last_brick = i_depth;
        }
        d_before = d_after;
        ion1.S = S_back;
        ion1.E = E_back;
        i_depth++;

        if(all_stop) {
#ifdef DEBUG
            fprintf(stderr, "All reactions have ceased by depth %g.\n", d_before.x/C_TFU);
#endif
            break;
        }
    }
    convolute_bricks(ws);
}

int process_reaction_files(simulation *sim, const jibal_isotope *jibal_isotopes, char * const *reaction_filenames, size_t n_reaction_filenames) {
    for(size_t i = 0; i < n_reaction_filenames; i++) {
        r33_file *rfile = r33_file_read(reaction_filenames[i]);
        if(!rfile) {
            return -1;
        }
        reaction *reaction_from_file = r33_file_to_reaction(jibal_isotopes, rfile);
        if(!reaction_from_file) {
            r33_file_free(rfile);
            return -1;
        }
        fprintf(stderr, "File: %s has a reaction with %s -> %s, product %s, theta %g deg\n", reaction_filenames[i],
                reaction_from_file->incident->name, reaction_from_file->target->name, reaction_from_file->product->name, reaction_from_file->theta/C_DEG);
        for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
            reaction *r = &sim->reactions[i_reaction];
            if(reaction_is_same(r, reaction_from_file)) {
                fprintf(stderr, "Replacing reaction %zu (%s with %s).\n", i_reaction, reaction_name(r), r->target->name);
                reaction_from_file->cs = r->cs; /* Adopt fallback cross-section from the reaction we are replacing */
                reaction_free(r);
                *r = *reaction_from_file;
                break;
            }
        }
    }
    return 0;
}

int assign_stopping(jibal_gsto *gsto, const simulation *sim) {
    struct sample *sample = sim->sample;
    if(!sample) {
        fprintf(stderr, "Could not assign stopping, because sample is not set!\n");
        return 1;
    }
    int Z2 = -1;
    for(size_t i = 0; i < sample->n_isotopes; i++) {
        if(Z2 == sample->isotopes[i]->Z) /* Z2 repeats, skip */
            continue;
        Z2 = sample->isotopes[i]->Z;
        if (!jibal_gsto_auto_assign(gsto, sim->beam_isotope->Z, Z2)) { /* This should handle RBS */
            fprintf(stderr, "Can not assign stopping.\n");
            return EXIT_FAILURE;
        }
        for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
            const reaction *r = &sim->reactions[i_reaction];
            if (!jibal_gsto_auto_assign(gsto, r->product->Z, Z2)) {
                fprintf(stderr, "Can not assign stopping.\n");
                return EXIT_SUCCESS;
            }
        }
    }
    return 0;
}

int print_spectra(const char *filename, const sim_workspace *ws, const gsl_histogram *exp) {
    char sep = ' ';
    if(!ws)
        return EXIT_FAILURE;
    FILE *f;
    if(!filename) {
        f = stdout; /* Another dirty hack, fopen_file_or_stream() should handle this */
    } else {
        f = fopen_file_or_stream(filename, "w");
    }
    if(!f)
        return EXIT_FAILURE;
    if(filename) {
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
            fprintf(f,"%c%.3lf", sep, exp->range[i]/C_KEV);
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
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

void fit_params_add(simulation *sim, const sample_model *sm, fit_params *params, const char *fit_vars) {
#ifdef DEBUG
    fprintf(stderr, "fitvars = %s\n", fit_vars);
#endif
    if(!fit_vars)
        return;
    char *token, *s, *s_orig;
    s_orig = s = strdup(fit_vars);
    assert(s != NULL);
    while ((token = strsep(&s, ",")) != NULL) { /* parse comma separated list of parameters to fit */
#ifdef DEBUG
        fprintf(stderr, "Thing to fit: \"%s\". Sim pointer is %p, sample model pointer is %p\n", token, (void *)sim, (void *)sm);
#endif
        if(sim) {
            for(size_t i_det = 0; i_det < sim->n_det; i_det++) { /* TODO: this adds parameters of all detectors, make this controllable for an individual detector */
                detector *det = sim->det[i_det];
                assert(det);
                if(strncmp(token, "calib", 5) == 0) {
                    fit_params_add_parameter(params, &det->slope); /* TODO: prevent adding already added things */
                    fit_params_add_parameter(params, &det->offset);
                    fit_params_add_parameter(params, &det->resolution);
                }
                if(strcmp(token, "slope") == 0) {
                    fit_params_add_parameter(params, &det->slope);
                }
                if(strcmp(token, "offset") == 0) {
                    fit_params_add_parameter(params, &det->offset);
                }
                if(strncmp(token, "reso", 4) == 0) {
                    fit_params_add_parameter(params, &det->resolution);
                }
                if(strncmp(token, "solid", 5) == 0) {
                    fit_params_add_parameter(params, &det->solid);
                }
            }
            if(strcmp(token, "fluence") == 0) {
                fit_params_add_parameter(params, &sim->fluence);
            }
            if(strcmp(token, "channeling_slope") == 0) {
                fit_params_add_parameter(params, &sim->channeling_slope);
            }
            if(strcmp(token, "channeling") == 0) {
                fit_params_add_parameter(params, &sim->channeling_offset);
            }

        }
        if(sm) {
            if(strncmp(token, "rough", 5) == 0 && strlen(token) > 5) {
                size_t i_layer = strtoul(token + 5, NULL, 10);
                if(i_layer >= 1 && i_layer <= sm->n_ranges) {
                    fit_params_add_parameter(params, &sm->ranges[i_layer - 1].rough.x);
                } else {
                    fprintf(stderr, "No layer %zu (parsed from \"%s\")\n", i_layer, token);
                }
            }
            if(strncmp(token, "thickness", 9) == 0 && strlen(token) > 9) {
                size_t i_layer = strtoul(token + 9, NULL, 10);
                if(i_layer >= 1 && i_layer <= sm->n_ranges) {
                    fit_params_add_parameter(params, &sm->ranges[i_layer - 1].x);
                } else {
                    fprintf(stderr, "No layer %zu (parsed from \"%s\")\n", i_layer, token);
                }
            }
            size_t i, j;
            if(sscanf(token, "conc%lu_%lu", &i, &j) == 2) {
                if(i >= 1 && i <= sm->n_ranges && j >= 1 && j <= sm->n_materials) {
                    fit_params_add_parameter(params, sample_model_conc_bin(sm, i - 1, j - 1));
                } else {
                    fprintf(stderr, "No element %lu in layer %lu\n", j, i);
                }
            }
        }
    }
    free(s_orig);
}


void print_bricks(const char *filename, const sim_workspace *ws) {
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return;
    for(size_t i = 0; i < ws->n_reactions; i++) {
        const sim_reaction *r = &ws->reactions[i];
        fprintf(f, "#%s %s\n", reaction_name(r->r), r->r->target->name);
        for(size_t j = 0; j <= r->last_brick; j++) {
            brick *b = &r->bricks[j];
            fprintf(f, "%2lu %2lu %8.3lf %8.3lf %8.3lf %8.3lf %12.3lf\n",
                    i, j, b->d.x/C_TFU, b->E_0/C_KEV, b->E/C_KEV, sqrt(b->S)/C_KEV, b->Q * ws->sim->fluence);
        }
        fprintf(f, "\n\n");
    }
    if(!(f == stdout || f == stderr))
        fclose(f);
}



void simulate_with_roughness(sim_workspace *ws) {
    double p_sr = ws->sim->fluence;
    size_t n_rl = 0; /* Number of rough layers */
    for(size_t i = 0; i < ws->sample->n_ranges; i++) {
        if(ws->sample->ranges[i].rough.model == ROUGHNESS_GAMMA)
            n_rl++;
    }
#ifdef DEBUG
    fprintf(stderr, "%zu rough layers\n", n_rl);
#endif
    if(!n_rl) {
        ion_set_angle(&ws->ion, 0.0*C_DEG, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        simulate(&ws->ion, depth_seek(ws->sample, 0.0*C_TFU), ws, ws->sample);
        return;
    }
    struct sample *sample_rough = sample_copy(ws->sample);
    size_t *index = malloc(sizeof(size_t) * n_rl);
    size_t *modulos = malloc(sizeof(size_t) * n_rl);
    size_t j = 0;
    thick_prob_dist **tpd = malloc(sizeof(thick_prob_dist *) * n_rl);
    for(size_t i = 0; i < ws->sample->n_ranges; i++) {
        if(ws->sample->ranges[i].rough.model == ROUGHNESS_GAMMA) {
#ifdef DEBUG
            fprintf(stderr, "Range %zu is rough (gamma), amount %g tfu, n = %zu spectra\n", i, ws->sample->ranges[i].rough.x/C_TFU, ws->sample->ranges[i].rough.n);
#endif
            assert(ws->sample->ranges[i].rough.n > 0 && ws->sample->ranges[i].rough.n < 1000);
            tpd[j] = thickness_probability_table_gen(ws->sample->ranges[i].x, ws->sample->ranges[i].rough.x, ws->sample->ranges[i].rough.n);
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
        for(size_t i_range = 0; i_range < ws->sample->n_ranges; i_range++) { /* Reset ranges for every iter */
            sample_rough->ranges[i_range].x = ws->sample->ranges[i_range].x;
        }
        for(size_t i = 0; i < n_rl; i++) {
            j = (i_iter / modulos[i]) % tpd[i]->n; /* "j"th roughness element */
            //fprintf(stderr, " %zu", j);

            size_t i_range = index[i];
            p *= tpd[i]->p[j].prob; /* Probability is multiplied by the "i"th roughness, element "j" */
            double x_diff = tpd[i]->p[j].x - ws->sample->ranges[i_range].x; /* Amount to change thickness of this and and all subsequent layers */
            for(; i_range < ws->sample->n_ranges; i_range++) {
                sample_rough->ranges[i_range].x += x_diff;
            }
#ifdef DEBUG
            fprintf(stderr, "Gamma roughness diff %g tfu (from %g tfu, index i_range=%zu), probability %.3lf%%)\n", x_diff/C_TFU, ws->sample->ranges[i_range].x/C_TFU, i_range, tpd[i]->p[j].prob*100.0);
            fprintf(stderr, "Gamma roughness, ranges");
            for(i_range = 0; i_range < ws->sample->n_ranges; i_range++) {
                fprintf(stderr, ", %zu: %g tfu ", i_range, sample_rough->ranges[i_range].x/C_TFU);
            }
            fprintf(stderr, "\n");
#endif
        }
        //fprintf(stderr, "\n");
        ws->fluence = p * p_sr;
        ion_set_angle(&ws->ion, 0.0, 0.0);
        ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
        simulate(&ws->ion, depth_seek(ws->sample, 0.0), ws, sample_rough);
    }
    for(size_t i = 0; i < n_rl; i++) {
        thickness_probability_table_free(tpd[i]);
    }
    sample_free(sample_rough);
    free(modulos);
    free(index);
    free(tpd);
}

void simulate_with_ds(sim_workspace *ws) {
    if(!ws) {
        fprintf(stderr, "No workspace, no simulation.\n");
        return;
    }
    double fluence = ws->fluence;
    simulate_with_roughness(ws);
    if(!ws->params.ds) {
        sim_workspace_calculate_sum_spectra(ws);
        return;
    }
    ion_set_angle(&ws->ion, 0.0, 0.0);
    ion_rotate(&ws->ion, ws->sim->sample_theta, ws->sim->sample_phi);
    ion ion1 = ws->ion;
    ion ion2 = ion1;
    depth d_before = depth_seek(ws->sample, 0.0);
    ws->params.rk4 = FALSE; /* This change is not reversed, nor reflected back to sim->params */
    ws->params.nucl_stop_accurate = FALSE;
    ws->params.mean_conc_and_energy = TRUE;
    fprintf(stderr, "\n");
    const jibal_isotope *incident = ws->sim->beam_isotope;
    while(1) {
        double E_front = ion1.E;
        if(E_front < ws->params.emin)
            break;
        depth d_after = stop_step(ws, &ion1, ws->sample, d_before, sqrt(ws->det->resolution+ion1.S)); /* TODO: step? */
        double thick_step = depth_diff(d_before, d_after);
        const depth d_halfdepth = {.x = (d_before.x + d_after.x)/2.0, .i = d_after.i}; /* Stop step performs all calculations in a single range (the one in output!). That is why d_after.i instead of d_before.i */
        double E_back = ion1.E;
        const double E_mean = (E_front + E_back) / 2.0;

        fprintf(stderr, "DS depth from %9.3lf tfu to %9.3lf tfu, E from %6.1lf keV to %6.1lf keV. p*sr = %g\n", d_before.x/C_TFU, d_after.x/C_TFU, E_front/C_KEV, E_back/C_KEV, fluence);
        double p_sum = 0.0;
        for(int i_polar = 0; i_polar < ws->params.ds_steps_polar; i_polar++) {
            const double ds_polar_min = 20.0*C_DEG;
            const double ds_polar_max = 180.0*C_DEG;
            double ds_polar_step = (ds_polar_max-ds_polar_min)/(ws->params.ds_steps_polar-1);
            double ds_polar = ds_polar_min + i_polar * ds_polar_step;
            double cs_sum = 0.0;
            for(size_t i = 0; i < ws->n_reactions; i++) {
                sim_reaction *r = &ws->reactions[i];
                double c = get_conc(ws->sample, d_halfdepth, r->i_isotope);
                if(c < ABUNDANCE_THRESHOLD)
                    continue;
                const jibal_isotope *target = r->r->target;
                if(incident->mass >= target->mass && ds_polar > asin(target->mass / incident->mass)) { /* Scattering not possible */
                    continue;
                }
#if 0
                double cs = 0.0;
                for(int polar_substep = 0; polar_substep < 9; polar_substep++) {
                    double ds_polar_sub = ds_polar_step*(1.0*(polar_substep-4)/9.0) + ds_polar;
                    cs += jibal_cross_section_rbs(incident, target, ds_polar_sub, E_mean, JIBAL_CS_RUTHERFORD) * sin(ds_polar_sub)/9.0;
                }
                cs_sum += c * cs;
#else
                /* TODO: integrals over cross sections? Maybe stepwise from ds_polar=PI. */
                double cs = jibal_cross_section_rbs(incident, target, ds_polar, E_mean, JIBAL_CS_ANDERSEN);
                cs *= sin(ds_polar); /* spherical coordinates and surface areas... */
                cs_sum += c * cs;
#endif
            }
            double fluence_tot = cs_sum * thick_step * (2.0 * C_PI) * ds_polar_step; /* TODO: check calculation after moving from p_sr to fluence!*/
            p_sum += fluence_tot;
            double fluence_azi = fluence_tot / (1.0 * (ws->params.ds_steps_azi));
            for(int i_azi = 0; i_azi < ws->params.ds_steps_azi; i_azi++) {
                ion2 = ion1;
                double ds_azi = 2.0 * M_PI * (1.0 * i_azi) / (ws->params.ds_steps_azi * 1.0);
                ion_rotate(&ion2, ds_polar, ds_azi); /* Dual scattering: first scattering to some angle (scattering angle: ds_polar). Note that this does not follow SimNRA conventions. */
                ws->fluence = fluence_azi * fluence;
#ifdef DEBUG
                if(d_before.x == 0.0) {
                    fprintf(stderr, "DS polar %.3lf, azi %.3lf, scatter %.3lf\n", ds_polar/C_DEG, ds_azi/C_DEG, scattering_angle(&ion2, ws)/C_DEG);
                }
#endif
                if(scattering_angle(&ion2, ws) > 30.0*C_DEG) {
                    simulate(&ion2, d_after, ws, ws->sample);
                }
            }
        }
        fluence -= p_sum * fluence;
        if(ws->sample->ranges[ws->sample->n_ranges-1].x - d_after.x < 0.01*C_TFU)
            break;
        d_before = d_after;
    }
    sim_workspace_calculate_sum_spectra(ws);
}
