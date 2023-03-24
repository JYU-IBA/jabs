#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <jibal.h>
#include <jibal_units.h>
#include <jibal_phys.h>
#include <jibal_kin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include "jabs_debug.h"
#include "scatint.h"

double scatint_sigma(scatint_params *p) {
    if(s_seek(p, p->theta_cm)) {
        DEBUGSTR("Did not work.");
        return 0.0;
    }
    if(cs(p)) {
        DEBUGSTR("Could not calculate cross section.");
        return 0.0;
    }
    return p->sigma_lab * p->ik_scaling;
}

double scatint_sigma_lab(scatint_params *p, double E_lab, double theta_lab) { /* This can be called after scatint_init(), returns lab cross section or zero on failure */
    if(!p) {
        return 0.0;
    }
    scatint_set_theta(p, theta_lab);
    scatint_set_energy(p, E_lab);
    return scatint_sigma(p);
}

int main_scatint(int argc, char **argv) {
    if(argc < 5) {
        fprintf(stderr, "Wrong number of arguments! Usage: scatint <incident> <target isotope> <lab angle> <energy>\n");
        return EXIT_FAILURE;
    }
    jibal *jibal = jibal_init(NULL);
    if(jibal->error) {
        fprintf(stderr, "Initializing JIBAL failed with error code: %i (%s)\n", jibal->error,
                jibal_error_string(jibal->error));
        return EXIT_FAILURE;
    }

    const jibal_isotope *incident = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    if(!incident) {
        fprintf(stderr, "There is no isotope %s in my database.\n", argv[1]);
        return -1;
    }
    const jibal_isotope *target = jibal_isotope_find(jibal->isotopes, argv[2], 0, 0);
    if(!target) {
        fprintf(stderr, "There is no isotope %s in my database.\n", argv[2]);
        return -1;
    }
    fprintf(stderr, "Z1 = %i\nZ2 = %i\n", incident->Z, target->Z);

    double theta = jibal_get_val(jibal->units, JIBAL_UNIT_TYPE_ANGLE, argv[3]);
    fprintf(stderr, "theta = %g deg\n", theta/C_DEG);
    if(theta < 0.0 || theta > C_PI) {
        fprintf(stderr, "Scattering angle must be between 0 and 180 degrees. You gave %g deg.\n", theta/C_DEG);
        return EXIT_FAILURE;
    }

    double E = jibal_get_val(jibal->units, JIBAL_UNIT_TYPE_ENERGY, argv[4]);
    fprintf(stderr, "E_lab = %g keV\n", E/C_KEV);
    if(E < 1.0*C_EV) {
        fprintf(stderr, "Energy must be above zero.\n");
        return EXIT_FAILURE;
    }
    if(E > 1.0e9 * C_EV) {
        fprintf(stderr, "Are you sure about your energy units? You gave %g MeV.\n", E/C_MEV);
        return EXIT_FAILURE;
    }

    potential_type pt = POTENTIAL_UNIVERSAL; /* Default */

    if(argc > 5) {
        if(strcmp(argv[5], "Andersen") == 0) {
            pt = POTENTIAL_ANDERSEN;
        } else if(strcmp(argv[5], "Rutherford") == 0) {
            pt = POTENTIAL_RUTHERFORD;
        } else if(strcmp(argv[5], "Bohr") == 0) {
            pt = POTENTIAL_BOHR;
        }
    }

    size_t n_params = 3;
    scatint_params **params = calloc(n_params, sizeof(scatint_params *));
    params[0] = scatint_init(REACTION_RBS, pt, incident, target);
    params[1] = scatint_init(REACTION_RBS_ALT, pt, incident, target);
    params[2] = scatint_init(REACTION_ERD, pt, incident, target);

    double theta_max = asin(target->mass/incident->mass);
    if(incident->mass >= target->mass) { /* Print maximum scattering angle for convenience. */
        fprintf(stderr, "theta_max = %g deg\n", theta_max/C_DEG);
    }

    for(size_t i_p = 0; i_p < n_params; i_p++) {
        scatint_params *p = params[i_p];
        fprintf(stderr, "\nReaction: %s\n",scatint_reaction_name(p->rt));
        scatint_set_energy(p, E);
        scatint_set_theta(p, theta);
        if(!scatint_has_solution(p)) {
            fprintf(stderr, "No solution.\n");
            continue;
        }
        print_cs(p, TRUE);
    }

    for(size_t i_p = 0; i_p < n_params; i_p++) {
        scatint_params_free(params[i_p]);
    }
    free(params);
    return EXIT_SUCCESS;
}

double potential_universal(double x) {
    return 0.1818*exp(-3.2*(x)) +  0.5099*exp(-0.9423*(x)) + 0.2802*exp(-0.4029*(x)) + 0.02817*exp(-0.2016*(x));
}

double potential_andersen(double x) {
    if(x < 1.0/1.586) {
        return 1.0 - 1.586 * x;
    } else {
        return 0.0;
    }
}

double potential_rutherford(double x) { /* Good old Coulomb */
    (void) x;
    return 1.0;
}

double potential_bohr(double x) {
    return exp(-(x));
}

double apsis_f (double x, void *p) {
    scatint_params * params = (scatint_params *)p;
    double s = (params->s);
    double E_rel = (params->E_rel);
    double out = 1.0 - (s*s)/(x*x) - params->potential(x) / x / E_rel;
//    fprintf(stderr, "apsis_f(s = %g, E_rel = %g, r = %g) = %g\n", s, E_rel, x, out); 
    return out;
}


int apsis(scatint_params *params) {
    int status;
    int iter = 0, max_iter = 100;
    gsl_function F;
    F.function = &apsis_f;
    F.params = params;
    double x_lo = IMPACT_MIN, x_hi = IMPACT_MAX;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *solver;
    T = gsl_root_fsolver_brent;
    solver = gsl_root_fsolver_alloc (T);
    double r;
    gsl_root_fsolver_set (solver, &F, x_lo, x_hi);
    double accuracy = 1e-9;
    do
    {
      iter++;
      status = gsl_root_fsolver_iterate(solver);
      r = gsl_root_fsolver_root(solver);
      x_lo = gsl_root_fsolver_x_lower(solver);
      x_hi = gsl_root_fsolver_x_upper(solver);
      status = gsl_root_test_interval (x_lo, x_hi, 0, accuracy);
#if 0
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
              iter, x_lo, x_hi,
              r, x_hi - x_lo);
#endif
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    if(iter == max_iter - 1) {
        fprintf(stderr, "Max iters reached!\n");
        params->R = -1.0;
        return EXIT_FAILURE;
    }
    gsl_root_fsolver_free(solver);
    params->R = r;
    return EXIT_SUCCESS;
}

double scat_theta_f(double x, void *params) {
    struct scatint_params *p = (struct scatint_params *) params;
    double f = sqrt(apsis_f(x, p));
    double out = 1.0/(x*x*f);
    return out;
}

int calc_scattering_angle(struct scatint_params *params) {
    double result, error;
    gsl_function F;
    F.function = &scat_theta_f;
    F.params = params;
    if(apsis(params)) {
        fprintf(stderr, "Could not calculate apsis.\n");
        return EXIT_FAILURE;
    }
    DEBUGVERBOSEMSG("Integration from %g to infinity, accuracy %g.", params->R, params->accuracy);
    gsl_set_error_handler_off();
    int status = gsl_integration_qagiu(&F, params->R, 0, params->accuracy, params->w->limit, params->w, &result, &error);
    DEBUGVERBOSEMSG("result          = % .18f", result);
    DEBUGVERBOSEMSG("estimated error = %.18f", error);
    DEBUGVERBOSEMSG("intervals       = %zu", params->w->size);
    params->theta_cm = C_PI - 2.0 * params->s * result;
    return status;
}

int cs(scatint_params *p) {
    DEBUGVERBOSEMSG("E_cm = %g keV", p->E_cm/C_KEV);
    DEBUGVERBOSEMSG("E_rel = %g", p->E_rel);
    DEBUGVERBOSEMSG("a = %g (%g angstrom)", p->a, p->a/ C_ANGSTROM);
    DEBUGVERBOSEMSG("c = %g (%g keV)", p->c, p->c/C_KEV);

    if(calc_scattering_angle(p)) {
        DEBUGSTR("Failed to calculate angle.");
        return EXIT_FAILURE;
    }
    scatint_params p2 = *p;
    p2.s = p->s * 1.00001 + 1e-9;
    if(calc_scattering_angle(&p2)) {
        DEBUGSTR("Failed to calculate angle (with difference).");
        return EXIT_FAILURE;
    }
    double cos_theta_cm_delta = cos(p2.theta_cm) - cos(p->theta_cm);
    double sigma_cm = 2.0 * C_PI * (pow2(p2.s)-pow2(p->s))/cos_theta_cm_delta; /* TODO: this is positive diff */
    sigma_cm *= p->a * p->a / (4.0 * C_PI);
    double theta_lab = atan2(sin(p->theta_cm), (cos(p->theta_cm) + p->m12));
    DEBUGVERBOSEMSG("Sigma conversion from C.M to lab is %g", p->sigma_lab_to_cm_ratio);
    double sigma_lab = fabs(p->sigma_lab_to_cm_ratio * sigma_cm);
    DEBUGVERBOSEMSG("theta_cm = %g (%.5lf deg)\n", p->theta_cm, p->theta_cm/C_DEG);
    DEBUGVERBOSEMSG("ratio %g", p->sigma_lab_to_cm_ratio);
    DEBUGVERBOSEMSG("Diff. scatt. cross-section = %g", sigma_cm);
    DEBUGVERBOSEMSG("In other units = %g mb/sr (C.M.)", sigma_cm / C_MB_SR);
    DEBUGVERBOSEMSG("Cross section (lab) = %g mb/sr", sigma_lab / C_MB_SR);
    DEBUGVERBOSEMSG("theta_lab = %g (%g deg)", theta_lab, 180.0*theta_lab/C_PI);
    p->theta_lab = theta_lab;
    p->sigma_lab = sigma_lab;
    p->sigma_cm = sigma_cm;
    return EXIT_SUCCESS;
}


scatint_params *scatint_init(reaction_type rt, potential_type pt, const jibal_isotope *incident, const jibal_isotope *target) {
    if(!incident || !target) {
        return NULL;
    }
    scatint_params *p = calloc(1, sizeof(scatint_params));
    p->rt = rt;
    p->pt = pt;
    if(rt == REACTION_ERD) { /* ERD will be done with inverse kinematics */
        p->incident = target;
        p->target = incident;
    } else if(rt == REACTION_RBS || rt == REACTION_RBS_ALT){
        p->incident = incident;
        p->target = target;
    } else {
        free(p);
        return NULL;
    }
    switch(p->pt) {
        case POTENTIAL_BOHR:
            p->a = screening_length_bohr(incident->Z, target->Z);
            p->potential = potential_bohr;
            break;
        case POTENTIAL_UNIVERSAL:
            p->a = screening_length_universal(incident->Z, target->Z);
            p->potential = potential_universal;
            break;
        case POTENTIAL_ANDERSEN:
            p->a = screening_length_andersen(incident->Z, target->Z);
            p->potential = potential_andersen;
            break;
        case POTENTIAL_RUTHERFORD:
            p->a = C_BOHR_RADIUS; /* This can be an arbitrary number (no screening length) but it is best to choose a number comparable to the other models */
            p->potential = potential_rutherford;
            break;
        case POTENTIAL_NONE:
        default:
            p->a = 0.0;
            p->potential = NULL;
            break;
    }
    if(!p->potential) {
        free(p);
        return NULL;
    }
    p->s = 0.0; /* Will be changed later */
    p->E_lab = 0.0;
    p->E_rel = 0.0;
    p->theta_cm = 0.0;
    p->theta_lab = 0.0;
    p->c = (1.0 / (4.0 * C_PI * C_EPSILON0)) * incident->Z * target->Z * C_E * C_E / p->a;
    p->Z1 = p->incident->Z;
    p->Z2 = p->target->Z;
    p->m1 = p->incident->mass;
    p->m2 = p->target->mass;
    p->m12 = p->incident->mass / p->target->mass;
    p->m21 = p->target->mass / p->incident->mass;
    p->theta_max = asin(target->mass/incident->mass); /* only valid when incident->mass >= target->mass */
    p->ik_scaling = 1.0; /* Will be changed for ERD when theta is changed */
    p->E_ik_ratio = target->mass / incident->mass; /* Remember, target and incident are already inverted */
    p->w = gsl_integration_workspace_alloc(INTEGRATION_WORKSPACE_N);
    p->accuracy = INTEGRATION_ACCURACY;
    return p;
}

void scatint_params_free(scatint_params *p) {
    if(!p) {
        return;
    }
    gsl_integration_workspace_free(p->w);
    free(p);
}

int scatint_set_energy(scatint_params *p, double E_lab) {
    if(p->rt == REACTION_ERD) { /* ERD: calculate inverse kinematics */
        p->E_lab = p->E_ik_ratio * E_lab;
    } else {
        p->E_lab = E_lab;
    }
    p->E_cm = p->E_lab * p->m21 / (1.0 + p->m21);
    p->E_rel = p->E_cm / p->c;
    return EXIT_SUCCESS;
}

int scatint_set_energy_cm(scatint_params *p, double E_cm) { /* alternative to scatint_set_energy() */
    if(p->rt == REACTION_ERD) {
        p->E_cm = p->E_ik_ratio * E_cm;
    } else {
        p->E_cm = E_cm;
    }
    p->E_lab = p->E_cm * (1.0 + p->m21) / p->m21;
    p->E_rel = p->E_cm / p->c;
    return EXIT_SUCCESS;
}

int scatint_set_theta(scatint_params *p, double theta_lab) {
    switch(p->rt) {
        case REACTION_RBS:
            p->theta_cm = theta_lab + asin(p->m12 * sin(theta_lab));
            p->theta_lab = theta_lab;
            p->ik_scaling = 1.0;
        break;
        case REACTION_RBS_ALT:
            p->theta_cm = C_PI - (asin(p->m12 * sin(theta_lab)) - theta_lab);
            p->theta_lab = theta_lab;
            p->ik_scaling = 1.0;
            break;
        case REACTION_ERD:
            p->theta_cm = C_PI - 2 * theta_lab; /* theta angle is now the lab angle of the recoil (usually called phi) */
            p->theta_lab = atan2(sin(p->theta_cm), (cos(p->theta_cm) + p->m12)); /* "RBS" angle, i.e. lab angle of the other guy (not the recoil). Note that m12 is for inverse kinematics. */
            p->ik_scaling = inverse_kinematics_erd_scale(p->theta_cm, p->theta_lab, theta_lab);
            break;
        default:
            return EXIT_FAILURE;
    }
    p->sigma_lab_to_cm_ratio = sigma_lab_to_cm(p->theta_cm, p->m12);
#ifdef DEBUG
    fprintf(stderr, "Angles calculated based on theta_lab = %g deg. Reaction type %i. theta_cm = %g deg, theta_lab = %g deg, ik_scaling = %g.\n", theta_lab/C_DEG, p->rt, p->theta_cm/C_DEG, p->theta_lab/C_DEG, p->ik_scaling);
#endif
    return EXIT_SUCCESS;
}

int s_seek(scatint_params *p, double theta_cm) {
    double s_low = 0.0;
    double s_high = IMPACT_MAX;
    const size_t n_iter_max = 100;
    double s = sqrt(IMPACT_MIN*IMPACT_MAX);
    double accuracy = IMPACT_FACTOR_ACCURACY; /* of theta_cm */
    if(!p->potential) {
        fprintf(stderr, "No potential has been set.\n");
        return EXIT_FAILURE;
    }
    for(size_t i = 0; i < n_iter_max; i++) { /* Find s for a given theta_cm. Small s => large theta_cm */
        p->s = s;
        if(calc_scattering_angle(p)) {
            DEBUGMSG("Fail on iter %zu. Could not calculate scattering angle with s = %g (%g m).\n", i+1, p->s, scatint_get_impact_parameter(p));
            return EXIT_FAILURE;
        }
        double result = p->theta_cm;
        DEBUGVERBOSEMSG("Iter %zu [%g, %g], tried %g got %g, diff %g", i+1, s_low, s_high, s, result, result - theta_cm);
        if(fabs(theta_cm - result) / theta_cm < accuracy) {
            return EXIT_SUCCESS;
        }
        if(result < theta_cm) {
            s_high = s;
        } else {
            s_low = s;
        }
        s = (s_low + s_high)/2.0;
    }
    p->s = -1.0; /* Invalidate */
    return EXIT_FAILURE;
}

double scatint_get_impact_parameter(const scatint_params *p) {
    return p->s * p->a;
}

int print_cs(scatint_params *p, int verbose) {
    const char *pot_name = scatint_potential_name(p->pt);
    if(s_seek(p, p->theta_cm)) {
        DEBUGSTR("Did not work.");
        return EXIT_FAILURE;
    }
    if(cs(p)) {
        DEBUGSTR("Could not calculate cross section.");
        return EXIT_FAILURE;
    }
    if(verbose) {
        if(p->rt == REACTION_ERD) {
            fprintf(stderr,
                    "This ERD cross section is calculated by inverse kinematics.\nERD cross section is %g times the RBS cross section for %g MeV %s scattering from %s to an angle of %g deg\n",
                    p->ik_scaling,
                    p->E_lab / C_MEV,
                    p->incident->name, /* Target and incident are named correctly for inverse kinematics */
                    p->target->name,
                    p->theta_lab / C_DEG
            );
        }
        fprintf(stderr, "Impact parameter is %g, i.e. %g meters.\n", p->s, scatint_get_impact_parameter(p));
        double aps = p->R * p->a;
        const double bohr_scaled = C_BOHR_RADIUS/p->Z2;
        const double nuclear_radius = 1.25e-15 * pow(p->target->A, 1.0/3.0);
        fprintf(stderr, "Apsis (%s potential) is %g meters, that is %g times the Bohr radius / Z2 and %g times nuclear radius, assuming it is 1.25e-15 * A^(1/3) m.\n", pot_name, aps, aps/bohr_scaled, aps/nuclear_radius);
#if 0
        double a_Coulomb = p->Z1 * p->Z2 * pow2(C_E) / p->E_lab / (4.0*C_PI*C_EPSILON0);
        fprintf(stderr, "Apsis (pure Coulomb) is %g meters, i.e. %g times the Bohr radius / Z2.\n", a_Coulomb, a_Coulomb/bohr_scaled);
#endif
    }
    double sigma_rutherford_cm = scatint_sigma_rutherford_cm(p);
    sigma_rutherford_cm *= p->ik_scaling;
    double sigma_rutherford_lab = sigma_rutherford_cm * p->sigma_lab_to_cm_ratio;
    double F_andersen = jibal_andersen_correction(p->Z1, p->Z2, p->E_cm, p->theta_cm);
    p->sigma_cm *= p->ik_scaling;
    p->sigma_lab *= p->ik_scaling;
    if(verbose) {
        fprintf(stderr, "E_lab = %g keV\n", p->E_lab / C_KEV);
        fprintf(stderr, "theta_cm = %g deg\n", p->theta_cm / C_DEG);
        fprintf(stderr, "E_cm = %g keV\n", p->E_cm / C_KEV);
        fprintf(stderr, "theta_lab = %lf deg (sanity check)\n", p->theta_lab / C_DEG);
        fprintf(stderr, "Sigma (Rutherford, analytical) = %g mb/sr (C.M.)\n", sigma_rutherford_cm / C_MB_SR);
        fprintf(stderr, "Sigma (Rutherford, analytical) = %g mb/sr (lab)\n", sigma_rutherford_lab / C_MB_SR);
        fprintf(stderr, "Sigma (Andersen, analytical) = %g mb/sr (C.M.)\n", F_andersen * sigma_rutherford_cm / C_MB_SR);
        fprintf(stderr, "Sigma (Andersen, analytical) = %g mb/sr (lab)\n", F_andersen * sigma_rutherford_lab / C_MB_SR);
        fprintf(stderr, "Sigma (%s) = %g mb/sr (C.M.)\n", pot_name, p->sigma_cm / C_MB_SR);
        fprintf(stderr, "Sigma (%s) = %g mb/sr (lab)\n", pot_name, p->sigma_lab / C_MB_SR);
        fprintf(stderr, "Screening factor (Andersen, analytical) = %lf\n", F_andersen);
        fprintf(stderr, "Screening factor (%s) = %lf\n", pot_name, p->sigma_cm / sigma_rutherford_cm);
    } else {
        printf(" %e %e %e", p->sigma_lab / C_MB_SR, F_andersen * sigma_rutherford_lab / C_MB_SR, sigma_rutherford_lab / C_MB_SR);
    }
    return EXIT_SUCCESS;
}

double scatint_sigma_rutherford_cm(const scatint_params *p) {
    return pow2((p->Z1*C_E*p->Z2*C_E)/(4.0*C_PI*C_EPSILON0))*pow2(1.0/(4.0*p->E_cm))*pow4(1.0/sin(p->theta_cm/2.0));
}

double screening_length_universal(int Z1, int Z2) {
    return 0.8853 * C_BOHR_RADIUS / (pow(Z1 * 1.0, 0.23) + pow(Z2 * 1.0, 0.23));
}
double screening_length_bohr(int Z1, int Z2) {
    return C_BOHR_RADIUS / sqrt(pow(Z1, 2.0/3.0) + pow(Z2, 2.0/3.0));
}

double screening_length_andersen(int Z1, int Z2) { /* Screening length used by Andersen */
    return 0.8853 * C_BOHR_RADIUS / sqrt(pow(Z1, 2.0/3.0) + pow(Z2, 2.0/3.0));
}

const char *scatint_potential_name(potential_type pt) {
    switch(pt) {
        case POTENTIAL_BOHR:
            return "Bohr";
        case POTENTIAL_UNIVERSAL:
            return "Universal";
        case POTENTIAL_ANDERSEN:
            return "Andersen";
        case POTENTIAL_RUTHERFORD:
            return "Rutherford";
        case POTENTIAL_NONE:
        default:
            return "None";
    }
}

const char *scatint_reaction_name(reaction_type rt) {
    switch(rt) {
        case REACTION_ERD:
            return "ERD";
        case REACTION_RBS:
            return "RBS";
        case REACTION_RBS_ALT:
            return "RBS alternative";
        case REACTION_NONE:
        default:
            return "None";
    }
}

int scatint_has_solution(scatint_params *p) {
    if(p->incident->mass >= p->target->mass && p->theta_lab > p->theta_max) {
        DEBUGMSG("No solution, %g u > %g u and %g deg > %g deg", p->incident->mass/C_U, p->target->mass/C_U, p->theta_lab/C_DEG, p->theta_max/C_DEG);
        return FALSE;
    }
    if(p->rt == REACTION_RBS_ALT) {
        if(p->target->mass < p->incident->mass) {
            return TRUE;
        } else {
            return FALSE;
        }
    }
    if(p->rt == REACTION_ERD && p->theta_cm < 0.0) { /* We can't rely on theta_lab due to inverse kinematics. phi > 90 deg (lab) results in negative theta_cm (quite unphysical..) so we use it instead */
        return FALSE;
    }
    return TRUE;
}

double inverse_kinematics_erd_scale(double theta_cm, double theta_lab, double phi_lab) { /* "RBS" angles theta_lab and theta_cm, phi is ERD in lab. */
    return fabs(4.0 * pow(sin(theta_lab), 2.0) * cos(theta_cm - theta_lab) * cos(phi_lab) / (pow(sin(theta_cm), 2.0)));
}
double sigma_lab_to_cm(double theta_cm, double m12) { /* m12 = m1/m2 */
    return fabs(pow((1.0 + pow2(m12) + 2.0 * m12 * cos(theta_cm)), 3.0 / 2.0) / (1.0 + m12 * cos(theta_cm)));
}
