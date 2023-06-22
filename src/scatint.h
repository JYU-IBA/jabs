#define IMPACT_MIN (1e-9) /* In units of screening length. Apsis can not be below this! */
#define IMPACT_MAX (100.0)
#define INTEGRATION_WORKSPACE_N (20000)
#define INTEGRATION_ACCURACY (1e-7)
#define IMPACT_FACTOR_ACCURACY (1e-9)

#include "reaction.h"

typedef enum potential_type {
    POTENTIAL_NONE = 0,
    POTENTIAL_BOHR = 1,
    POTENTIAL_UNIVERSAL = 2,
    POTENTIAL_ANDERSEN = 3,
    POTENTIAL_RUTHERFORD = 4,
    POTENTIAL_TF_SOMMERFELD = 5,
    POTENTIAL_TEST = 6
} potential_type;

typedef struct scatint_params {
    reaction_type rt;
    potential_type pt;
    const jibal_isotope *incident;
    const jibal_isotope *target;
    double s; /* Impact parameter, in units of a */
    double E_rel; /* Energy (C.M), in units of c */
    double a; /* Screening length */
    double c; /* Energy parameter */
    double R; /* Apsis */
    int Z1;
    int Z2;
    double m1;
    double m2;
    double m12; /* incident/target mass ratio */
    double m21; /* target/incident mass ratio */
    double E_cm;
    double E_lab;
    double sigma_cm;
    double sigma_lab;
    double theta_cm;
    double theta_lab;
    double theta_max;
    double ik_scaling;
    double E_ik_ratio;
    double sigma_lab_to_cm_ratio;
    double (*potential)(double);
    gsl_integration_workspace *w;
    double accuracy;
} scatint_params;

double scatint_sigma_lab(scatint_params *p, double E_lab, double theta_lab); /* This is the easiest one to use */
double scatint_sigma(scatint_params *p);
double potential_universal(double x);
double potential_andersen(double x);
double potential_rutherford(double x);
double screening_length_universal(int Z1, int Z2);
double screening_length_test(int Z1, int Z2, double em);
double screening_length_andersen(int Z1, int Z2);
double screening_length_bohr(int Z1, int Z2);
double screening_length_thomas_fermi(int Z2);
double potential_bohr(double x);
double apsis_f (double x, void *p);
int apsis(scatint_params *params);
double scat_theta_f(double x, void *params);
int calc_scattering_angle(struct scatint_params *params);
int cs(scatint_params *p);
scatint_params *scatint_init(reaction_type rt, potential_type pt, const jibal_isotope *incident, const jibal_isotope *target);
void scatint_params_free(scatint_params *p);
int scatint_set_energy(scatint_params *p, double E_lab);
int scatint_set_energy_cm(scatint_params *p, double E_cm);
int scatint_set_theta(scatint_params *p, double theta_lab);
double scatint_get_impact_parameter(const scatint_params *p); /* Returns impact parameter in SI units (m) */
int s_seek(scatint_params *p, double theta_cm);
int print_cs(scatint_params *p, int verbose);
double scatint_sigma_rutherford_cm(const scatint_params *p);
const char *scatint_potential_name(potential_type pt);
const char *scatint_reaction_name(reaction_type rt);
int scatint_has_solution(scatint_params *p);
double inverse_kinematics_erd_scale(double theta_cm, double theta_lab, double phi_lab);
double sigma_lab_to_cm(double theta_cm, double m12);
