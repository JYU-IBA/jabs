#include <stdlib.h>
#include <jibal_cs.h>
#include "testplugin.h"
#include "testplugin_version.h"

const char *name(void) {
    return "JaBS test plugin";
}

const char *version(void) {
    return testplugin_VERSION;
}

jabs_plugin_type type(void) {
    return JABS_PLUGIN_CS;
}

jabs_plugin_reaction *reaction_init(const jibal_isotope *isotopes, const jibal_isotope *incident, const jibal_isotope *target, int *argc, char * const **argv) {
    (void) isotopes;
    (void) argc;
    (void) argv;
    jabs_plugin_reaction *r = malloc(sizeof(jabs_plugin_reaction));
    r->incident = incident;
    r->target = target;
    r->product = incident; /* RBS, EBS */
    r->product_heavy = target;
    r->cs = testplugin_cs;
    r->E_min = 0.0;
    r->E_max = 1000.0 * C_MEV;
    r->yield = 1.0;
    return r;
}

void reaction_free(jabs_plugin_reaction *r) {
    free(r);
}

double testplugin_cs(const jabs_plugin_reaction *r, double theta, double E) {
    return jibal_cross_section_rbs(r->incident, r->target, theta, E, JIBAL_CS_RUTHERFORD); /* just a simple test */
}
