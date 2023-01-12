#include <jibal_masses.h>

typedef enum jabs_plugin_type {
    PLUGIN_NONE = 0,
    PLUGIN_CS = 1, /* Cross section evaluation plugin */
    PLUGIN_SPECTRUM_READER = 2 /* Spectrum reader plugin */
} jabs_plugin_type;

typedef struct jabs_plugin_reaction {
    const jibal_isotope *incident;
    const jibal_isotope *target;
    const jibal_isotope *product;
    const jibal_isotope *product_heavy;
} jabs_plugin_reaction;
