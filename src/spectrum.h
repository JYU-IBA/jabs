/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_SPECTRUM_H
#define JABS_SPECTRUM_H

#include "detector.h"
#include "simulation.h"
#include "histogram.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RESULT_SPECTRA_EXPERIMENTAL (0)
#define RESULT_SPECTRA_SIMULATED (1)
#define RESULT_SPECTRA_UNCERTAINTY_NEGATIVE (2) /* Simulated spectrum minus uncertainty (fit error?) */
#define RESULT_SPECTRA_UNCERTAINTY_POSITIVE (3) /* Simulated spectrum plus uncertainty (fit error?) */
#define RESULT_SPECTRA_N_FIXED (4)
#define RESULT_SPECTRA_REACTION_SPECTRUM(x) ((x) + RESULT_SPECTRA_N_FIXED)

typedef struct {
    jabs_histogram *histo;
    char *name;
    const jibal_isotope *target_isotope;
    reaction_type type;
} result_spectrum;

typedef struct result_spectra {
    result_spectrum *s; /* Use macros to determine index (see above), starts with experimental and simulated spectra, then there are others, total n_spectra (below) */
    size_t n_spectra; /*  RESULT_SPECTRA_N_FIXED + number of reaction spectra */
} result_spectra;

result_spectra *result_spectra_alloc(size_t n);
void result_spectra_free(result_spectra *spectra);
int result_spectra_copy(result_spectra *dest, const result_spectra *src);
int result_spectrum_copy(result_spectrum *dest, const result_spectrum *src);
int result_spectrum_set(result_spectrum *dest, const jabs_histogram *h, const char *name, const jibal_isotope *target_isotope, reaction_type type);
size_t result_spectra_n_ch(const result_spectra *spectra);
jabs_histogram *result_spectrum_histo(const result_spectra *spectra, size_t i_spectrum);
jabs_histogram *result_spectra_reaction_histo(const result_spectra *spectra, size_t i_reaction);
jabs_histogram *result_spectra_simulated_histo(const result_spectra *spectra);
jabs_histogram *result_spectra_experimental_histo(const result_spectra *spectra);
jabs_histogram *result_spectra_uncertainty_histo_negative(const result_spectra *spectra);
jabs_histogram *result_spectra_uncertainty_histo_positive(const result_spectra *spectra);

jabs_histogram *spectrum_read(const char *filename, size_t skip, size_t channels_max, size_t column, size_t compress);
jabs_histogram *spectrum_read_detector(const char *filename, const detector *det);
#ifdef __cplusplus
}
#endif
#endif // JABS_SPECTRUM_H
