/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_NUCLEAR_STOPPING_H
#define JABS_NUCLEAR_STOPPING_H
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <jibal_masses.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double k;
    double eps0;
#ifdef DEBUG
    const jibal_isotope *target; /* target is not needed, but it is here for debugging purposes */
#endif
} nucl_stop_pair;

/* Universal nuclear stopping (reduced, as a function of reduced energy epsilon) is tabulated for epsilon in [2^NUCLEAR_STOPPING_TABLE_EXP_MIN, 2^NUCLEAR_STOPPING_TABLE_EXP_MAX)
 * with 2^NUCLEAR_STOPPING_TABLE_BITS segments per octave. Each segment is a cubic Hermite polynomial. Maximum relative error compared to the formula is about 5e-8. */
#define NUCLEAR_STOPPING_TABLE_BITS (5)
#define NUCLEAR_STOPPING_TABLE_EXP_MIN (-20)
#define NUCLEAR_STOPPING_TABLE_EXP_MAX (24)
#define NUCLEAR_STOPPING_TABLE_SEGMENTS ((NUCLEAR_STOPPING_TABLE_EXP_MAX - NUCLEAR_STOPPING_TABLE_EXP_MIN) << NUCLEAR_STOPPING_TABLE_BITS)

typedef struct {
    nucl_stop_pair *t; /* table, indexing refers to ith element in jibal_isotopes (as given to nuclear_stopping_new()) */
    size_t n_isotopes;
    int refcount; /* nuclear_stopping_shared_copy() increases refcount, free() decreases. Actual free on zero. */
    const jibal_isotope *incident;
    double *universal; /* Polynomial coefficients (4 per segment) of tabulated universal nuclear stopping */
} nuclear_stopping;

nuclear_stopping *nuclear_stopping_new(const jibal_isotope *incident, const jibal_isotope *isotopes);
nuclear_stopping *nuclear_stopping_shared_copy(nuclear_stopping *ns);
void nuclear_stopping_free(nuclear_stopping *ns);
double nuclear_stopping_universal_exact(double epsilon); /* Reduced nuclear stopping (ZBL universal), calculated from the formula */
inline double nuclear_stopping_universal(const nuclear_stopping *ns, double epsilon) { /* Same as nuclear_stopping_universal_exact(), but from table */
    uint64_t bits;
    memcpy(&bits, &epsilon, sizeof(bits));
    int e = (int)(bits >> 52) - 1023; /* Binary exponent, epsilon is positive */
    if(e < NUCLEAR_STOPPING_TABLE_EXP_MIN || e >= NUCLEAR_STOPPING_TABLE_EXP_MAX) {
        return nuclear_stopping_universal_exact(epsilon);
    }
    uint64_t m_bits = (bits & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL; /* Mantissa as a number in [1, 2) */
    double m;
    memcpy(&m, &m_bits, sizeof(m));
    double u = (m - 1.0) * (double)(1 << NUCLEAR_STOPPING_TABLE_BITS);
    int j = (int)u; /* Segment within octave */
    double t = u - j; /* Position within segment, [0, 1) */
    const double *c = ns->universal + 4 * (((e - NUCLEAR_STOPPING_TABLE_EXP_MIN) << NUCLEAR_STOPPING_TABLE_BITS) + j);
    return c[0] + t * (c[1] + t * (c[2] + t * c[3]));
}
#ifdef __cplusplus
}
#endif
#endif //JABS_NUCLEAR_STOPPING_H
