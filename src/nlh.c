/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <assert.h>
#include "nlh.h"

const nlhcoeff *nlh_coeffs(int Z1, int Z2) { /* No checking of boundaries is performed, except in asserts. */
    assert(Z1 > 0 && Z1 <= NLH_MAX_Z);
    assert(Z2 > 0 && Z1 <= NLH_MAX_Z);
    int i, j;
    if(Z1 < Z2) {
        i = Z1 - 1;
        j = Z2 - 1;
    } else { /* Potentials are symmetric, switch Z1 and Z2 */
        i = Z2 - 1;
        j = Z1 - 1;
    }

    int index = i * (NLH_MAX_Z + (NLH_MAX_Z - i) - 1)/2 + j;
    const nlhcoeff *coeff = &nlhcoeffs[index];
    assert(coeff->Z1 == Z1 && coeff->Z2 == Z2);
    return coeff;
}
