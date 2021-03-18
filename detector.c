#include <stdio.h>
#include "detector.h"
extern inline double detector_calibrated(const detector *det, size_t ch);

int detector_sanity_check(const detector *det) {
    if (det->resolution <= 0.0) {
        fprintf(stderr, "Warning: detector resolution (%g) is negative.\n", det->resolution);
        return -1;
    }
    if (det->slope <= 0.0) {
        fprintf(stderr, "Warning: detector slope (%g) is negative.\n", det->slope);
        return -1;
    }
    return 0;
}
