#include <string.h>
#include <stdlib.h>
#include "histogram.h"

jabs_histogram *jabs_histogram_clone(const jabs_histogram *h_orig) {
    if(!h_orig) {
        return NULL;
    }
    jabs_histogram *h = malloc(sizeof(jabs_histogram));
    size_t n = h_orig->n;
    h->n = n;
    size_t range_bytes = (n + 1) * sizeof(double);
    size_t bin_bytes = n * sizeof(double);
    h->range = malloc(range_bytes);
    h->bin =  malloc(bin_bytes);
    if(!h->range || !h->bin) {
        return NULL;
    }
    memcpy(h->range, h_orig->range, range_bytes);
    memcpy(h->bin, h_orig->bin, bin_bytes);
    return h;
}
