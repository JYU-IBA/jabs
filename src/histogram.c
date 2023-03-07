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

void jabs_histogram_free(jabs_histogram *h) {
    if(!h) {
        return;
    }
    free(h->bin);
    free(h->range);
    free(h);
}

jabs_histogram *jabs_histogram_alloc(size_t n) {
    if(n == 0) {
        return NULL;
    }
    jabs_histogram *h = malloc(sizeof(jabs_histogram));
    h->n = n;
    h->range = calloc(n + 1, sizeof (double));
    h->bin = calloc(n, sizeof (double));
    if(!h->range || !h->bin) {
        jabs_histogram_free(h);
        return NULL;
    }
    return h;
}

void jabs_histogram_reset(jabs_histogram *h) {
    if(!h) {
        return;
    }
    memset(h->bin, 0, h->n * sizeof(double));
}

void jabs_histogram_scale(jabs_histogram *h, double scale) {
    if(!h) {
        return;
    }
    for(size_t i = 0; i < h->n; i++) {
        h->bin[i] *= scale;
     }
}
