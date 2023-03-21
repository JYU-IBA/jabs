#include <string.h>
#include <stdlib.h>
#include <jibal_phys.h>
#include "histogram.h"

extern inline double jabs_histogram_get(const jabs_histogram *h, size_t i);

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

int jabs_histogram_compare(const jabs_histogram *h1, const jabs_histogram *h2, size_t low, size_t high, double *out) {
    if(!h1 || !h2)
        return EXIT_FAILURE;
    if(h1->n == 0 || h2->n == 0)
        return EXIT_FAILURE;
    if(high >= h1->n || high >= h2->n)
        return EXIT_FAILURE;
    if(low >= high)
        return EXIT_FAILURE;
    double sum = 0.0;
    size_t n = 0;
    for(size_t i = low; i <= high; i++) {
        if(h2->bin[i] == 0.0) {
            continue;
        }
        sum += pow2((h1->bin[i] - h2->bin[i]))/h2->bin[i];
        n++;
    }
    *out = sqrt(sum)/(1.0*n);
    return EXIT_SUCCESS;
}


double jabs_histogram_roi(const jabs_histogram *h, size_t low, size_t high) {
    if(!h || h->n == 0)
        return 0.0;
    if(low >= h->n)
        return 0.0;
    if(high >= h->n)
        high = h->n - 1;
    double sum = 0.0;
    for(size_t i = low; i <= high; i++) {
        sum += h->bin[i];
    }
    return sum;
}

size_t jabs_histogram_channels_in_range(const jabs_histogram *h, size_t low, size_t high) {
    if(!h || h->n == 0) {
        return 0;
    }
    if(high < low) {
        return 0;
    }
    if(low >= h->n) {
        return 0;
    }
    if(high >= h->n) {
        high = h->n - 1;
    }
    return high - low + 1;
}
