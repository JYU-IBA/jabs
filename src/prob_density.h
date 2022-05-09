#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef struct {
    double x;
    double p;
} prob_point;

typedef struct {
    size_t n;
    prob_point *points;
} prob_dist; /* Discrete probability distribution */


prob_dist *prob_dist_alloc(size_t n);
void prob_dist_free(prob_dist *pd);
prob_dist *prob_dist_gaussian(size_t n);
