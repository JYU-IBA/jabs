#include "gsl_inline.h"

#ifndef DEBUG
extern inline void jabs_gsl_vector_set(gsl_vector *v, size_t i, double x);
extern inline double jabs_gsl_vector_get(const gsl_vector *v, size_t i);
extern inline void jabs_gsl_matrix_set(gsl_matrix *m, size_t i, size_t j, double x);
extern inline double jabs_gsl_matrix_get(const gsl_matrix *m, size_t i, size_t j);
#endif
