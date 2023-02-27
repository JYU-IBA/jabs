#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* This file provides replacements for some commonly used GSL vector/matrix inline functions. This should guarantee use of inlines no matter how GSL library was compiled. No range checking is performed. */

#ifdef DEBUG
#define jabs_gsl_vector_set(v, i, x) gsl_vector_set(v, i, x)
#define jabs_gsl_vector_get(v, i) gsl_vector_get(v, i)
#define jabs_gsl_matrix_set(m, i, j, x) gsl_matrix_set(m, i, j, x)
#define jabs_gsl_matrix_get(m, i, j) gsl_matrix_get(m, i, j)
#else
inline void jabs_gsl_vector_set(gsl_vector *v, const size_t i, const double x) {v->data[i * v->stride] = x;}
inline double jabs_gsl_vector_get(const gsl_vector *v, const size_t i) {return v->data[i * v->stride];}
inline void jabs_gsl_matrix_set(gsl_matrix *m, const size_t i, const size_t j, const double x) {m->data[i * m->tda + j] = x;}
inline double jabs_gsl_matrix_get(const gsl_matrix *m, const size_t i, const size_t j) {return m->data[i * m->tda + j];}
#endif
