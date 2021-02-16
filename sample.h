#ifndef _SAMPLE_H_
#define _SAMPLE_H_

#include <stdio.h>
#include <jibal_masses.h>
#include <jibal_material.h>
#include <jibal_layer.h>

typedef struct {
    int n_isotopes;
    int n_ranges;
    const jibal_isotope **isotopes; /* table, size is n_isotopes */
    double *cranges;
    double *cbins; /* 2D-table: size is n_isotopes * n_ranges  */
    int i_range_accel;
} sample;

double get_conc(sample *s, double x, int i);
int get_concs(sample *s, double x, double *out);
int get_range_bin(sample *s, double x);

sample sample_from_layers(jibal_layer **layers, int n_layers);
void print_sample(FILE *f, sample *sample);
#endif /* _SAMPLE_H_ */
