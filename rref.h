#include <stdint.h>

#ifndef RREF_H
#define RREF_H

int is_zero(double val);

int is_double_eq(double v1, double v2);

void rref(double *matrix, uint32_t width, uint32_t height);

void print_matrix(double *matrix, uint32_t width, uint32_t height);

#endif
