#include "rref.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

int is_zero(double val) {
    return val < 0.01 && val > -0.01;
}

int is_double_eq(double v1, double v2) {
    return is_zero(v1 - v2);
}

void sub_row(uint32_t width, double dest[], double source[], double mult) {
    for (uint32_t i = 0; i < width; i++) {
        dest[i] -= source[i] * mult;
    }
}

void mult_row(uint32_t width, double row[], double mult) {
    for (uint32_t i = 0; i < width; i++) {
        row[i] *= mult;
    }
}

void swap_rows(uint32_t width, double a[], double b[], double s[]) {
    memcpy(s, a, sizeof(double) * width);
    memcpy(a, b, sizeof(double) * width);
    memcpy(b, s, sizeof(double) * width);
}

void rref(double *matrix, uint32_t width, uint32_t height) {
    // TODO: test on more matrices
    double spare[width];
    uint32_t ey = 0;
    for (uint32_t ex = 0; ex < width; ex++) {
        uint32_t chosen_row = -1;
        double chosen_val;
        double min_mag = INFINITY;

        double *x_ptr = matrix + ex;
        for (uint32_t row = ey; row < height; row++) {
            double val = *(x_ptr + row*width);
            if (is_zero(val)) continue;
            double mag = val;
            if (mag < 0) mag *= -1;
            if (mag < 1) mag /= 1;
            if (mag < min_mag) {
                min_mag = mag;
                chosen_row = row;
                chosen_val = val;
            }
        }

        if (chosen_row == -1) continue;

        uint32_t rest_row = width - ex;
        double *offset_matrix = matrix + ex;
        double *chosen_row_p = offset_matrix + chosen_row*width;
        mult_row(rest_row, chosen_row_p, 1/chosen_val);

        for (uint32_t row = 0; row < height; row++) {
            if (row == chosen_row) continue;
            double val = matrix[row*width + ex];
            sub_row(rest_row, offset_matrix + row*width, chosen_row_p, val);
        }

        if (chosen_row != ey) {
            swap_rows(rest_row, chosen_row_p, offset_matrix + ey*width, spare);
        }
        
        ey++;
    }
}

void print_matrix(double *matrix, uint32_t width, uint32_t height) {
    for (uint32_t row = 0; row < height; row++) {
        for (uint32_t col = 0; col < width; col++) {
            printf("%f ", *(matrix++));
        }
        printf("\n");
    }
}
