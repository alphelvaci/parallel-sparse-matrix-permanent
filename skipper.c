#pragma ide diagnostic ignored "NullDereference"
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define COLS 35
#define ROWS COLS

#include "skipper.h"

int pow_neg1(long long unsigned a) {
    return 1 - 2 * ((int)a & 0b1);
}

long long unsigned pow2(long long unsigned a) {
    return (long long unsigned)0b1 << a;
}

unsigned next(unsigned g, unsigned j) {
    //TODO
    return 0;
}

double skip_per(crs_t crs, ccs_t ccs) {
    double x[ROWS];

    double p = 1;
    for (int i=0; i < ROWS; i++) {
        double sum = 0;
        for (int ptr = crs.row_pointers[i]; ptr < crs.row_pointers[i+1]; ptr++) {
            sum += crs.row_values[ptr];
        }
        x[i] = crs.row_values[crs.row_pointers[i+1] - 1] - sum/2;
        p *= x[i];
    }

    long long unsigned prev_g = 0;
    long long unsigned curr_g;

    long long unsigned i = 1;
    while (i < pow2(ROWS-1)) {
        curr_g = i ^ (i >> 1);
        long long unsigned g_diff = prev_g ^ curr_g;
        long long unsigned mask = 0b1;
        unsigned bit_index = 0;
        while (g_diff > 0) {
            if (g_diff & mask) {
                g_diff = g_diff & ~mask;

                int curr_g_j = (curr_g & mask) != 0;

                double s = 2 * curr_g_j - 1;
                for (int ptr = ccs.column_pointers[bit_index]; ptr < ccs.column_pointers[bit_index+1]; ptr++) {
                    int row = ccs.rows[ptr];
                    double value = ccs.column_values[ptr];
                    x[row] += s * value;
                }
            }
            bit_index += 1;
            mask = mask << 1;
        }

        double prod = 1;
        for (int j=0; j < ROWS; j++) {
            prod *= x[j];
        }
        p += pow_neg1(i) * prod;

        prev_g = curr_g;

        //TODO calculate i using next
        i += 1;
    }

    return p * (4 * (ROWS % 2) - 2);
}

int main() {
    FILE* file = fopen("ey35_02.mat", "r");

    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    int cols, non_zeros;
    if (fscanf(file, "%d %d", &cols, &non_zeros) != 2 || cols != COLS) {
        perror("Invalid file");
        return EXIT_FAILURE;
    }

    double matrix[ROWS][COLS] = {0.0};
    for (int i=0; i < non_zeros; i++) {
        int row, col;
        double value;
        if (fscanf(file, "%d %d %lf", &row, &col, &value) != 3) {
            perror("Invalid file");
            return EXIT_FAILURE;
        }
        matrix[row][col] = value;
    }

    fclose(file);

    int rows = cols;
    printf("Matrix read. Size %dx%d with %d non-zeros.\n", rows, cols, non_zeros);

    crs_t matrix_crs = {
            malloc((rows+ 1) * sizeof(int)),
            malloc(non_zeros * sizeof(int)),
            malloc(non_zeros * sizeof(double))
    };
    matrix_crs.row_pointers[0] = 0;
    for (int i=0; i < rows; i++) {
        matrix_crs.row_pointers[i+1] = matrix_crs.row_pointers[i];
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] != 0) {
                matrix_crs.columns[matrix_crs.row_pointers[i+1]] = j;
                matrix_crs.row_values[matrix_crs.row_pointers[i+1]] = matrix[i][j];
                matrix_crs.row_pointers[i+1] += 1;
            }
        }
    }

    ccs_t matrix_ccs = {
            malloc((cols + 1) * sizeof(int)),
            malloc(non_zeros * sizeof(int)),
            malloc(non_zeros * sizeof(double))
    };
    matrix_ccs.column_pointers[0] = 0;
    for (int j = 0; j < cols; j++) {
        matrix_ccs.column_pointers[j+1] = matrix_ccs.column_pointers[j];
        for (int i = 0; i < rows; i++) {
            if (matrix[i][j] != 0) {
                matrix_ccs.rows[matrix_ccs.column_pointers[j+1]] = i;
                matrix_ccs.column_values[matrix_ccs.column_pointers[j+1]] = matrix[i][j];
                matrix_ccs.column_pointers[j+1] += 1;
            }
        }
    }

    if (matrix_crs.row_pointers[rows] != non_zeros || matrix_ccs.column_pointers[cols] != non_zeros) {
        perror("CRS/CCS error");

        free(matrix_crs.row_pointers);
        free(matrix_crs.columns);
        free(matrix_crs.row_values);

        free(matrix_ccs.column_pointers);
        free(matrix_ccs.rows);
        free(matrix_ccs.column_values);

        return EXIT_FAILURE;
    }

    double permanent = skip_per(matrix_crs, matrix_ccs);
    printf("The permanent is %lf.", permanent);

    free(matrix_crs.row_pointers);
    free(matrix_crs.columns);
    free(matrix_crs.row_values);

    free(matrix_ccs.column_pointers);
    free(matrix_ccs.rows);
    free(matrix_ccs.column_values);

    return EXIT_SUCCESS;
}
