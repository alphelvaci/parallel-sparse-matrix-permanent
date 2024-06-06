#pragma ide diagnostic ignored "NullDereference"
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <stdlib.h>

#define COLS 7
#define ROWS COLS

#include "skipper.h"

long long unsigned pow2(unsigned exp) {
    return (long long unsigned)0b1 << exp;
}

long long unsigned next(long long unsigned g, unsigned j) {
    if (g < pow2(j)) {
        return pow2(j);
    } else {
        return g + pow2(j + 1) - ((g - pow2(j) % pow2(j + 1)));
    }
}

long long unsigned min_next(long long unsigned g, unsigned i, crs_t crs) {
    long long unsigned min = UINT64_MAX;
    for (int ptr = crs.row_pointers[i]; ptr < crs.row_pointers[i+1]; ptr++) {
        if (next(g, crs.columns[ptr]) < min) {
            min = next(g, crs.columns[ptr]);
        }
    }
    return min;
}

long long unsigned next_max(long long unsigned g, crs_t crs) {
    long long unsigned max = g + 1;
    for (int i = 0; i < ROWS; i++) {
        if(min_next(g, i, crs) > max) {
            max = min_next(g, i, crs);
        }
    }
    return max;
}

double skip_per(crs_t crs, ccs_t ccs) {
    double x[ROWS] = {0.0};

    for (int ptr = ccs.column_pointers[ROWS-1]; ptr < ccs.column_pointers[ROWS]; ptr++) {
        x[ccs.rows[ptr]] = ccs.column_values[ptr];
    }

    for (int i=0; i < ROWS; i++) {
        double sum = 0;
        for (int ptr = crs.row_pointers[i]; ptr < crs.row_pointers[i+1]; ptr++) {
            sum += crs.row_values[ptr];
        }
        x[i] -= sum / 2;
    }

    double p = 1;
    for (int i=0; i < ROWS; i++) {
        p *= x[i];
        if (p == 0) {
            break;
        }
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
                if (curr_g_j) {
                    for (int ptr = ccs.column_pointers[bit_index]; ptr < ccs.column_pointers[bit_index+1]; ptr++) {
                        x[ccs.rows[ptr]] += ccs.column_values[ptr];
                    }
                } else {
                    for (int ptr = ccs.column_pointers[bit_index]; ptr < ccs.column_pointers[bit_index+1]; ptr++) {
                        x[ccs.rows[ptr]] -= ccs.column_values[ptr];
                    }
                }
            }
            bit_index += 1;
            mask = mask << 1;
        }

        double prod = 1;
        for (int j=0; j < ROWS; j++) {
            prod *= x[j];
        }

        if ((int)prod != 0) {
            if (i & 0b1) {
                p -= prod;
            } else {
                p += prod;
            }
            i += 1;
        } else {
            i = next_max(i, crs);
        }

        prev_g = curr_g;
    }

    return p * (4 * (ROWS & 0b1) - 2);
}

int main() {
    FILE* file = fopen("test7.mat", "r");

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
