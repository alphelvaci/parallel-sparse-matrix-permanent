#pragma ide diagnostic ignored "NullDereference"
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <stdlib.h>

#define COLS 35
#define ROWS COLS

#include "skipper.h"

double skip_per(crs_t crs, ccs_t ccs) {
    return 0.0;
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
