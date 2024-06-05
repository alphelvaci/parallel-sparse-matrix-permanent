#define COLS 5
#define ROWS COLS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "skipper.h"

int gray_code(int n) {
    return n ^ (n >> 1);
}

int gray_diff(int g1, int g2) {
    return g1 ^ g2;
}

int next_g(int g, int n, double* x, crs_t crs) {
    int max_next = g + 1;
    for (int i = 0; i < n; i++) {
        if (x[i] == 0) {
            for (int ptr = crs.row_pointers[i]; ptr < crs.row_pointers[i + 1]; ptr++) {
                int j = crs.columns[ptr];
                int next = (g < (1 << j)) ? (1 << j) : (g + (1 << (j + 1)) - ((g - (1 << j)) % (1 << (j + 1))));
                if (next > max_next) {
                    max_next = next;
                }
            }
        }
    }
    return max_next;
}

double skip_per(crs_t crs, ccs_t ccs) {
    int n = ROWS;
    double x[ROWS];
    double p = 1.0;

    // Initial computation for x
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int ptr = crs.row_pointers[i]; ptr < crs.row_pointers[i + 1]; ptr++) {
            sum += crs.row_values[ptr];
        }
        x[i] = crs.row_values[crs.row_pointers[i + 1] - 1] - sum / 2.0;
    }

    for(int i=0; i<n; i++){
        p *= x[i];
    }

    int g_prev = 0;
    int g = 1;
    while (g < (1 << (n - 1))) {
        int grdiff = gray_diff(g, g_prev);
        for (int j = 0; j < n; j++) {
            if ((grdiff & (1 << j)) != 0) {
                grdiff = grdiff & ~(1 << j);
                int s = 2 * ((gray_code(g) & (1 << j)) != 0) - 1;
                for (int ptr = ccs.column_pointers[j]; ptr < ccs.column_pointers[j + 1]; ptr++) {
                    int row = ccs.rows[ptr];
                    double val = ccs.column_values[ptr];
                    x[row] += s * val;
                }
            }
        }

        double prod = 1.0;
        for (int i = 0; i < n; i++) {
            prod *= x[i];
        }
        p += ((g % 2 == 0) ? 1 : -1) * prod; 

        g_prev = g;
        g = next_g(g, n, x, crs);
    }

    return p * (4 * (n % 2) - 2);
}

// Function to read the matrix from a file and store it in CRS and CCS formats
void read_matrix(const char *filename, crs_t *crs, ccs_t *ccs) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    int rows, non_zeros;
    int cols = COLS;
    if (fscanf(file, "%d %d", &rows, &non_zeros) != 2 || rows != ROWS || cols != COLS) {
        perror("Invalid file format");
        fclose(file);
        exit(EXIT_FAILURE);
    }


    crs->row_pointers = (int *)malloc((ROWS + 1) * sizeof(int));
    crs->columns = (int *)malloc(non_zeros * sizeof(int));
    crs->row_values = (double *)malloc(non_zeros * sizeof(double));

    ccs->column_pointers = (int *)malloc((COLS + 1) * sizeof(int));
    ccs->rows = (int *)malloc(non_zeros * sizeof(int));
    ccs->column_values = (double *)malloc(non_zeros * sizeof(double));


    for (int i = 0; i <= ROWS; i++) crs->row_pointers[i] = 0;
    for (int j = 0; j <= COLS; j++) ccs->column_pointers[j] = 0;

    double matrix[ROWS][COLS] = {0.0};
    for (int i = 0; i < non_zeros; i++) {
        int row, col;
        double value;
        if (fscanf(file, "%d %d %lf", &row, &col, &value) != 3) {
            perror("Invalid file format");
            fclose(file);
            exit(EXIT_FAILURE);
        }
        matrix[row][col] = value;
    }
    fclose(file);

    // Convert to CRS format
    int nz = 0;
    for (int i = 0; i < ROWS; i++) {
        crs->row_pointers[i] = nz;
        for (int j = 0; j < COLS; j++) {
            if (matrix[i][j] != 0.0) {
                crs->columns[nz] = j;
                crs->row_values[nz] = matrix[i][j];
                nz++;
            }
        }
    }
    crs->row_pointers[ROWS] = nz;

    // Convert to CCS format
    nz = 0;
    for (int j = 0; j < COLS; j++) {
        ccs->column_pointers[j] = nz;
        for (int i = 0; i < ROWS; i++) {
            if (matrix[i][j] != 0.0) {
                ccs->rows[nz] = i;
                ccs->column_values[nz] = matrix[i][j];
                nz++;
            }
        }
    }
    ccs->column_pointers[COLS] = nz;
}

int main(int argc, char *argv[]) {
    const char *filename = "test5.mat";
    crs_t crs;
    ccs_t ccs;

    read_matrix(filename, &crs, &ccs);

    double permanent = skip_per(crs, ccs);
    printf("SkipPer Permanent: %lf\n", permanent);

    // Free allocated memory
    free(crs.row_pointers);
    free(crs.columns);
    free(crs.row_values);
    free(ccs.column_pointers);
    free(ccs.rows);
    free(ccs.column_values);

    return 0;
}
