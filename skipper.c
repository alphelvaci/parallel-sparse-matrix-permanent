#pragma ide diagnostic ignored "NullDereference"
#pragma ide diagnostic ignored "cert-err34-c"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdbool.h>

// gcc -o skipper skipper.c -O3 -std=c11 -fopenmp -lm -DCOLS=35
// ./skipper sample_matrices/football.mat 16 1

#include "skipper.h"

long long unsigned pow2(unsigned exp) {
    return (long long unsigned)0b1 << exp;
}

// unused
long long unsigned next(long long unsigned g, unsigned j) {
    if (g < pow2(j)) {
        return pow2(j);
    } else {
        return g + pow2(j + 1) - ((g - pow2(j) % pow2(j + 1)));
    }
}

// unused
long long unsigned min_next(long long unsigned g, unsigned i, crs_t crs) {
    long long unsigned min = UINT64_MAX;
    for (int ptr = crs.row_pointers[i]; ptr < crs.row_pointers[i+1]; ptr++) {
        if (next(g, crs.columns[ptr]) < min) {
            min = next(g, crs.columns[ptr]);
        }
    }
    return min;
}

// unused
long long unsigned next_max(long long unsigned g, crs_t crs) {
    long long unsigned max = g + 1;
    for (int i = 0; i < ROWS; i++) {
        if(min_next(g, i, crs) > max) {
            max = min_next(g, i, crs);
        }
    }
    return max;
}

// unused
void skipOrd(int *colPerm, int *rowPerm, int *degs, int n, const crs_t *crs) {
    bool *rowVisited = (bool *)malloc(n * sizeof(bool));
    memset(rowVisited, false, n * sizeof(bool));
    memset(degs, 0, n * sizeof(int));
    
    for (int j = 0; j < n; j++) {
        for (int ptr = crs->row_pointers[j]; ptr < crs->row_pointers[j + 1]; ptr++) {
            degs[crs->columns[ptr]]++;
        }
    }

    for (int j = 0; j < n; j++) {
        int minDegIndex = -1;
        for (int k = 0; k < n; k++) {
            if (degs[k] != -1 && (minDegIndex == -1 || degs[k] < degs[minDegIndex])) {
                minDegIndex = k;
            }
        }
        colPerm[j] = minDegIndex;
        degs[minDegIndex] = -1;

        for (int ptr = crs->row_pointers[minDegIndex]; ptr < crs->row_pointers[minDegIndex + 1]; ptr++) {
            int row = crs->columns[ptr];
            if (!rowVisited[row]) {
                rowVisited[row] = true;
                rowPerm[j] = row;
                for (int k = 0; k < n; k++) {
                    if (crs->columns[k] == row) {
                        degs[k]--;
                    }
                }
            }
        }
    }

    free(rowVisited);
}

double skip_per(crs_t crs, ccs_t ccs, int threads) {
    // int *colPerm = (int *)malloc(ROWS * sizeof(int));
    // int *rowPerm = (int *)malloc(ROWS * sizeof(int));
    // int *degs = (int *)malloc(ROWS * sizeof(int));

    // Apply SkipOrd to reorder the matrix
    // skipOrd(colPerm, rowPerm, degs, ROWS, &crs);

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

    double p, init_prod = 1.0;
    for (int i=0; i < ROWS; i++) {
        init_prod *= x[i];
    }
    p = init_prod;

    unsigned long long start, end, chunk_size;
    start = 1;
    end = pow2(ROWS-1);
    chunk_size = (end - start + 1) / threads + 1;

    #pragma omp parallel num_threads(threads)
    {
        double local_x[ROWS];
        memcpy(local_x, x, sizeof(double) * ROWS);

        int tid = omp_get_thread_num();
        unsigned long long t_start = start + tid * chunk_size;
        unsigned long long t_end = fmin(start + ((tid+1) * chunk_size), end);
        double t_p = 0;

        unsigned long long t_gray;
        unsigned long long t_prev_gray = 0;

        
        for(long long unsigned i = t_start; i < t_end;) {
            t_gray = i ^ (i >> 1);
            unsigned bit_index = 0;
            long long unsigned mask = 0b1;
            long long unsigned g_diff = t_prev_gray ^ t_gray;
            while ( g_diff > 0) {
                if (g_diff & mask) {
                    g_diff = g_diff & ~mask;

                    int curr_g_j = (t_gray & mask) != 0;
                    if (curr_g_j) {
                        for (int ptr = ccs.column_pointers[bit_index]; ptr < ccs.column_pointers[bit_index+1]; ptr++) {
                            local_x[ccs.rows[ptr]] += ccs.column_values[ptr];
                        }
                    } else {
                        for (int ptr = ccs.column_pointers[bit_index]; ptr < ccs.column_pointers[bit_index+1]; ptr++) {
                            local_x[ccs.rows[ptr]] -= ccs.column_values[ptr];
                        }
                    }
                }
                bit_index += 1;
                mask = mask << 1;
            }

            double prod = 1;
            for (int j=0; j < ROWS; j++) {
                prod *= local_x[j];
            }

        if ((int)prod != 0) {
            t_p += ((i & 1ULL) ? -1.0 : 1.0) * prod;
            i++;
        } else {
            // i = next_max(i, crs);
            i++;
        }
        t_prev_gray = t_gray;

        }
        #pragma omp critical
        p += t_p;
    }

    return p * (4 * (ROWS & 0b1) - 2);
}

// Function to read the matrix from a file and store it in CRS and CCS formats
void read_matrix(FILE *file, crs_t *crs, ccs_t *ccs) {

    int rows, non_zeros;
    int cols = COLS;
    if (fscanf(file, "%d %d", &rows, &non_zeros) != 2 || rows != ROWS) {
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
            perror("Invalid file format of values");
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
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_file> <thread_num> <num_iterations>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *filename = argv[1];
    int threads = atoi(argv[2]);
    int iterations = atoi(argv[3]);
    crs_t crs;
    ccs_t ccs;

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    read_matrix(file, &crs, &ccs);
    printf("Matrix are read\n");
    double t_time = 0.0;
    double permanent;
    for(int i=0; i<iterations; i++){
        double start = omp_get_wtime();
        permanent = skip_per(crs, ccs, threads);
        double end = omp_get_wtime();
        t_time += (end - start);
    }
    printf("SkipPer Permanent: %.2e\n", permanent);
    printf("Average time for %d iterations: %f\n", iterations, t_time/iterations);

    // Free allocated memory
    free(crs.row_pointers);
    free(crs.columns);
    free(crs.row_values);
    free(ccs.column_pointers);
    free(ccs.rows);
    free(ccs.column_values);

    return 0;
}
