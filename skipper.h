#ifndef COLS
#define COLS 32
#endif

#ifndef ROWS
#define ROWS COLS
#endif

#ifndef SKIPPER
#define SKIPPER
typedef struct crs_s {
    int* row_pointers;
    int* columns;
    double* row_values;
} crs_t;

typedef struct ccs_s {
    int* column_pointers;
    int* rows;
    double* column_values;
} ccs_t;

void skip_ord(int* matrix[ROWS][COLS]);
double skip_per(crs_t crs, ccs_t ccs, int threads);
#endif
