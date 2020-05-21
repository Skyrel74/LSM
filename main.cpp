#include <stdio.h> 

#include <stdlib.h> 



int slau(int n, const double** a, const double* b, double* x);



int main(int argc, char* argv[]) {

    int n, i, j;
    double check;
    FILE* file;

    file = fopen("first.txt", "r");
    fscanf(file, "%i", &n);

    double** a = (double**)malloc(n * sizeof(double*));
    double* b = (double*)malloc(n * sizeof(double));
    double* x = (double*)malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        a[i] = (double*)malloc(n * sizeof(double));
        for (j = 0; j < n; j++) {
            fscanf(file, "%lf", &a[i][j]);
            printf("%12.2fx%i", a[i][j], j + 1);
        }
        fscanf(file, "%lf", &b[i]);
        printf("%12.2f \n", b[i]);
    }

    putchar('\n');

    if (slau(n, a, b, x)) {
        printf("\nThe roots of the system:\n");
        for (i = 0; i < n; i++)
            printf("  x%i = %.3lf \n", i + 1, x[i]);
        printf("\nCheck:\n");
        for (i = 0; i < n; i++) {
            check = 0.;
            for (j = 0; j < n; j++) {
                check += a[i][j] * x[j];
                printf("%7.2lf*%.2lf ", a[i][j], x[j]);
                if (j < (n - 1))
                    printf(" +");
                else
                    printf(" =");
            }
            printf("  %.2lf\n", check);
        }
    } else printf("The system cannot be solved by this method!\n");

    for (i = 0; i < n; i++)
        free(a[i]);

    free(a);
    free(b);
    free(x);

    return 0;
}



int slau(int n, const double** a, const double* b, double* x) {

    int i, j, k;
    double tmp;
    double** LU, * y;

    LU = (double**)malloc(n * sizeof(double*));

    for (i = 0; i < n; i++)
        LU[i] = (double*)malloc(n * sizeof(double));

    y = (double*)malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            tmp = 0.;
            if (j > i) {
                for (k = 0; k < i; k++)
                    tmp += LU[i][k] * LU[k][j];
                LU[i][j] = (a[i][j] - tmp) / LU[i][i];
            }
            else {
                for (k = 0; k < j; k++)
                    tmp += LU[i][k] * LU[k][j];
                LU[i][j] = a[i][j] - tmp;

                if (LU[i][j] == 0. && i == j) {
                    for (k = 0; k < n; k++) {
                        free(LU[k]);
                    }

                    free(LU);
                    free(y);

                    return 0;
                }
            }
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%12.4lf\t", LU[i][j]);
        putchar('\n');
    }

    for (i = 0; i < n; i++) {
        tmp = 0.;
        for (k = 0; k < i; k++)
            tmp += LU[i][k] * y[k];
        y[i] = (b[i] - tmp) / LU[i][i];
    }

    for (i = n - 1; i >= 0; i--) {
        tmp = 0.;
        for (k = i + 1; k < n; k++)
            tmp += LU[i][k] * x[k];
        x[i] = y[i] - tmp;
    }

    for (i = 0; i < n; i++) free(LU[i]);

    free(LU);
    free(y);

    return 1;
}