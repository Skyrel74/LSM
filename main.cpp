#include <stdio.h> 

#include <stdlib.h> 

#include <math.h> 



int slau(int n, double** a, double* x);

double pol(int n, double x, double* c);

int lms(int n, int p, double* x, double* y, double* c);



int main(int argc, char* argv[]) {

    int i, n, p;

    double dx;

    FILE* file;

    file = fopen("vector.txt", "r");

    fscanf(file, "%i", &n);

    fscanf(file, "%i", &p);

    double* x = (double*)malloc(n * sizeof(double));

    double* y = (double*)malloc(n * sizeof(double));

    double* c = (double*)malloc((p + 1) * sizeof(double));

    printf("x ");

    for (i = 0; i < n; i++) {

        fscanf(file, "%lf", &x[i]);

        printf("%5.2lf ", x[i]);

    }

    printf("\ny ");

    for (i = 0; i < n; i++) {

        fscanf(file, "%lf", &y[i]);

        printf("%5.2lf ", y[i]);

    }



    if (lms(n, p, x, y, c)) {

        dx = x[0];

        while (dx <= x[n - 1]) {

            printf("\n%1.1lf %lf", dx, pol(n, dx, c));

            dx += 0.1;

        }

    }

    free(x);

    free(y);

    return 0;

}



int slau(int n, double** a, double* x) {

    int i, j, k;

    double** A, tmp;

    A = (double**)malloc((n + 1) * sizeof(double*));

    for (i = 0; i < n; i++) {

        A[i] = (double*)malloc(n * sizeof(double));

        for (j = 0; j <= n; j++) {

            A[i][j] = a[i][j];

        }

    }

    for (i = 0; i < n; i++) {

        tmp = A[i][i];

        if (fabs(tmp) < 1e-6) {

            for (j = 0; j <= n; j++)

                free(A[j]);

            free(A);

            return 0;

        }

        for (j = i; j <= n; j++) {

            A[i][j] /= tmp;

        }

        for (j = 0; j < n; j++) {

            if (j == i)

                continue;

            tmp = A[j][i];

            for (k = 0; k <= n; k++)

                A[j][k] -= tmp * A[i][k];

        }

    }

    for (i = 0; i < n; i++) {

        x[i] = A[i][n];

        free(A[i]);

    }

    free(A);

    return 1;

}



double pol(int n, double x, double* c) {

    double y;

    int i;

    y = c[n];

    for (i = n; i > 0; i--)

        y = y * x + c[i - 1];

    return y;

}



int lms(int n, int p, double* x, double* y, double* c) {

    double** a;

    int i, j, k, iret;

    a = (double**)malloc((p + 2) * sizeof(double*));

    for (i = 0; i <= p; i++)

        a[i] = (double*)malloc((p + 1) * sizeof(double));

    for (i = 0; i <= p; i++) {

        a[0][i] = 0.;

        for (k = 0; k < n; k++)

            a[0][i] += pow(x[k], i);

        k = 0;

        for (j = i - 1; j >= 0; j--) {

            k++;

            a[k][j] = a[0][i];

        }

    }

    for (j = 1; j <= p; j++) {

        a[j][p] = 0.;

        for (k = 0; k < n; k++)

            a[j][p] += pow(x[k], p + j);

        i = p;

        for (k = j + 1; k <= p; k++) {

            i--;

            a[k][i] = a[j][p];

        }

    }

    for (i = 0; i <= p; i++) {

        a[i][p + 1] = 0.;

        for (j = 0; j < n; j++)

            a[i][p + 1] += y[j] * pow(x[j], i);

    }

    iret = slau(p + 1, a, c);

    for (i = 0; i <= p; i++)

        printf("\nC[%i] = %lf", i, c[i]);

    for (i = 0; i <= p; i++)

        free(a[i]);

    free(a);

    return iret;

}