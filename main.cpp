#include <iostream>
#include <math.h>

// vairan 2 v 1

double ds(double x, double s) {
    return -pow(s, 2) - cos(x) * s - sin(x);
}

double dt(double x, double *s, int i, double t) {
    return 1 - sin(x) - s[i] * t - t * cos(x);
}

double y(double x) {
    return sin(x);
}

double dy(double *s, double *y, double *t, int i) {
    return s[i] * y[i] + t[i];
}

double eulerMethod(double x, double y1, double h) {
    double y2, ys;
    ys = y1 + h / 2 * ds(x, y1);
    y2 = y1 + h * ds(x + h / 2, ys);
    return y2;
}

double eulerMethod2(double x, double y1, double h, int i, double *s) {
    double y2, ys;
    ys = y1 + h / 2 * dt(x, s, i, y1);
    y2 = y1 + h * dt(x + h / 2, s, i, ys);
    return y2;
}

double* EulerProcess(double x_a, double y_a, double x_b, int n)
{
    double h = (x_b-x_a)/n;
    double* sol = (double*)malloc((n+1)*sizeof(double));
    sol[0] = y_a;
    for (int i = 1; i <= n; i++) {
        sol[i] = eulerMethod(x_a, sol[i-1], h);
        x_a += h;
    }
    return sol;
}

double* EulerProcess2(double x_a, double y_a, double x_b, int n, double *s)
{
    double h = (x_b-x_a)/n;
    double* sol = (double*)malloc((n+1)*sizeof(double));
    sol[0] = y_a;
    for (int i = 1; i <= n; i++) {
        sol[i] = eulerMethod2(x_a, sol[i-1], h, i - 1, s);
        x_a += h;
    }
    return sol;
}

int main() {
    double a = 0, b = M_PI/2;
    double p = 8;
    int n = pow(2, p);
    double h = (b-a)/n, x_a = a;

    double* sol_s = EulerProcess(a, 75, b, n);

    double* sol_x = (double*)malloc((n+1)*sizeof(double));
    for (int i = 0; i <= n; i++) {
        sol_x[i] = x_a;
        x_a += h;
    }

    double* sol_t = EulerProcess2(a, 0, b, n, sol_s);

    double* sol_y = (double*)malloc((n+1)*sizeof(double));
    for (int i = 0; i <= n; i ++) {
        sol_y[i] = y(sol_x[i]);
    }

    double* sol_dt = (double*)malloc((n+1)*sizeof(double));
    for (int i = 0; i <= n; i ++) {
        sol_dt[i] = dy(sol_s, sol_y, sol_t, i);
    }

    double* sol_y2 = (double*)malloc((n+1)*sizeof(double));
    sol_y2[n] = y(b);

    for (int i = n - 1; i >= 0; i --) {
        sol_y2[i] = sol_y2[i+1] - h * (sol_s[i+1] * sol_y[i+1] + sol_t[i+1]);
    }

    FILE *output;
    output = fopen("data.txt", "w");

    FILE *factError;
    factError = fopen("factError.txt", "w");

    for (int i = 0; i <= n; i ++) {
        fprintf(output, "%.16f %.16f %.16f\n", sol_x[i], sol_y2[i], sol_y[i]);
        double error = abs(sol_y2[i] - sol_y[i]);
        fprintf(factError, "%.16f %.16f\n", sol_x[i], error);
    }    
}