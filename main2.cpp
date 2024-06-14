#include <iostream>
#include <math.h>

// vairan 2 v 1

double ds(double x, double s) {
    return -pow(s, 2) - cos(x) * s - sin(x);
}

double eulerMethod(double x, double y1, double h) {
    double y2, ys;
    ys = y1 + h / 2 * ds(x, y1);
    y2 = y1 + h * ds(x + h / 2, ys);
    return y2;
}

double* EulerProcess(double x_a, double y_a, double x_b, int div)
{
    double h = (x_b-x_a)/div;
    double* sol = (double*)malloc((div+1)*sizeof(double));
    sol[0] = y_a;
    for (int i = 1; i <= div; i++) {
        sol[i] = eulerMethod(x_a, sol[i-1], h);
        x_a += h;
    }
    return sol;
}


int main() {
    double a = 0, b = M_PI/2;
    double p = 5;
    double n = pow(2, p);
    double h = (b-a)/n, x_a = a;

    double* sol_s = EulerProcess(a, 2, b, 4 * n);

    double* sol_x = (double*)malloc((n+1)*sizeof(double));
    for (int i = 1; i <= n; i++) {
        sol_x[i] = x_a;
        x_a += h;
    }


}