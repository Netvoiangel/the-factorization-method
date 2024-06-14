# define size 129
#include <stdio.h>
#include <math.h>
double f(double x) {
    return sin(x);
}
double f_der(double x) {
    return exp(x);
}
double s_d(double x, double s) {
    return -(-exp(x) - s + (exp(x) + 1) * s * s) / (exp(x) + 1);
}
double t_d(double x, double t, double s) {
    return (exp(x) - 1 + t - (exp(x) + 1) * s * t) / (exp(x) + 1);
}
double fst(double y, double s, double t) {
    return s * y + t;
}
double euler_y(double y1, double h, double s1, double s2, double t1, double t2) {
    double y2, ys;
    ys = y1 + h / 2 * fst(y1, s1, t1);
    y2 = y1 + h * fst(ys, s2, t2);
    return y2;
}
double euler_s(double x, double y1, double h) {
    double y2, ys;
    ys = y1 + h / 2 * s_d(x, y1);
    y2 = y1 + h * s_d(x + h / 2, ys);
    return y2;
}
double euler_t(double x, double y1, double h, double s1, double s2) {
    double y2, ys;
    ys = y1 + h / 2 * t_d(x, y1, s1);
    y2 = y1 + h * t_d(x + h / 2, ys, s2);
    return y2;
}

int main() {
    // (exp(x) + 1) * y'' - y' - y * exp(x) = exp(x)
    // y = exp(x) - 1
    double a = 0, b = 1;
    double A, B, beta = 0.8, p = 5, h, x = b;
    A = f(a);
    B = f(b) + beta * f_der(b);
    double s1 = -1 / 0.3;
    double t1 = B / beta;
    double y1 = f(b);
    double S[size] = { 0 };
    double T[size] = { 0 };
    double Y[size] = { 0 };
    int num = 1;
    h = (b - a) / pow(2, p + 2);
    S[0] = -1 / beta;
    x = b - h;
    while (x >= a) {
        S[num] = euler_s(x, S[num - 1], -h);
        num++; x -= h;
    }
    h *= 2;
    x = b - h;
    num = 1;
    T[0] = A / beta;
    while (x >= a) {
        T[num] = euler_t(x, T[num - 1], -h, S[2 * num], S[2 * num - 1]);
        num++; x -= h;
    }
    h *= 2;
    x = a + h;
    Y[0] = A;
    num = pow(2, p + 2);
    int n = 1;
    while (x <= b) {
        Y[n] = euler_y(Y[n - 1], h, S[num], S[num - 2], T[num / 2], T[num / 2 - 1]);
        num -= 4; x += h; n++;
    }
    x = a;
    for (int i = 0; i < 34; i++) printf("%lf %lf\n", Y[i], f(x+i*h));
}