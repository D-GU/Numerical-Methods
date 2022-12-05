#include "library.h"

/*
double f(double x, double u, double v) {
    return u + v + exp(x) * (1 - x * x);
}

double g(double x, double u, double v) {
    return 2 * u * (1 + x) - v;
}

double u(double x) {
    return x * exp(x);
}

double v(double x) {
    return x * x * exp(x);
}

double F(double u, double v, double alpha, double betta, double gamma) {
    return double(alpha * u + betta * v - gamma);
}

void rungeTest() {
    int counter = 0;

    double start = 0;
    double end = 1;
    double h0 = 0.3;
    double eps = 1e-4;
    double y0 = 0;
    double solution;

    std::vector<double> result(3);
    NumericalMethods::Runge_Kutta runge(start, y0, h0, eps);

    while (runge.start_ < end) {
        counter++;
        result = runge.jump(f);

        runge.start_ = result[0];
        runge.y0_ = result[1];

        solution = fabs((runge.y0_ - runge.start_) * exp(runge.start_));
        std::cout << "Шаг " << counter << ": Значение узла: " << result[2] << " |Приближенное решение: " << solution
                  << "\n";

        if (end - runge.start_ < runge.h_) {
            runge.h_ = end - runge.start_;
        }
    }
}

void shotsTest() {
    double start = 0.0;
    double end = 1.0;

    double alpha0 = 0.0;
    double betta0 = 1.0;
    double gamma0 = 0.0;

    double alpha1 = 1.0;
    double betta1 = -1.0;
    double gamma1 = 0.0;

    double ksi0 = 0;
    double ksi1 = 0;

    double h = 0.3;

    //double `eps = 1e-7
    double eps = 1e-4;

    NumericalMethods::ShootingMethod shots(f,
                                           g,
                                           F,
                                           start,
                                           end,
                                           alpha0,
                                           betta0,
                                           gamma0,
                                           alpha1,
                                           betta1,
                                           gamma1,
                                           h,
                                           eps,
                                           ksi0,
                                           ksi1);

    double sol = shots.shootingSolution();

}
*/

double f(double x) {
    return cos(x);
}

std::vector<double> g(std::vector<double> &vec, double p, double q) {
    std::vector<double> result(vec.size(), 0);
    for (int i = 0; i < result.size(); i++) {
        result[i] = -cos(vec[i]) - (p * sin(vec[i])) + (q * cos(vec[i]));
    }

    /*for (double &i: vec) {
        i = -cos(i) - (p * sin(i)) + (q * cos(i));
    }*/

    return result;
}

void gridTest() {
    double start = 0;
    double end = M_PI;

    int n = 32;
    int lin_space_num = n + 1;

    double p = 3;
    double q = -3;

    std::vector<double> alpha{1, 1};
    std::vector<double> betta{0, 0};
    std::vector<double> gamma{cos(start), cos(end)};

    std::vector<double> a(n + 1, 0);
    std::vector<double> b(n + 1, 0);
    std::vector<double> c(n + 1, 0);
    std::vector<double> d(n + 1, 0);

    NumericalMethods::GridMethod gridMethod(
            f, g, n, start, end, p, q, lin_space_num, alpha, betta, gamma, a, b, c, d
    );

    std::vector<double> s = gridMethod.solveGrid();
    std::vector<double> gr = gridMethod.getLinSpace();

    for (int i = 0; i < n + 1; i++) {
        std::cout << "x[" << i << "] = " << gr[i] << "\t y[" << i << "] = " << s[i] << " " << cos(gr[i]) << "\n";
    }
}

int main() {
    gridTest();
    return 0;
}