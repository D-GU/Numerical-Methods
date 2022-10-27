#include <iostream>
#include <vector>
#include "library.h"

double f(double x) {
    /// Function
    return exp(x) * x;
}

double fDerivative(double x, double y) {
    /// Derivative of a function
    return exp(x) + y;
}

void rungeTest() {
    double start = 0;
    double end = 1;
    double h = 0.3;
    double eps = 1e-4;
    double y0 = 0;

    std::vector<double> solution{0, 0};

    Runge_Kutta runge(start, end, y0, h, eps);
    solution = runge.rungeKuttaSolution(fDerivative, f);

}

int main() {
    rungeTest();
}