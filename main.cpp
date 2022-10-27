#include <iostream>
#include <vector>
#include "library.h"

double f(double x) {
    return exp(x) * x;
}

double fDerivative(double x, double y) {
    /// Derivative of a function
    return exp(x) + y;
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
    Runge_Kutta runge(start, y0, h0, eps);

    while (runge.start_ < end) {
        counter++;
        result = runge.jump(fDerivative);

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

int main() {
    rungeTest();
}