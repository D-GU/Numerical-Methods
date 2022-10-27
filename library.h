#ifndef NUMERICALMETHODSLIB_LIBRARY_H
#define NUMERICALMETHODSLIB_LIBRARY_H

#include <iostream>
#include <vector>
#include <string>
#include <random>

// Function, that returns random double number within the diapason
double randomDouble(double start, double end);

/*
 Dichotomy method. Method that finds the root of the function on specific section.

 Example:
        Dichotomy dichotomy(a, b, epsilon);
        dichotomy.calculateRoot(function);
*/

class Dichotomy {
public:
    explicit Dichotomy(double start, double end, double epsilon);

    virtual ~Dichotomy() = default;

public:
    virtual double calculateRoot(double function(double x), int *div_counter);

private:
    double start_; /// start of the section
    double end_; /// end of the section
    double pivot_; /// middle of the section (end of the new section)
    double eps_; /// accuracy

};

/*Class Secant. A derivative class of Dichotomy.
To initialize the object Secant uses the same constructor as Dichotomy class.
To calculate root of the function Secant uses the same method as Dichotomy class.

If in Dichotomy the end of the section calculates by the formula: (a + b) / 2
in secant the formula is ((-f(a) * b) + (a * f(b)) / (-f(a) * f(b)).
*/
class Secant : protected Dichotomy {
public:

    explicit Secant(double start,
                    double end,
                    double epsilon) : Dichotomy(start, end, epsilon) {};

    ~Secant() override = default;

public:
    double calculateRoot(double function(double x), int *div_counter) override;

private:
    double start_; /// start of the section
    double end_; /// end of the section
    double pivot_; /// middle of the section (end of the new section)
    double eps_; /// accuracy

};

/*Gaussian method of solving linear equations.

The object could be initialized in two different ways.
First one being reading the data from the files
and the second one is manually passing the data in to the constructor.

The data is the matrix A and vector b. In this project
matrix is banded-matrix.

Example:
       1) Gauss gauss(size,
                      MATRIX_PATH,
                      VECTOR_PATH,
                      true);
          solution_vec = gauss.findSolutionGauss();

       2) std::vector<double> matrix = {x1, x2, ..., xn};
          std::vector<double> vector = {y1, y2, ..., yn};
          Gauss gauss(size,
                      matrix,
                      vector,
                      true);
          std::vector<double> solution_vec = gauss.findSolutionGauss();
*/
class Gauss {
public:
    explicit Gauss(short size,
                   const std::string &matrix_filepath,
                   const std::string &vector_filepath,
                   bool create = true);

    explicit Gauss(short power,
                   std::vector<double> &matrix,
                   std::vector<double> &free_vector_);

    virtual ~Gauss() = default;

public:
    std::vector<double> findSolutionGauss();

private:
    std::vector<double> matrix_; /// matrix A
    std::vector<double> vec_; /// vector b
    std::vector<double> result_vec_; /// vector x(solution of linear system)

    short int size_; /// dim(A)
};

/*Newton Interpolation.

 Newton's method of function interpolation in
 specific point

 Example:
          std::vector<double> nodes = {1.0, 2.0, 3.0};
          double exp_point = 3.14;
          short int power = 3;
          Newton newton(f, nodes, exp_point, power);
          double result = newton.calculateNewton();
*/
class Newton {
public:

    explicit Newton(double function(double x),
                    const std::vector<double> &nodes,
                    double exp_point,
                    short power);

    virtual ~Newton() = default;

public:
    double calculateNewton();

private:
    std::vector<double> nodes_; /// nodes
    std::vector<double> fvalue_; /// value in nodes
    std::vector<double> difference_; /// table of divided difference
    std::vector<double> diag_vec_; /// diagonal vector of difference_ table

    short power_;

    double exp_point_;
};

/*Numeric method of integration

 Simpson's formula is the one of the quadratic
 integration formulas of Newton-Cotes.

 Example:
        IntegralMethod integral_method(start, end, eps, step);
        double result =  integral_method.calculateSimpson(f);
*/
class IntegralMethod {
public:
    explicit IntegralMethod(double start, double end, double eps, int step);

    virtual ~IntegralMethod() = default;

public:
    virtual double calculateSimpson(double function(double x));

private:
    double start_;
    double end_;
    double width_;
    double eps_;

    int step_;

};

/* Cubical Spline function interpolation in specific point
 Sweep method is method of three diagonal matrix where
 we find p,q and m coefficients

 Example:
        CubeSpline spline(short power,
                          double function(double x),
                          std::vector<double> &nodes,
                          double alpha0,
                          double alpha1,
                          double beta0,
                          double beta1,
                          double gamma0,
                          double gamma1,
                          double exp_point);
        double result  = spline.calculateCubeSpline();
*/
class CubeSpline {
public:
    explicit CubeSpline(short power,
                        double function(double x),
                        std::vector<double> &nodes,
                        double alpha0,
                        double alpha1,
                        double beta0,
                        double beta1,
                        double gamma0,
                        double gamma1,
                        double exp_point);

    virtual ~CubeSpline() = default;

public:
    double calculateCubeSpline();

private:
    void sweepMethod();

private:
    std::vector<double> nodes_;
    std::vector<double> fvalue_;;
    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
    std::vector<double> m_;
    std::vector<double> q_;
    std::vector<double> p_;

    double exp_point_;
    double alpha0_; /// Boundary condition coefficient a0
    double alpha1_; /// Boundary condition coefficient a1
    double beta0_; /// Boundary condition coefficient b0
    double beta1_; /// Boundary condition coefficient b1
    double gamma0_; /// Boundary condition a0 * g(x)' + a1 * g(x)" = g0
    double gamma1_; /// Boundary condition b0 * g(x)' + b1 * g(x)" = g1

    short power_;

};

/*MonteCarlo method of calculating integral
 We got two methods here

 First one is the method where we calculate
 random numbers within the square (1 x 1)

 Second one is the method where we calculate
 random numbers within the cube (1 x 1)

 Example:
        MonteCarlo monte_carlo(power);
        double result1 = monte_carlo.calculateMonteCarlo1(f);
       double result2 = monte_carlo.calculateMonteCarlo1(f);
*/
class MonteCarlo {
public:
    explicit MonteCarlo(long power);

    virtual ~MonteCarlo() = default;

public:
    double calculateMonteCarlo1(double function(double x, double y));

    double calculateMonteCarlo2(double function(double x, double y));

private:
    long power_;

    double random_x_ = 0;
    double random_y_ = 0;
    double random_z_ = 0;
    double summary_ = 0;

};

class Runge_Kutta {
public:
    Runge_Kutta(double start, double end, double y0, double h, double eps);

    ~Runge_Kutta() = default;

public:
    std::vector<double> jump(double function(double x, double y));

    double step(double function(double x, double y));

    std::vector<double> rungeKuttaSolution(double fDerivative(double x, double y), double f(double x));

public:
    double start_;
    double y0_;
    double h_;
    double end_;

private:
    double eps_;

    std::vector<double> result_{0, 0, 0};


};

#endif //NUMERICALMETHODSLIB_LIBRARY_H
