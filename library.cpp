#include "library.h"
#include "cmath"
#include <fstream>
#include <utility>


double randomDouble(double start, double end) {
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distribution(start, end);
    return distribution(eng);
}

NumericalMethods::Dichotomy::Dichotomy(double start,
                                       double end,
                                       double epsilon) :
        start_(start),
        end_(end),
        eps_(epsilon) {

}

double NumericalMethods::Dichotomy::calculateRoot(double (*function)(double),
                                                  int *div_counter) {
    /*
     * First we calculate function values on the start
     * and the end of the section.
     * If function has the same sign there
     * then it has no roots - program is done
     *
     * Otherwise, while absolute of (end - start) is
     * greater than our epsilon, we continue doing our actions
     * Finding the pivot and counting function value there and then
     * checking sign of the function values on start and on pivot gives
     * us the information of the function root.
     */

    double function_a_ = function(start_);
    double function_b_ = function(end_);
    double function_c_ = 0;

    if (function_a_ * function_b_ > 0) {
        std::cout << "Function has no root on this section." << "\n";
        return 0;
    }

    while (fabs(end_ - start_) >= eps_) {
        pivot_ = (start_ + end_) * 0.5;
        function_c_ = function(pivot_);
        *div_counter += 1;

        if (function_a_ * function_c_ <= 0) {
            end_ = pivot_;
            function_b_ = function_c_;
        } else {
            start_ = pivot_;
            function_a_ = function_c_;
        }

        if (fabs(end_ - start_) < eps_) {
            return pivot_;
        }
    }

    return 0;
}

double NumericalMethods::Secant::calculateRoot(double (*function)(double), int *div_counter) {
    /* In the secant method we do the same thing as we did
     * in the dichotomy method but the pivot calculation is different
     *
     * pivot = ((-f(a) * b) + (a * f(b)) / (-f(a) * f(b))
     * */

    double temp;

    double function_a_ = function(start_);
    double function_b_ = function(end_);
    double function_c_ = 0;

    while (function_a_ * function_b_ < 0) {
        pivot_ = (-function_a_ * end_ + function_b_ * start_) / (-function_a_ + function_b_);
        function_c_ = function(pivot_);
        *div_counter += 1;

        if (function_a_ * function_c_ <= 0) {
            end_ = pivot_;
            function_b_ = function_c_;
        } else {
            start_ = pivot_;
            function_a_ = function_c_;
        }

        if (fabs(function_c_) < eps_ && fabs(pivot_ - temp) < eps_) {
            return pivot_;
        }

        temp = pivot_;
    }

    return 0;
}

NumericalMethods::Gauss::Gauss(short size,
                               const std::string &matrix_filepath,
                               const std::string &vector_filepath,
                               bool create) : size_(size) {

    std::ifstream numbers;
    numbers.open(matrix_filepath);

    matrix_.resize(size_ * size_);
    vec_.resize(size_);

    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            numbers >> matrix_[i * size_ + j];
        }
    }
    numbers.close();

    std::ifstream free_row;
    free_row.open(vector_filepath);

    if (create) {
        for (int i = 0; i < size_; i++) {
            free_row >> vec_[i];
        }
    }

    free_row.close();
}

NumericalMethods::Gauss::Gauss(short size,
                               std::vector<double> &matrix,
                               std::vector<double> &free_vector) :
        size_(size),
        matrix_(matrix),
        vec_(free_vector) {
}

std::vector<double> NumericalMethods::Gauss::findSolutionGauss() {
    /* The Gaussian method of finding the algebraic system
     * solution.
     *
     * First we find the maximum element in each column
     * If max element is less the epsilon then it equals to 0
     * and matrix is invertible.
     *
     * We swap the rows until the max element is not
     * on diagonal
     *
     * Then we subtract the rows, and we eliminated matrix, so we could find
     * the coefficients
     *
     * We found coefficients via Reverse Gaussian algorithm.
     **/
    const double eps_ = 1e-5;
    int index_;

    double maximum_;
    double temp_;


    result_vec_.resize(size_);

    result_vec_ = vec_;

    for (int rows = 0; rows < size_; rows++) {
        maximum_ = fabs(matrix_[rows * size_ + rows]);
        index_ = rows;

        for (int i = rows + 1; i < size_; i++) {
            if (fabs(matrix_[rows * size_ + rows]) > maximum_) {
                maximum_ = fabs(matrix_[rows * size_ + rows]);
                index_ = i;
            }
        }

        if (maximum_ < eps_) {
            std::cout << "Matrix has no solution " << "\n";
            return {};
        }

        for (int j = 0; j < size_; j++) {
            std::swap(matrix_[rows * size_ + j], matrix_[index_ * size_ + j]);
        }

        std::swap(vec_[rows], vec_[index_]);

        for (int i = rows + 1; i < size_; i++) {
            temp_ = matrix_[i * size_ + rows] / matrix_[rows * size_ + rows];
            for (int j = rows; j < size_; j++) {
                matrix_[i * size_ + j] -= temp_ * matrix_[rows * size_ + j];
            }
            vec_[i] -= temp_ * vec_[rows];
        }
    }

    temp_ = 0.0;

    if (matrix_[(size_ - 1) * size_ + (size_ - 1)] == 0) {
        if (vec_[size_ - 1] == 0) {
            std::cout << "Matrix has infinite solutions." << "\n";
            return {};
        } else {
            std::cout << "Matrix has no solution." << "\n";
            return {};
        }
    } else {
        for (int i = size_ - 1; i >= 0; i--) {
            temp_ = 0.0;
            for (int j = i + 1; j < size_; j++) {
                temp_ += matrix_[i * size_ + j] * result_vec_[j];
            }
            result_vec_[i] = (vec_[i] - temp_) / matrix_[i * size_ + i];
        };
    }

    return result_vec_;
}

std::vector<double> NumericalMethods::Gauss::getReversedMatrix() {
    reversed_matrix_.resize((2 * size_) * size_);

    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            reversed_matrix_[i * size_ + j] = matrix_[i * size_ + j];
        }
    }

    for (int i = 0; i < 2 * size_; i++) {
        for (int j = size_; j < 2 * size_; j++) {
            (i + size_ == j) ? reversed_matrix_[i * size_ + j] = 1 : reversed_matrix_[i * size_ + j] = 0;
        }
    }

    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < 2 * size_; j++) {
            std::cout << reversed_matrix_[i * size_ + j] << " ";
        }
        std::cout << "\n";
    }

    const double eps_ = 1e-5;
    int index_;

    double maximum_;
    double temp_;

    for (int rows = 0; rows < size_; rows++) {
        maximum_ = fabs(matrix_[rows * size_ + rows]);
        index_ = rows;

        for (int i = rows + 1; i < size_; i++) {
            if (fabs(matrix_[rows * size_ + rows]) > maximum_) {
                maximum_ = fabs(matrix_[rows * size_ + rows]);
                index_ = i;
            }
        }

        if (maximum_ < eps_) {
            std::cout << "Matrix has no solution " << "\n";
            return {};
        }

        for (int j = 0; j < size_; j++) {
            std::swap(matrix_[rows * size_ + j], matrix_[index_ * size_ + j]);
        }

        std::swap(vec_[rows], vec_[index_]);

        for (int i = rows + 1; i < size_; i++) {
            temp_ = matrix_[i * size_ + rows] / matrix_[rows * size_ + rows];
            for (int j = rows; j < size_; j++) {
                matrix_[i * size_ + j] -= temp_ * matrix_[rows * size_ + j];
            }
            vec_[i] -= temp_ * vec_[rows];
        }
    }

    return reversed_matrix_;
}

NumericalMethods::Newton::Newton(double (*function)(double),
                                 const std::vector<double> &nodes,
                                 double exp_point,
                                 short power) :
        exp_point_(exp_point),
        power_(power + 1) {

    fvalue_.resize(power_);

    nodes_.resize(power_);
    nodes_ = nodes;

    for (int i = 0; i <= nodes_.size(); i++) {
        fvalue_[i] = function(nodes[i]);
    }

    difference_.resize(power_ * power_);

    for (int i = 0; i < power_; i++) {
        difference_[i * power_ + 0] = fvalue_[i];
    }

}

double NumericalMethods::Newton::calculateNewton() {
/* First we calculate the divided difference table
 * by the formula: (f_i(x) - f_j(x)) / (x_i - x_j)
 * We calculate this first, so we don't have to calculate
 * derivative everytime we have to.
 * Needed derivative will be on the table diagonal
 *
 * After that we just use Newton's formula to count the
 * value of the function in our point.
*/
    double mul = 1.0;
    double result = 0.0;

    diag_vec_.resize(power_ * power_);

    for (int i = 1; i < power_; i++) {
        for (int j = 1; j < power_; j++) {
            if (i < j) {
                difference_[i * power_ + j] = 0;
            } else {
                difference_[i * power_ + j] =
                        (difference_[i * power_ + (j - 1)] -
                         difference_[(j - 1) * power_ + (j - 1)]) /
                        (nodes_[i] - nodes_[j - 1]);
            }
        }
    }

    for (int i = 0; i < power_; i++) {
        diag_vec_[i] = difference_[i * power_ + i];
    }

    for (int i = 0; i < power_; i++) {
        mul = 1.0;
        for (int j = 0; j < i; j++) {
            mul *= (exp_point_ - nodes_[j]);
        }
        result += diag_vec_[i] * mul;
        std::cout << "Interval value " << i << ": " << result << "\n";
    }

    return result;
}

NumericalMethods::IntegralMethod::IntegralMethod(double start, double end, double eps, int step) :
        start_(start),
        end_(end),
        width_((end - start) / step),
        eps_(eps),
        step_(step) {

}

double NumericalMethods::IntegralMethod::calculateSimpson(double function(double x)) {
    /* For start, we count sum of f(0) + f(1)
     * s1 - sum of odd nodes
     * s2 - sum of even nodes
     *
     * After we find our we can count I1 via formula:
     * (b - a) / (3 * n) * (s0 + 4 * s1 + 2 * s2)
     *
     * Then we use Runge's rule which is used to evaluate the error
     * We basically divide each section by two, so we can calculate s0 and s1 again
     * and if I1 and I2 is less than epsilon we return I2 which is the value of integral.
     * */
    const double temp = 15. / 16;

    double integral1;
    double integral2;

    double s0 = function(start_) + function(end_);
    double s1 = 0.0;
    double s2 = 0.0;

    eps_ *= temp;

    for (int i = 1; i < step_; i += 2) {
        s1 += function(start_ + (i * width_));
    }

    for (int i = 0; i < step_; i += 2) {
        s2 += function(start_ + (i * width_));
    }

    integral1 = ((end_ - start_) / (3 * step_)) * (s0 + 4 * s1 + 2 * s2);

    step_ *= 2;
    width_ /= 2;

    s2 += s1;
    s1 = 0;

    for (int i = 1; i < step_; i += 2) {
        s1 += function(start_ + (i * width_));
    }

    integral2 = ((end_ - start_) / (3 * step_)) * (s0 + 4 * s1 + 2 * s2);

    while (fabs(integral2 - integral1) >= eps_) {
        step_ *= 2;
        width_ /= 2;

        integral1 = integral2;
        s2 += s1;
        s1 = 0;

        for (int i = 1; i < step_; i += 2) {
            s1 += function(start_ + (i * width_));
        }

        integral2 = ((end_ - start_) / (3 * step_)) * (s0 + 4 * s1 + 2 * s2);
    }

    return integral2;
}

NumericalMethods::CubeSpline::CubeSpline(short power,
                                         double (*function)(double),
                                         std::vector<double> &nodes,
                                         double alpha0,
                                         double alpha1,
                                         double beta0,
                                         double beta1,
                                         double gamma0,
                                         double gamma1,
                                         double exp_point) :
        alpha0_(alpha0),
        alpha1_(alpha1),
        beta0_(beta0),
        beta1_(beta1),
        gamma0_(gamma0),
        gamma1_(gamma1),
        exp_point_(exp_point),
        nodes_(nodes),
        power_(power - 1) {

    a_.resize(power_ + 1);
    b_.resize(power_ + 1);
    c_.resize(power_ + 1);
    d_.resize(power_ + 1);

    p_.resize(power_);
    q_.resize(power_);
    m_.resize(power_ + 1);

    fvalue_.resize(power_ + 1);

    for (int i = 0; i < power_ + 1; i++) {
        fvalue_[i] = function(nodes_[i]);
    }

}

void NumericalMethods::CubeSpline::sweepMethod() {
    p_[0] = -c_[0] / b_[0];
    q_[0] = d_[0] / b_[0];

    for (int i = 1; i < power_; i++) {
        p_[i] = -c_[i] / (a_[i] * p_[i - 1] + b_[i]);
        q_[i] = (d_[i] - a_[i] * q_[i - 1]) / (a_[i] * p_[i - 1] + b_[i]);
    }

    m_[power_] = (d_[power_] - a_[power_] * q_[power_ - 1]) / (a_[power_] * p_[power_ - 1] + b_[power_]);

    for (int i = power_; i > 0; i--) {
        m_[i - 1] = p_[i - 1] * m_[i] + q_[i - 1];
    }

}

double NumericalMethods::CubeSpline::calculateCubeSpline() {
    double result_ = 0;

    int index = 0;

    a_[0] = 0;
    b_[0] = -alpha0_ * (nodes_[1] - nodes_[0]) / 3 + beta0_;
    c_[0] = alpha0_ * (nodes_[1] - nodes_[0]) / 6;
    d_[0] = gamma0_ - alpha0_ * (fvalue_[1] - fvalue_[0]) / (nodes_[1] - nodes_[0]);

    for (int i = 1; i < power_; i++) {
        a_[i] = (nodes_[i] - nodes_[i - 1]) / 6;
        b_[i] = (nodes_[i + 1] - nodes_[i - 1]) / 3;
        c_[i] = (nodes_[i + 1] - nodes_[i]) / 6;
        d_[i] = (fvalue_[i + 1] - fvalue_[i]) / (nodes_[i + 1] - nodes_[i]) -
                (fvalue_[i] - fvalue_[i - 1]) / (nodes_[i] - nodes_[i - 1]);
    }

    a_[power_] = alpha1_ * (nodes_[power_] - nodes_[power_ - 1]) / 6;
    b_[power_] = alpha1_ * (nodes_[power_] - nodes_[power_ - 1]) / 3 + beta1_;
    c_[power_] = 0;
    d_[power_] = gamma1_ - alpha1_ * (fvalue_[power_] - fvalue_[power_ - 1])
                           / (nodes_[power_] - nodes_[power_ - 1]);

    sweepMethod();

    for (int i = 1; i <= power_; i++) {
        if (exp_point_ >= nodes_[i - 1] && exp_point_ <= nodes_[i]) {
            index = i;
            break;
        }
    }

    result_ = m_[index - 1] * pow(nodes_[index] - exp_point_, 3) / (6 * (nodes_[index] - nodes_[index - 1])) +
              m_[index] * pow(exp_point_ - nodes_[index - 1], 3) / (6 * (nodes_[index] - nodes_[index - 1])) +
              (fvalue_[index - 1] - (m_[index - 1] * pow(nodes_[index] - nodes_[index - 1], 2)) / 6) *
              (nodes_[index] - exp_point_) / (nodes_[index] - nodes_[index - 1]) +
              (fvalue_[index] - (m_[index] * pow(nodes_[index] - nodes_[index - 1], 2)) / 6) *
              (exp_point_ - nodes_[index - 1]) / (nodes_[index] - nodes_[index - 1]);

    return result_;
}

NumericalMethods::MonteCarlo::MonteCarlo(long power) :
        power_(power),
        random_x_(0),
        random_y_(0),
        random_z_(0) {

}

double NumericalMethods::MonteCarlo::calculateMonteCarlo1(double function(double x, double y)) {
    for (int i = 0; i <= power_; i++) {
        random_x_ = randomDouble(0, 1);
        random_y_ = randomDouble(0, 1);

        summary_ += function(random_x_, random_y_);
    }

    return summary_ / power_;
}

double NumericalMethods::MonteCarlo::calculateMonteCarlo2(double (*function)(double, double)) {
    double counter = 0.0;

    for (int i = 0; i <= power_; i++) {
        random_x_ = randomDouble(0, 1);
        random_y_ = randomDouble(0, 1);
        random_z_ = randomDouble(0, 1);

        if (function(random_x_, random_y_) >= random_z_) {
            counter += 1.0;
        }
    }

    return double(counter) / double(power_);
}

NumericalMethods::Runge_Kutta::Runge_Kutta(double start,
                                           double y0,
                                           double h,
                                           double eps)
        : start_(start),
          y0_(y0),
          h_(h),
          eps_(eps) {
}

double NumericalMethods::Runge_Kutta::step(double function(double, double)) {
    double phi0 = h_ * function(start_, y0_);
    double phi1 = h_ * function(start_ + h_ / 2, y0_ + phi0 / 2);
    double phi2 = h_ * function(start_ + h_ / 2, y0_ + phi1 / 2);
    double phi3 = h_ * function(start_ + h_, y0_ + phi2);

    return y0_ + ((phi0 + 2 * phi1 + 2 * phi2 + phi3) / 6);
}


std::vector<double> NumericalMethods::Runge_Kutta::jump(double function(double, double)) {
    double yh = step(function);
    double yh2;
    double y2h2;
    double temp = y0_;

    while (true) {
        h_ /= 2;
        yh2 = step(function);
        start_ += h_;
        y0_ = yh2;

        y2h2 = step(function);

        start_ -= h_;
        h_ *= 2;
        y0_ = temp;

        if (fabs(yh - y2h2) < eps_) {
            return std::vector<double>{start_ + h_, yh, h_};
        } else {
            h_ /= 2;
            yh = yh2;
        }
    }
}

NumericalMethods::ShootingMethod::ShootingMethod(double f(double, double, double),
                                                 double g(double, double, double),
                                                 double F(double, double, double, double, double),
                                                 double start,
                                                 double end,
                                                 double alpha0,
                                                 double betta0,
                                                 double gamma0,
                                                 double alpha1,
                                                 double betta1,
                                                 double gamma1,
                                                 double h,
                                                 double eps,
                                                 double ksi0,
                                                 double ksi1)
        : f_(f),
          g_(g),
          F_(F),
          start_(start),
          end_(end),
          alpha0_(alpha0),
          betta0_(betta0),
          gamma0_(gamma0),
          alpha1_(alpha1),
          betta1_(betta1),
          gamma1_(gamma1),
          h_(h),
          eps_(eps),
          ksi0_(ksi0),
          ksi1_(ksi1) {

}

std::vector<double> NumericalMethods::ShootingMethod::step() {
    double phi_0 = h_ * f_(start_, u0_, v0_);
    double psi_0 = h_ * g_(start_, u0_, v0_);
    double phi_1 = h_ * f_(start_ + h_ / 2, u0_ + phi_0 / 2, v0_ + psi_0 / 2);
    double psi_1 = h_ * g_(start_ + h_ / 2, u0_ + phi_0 / 2, v0_ + psi_0 / 2);
    double phi_2 = h_ * f_(start_ + h_ / 2, u0_ + phi_1 / 2, v0_ + psi_1 / 2);
    double psi_2 = h_ * g_(start_ + h_ / 2, u0_ + phi_1 / 2, v0_ + psi_1 / 2);
    double phi_3 = h_ * f_(start_ + h_, u0_ + phi_2, v0_ + psi_2);
    double psi_3 = h_ * g_(start_ + h_, u0_ + phi_2, v0_ + psi_2);
    double y1 = u0_ + (phi_0 + 2 * phi_1 + 2 * phi_2 + phi_3) / 6;
    double y2 = v0_ + (psi_0 + 2 * psi_1 + 2 * psi_2 + psi_3) / 6;

    return {y1, y2};
}

void NumericalMethods::ShootingMethod::jump() {
    std::vector<double> s1 = step();
    std::vector<double> s_half_(0, 2);
    std::vector<double> s2(0, 2);

    double y1_ = 0;
    double y2_ = 0;

    double temp1 = u0_;
    double temp2 = v0_;

    while (true) {
        h_ /= 2;
        s_half_ = step();

        start_ += h_;

        y1_ = s_half_[0];
        y2_ = s_half_[1];

        u0_ = y1_;
        v0_ = y2_;

        s2 = step();

        start_ -= h_;
        h_ *= 2;
        u0_ = temp1;
        v0_ = temp2;

        if (fabs(s1[0] - s2[0]) < eps_ && fabs(s1[1] - s2[1]) < eps_) {
            break;
        } else {
            h_ /= 2;
            s1 = s_half_;
        }
    }

    start_ += h_;

    double sol1 = s2[0];
    double sol2 = s2[1];

    u0_ = sol1;
    v0_ = sol2;

}


double NumericalMethods::ShootingMethod::shootingSolution() {
    double temp_h = h_;
    double temp_start = start_;
    double temp_end = end_;

    double temp;

    std::vector<double> jump_solution(0, 3);

    std::cout << "?????????????? ksi0: ";
    std::cin >> ksi0_;

    u0_ = ksi0_;
    v0_ = (gamma0_ - alpha0_ * ksi0_) / betta0_;

    while (start_ < end_) {
        jump();

        if (end_ - start_ < h_) {
            h_ = end_ - start_;
        }
    }

    double F0 = F_(u0_, v0_, alpha1_, betta1_, gamma1_);
    std::cout << "F{ksi}: " << F0 << "\n";

    double temp_u0 = u0_;
    double temp_v0 = v0_;

    solution_ = 1;

    double F1;
    double temp_u1;
    double temp_v1;
    double ksi_counter = 0;

    while (true) {
        solution_++;
        start_ = temp_start;
        h_ = temp_h;

        ksi_counter++;
        std::cout << "?????????????? ksi{" << ksi_counter << "}: ";
        std::cin >> ksi1_;

        u0_ = ksi1_;
        v0_ = (gamma0_ - alpha0_ * ksi1_) / betta0_;

        while (start_ < end_) {
            jump();

            if (end_ - start_ < h_) {
                h_ = end_ - start_;
            }
        }

        F1 = F_(u0_, v0_, alpha1_, betta1_, gamma1_);
        std::cout << "F{ksi}: " << F1 << "\n";

        temp_u1 = u0_;
        temp_v1 = v0_;

        if (F0 * F1 < 0) {
            break;
        } else {
            F0 = F1;
            temp_u0 = temp_u1;
            temp_v0 = temp_v1;
            ksi0_ = ksi1_;
        }
    }

    std::cout << "???????????????? F ???????????????? ????????" << "\n";
    std::cout << "\n ???????????? ???????????????????????????? ????????????????" << "\n";

    double F_half;
    double u_half;
    double v_half;

    /// ???????????????????????????? ????????????????
    while (true) {
        shots_fired_++;
        start_ = temp_start;
        h_ = temp_h;

        ksi_half_ = (ksi0_ + ksi1_) * 0.5;

        u0_ = ksi_half_;
        v0_ = (gamma0_ - alpha0_ * ksi_half_) / betta0_;

        while (start_ < end_) {
            jump();

            if (end_ - start_ < h_) {
                h_ = end_ - start_;
            }
        }

        F_half = F_(u0_, v0_, alpha1_, betta1_, gamma1_);
        u_half = u0_;
        v_half = v0_;

        if (F0 * F_half < 0) {
            F1 = F_half;
            temp_u1 = u_half;
            temp_v1 = v_half;
            ksi1_ = ksi_half_;
        } else {
            F0 = F_half;
            temp_u0 = u_half;
            temp_v0 = v_half;
            ksi0_ = ksi_half_;
        }

        if (fabs(F0 - F1) < eps_ * 10 || fabs(F0) < eps_ * 10) {
            break;
        }
    }

    std::cout << "???????????????????? ???????????????? ??????????: " << solution_ + shots_fired_ << "\n";
    std::cout << "???????????????????? ???????????????????????????? ??????????????????: " << shots_fired_ << "\n";
    std::cout << "???????????????? ?????????? F(ksi) = 0: " << F0 << "\n";
    std::cout << "U: " << temp_u0 << " | V: " << temp_v0 << "\n";

    return F0;
}

NumericalMethods::GridMethod::GridMethod(double (*f)(double),
                                         std::vector<double> (*g)(std::vector<double> &, double, double),
                                         int n,
                                         double start,
                                         double end,
                                         double p,
                                         double q,
                                         int lin_space_num,
                                         std::vector<double> alpha,
                                         std::vector<double> betta,
                                         std::vector<double> gamma,
                                         std::vector<double> a,
                                         std::vector<double> b,
                                         std::vector<double> c,
                                         std::vector<double> d) :
        f_(f),
        g_(g),
        start_(start),
        end_(end),
        h_(start + end / n),
        n_(n),
        p_(p),
        q_(q),
        lin_space_num_(lin_space_num),
        alpha_(std::move(alpha)),
        betta_(std::move(betta)),
        gamma_(std::move(gamma)),
        a_(std::move(a)),
        b_(std::move(b)),
        c_(std::move(c)),
        d_(std::move(d)) {

    solution_.resize(n_);
}

std::vector<double> NumericalMethods::GridMethod::getLinSpace() {
    std::vector<double> linspaced;

    double delta = (end_ - start_) / (lin_space_num_) - 1;

    for (int i = 0; i < lin_space_num_ - 1; i++) {
        linspaced.push_back(start_ + delta * i);
    }

    linspaced.push_back(end_);

    return linspaced;
}

std::vector<double> NumericalMethods::GridMethod::solveGrid() {
    d_ = getLinSpace();
    d_ = g_(d_, p_, q_);

    for (int i = 1; i < n_; i++) {
        a_[i] = (1 - (h_ * 0.5) * p_) / (h_ * h_);
        b_[i] = (-2 + (h_ * h_) * q_) / (h_ * h_);
        c_[i] = (1 + (h_ * 0.5) * p_) / (h_ * h_);
    }

    b_[0] = alpha_[0] - betta_[0] / h_;
    c_[0] = betta_[0] / h_;
    d_[0] = gamma_[0];

    b_[n_] = alpha_[1] - betta_[1] / h_;
    c_[n_] = betta_[1] / h_;
    d_[n_] = gamma_[1];

    return sweepMethod();
}

std::vector<double> NumericalMethods::GridMethod::sweepMethod() {
    int size = a_.size();
    std::vector<double> P(size - 1, 0);
    std::vector<double> Q(size - 1, 0);
    std::vector<double> m(a_.size(), 0);

    P[0] = -c_[0] / b_[0];
    Q[0] = d_[0] / b_[0];

    for (int i = 1; i < size - 1; i++) {
        P[i] = -c_[i] / (a_[i] * P[i - 1] + b_[i]);
        Q[i] = (d_[i] - a_[i] * Q[i - 1]) / (a_[i] * P[i - 1] + b_[i]);
    }

    m[size - 1] = (d_[size - 1] - a_[size - 1] * Q[size - 2]) / (a_[size - 1] * P[size - 2] + b_[size - 1]);

    for (int i = size - 2; i > -1; i--) {
        m[i] = P[i] * m[i + 1] + Q[i];
    }

    return m;
}


