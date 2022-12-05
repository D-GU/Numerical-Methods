//
// Created by Данил  Казаев on 15.11.2022.
//

#ifndef NUMERICAL_METHODS_STRUCTURES_H
#define NUMERICAL_METHODS_STRUCTURES_H

#include <vector>
#include <string>
#include <unordered_map>

namespace Structures {
    class Matrix {
    public:
        Matrix(int rows, int cols);

        ~Matrix() = default;

    public:
        Matrix &operator*(const Matrix & m);

        Matrix &operator+(const Matrix &);

    private:
        int rows_;
        int cols_;

        std::vector<double> matrix_;
    };
}


#endif //NUMERICAL_METHODS_STRUCTURES_H
