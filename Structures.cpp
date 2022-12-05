//
// Created by Данил  Казаев on 15.11.2022.
//

#include "Structures.h"

Structures::Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
    matrix_.resize(rows_ * cols_, 0);
}

Structures::Matrix &Structures::Matrix::operator*(const Structures::Matrix & m) {
    Matrix temp(rows_ * cols_, 0);

    for (int i = 0; i < temp.rows_; i++) {
        for (int j = 0; j < temp.cols_; j++) {
            for (int k = 0; k < temp.rows_; k++) {
            }
        }
    }
    return *this;
}
