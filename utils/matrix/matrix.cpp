#include "./matrix.h"

Matrix::Matrix(size_t n, double fillValue): n_(n), m_(n) {
    matrix_ = std::vector<std::vector<double>>(n, std::vector<double>(n, fillValue));
}

Matrix::Matrix(size_t n, size_t m, double fillValue): n_(n), m_(m) {
    matrix_ = std::vector<std::vector<double>>(m, std::vector<double>(n, fillValue));
}

Matrix::Matrix(const std::initializer_list<std::vector<double>>& list) {
    n_ = list.begin()->size();
    m_ = list.size();

    for (const auto& row: list) {
        matrix_.push_back(row);
    }
}

std::vector<double>& Matrix::operator[](size_t index) noexcept {
    return matrix_[index];
}

const std::vector<double>& Matrix::operator[](size_t index) const noexcept {
    return matrix_[index];
}

size_t Matrix::getN() const noexcept {
    return n_;
}

size_t Matrix::getM() const noexcept {
    return m_;
}

const std::vector<unsigned int>& Matrix::getRowPermutations() const noexcept {
    return rowPermutations_;
}

bool Matrix::isLUMatrix() const noexcept {
    return isLUMatrix_;
}

std::pair<unsigned int, std::vector<unsigned int>> Matrix::getPermutationsInfo() const noexcept {
    return { permutationsCounter_, rowPermutations_ };
}

Matrix Matrix::operator*(const Matrix& rhs) const {
    if (n_ != rhs.m_) {
        throw NonConformable();
    }

    Matrix result(rhs.n_, m_);

    for (size_t row = 0; row != m_; ++row) {
        for (size_t column = 0; column != rhs.n_; ++column) {
            for (size_t columnElem = 0; columnElem != rhs.m_; ++columnElem) {
                result[row][column] += matrix_[row][columnElem] * rhs[columnElem][column];
            }
        }
    }

    return result;
}

bool Matrix::operator==(const Matrix& rhs) const noexcept {
    if (n_ != rhs.n_ || m_ != rhs.m_) return false;

    for (size_t row = 0; row != m_; ++row) {
        for (size_t column = 0; column != n_; ++column) {
            if (std::abs(matrix_[row][column] - rhs[row][column]) > std::numeric_limits<float>::epsilon()) {
                return false;
            }
        }
    }

    return true;
}

double Matrix::determinant() const {
    if (n_ != m_) {
        throw NotSquareMatrix();
    } else if (!isLUMatrix_) {
        throw NotSuitableMatrixForThisOperation();
    }

    double det = (permutationsCounter_ % 2 == 0) ? 1 : -1; 
    for (size_t i = 0; i != n_; ++i) {
        det *= matrix_[i][i];
    }

    return det;
}

Matrix Matrix::LUDecomposition() const {
    Matrix LUMatrix = *this;
    LUMatrix.LUDecompositionInplace();

    return LUMatrix;
}

void Matrix::LUDecompositionInplace() {
    if (n_ != m_) {
        throw NotSquareMatrix();
    }

    initializeRowPermutationsVector();

    for (size_t column = 0; column != n_ - 1; ++column) {
        size_t pivotElementIndex = column;

        for (size_t columnElement = column + 1; columnElement != n_; ++columnElement) {
            if (std::abs(matrix_[columnElement][column]) > std::abs(matrix_[pivotElementIndex][column])) {
                pivotElementIndex = columnElement;
            }
        }

        if (pivotElementIndex != column) {
            ++permutationsCounter_;
            std::swap(matrix_[pivotElementIndex], matrix_[column]);
            std::swap(rowPermutations_[pivotElementIndex], rowPermutations_[column]);
        }

        for (size_t element = column + 1; element != n_; ++element) {
            double L = matrix_[element][column] / matrix_[column][column];

            for (size_t num = column; num != n_; ++num) {
                matrix_[element][num] = matrix_[element][num] - L * matrix_[column][num];
            }

            matrix_[element][column] = L;
        }
    }

    isLUMatrix_ = true;
}

void Matrix::initializeRowPermutationsVector() {
    rowPermutations_.resize(n_);

    for (size_t i = 0; i != n_; ++i) {
        rowPermutations_[i] = i;
    }
}