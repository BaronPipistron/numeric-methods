#include <iostream>

#include "../../utils/matrix/matrix.h"
#include "../../utils/slae/SLAE.h"

[[nodiscard]]
Matrix findInverseMatrix(const Matrix& LUMatrix) {
    size_t matrixSize = LUMatrix.getN();
    Matrix inverseMatrix(matrixSize);
    std::vector<double> unitMatrixColumn(matrixSize, 0);
    SLAE slae(LUMatrix, unitMatrixColumn);

    for (size_t i = 0; i != matrixSize; ++i) {
        unitMatrixColumn[i] = 1;

        std::vector<double> inverseMatrixColumn = slae.solveSystemLU(unitMatrixColumn);

        for (size_t j = 0; j != matrixSize; ++j) {
            inverseMatrix[j][i] = inverseMatrixColumn[j];
        }

        unitMatrixColumn[i] = 0;
    }

    return inverseMatrix;
}

bool solutionCheck(
    const Matrix& inputMatrix, 
    const std::vector<double>& bColumn, 
    const std::vector<double>& xColumn) 
{
    size_t matrixSize = inputMatrix.getN();

    for (size_t row = 0; row != matrixSize; ++row) {
        double getNum = 0;
        for (size_t column = 0; column != matrixSize; ++column) {
            getNum += inputMatrix[row][column] * xColumn[column];
        }

        if (std::abs(getNum - bColumn[row]) > std::numeric_limits<float>::epsilon()) {
            return false;
        }
    }

    return true;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (size_t i = 0; i != m.getM(); ++i) {
        for (size_t j = 0; j != m.getN(); ++j) {
            os << m[i][j] << ' ';
        }
        os << '\n';
    }

    return os;
}

int main() {
    size_t matrixSize;
    std::cout << "Input matrix size (one positive number): ";
    std::cin >> matrixSize;

    Matrix matrix(matrixSize);
    std::cout << "Input system matrix:\n";
    for (size_t i = 0; i != matrixSize; ++i) {
        for (size_t j = 0; j != matrixSize; ++j) {
            std::cin >> matrix[i][j];
        }
    }
    Matrix toCheckMatrix = matrix;

    std::vector<double> b(matrixSize);
    std::cout << "Input b column:\n";
    for (size_t i = 0; i != matrixSize; ++i) {
        std::cin >> b[i];
    }
    std::cout << '\n' << std::endl;

    matrix.LUDecompositionInplace();
    SLAE slae(matrix, b);

    std::vector<double> xVector = slae.solveSystemLU();
    Matrix inverseMatrix = findInverseMatrix(matrix);

    std::cout << "System Matrix: " << std::endl;
    std::cout << toCheckMatrix << std::endl;

    std::cout << "LU matrix:\n" << matrix << std::endl;
    std::cout << "Determinant: " << matrix.determinant() << std::endl;

    std::cout << "Answer: " << std::endl;
    for (size_t i = 0; i != xVector.size(); ++i) { 
        std::cout << "x_" << i + 1 << " = " << xVector[i] << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Inverse Matrix:\n" << inverseMatrix << std::endl;

    if (solutionCheck(toCheckMatrix, b, xVector)) {
        std::cout << "Solution find correctly" << std::endl;
    } else {
        std::cout << "ERROR: Solution aren't correctly" << std::endl;
    }

    Matrix identityMatrix(toCheckMatrix.getN());
    for (size_t i = 0; i != toCheckMatrix.getN(); ++i) {
        identityMatrix[i][i] = 1;
    }

    if (toCheckMatrix * inverseMatrix == identityMatrix) {
        std::cout << "Inverse Matrix find correctly" << std::endl;
    } else {
        std::cout << "ERROR: Inverse Matrix aren't correctly" << std::endl;
    }

    return 0;
}