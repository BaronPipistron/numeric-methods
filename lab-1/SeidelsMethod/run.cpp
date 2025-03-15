#include <iostream>

#include "../../utils/matrix/matrix.h"
#include "../../utils/slae/SLAE.h"

bool solutionCheck(
    const Matrix& inputMatrix, 
    const std::vector<double>& bColumn, 
    const std::vector<double>& xColumn,
    double eps) 
{
    size_t matrixSize = inputMatrix.getN();

    for (size_t row = 0; row != matrixSize; ++row) {
        double getNum = 0;
        for (size_t column = 0; column != matrixSize; ++column) {
            getNum += inputMatrix[row][column] * xColumn[column];
        }

        if (std::abs(getNum - bColumn[row]) > eps) {
            return false;
        }
    }

    return true;
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

    std::vector<double> bColumn(matrixSize);
    std::cout << "Input b column:\n";
    for (size_t i = 0; i != matrixSize; ++i) {
        std::cin >> bColumn[i];
    }
    std::cout << '\n';

    double eps;
    std::cout << "Input error rate: ";
    std::cin >> eps;

    SLAE slae(matrix, bColumn);

    double matrixNorm = slae.calculateMatrixAlphaNorm();
    if (matrixNorm >= 1) {
        std::cout << "ERROR: Seidels Method not applicable - matrix norm is " \
                  << matrixNorm << std::endl;
        
        return 0;
    } else {
        std::cout << "Seidels Method applicable" <<std::endl;
    }

    std::pair<unsigned int, std::vector<double>> seidelsMethodResult = slae.solveSystemSeidelsMethod(eps);
    unsigned int amountOfIterations = seidelsMethodResult.first;
    std::vector<double> xColumn = std::move(seidelsMethodResult.second);

    std::cout << "\nMatrix Norm: " << matrixNorm << std::endl;
    std::cout << "\nAmount of iterations: " << amountOfIterations << std::endl;

    std::cout << "\nAnswer:" << std::endl;
    for (size_t i = 0; i != matrixSize; ++i) {
        std::cout << "x_" << i + 1 << " = " << xColumn[i] << std::endl;
    }
    std::cout << std::endl;

    if (solutionCheck(matrix, bColumn, xColumn, eps)) {
        std::cout << "Solution find correctly" << std::endl;
    } else {
        std::cout << "ERROR: Solution aren't correctly" << std::endl;
    }

    return 0;
}