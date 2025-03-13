#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "./include/TDM.h"

bool checkApplicabilityTMA(const TDM& tdm) {
    size_t matrixSize = tdm.getN();

    for (size_t i = 1; i != matrixSize - 1; ++i) {
        if (std::abs(tdm.getElement(i, i)) < std::abs(tdm.getElement(i, i - 1)) + std::abs(tdm.getElement(i, i + 1))) {
            return false;
        }
    }

    if (std::abs(tdm.getElement(0, 0)) < std::abs(tdm.getElement(0, 1))) {
        return false;
    } else if (std::abs(tdm.getElement(matrixSize - 1, matrixSize - 1)) < \
               std::abs(tdm.getElement(matrixSize - 1, matrixSize - 2))) {
        return false;
    }

    return true;
}

std::vector<double> tridiagonalMatrixAlgorithm(
    const TDM& tdm,
    const std::vector<double>& bColumn) 
{
    size_t matrixSize = tdm.getN();
    std::vector<double> xColumn(matrixSize);
    std::vector<double> P(matrixSize);
    std::vector<double> Q(matrixSize);

    P[0] = -tdm.getElement(0, 1) / tdm.getElement(0, 0);
    Q[0] = bColumn[0] / tdm.getElement(0, 0);

    for (size_t i = 1; i != matrixSize - 1; ++i) {
        P[i] = -tdm.getElement(i, i + 1) / (tdm.getElement(i, i - 1) * P[i - 1] + tdm.getElement(i, i));
        Q[i] = (bColumn[i] - tdm.getElement(i, i - 1) * Q[i - 1]) / \
               (tdm.getElement(i, i - 1) * P[i - 1] + tdm.getElement(i, i));
    }

    xColumn[matrixSize - 1] = (bColumn[matrixSize - 1] - tdm.getElement(matrixSize - 1, matrixSize - 2) * Q[matrixSize - 2]) / \
                              (tdm.getElement(matrixSize - 1, matrixSize - 2) * P[matrixSize - 2] + \
                               tdm.getElement(matrixSize - 1, matrixSize - 1));
    for (int i = matrixSize - 2; i >= 0; --i) {
        xColumn[i] = P[i] * xColumn[i + 1] + Q[i];
    }

    return xColumn;
}

bool solutionCheck(
    const TDM& tdm,
    const std::vector<double>& bColumn,
    const std::vector<double>& xColumn
) 
{
    size_t matrixSize = tdm.getN();

    double firstEquation = tdm.getElement(0, 0) * xColumn[0] + tdm.getElement(0, 1) * xColumn[1];
    
    if (std::abs(firstEquation - bColumn[0]) > std::numeric_limits<float>::epsilon()) {
        return false;
    }

    for (size_t i = 1; i != matrixSize - 1; ++i) {
        double getNum = tdm.getElement(i, i - 1) * xColumn[i - 1] + \
                        tdm.getElement(i, i) * xColumn[i] + \
                        tdm.getElement(i, i + 1) * xColumn[i + 1];

        if (std::abs(getNum - bColumn[i]) > std::numeric_limits<float>::epsilon()) {
            return false;
        }
    }

    double lastEquation = tdm.getElement(matrixSize - 1, matrixSize - 2) * xColumn[matrixSize - 2] + \
                          tdm.getElement(matrixSize - 1, matrixSize - 1) * xColumn[matrixSize - 1];

    if (std::abs(lastEquation - bColumn[matrixSize - 1]) > std::numeric_limits<float>::epsilon()) {
        return false;
    }
                                       
    return true;
}

int main() {
    size_t matrixSize;
    std::cout << "Input Matrix Size: ";
    std::cin >> matrixSize;

    TDM tdm(matrixSize);
    std::cout << "Input Matrix:\n";
    std::cin >> tdm;

    std::vector<double> bColumn(matrixSize);
    std::cout << "Input b column:\n";
    for (size_t i = 0; i != matrixSize; ++i) {
        std::cin >> bColumn[i];
    }

    std::cout << '\n';
    if (!checkApplicabilityTMA(tdm)) {
        std::cout << "ERROR: Tridiagonal Matrix Algorithm not applicable" << std::endl;
        
        return 0;
    } else {
        std::cout << "Tridiagonal Matrix Algorithm applicable" <<std::endl;
    }

    std::vector<double> xColumn = tridiagonalMatrixAlgorithm(tdm, bColumn);

    std::cout << "Answer:" << std::endl;
    for (size_t i = 0; i != matrixSize; ++i) {
        std::cout << "x_" << i + 1 << " = " << xColumn[i] << std::endl;
    }
    std::cout << std::endl;

    if (solutionCheck(tdm, bColumn, xColumn)) {
        std::cout << "Solution find correctly" << std::endl;
    } else {
        std::cout << "ERROR: Solution aren't correctly" << std::endl;
    }

    return 0;
}