#include "./SLAE.h"

SLAE::SLAE(const Matrix& matrixOfSystem, const std::vector<double>& bColumn):
    matrixOfSystem_(matrixOfSystem), bColumn_(bColumn) {}

std::vector<double> SLAE::solveSystemLU() {
    if (!matrixOfSystem_.isLUMatrix()) {
        matrixOfSystem_.LUDecompositionInplace();
    }

    if (!bColumnPermutated_) {
        bColumnPermutation();
    }

    size_t matrixSize = matrixOfSystem_.getN();
    std::vector<double> z(matrixSize);

    for (size_t i = 0; i != matrixSize; ++i) {
        z[i] = bColumn_[i];

        for (size_t j = 0; j != i; ++j) {
            z[i] -= z[j] * matrixOfSystem_[i][j];
        }
    }

    std::vector<double> xVector(matrixSize);
    for (int64_t i = matrixSize - 1; i >= 0; --i) {
        xVector[i] = z[i];

        for (size_t j = i + 1; j != matrixSize; ++j) {
            xVector[i] -= xVector[j] * matrixOfSystem_[i][j];
        }

        xVector[i] /= matrixOfSystem_[i][i];
    }

    return xVector;
}

std::vector<double> SLAE::solveSystemLU(const std::vector<double>& otherBColumn) {
    bColumn_ = otherBColumn;
    bColumnPermutated_ = false;

    return solveSystemLU();
}

std::pair<unsigned int, std::vector<double>> SLAE::solveSystemSimpleIterationsMethod(double eps) {
    size_t matrixSize = matrixOfSystem_.getN();
    double matrixAlphaNorm = calculateMatrixAlphaNorm();
    std::vector<double> xColumn(matrixSize);

    for (size_t i = 0; i != matrixSize; ++i) {
        xColumn[i] = bColumn_[i] / matrixOfSystem_[i][i];
    }

    double currentEps = 1;
    unsigned int iterationsCounter = 0;

    while (eps < currentEps) {
        std::vector<double> tmpXColumn(matrixSize);

        for (size_t i = 0; i != matrixSize; ++i) {
            tmpXColumn[i] = bColumn_[i] / matrixOfSystem_[i][i];
            for (size_t j = 0; j != matrixSize; ++j) {  
                if (i == j) {
                    continue;
                }

                tmpXColumn[i] -= xColumn[j] * matrixOfSystem_[i][j] / matrixOfSystem_[i][i];
            }
        }

        currentEps = calculateCurrentEpsIterations(tmpXColumn, xColumn, matrixAlphaNorm);
        xColumn = std::move(tmpXColumn);

        ++iterationsCounter;
    }

    return {iterationsCounter, xColumn};
}

std::pair<unsigned int, std::vector<double>> SLAE::solveSystemSimpleIterationsMethod(
    const std::vector<double>& otherBColumn,
    double eps)
{
    bColumn_ = otherBColumn;
    bColumnPermutated_ = false;

    return solveSystemSimpleIterationsMethod(eps);
}

std::pair<unsigned int, std::vector<double>> SLAE::solveSystemSeidelsMethod(double eps) {
    size_t matrixSize = matrixOfSystem_.getN();
    double matrixAlphaNorm = calculateMatrixAlphaNorm();
    double matrixCNorm = calculateMatrixCNorm();
    std::vector<double> xColumn(matrixSize);

    for (size_t i = 0; i != matrixSize; ++i) {
        xColumn[i] = bColumn_[i] / matrixOfSystem_[i][i];
    }

    double currentEps = 1;
    unsigned int iterationsCounter = 0;

    while (eps < currentEps) {
        std::vector<double> tmpXColumn(matrixSize);

        for (size_t i = 0; i != matrixSize; ++i) {
            tmpXColumn[i] = bColumn_[i] / matrixOfSystem_[i][i];
            for (size_t j = 0; j != matrixSize; ++j) {  
                if (i == j) {
                    continue;
                }

                double x_k = (j < i) ? tmpXColumn[j] : xColumn[j];
                tmpXColumn[i] -= x_k * matrixOfSystem_[i][j] / matrixOfSystem_[i][i];
            }
        }

        currentEps = calculateCurrentEpsSeidels(
            tmpXColumn, 
            xColumn, 
            matrixAlphaNorm, 
            matrixCNorm
        );
        xColumn = std::move(tmpXColumn);

        ++iterationsCounter;
    }

    return {iterationsCounter, xColumn};
}

std::pair<unsigned int, std::vector<double>> SLAE::solveSystemSeidelsMethod(
    const std::vector<double>& otherBColumn,
    double eps)
{
    bColumn_ = otherBColumn;
    bColumnPermutated_ = false;

    return solveSystemSeidelsMethod(eps);
}

void SLAE::bColumnPermutation() {
    size_t bSize = bColumn_.size();
    std::vector<double> newB(bSize);
    const std::vector<unsigned int>& rowPermutations = matrixOfSystem_.getRowPermutations();

    for (size_t i = 0; i != bSize; ++i) {
        newB[i] = bColumn_[rowPermutations[i]];
    }

    bColumn_ = std::move(newB);
    bColumnPermutated_ = true;
}

double SLAE::calculateMatrixAlphaNorm() const {
    size_t matrixSize = matrixOfSystem_.getN();
    double matrixNorm = 0;

    for (size_t i = 0; i != matrixSize; ++i) {
        for (size_t j = 0; j != matrixSize; ++j) {
            if (i == j) {
                continue;
            }

            double alpha_ij = -matrixOfSystem_[i][j] / matrixOfSystem_[i][i]; 
            matrixNorm += alpha_ij * alpha_ij;
        }
    }

    return std::sqrt(matrixNorm);
}

double SLAE::calculateCurrentEpsIterations(
    const std::vector<double>& tmpXColumn,
    const std::vector<double>& prevXColumn,
    double matrixAlphaNorm) const
{
    size_t xColumnSize = tmpXColumn.size();
    double vectorNorm = 0;

    for (size_t i = 0; i != xColumnSize; ++i) {
        vectorNorm += (tmpXColumn[i] - prevXColumn[i]) * (tmpXColumn[i] - prevXColumn[i]);
    }

    vectorNorm = std::sqrt(vectorNorm);
    double currentEps = vectorNorm * matrixAlphaNorm / (1 - matrixAlphaNorm);

    return currentEps;
}

double SLAE::calculateMatrixCNorm() const {
    size_t matrixSize = matrixOfSystem_.getN();
    double matrixNorm = 0;

    for (size_t i = 0; i != matrixSize; ++i) {
        for (size_t j = i + 1; j != matrixSize; ++j) {
            double alpha_ij = -matrixOfSystem_[i][j] / matrixOfSystem_[i][i];
            matrixNorm += alpha_ij * alpha_ij;
        }
    }

    return std::sqrt(matrixNorm);
}

double SLAE::calculateCurrentEpsSeidels(
    const std::vector<double>& tmpXColumn,
    const std::vector<double>& prevXColumn,
    double matrixAlphaNorm,
    double matrixCNorm) const
{
    size_t xColumnSize = tmpXColumn.size();
    double vectorNorm = 0;

    for (size_t i = 0; i != xColumnSize; ++i) {
        vectorNorm += (tmpXColumn[i] - prevXColumn[i]) * (tmpXColumn[i] - prevXColumn[i]);
    }

    vectorNorm = std::sqrt(vectorNorm);
    double currentEps = vectorNorm * matrixCNorm / (1 - matrixAlphaNorm);

    return currentEps;
}
