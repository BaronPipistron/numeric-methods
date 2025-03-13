#include "../../include/SLAE.h"

SLAE::SLAE(const Matrix& matrixOfSystem, const std::vector<double>& bColumn):
    matrixOfSystem_(matrixOfSystem), bColumn_(bColumn) {}

std::vector<double> SLAE::solveSystem() {
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

std::vector<double> SLAE::solveSystem(const std::vector<double>& otherBColumn) {
    bColumn_ = otherBColumn;
    bColumnPermutated_ = false;

    return solveSystem();
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