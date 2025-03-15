#ifndef _SLAE_H_
#define _SLAE_H_

#include <cmath>
#include <vector>

#include "../matrix/matrix.h"

/*
SLAE - System of Linear Algebraic Equations
*/

class SLAE {
  public:
    SLAE(const Matrix& matrixOfSystem, const std::vector<double>& bColumn);

    [[nodiscard]]
    std::vector<double> solveSystemLU();

    [[nodiscard]]
    std::vector<double> solveSystemLU(const std::vector<double>& otherBColumn);

    double calculateMatrixAlphaNorm() const;

    [[nodiscard]]
    std::pair<unsigned int, std::vector<double>> solveSystemSimpleIterationsMethod(double eps);

    [[nodiscard]]
    std::pair<unsigned int, std::vector<double>> solveSystemSimpleIterationsMethod(
        const std::vector<double>& otherBColumn,
        double eps
    );

    [[nodiscard]]
    std::pair<unsigned int, std::vector<double>> solveSystemSeidelsMethod(double eps);

    [[nodiscard]]
    std::pair<unsigned int, std::vector<double>> solveSystemSeidelsMethod(
        const std::vector<double>& otherBColumn,
        double eps
    );

  private:
    void bColumnPermutation();

    double calculateCurrentEpsIterations(
      const std::vector<double>& tmpXColumn,
      const std::vector<double>& prevXColumn,
      double matrixAlphaNorm
    ) const;

    double calculateMatrixCNorm() const;
    double calculateCurrentEpsSeidels(
      const std::vector<double>& tmpXColumn,
      const std::vector<double>& prevXColumn,
      double matrixAlphaNorm,
      double matrixCNorm
    ) const;

  private:
    Matrix matrixOfSystem_;
    std::vector<double> bColumn_;
    bool bColumnPermutated_ = false;
};

#endif // _SLAE_H_