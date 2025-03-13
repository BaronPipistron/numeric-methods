#ifndef _SLAE_H_
#define _SLAE_H_

#include <vector>

#include "./matrix.h"

/*
SLAE - System of Linear Algebraic Equations
*/

class SLAE {
  public:
    SLAE(const Matrix& matrixOfSystem, const std::vector<double>& bColumn);

    [[nodiscard]]
    std::vector<double> solveSystem();

    [[nodiscard]]
    std::vector<double> solveSystem(const std::vector<double>& otherBColumn);

  private:
    void bColumnPermutation();

  private:
    Matrix matrixOfSystem_;
    std::vector<double> bColumn_;
    bool bColumnPermutated_ = false;
};

#endif // _SLAE_H_