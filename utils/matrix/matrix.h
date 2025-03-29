#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cmath>
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

class Matrix {
  public: // methods
    Matrix(size_t n, double fillValue = 0.);
    Matrix(size_t n, size_t m, double fillValue = 0.);
    Matrix(const std::initializer_list<std::vector<double>>& list);

    ~Matrix() = default;

    std::vector<double>& operator[](size_t index) noexcept;
    const std::vector<double>& operator[](size_t index) const noexcept;

    [[nodiscard]]
    size_t getN() const noexcept;

    [[nodiscard]]
    size_t getM() const noexcept;

    [[nodiscard]]
    const std::vector<unsigned int>& getRowPermutations() const noexcept;

    [[nodiscard]]
    bool isLUMatrix() const noexcept;

    [[nodiscard]]
    std::pair<unsigned int, std::vector<unsigned int>> getPermutationsInfo() const noexcept;

    [[nodiscard]]
    Matrix operator*(const Matrix& rhs) const;

    [[nodiscard]]
    bool operator==(const Matrix& rhs) const noexcept;

    double determinant() const;

    [[nodiscard]]
    Matrix LUDecomposition() const;

    void LUDecompositionInplace();

  public: // exceptions
    struct NonConformable { };
    struct NotSquareMatrix { };
    struct NotSuitableMatrixForThisOperation { };

  private: // methods
    void initializeRowPermutationsVector();

  private: // fields
    size_t n_;
    size_t m_;
    std::vector<std::vector<double>> matrix_;

    unsigned int permutationsCounter_ = 0;
    std::vector<unsigned int> rowPermutations_;

    bool isLUMatrix_ = false;
};

#endif // _MATRIX_H_