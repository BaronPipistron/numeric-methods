#ifndef _TDM_H_
#define _TDM_H_

#include <iostream>
#include <vector>

/*
TDM - Tridiagonal Matrix
*/

class TDM {
  public:
    TDM(size_t n);

    double getElement(size_t i, size_t j) const noexcept;
    bool setElement(double value, size_t i, size_t j) noexcept;

    [[nodiscard]]
    size_t getN() const noexcept;

    [[nodiscard]]
    const std::vector<double>& getUpperDiagonal() const noexcept;

    [[nodiscard]]
    const std::vector<double>& getMainDiagonal() const noexcept;

    [[nodiscard]]
    const std::vector<double>& getLowerDiagonal() const noexcept;

    friend std::istream& operator>>(std::istream& is, TDM& tdm);
    friend std::ostream& operator<<(std::ostream& os, const TDM& tdm);

  private:
    size_t n_;

    std::vector<double> upperDiagonal_;
    std::vector<double> mainDiagonal_;
    std::vector<double> lowerDiagonal_;
};

#endif // _TDM_H_S