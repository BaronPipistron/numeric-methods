#include "../../include/TDM.h"

TDM::TDM(size_t n): n_(n) {
    upperDiagonal_.resize(n - 1, 0.);
    mainDiagonal_.resize(n, 0.);
    lowerDiagonal_.resize(n, 0.);
}

double TDM::getElement(size_t i, size_t j) const noexcept {
    if (i == j) {
        return mainDiagonal_[i];
    } else if (i == j + 1) {
        return lowerDiagonal_[j];
    } else if (i + 1 == j) {
        return upperDiagonal_[i];
    }

    return 0.;
}

bool TDM::setElement(double value, size_t i, size_t j) noexcept {
    if (i == j) {
        mainDiagonal_[i] = value;
        return true;
    } else if (i == j + 1) {
        lowerDiagonal_[j] = value;
        return true;
    } else if (i + 1 == j) {
        upperDiagonal_[i] = value;
        return true;
    }

    return false;
}

size_t TDM::getN() const noexcept {
    return n_;
}

const std::vector<double>& TDM::getUpperDiagonal() const noexcept {
    return upperDiagonal_;
}

const std::vector<double>& TDM::getMainDiagonal() const noexcept {
    return mainDiagonal_;
}

const std::vector<double>& TDM::getLowerDiagonal() const noexcept {
    return lowerDiagonal_;
}

std::istream& operator>>(std::istream& is, TDM& tdm) {
    for (size_t i = 0; i != tdm.n_; ++i) {
        for (size_t j = 0; j != tdm.n_; ++j) {
            double num;
            is >> num;
            
            if (i == j) {
                tdm.mainDiagonal_[i] = num;
            } else if (i == j + 1) {
                tdm.lowerDiagonal_[j] = num;
            } else if (i + 1 == j) {
                tdm.upperDiagonal_[i] = num;
            }
        }
    }

    return is;
}

std::ostream& operator<<(std::ostream& os, const TDM& tdm) {
    for (size_t i = 0; i != tdm.n_; ++i) {
        for (size_t j = 0; j != tdm.n_; ++j) {
            if (i == j) {
                os << tdm.mainDiagonal_[i];
            } else if (i == j + 1) {
                os << tdm.lowerDiagonal_[j];
            } else if (i + 1 == j) {
                os << tdm.upperDiagonal_[i];
            } else {
                os << 0;
            }
            os << ' ';
        }
    }

    return os;
}