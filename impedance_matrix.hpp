#ifndef IMPEDANCE_MATRIX_HPP
#define IMPEDANCE_MATRIX_HPP

#include <vector>
#include <complex>

class ImpedanceMatrix {
public:
    using Complex = std::complex<double>;

    explicit ImpedanceMatrix(int n = 0)
        : n_(n), data_(n*n, Complex(0.0, 0.0)) {}

    int size() const { return n_; }

    Complex& operator()(int i, int j) {
        return data_[i*n_ + j];
    }

    const Complex& operator()(int i, int j) const {
        return data_[i*n_ + j];
    }

private:
    int n_;
    std::vector<Complex> data_;
};

#endif
