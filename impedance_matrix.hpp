#ifndef IMPEDANCE_MATRIX_HPP
#define IMPEDANCE_MATRIX_HPP

#include <complex>
#include <vector>
//#include <stdexcept>

class ImpedanceMatrix {
public:
    using Complex = std::complex<double>;

    explicit ImpedanceMatrix(int n)
        : n_(n), data_(n * n, Complex(0.0, 0.0)) {}

    int size() const {
        return n_;
    }

    Complex& operator()(int i, int j) {
        return data_.at(i * n_ + j);
    }

    const Complex& operator()(int i, int j) const {
        return data_.at(i * n_ + j);
    }

private:
    int n_;
    std::vector<Complex> data_;
};

#endif
