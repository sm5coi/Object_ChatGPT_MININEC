#ifndef MOM_SOLVER_HPP
#define MOM_SOLVER_HPP

#include <vector>
#include <complex>
#include "impedance_matrix.hpp"

class MoMSolver {
public:
    using Complex = std::complex<double>;

    // Lös Z * I = V
    std::vector<Complex>
    solveCurrents(
        const ImpedanceMatrix& Z,
        int feedSegment
        ) const;

    // Beräkna ingångsimpedans
    Complex
    inputImpedance(
        const ImpedanceMatrix& Z,
        int feedSegment
        ) const;

private:
    std::vector<Complex>
    gaussianElimination(
        std::vector<std::vector<Complex>> A,
        std::vector<Complex> b
        ) const;
};

#endif
