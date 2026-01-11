#ifndef MOM_SOLVER_HPP
#define MOM_SOLVER_HPP

#include <vector>
#include <complex>
#include "impedance_matrix.hpp"

class ImpedanceMatrix;
struct Excitation;

class MoMSolver {
public:
    using Complex = std::complex<double>;

    // Bygg V-vektorn (MININEC-style EX)
    std::vector<Complex>
    buildVoltageVector(int N, const Excitation& ex) const;

    // Lös Z·I = V
    std::vector<Complex>
    solveCurrents(
        const ImpedanceMatrix& Z,
        const std::vector<Complex>& V
        ) const;

    std::vector<MoMSolver::Complex>
    gaussianElimination(
        std::vector<std::vector<Complex>> A,
        std::vector<Complex> b
        ) const;

    MoMSolver::Complex
    inputImpedance(
        const ImpedanceMatrix& Z,
        int feedSegment
        ) const;

    std::vector<MoMSolver::Complex>
    solveCurrents(
        const ImpedanceMatrix& Z,
        int feedSegment
        ) const;
};


#endif
