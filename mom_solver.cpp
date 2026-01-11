#include "mom_solver.hpp"
#include "impedance_matrix.hpp"
#include "excitation.hpp"
#include <cmath>

std::vector<MoMSolver::Complex>
MoMSolver::buildVoltageVector(int N, const Excitation& ex) const
{
    std::vector<Complex> V(N, {0.0, 0.0});

    V[ex.segment] = {
        ex.voltage * std::cos(ex.phase),
        ex.voltage * std::sin(ex.phase)
    };

    return V;
}

// ------------------------------------------------------------
// Public API
// ------------------------------------------------------------
std::vector<MoMSolver::Complex>
MoMSolver::solveCurrents(
    const ImpedanceMatrix& Z,
    int feedSegment
    ) const
{
    int N = Z.size();
    if (feedSegment < 0 || feedSegment >= N) {
        throw std::runtime_error("Invalid feed segment index");
    }

    // Build dense matrix A and RHS b
    std::vector<std::vector<Complex>> A(N, std::vector<Complex>(N));
    std::vector<Complex> b(N, Complex(0.0, 0.0));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = Z(i,j);
        }
    }

    // Delta-gap voltage source: V = 1V at feed segment
    b[feedSegment] = Complex(1.0, 0.0);

    return gaussianElimination(A, b);
}

MoMSolver::Complex
MoMSolver::inputImpedance(
    const ImpedanceMatrix& Z,
    int feedSegment
    ) const
{
    auto I = solveCurrents(Z, feedSegment);

    // Zin = V / I_feed, with V = 1V
    return Complex(1.0, 0.0) / I[feedSegment];
}

std::vector<MoMSolver::Complex>
MoMSolver::gaussianElimination(
    std::vector<std::vector<Complex>> A,
    std::vector<Complex> b
    ) const
{
    int N = static_cast<int>(b.size());

    // Forward elimination
    for (int k = 0; k < N; ++k) {

        if (std::abs(A[k][k]) < 1e-12) {
            throw std::runtime_error("Singular matrix in MoMSolver");
        }

        for (int i = k + 1; i < N; ++i) {
            Complex factor = A[i][k] / A[k][k];

            for (int j = k; j < N; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    std::vector<Complex> x(N);
    for (int i = N - 1; i >= 0; --i) {
        Complex sum = b[i];
        for (int j = i + 1; j < N; ++j) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }

    return x;
}
