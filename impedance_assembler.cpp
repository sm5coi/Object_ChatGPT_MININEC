#include "impedance_assembler.hpp"
#include "impedance_matrix.hpp"
#include "mininec_kernel.hpp"
#include "geometry.hpp"
#include "Vec3.hpp"
#include <iostream>

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
ImpedanceAssembler::ImpedanceAssembler(const MininecKernel& kernel)
    : kernel_(kernel)
{
}

ImpedanceMatrix
ImpedanceAssembler::build(const Geometry& geom) const
{
    const int N = static_cast<int>(geom.segments.size());
    ImpedanceMatrix Z(N);

    for (int m = 0; m < N; ++m)
    {
        // diagonal
        Z(m,m) = kernel_.selfImpedanceAnalytic(geom, m);

        for (int n = m+1; n < N; ++n)
        {
            Complex Zmn(0.0,0.0);

            const auto& segM = geom.segments[m];
            const auto& segN = geom.segments[n];

            if (segM.wire != segN.wire)
            {
                Zmn = Complex(0.0,0.0);
            }
            else
            {
                // --- vector ---
                Complex psiV =
                    0.5 * (
                        kernel_.vectorPotentialContributionMININEC(geom, m, n, +0.5) +
                        kernel_.vectorPotentialContributionMININEC(geom, m, n, -0.5)
                        );

                Vec3 Tm = segM.unitDir(geom);
                Vec3 An = segN.unitDir(geom);
                double dotTA = Tm.dot(An);

                Complex Zvec = dotTA * psiV;

                // --- gradPhi ---
                Complex Zgrad(0.0,0.0);
                if (!geom.areAdjacentOnSameWire(m, n))
                    Zgrad = kernel_.gradPhiContributionMININEC(geom, m, n);

                Zmn = Zvec + Zgrad;
            }

            // force reciprocity
            Z(m,n) = Zmn;
            Z(n,m) = Zmn;
        }
    }

    return Z;
}

/*
ImpedanceMatrix
ImpedanceAssembler::build(const Geometry& geom) const
{
    int N = static_cast<int>(geom.segments.size());
    ImpedanceMatrix Z(N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // 1) SELF först och alltid
            if (i == j)
            {
                Z(i,j) = kernel_.selfImpedanceAnalytic(geom, i);
                continue;
            }

            const auto& segI = geom.segments[i];
            const auto& segJ = geom.segments[j];

            Complex psiS(0.0, 0.0);
            Complex psiV(0.0, 0.0);
            Complex grad(0.0, 0.0);

            // 2) Same-wire mutual (i != j)
            if (segI.wire == segJ.wire)
            {
                // psiS = kernel_.scalarPotentialMININEC(geom, i, j, +0.5, +1);
                // psiV = kernel_.vectorPotentialMININEC(geom, i, j, +0.5, +1);
                // grad = 0
                double obsOffset = (i >= j) ? +0.5 : -0.5;
                int    srcMode   = (i >= j) ? +1   : -1;

                psiS = kernel_.scalarPotentialMININEC(geom, i, j, obsOffset, srcMode);
                psiV = kernel_.vectorPotentialMININEC(geom, i, j, obsOffset, srcMode);

            }
            else
            {
                psiS = kernel_.scalarPotential(geom, i, j);
                psiV = kernel_.vectorPotential(geom, i, j);
            }

            if (!geom.areAdjacentOnSameWire(i, j)) {
                    grad = kernel_.gradPhiContributionMININEC(geom, i, j);
            }

            if (i == 5 && j == 6) {
                std::cout << "gradPhi(5,6) = " << grad << "\n";
            }

            Vec3 Ti = segI.unitDir(geom);
            Vec3 Aj = segJ.unitDir(geom);
            double dotTA = Ti.dot(Aj);

            Complex vecTerm = dotTA * psiV;

            if (i == 0 && j == 0) {
                std::cout << "Z(" << i << "," << j << ")\n";
                std::cout << "  psiS = " << psiS << "\n";
                std::cout << "  psiV = " << psiV << "\n";
                std::cout << "  grad = " << grad << "\n";
            }

            if (i == 2 && j == 8) {
                std::cout << "Z(" << i << "," << j << ")\n";
                std::cout << "  psiS = " << psiS << "\n";
                std::cout << "  psiV = " << psiV << "\n";
                std::cout << "  grad = " << grad << "\n";
            }

            if (i == 5 && j == 6) {
                std::cout << "Z(" << i << "," << j << ")\n";
                std::cout << "  psiS = " << psiS << "\n";
                std::cout << "  psiV = " << psiV << "\n";
                std::cout << "  grad = " << grad << "\n";
            }

            Z(i,j) = psiS + vecTerm + grad;
        }
    }
    return Z;
}
*/

std::complex<double>
ImpedanceAssembler::gradPhiContributionMININEC(
    const Geometry& geom,
    int I,
    int J
    ) const
{
    using Cx = std::complex<double>;

    const double LJ = geom.segmentLength(J);
    if (LJ < 1e-12) return Cx(0.0, 0.0);

    // ψ(M+½, N, N+1)
    Cx psi_p_np = kernel_.scalarPotentialMININEC(geom, I, J, +0.5, +1);

    // ψ(M−½, N, N+1)
    Cx psi_m_np = kernel_.scalarPotentialMININEC(geom, I, J, -0.5, +1);

    // ψ(M+½, N−1, N)
    Cx psi_p_nm = kernel_.scalarPotentialMININEC(geom, I, J, +0.5, -1);

    // ψ(M−½, N−1, N)
    Cx psi_m_nm = kernel_.scalarPotentialMININEC(geom, I, J, -0.5, -1);

    // gradPhi ≈ [ (ψp-ψm)/L + (ψp2-ψm2)/L ]
    Cx grad = (psi_p_np - psi_m_np) / LJ
              + (psi_p_nm - psi_m_nm) / LJ;

    // Teckenfaktorer F4/F5 kommer senare (koppling/strömorientering)
    return grad;
}
