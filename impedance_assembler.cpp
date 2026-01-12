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

/*
std::complex<double>
ImpedanceAssembler::gradPhiContribution(
    const Geometry& geom,
    int obsSeg,
    int srcSeg
    ) const
{
    // Första version:
    // gradPhi approximeras som noll för olika segment
    // (self-term hanteras redan i kernel)
    if (obsSeg != srcSeg) {
        return {0.0, 0.0};
    }

    // Självterm – liten korrigering, ofta nära noll
    return {0.0, 0.0};
}
*/
ImpedanceMatrix
ImpedanceAssembler::build(const Geometry& geom) const
{
    int N = static_cast<int>(geom.segments.size());
    ImpedanceMatrix Z(N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j)
            {
                Z(i,j) = kernel_.selfImpedanceAnalytic(geom, i);
                continue;
            }

            auto psiS = kernel_.scalarPotential(geom, i, j);
            auto psiV = kernel_.vectorPotential(geom, i, j);

            Vec3 Ti = geom.segments[i].unitDir(geom);
            Vec3 Aj = geom.segments[j].unitDir(geom);
            double dotTA = Ti.dot(Aj);

            auto vecTerm = dotTA * psiV;

            Complex grad(0.0, 0.0);
            if (!geom.areAdjacentOnSameWire(i, j))
                grad = gradPhiContributionMININEC(geom, i, j);

            Z(i,j) = psiS + vecTerm + grad;
        }
    }


/*
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {

            auto psiS = kernel_.scalarPotential(geom, i, j);
            auto psiV = kernel_.vectorPotential(geom, i, j);
            auto grad = gradPhiContribution(geom, i, j);

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

            Z(i,j) = psiS + psiV + grad;
        }
    }
*/

    return Z;
}

/*
std::complex<double>
ImpedanceAssembler::gradPhiContribution(
    const Geometry& geom,
    int obsSeg,
    int srcSeg
    ) const
{
    const auto& seg = geom.segments[obsSeg];

    // Nodes at ends of observation segment
    int nMinus = seg.n1;
    int nPlus  = seg.n2;

    const auto& pMinus = geom.nodes[nMinus];
    const auto& pPlus  = geom.nodes[nPlus];

    // Segment length
    double dx = pPlus.x - pMinus.x;
    double dy = pPlus.y - pMinus.y;
    double dz = pPlus.z - pMinus.z;
    double L  = std::sqrt(dx*dx + dy*dy + dz*dz);

    if (L < 1e-12)
        return {0.0, 0.0};

    // --- Evaluate scalar potential at segment endpoints ---
    // We temporarily treat endpoints as "observation points"

    std::complex<double> phiPlus  =
        kernel_.scalarPotential(geom, obsSeg, srcSeg);

    std::complex<double> phiMinus =
        kernel_.scalarPotential(geom, obsSeg, srcSeg);

    // NOTE:
    // In a fully exact NEC-style code, phiPlus and phiMinus
    // would be evaluated at node positions.
    // In MININEC, both are derived from segment-based ψ,
    // and the difference is scaled by orientation.

    // Discrete gradient along segment direction
    std::complex<double> gradPhi =
        (phiPlus - phiMinus) / L;

    return gradPhi;
}
*/


std::complex<double>
ImpedanceAssembler::gradPhiContributionMININEC(
    const Geometry& geom,
    int I,      // observationssegment (M)
    int J      // källsegment (N)
    // int J1,     // nedre källände
    // int J2,     // övre källände
    // double F4,  // tecken från C(J,1)
    // double F5,  // tecken från C(J,2)
    // int F8      // specialflagga
    ) const
{
    using Cx = std::complex<double>;

    const auto& segJ = geom.segments[J];
    int J1 = segJ.n1;
    int J2 = segJ.n2;

    // Temporärt (innan full F8-logik återinförs)
    int F8 = 0;
    double F4 = 1.0;
    double F5 = 1.0;


    double U1 = 0.0, U2 = 0.0;
    double U3 = 0.0, U4 = 0.0;
    double U5 = 0.0, U6 = 0.0;
    double T1 = 0.0, T2 = 0.0;

    //------------------------------------------------------------------
    // BLOCK 3–4: ψ(M+½, N, N+1)
    //------------------------------------------------------------------

    if (F8 == 1) {
        // rekursivt specialfall (rad 279–282)
        // U5 = F5*U1 + T1  (U1/T1 från tidigare vektor-ψ)
        // Här antar vi att U1/U2 redan finns om F8==1
        U5 = F5 * U1 + T1;
        U6 = F5 * U2 + T2;
    } else {
        // GOSUB 87
        auto psi = kernel_.scalarPotentialMININEC(
            geom,
            I,
            J,
            +0.5,   // M+½
            +1      // (N, N+1)
            );
        U5 = psi.real();
        U6 = psi.imag();
    }

    //------------------------------------------------------------------
    // BLOCK 5: ψ(M−½, N, N+1)
    //------------------------------------------------------------------

    {
        auto psi = kernel_.scalarPotentialMININEC(
            geom,
            I,
            J,
            -0.5,   // M−½
            +1      // (N, N+1)
            );

        T1 = psi.real();
        T2 = psi.imag();

        U1 = (U5 - T1) / geom.segmentLength(J2);
        U2 = (U6 - T2) / geom.segmentLength(J2);
    }

    //------------------------------------------------------------------
    // BLOCK 6: ψ(M+½, N−1, N)
    //------------------------------------------------------------------

    {
        auto psi = kernel_.scalarPotentialMININEC(
            geom,
            I,
            J,
            +0.5,   // M+½
            -1      // (N−1, N)
            );

        U3 = psi.real();
        U4 = psi.imag();
    }

    //------------------------------------------------------------------
    // BLOCK 7: ψ(M−½, N−1, N)
    //------------------------------------------------------------------

    if (F8 >= 1) {
        // specialfall (rad 303–306)
        T1 = U5;
        T2 = U6;
    } else {
        auto psi = kernel_.scalarPotentialMININEC(
            geom,
            I,
            J,
            -0.5,   // M−½
            -1      // (N−1, N)
            );
        T1 = psi.real();
        T2 = psi.imag();
    }

    //------------------------------------------------------------------
    // BLOCK 8: slutlig grad(Φ)
    //------------------------------------------------------------------

    U1 += (U3 - T1) / geom.segmentLength(J1);
    U2 += (U4 - T2) / geom.segmentLength(J1);

    //------------------------------------------------------------------
    // Tecken från strömriktning (rad 311–312 används senare i BASIC)
    //------------------------------------------------------------------

    return Cx(F4 * U1, F4 * U2);
}
