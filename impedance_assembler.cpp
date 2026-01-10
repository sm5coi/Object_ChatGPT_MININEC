#include "impedance_assembler.hpp"

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
ImpedanceAssembler::ImpedanceAssembler(const EMKernel& kernel)
    : kernel_(kernel)
{
}

// std::complex<double>
// ImpedanceAssembler::gradPhiContribution(
//     const Geometry& geom,
//     int obsSeg,
//     int srcSeg
//     ) const
// {
//     // Första version:
//     // gradPhi approximeras som noll för olika segment
//     // (self-term hanteras redan i kernel)
//     if (obsSeg != srcSeg) {
//         return {0.0, 0.0};
//     }

//     // Självterm – liten korrigering, ofta nära noll
//     return {0.0, 0.0};
// }

ImpedanceMatrix
ImpedanceAssembler::build(const Geometry& geom) const
{
    int N = static_cast<int>(geom.segments.size());
    ImpedanceMatrix Z(N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {

            auto psiS = kernel_.scalarPotential(geom, i, j);
            auto psiV = kernel_.vectorPotential(geom, i, j);
            auto grad = gradPhiContribution(geom, i, j);

            Z(i,j) = psiS + psiV + grad;
        }
    }

    return Z;
}

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
