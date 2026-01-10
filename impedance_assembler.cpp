#include "impedance_assembler.hpp"

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
ImpedanceAssembler::ImpedanceAssembler(const EMKernel& kernel)
    : kernel_(kernel)
{
}

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
