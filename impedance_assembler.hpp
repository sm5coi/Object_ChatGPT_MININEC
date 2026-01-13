#ifndef IMPEDANCE_ASSEMBLER_HPP
#define IMPEDANCE_ASSEMBLER_HPP

#include "geometry.hpp"
#include "mininec_kernel.hpp"

#include <complex>

class Geometry;
class ImpedanceMatrix;
class MininecKernel;

class ImpedanceAssembler {
public:
    explicit ImpedanceAssembler(const MininecKernel& kernel);

    ImpedanceMatrix build(const Geometry& geom) const;

private:
    using Complex = std::complex<double>;

    const MininecKernel& kernel_;

    // -------------------------------------------------
    // Full MININEC-grad(Φ), rader 247–312
    // -------------------------------------------------
    Complex gradPhiContributionMININEC(
        const Geometry& geom,
        int obsSeg,   // M
        int srcSeg    // N
        ) const;
};

#endif
