#ifndef IMPEDANCE_ASSEMBLER_HPP
#define IMPEDANCE_ASSEMBLER_HPP

#include "geometry.hpp"
#include "kernel.hpp"
#include "impedance_matrix.hpp"

class ImpedanceAssembler {
public:
    explicit ImpedanceAssembler(const EMKernel& kernel);

    ImpedanceMatrix build(const Geometry& geom) const;

private:
    const EMKernel& kernel_;

    // grad(Φ)-delen – separat för tydlighet
    std::complex<double>
    gradPhiContribution(
        const Geometry& geom,
        int obsSeg,
        int srcSeg
        ) const;
};

#endif
