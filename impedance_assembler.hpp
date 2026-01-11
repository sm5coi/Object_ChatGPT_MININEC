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




// #include "geometry.hpp"
// #include "kernel.hpp"
// #include "impedance_matrix.hpp"

// class ImpedanceAssembler {
// public:
//     explicit ImpedanceAssembler(const EMKernel& kernel);

//     ImpedanceMatrix build(const Geometry& geom) const;

// private:
//     const EMKernel& kernel_;

//     // grad(Φ)-delen – separat för tydlighet
//     std::complex<double>
//     gradPhiContributionMININEC(
//         const Geometry& geom,
//         int I,      // observationssegment (M)
//         int J,      // källsegment (N)
//         int J1,     // nedre källände
//         int J2,     // övre källände
//         double F4,  // tecken från C(J,1)
//         double F5,  // tecken från C(J,2)
//         int F8      // specialflagga
//         ) const;
// };

#endif
