#ifndef KERNEL_MININEC_HPP
#define KERNEL_MININEC_HPP

#include <complex>
#include "kernel.hpp"
#include "geometry.hpp"

class MininecKernel : public EMKernel {
public:
    explicit MininecKernel(double frequency);

    std::complex<double>
    scalarPotential(
        const Geometry& geom,
        int obsSeg,
        int srcSeg
        ) const override;

    std::complex<double>
    vectorPotential(
        const Geometry& geom,
        int obsSeg,
        int srcSeg
        ) const override;

private:
    // Physical constants
    double k_;     // wave number
    double w2_;    // MININEC normalization (~15)

    // ---- internal helpers ----
    std::complex<double>
    psiScalarSelf(const Geometry& geom, int seg) const;

    std::complex<double>
    psiVectorSelf(const Geometry& geom, int seg) const;

    std::complex<double>
    psiGauss(
        const Geometry& geom,
        int obsSeg,
        int srcSeg,
        bool vectorCase
        ) const;

    void kernel28(
        double& re, double& im,
        double dx, double dy, double dz
        ) const;
};

#endif
