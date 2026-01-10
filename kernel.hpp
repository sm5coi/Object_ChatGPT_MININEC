#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <complex>
#include "geometry.hpp"

class EMKernel {
public:
    virtual std::complex<double>
    scalarPotential(
        const Geometry& geom,
        int obsSeg,
        int srcSeg
        ) const = 0;

    virtual std::complex<double>
    vectorPotential(
        const Geometry& geom,
        int obsSeg,
        int srcSeg
        ) const = 0;

    virtual ~EMKernel() = default;
};

#endif
