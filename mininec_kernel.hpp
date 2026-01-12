#ifndef MININEC_KERNEL_H
#define MININEC_KERNEL_H

#include "geometry.hpp"
#include "types.hpp"
#include <complex>

class Geometry;

class MininecKernel {
public:
    //using Complex = std::complex<double>;

    explicit MininecKernel(double k, double srm);

    Complex scalarPotential(const Geometry&, int obsSeg, int srcSeg) const;

    Complex gradPhiContributionMININEC(const Geometry&, int obsSeg, int srcSeg) const;

    // TEMPORÄR STUB
    Complex vectorPotential(const Geometry&, int obsSeg, int srcSeg) const;
    Complex selfImpedanceAnalytic(const Geometry&, int seg) const;
    // -------------------------------------------------
    // Standard scalar ψ
    // Motsvarar: ψ(M, N, N+½)
    // Används av vektorpotentialdelen
    // -------------------------------------------------
    // Complex scalarPotential(
    //     const Geometry& geom,
    //     int obsSeg,
    //     int srcSeg
    //     ) const;

    // -------------------------------------------------
    // MININEC-exakt scalar ψ (GOSUB 87)
    // Används ENDAST av grad(Φ)
    // -------------------------------------------------
    Complex scalarPotentialMININEC(
        const Geometry& geom,
        int obsSeg,          // M
        int srcSeg,          // N
        double obsOffset,    // +0.5 eller −0.5
        int srcMode          // +1 = (N,N+1), −1 = (N−1,N)
        ) const;

    // -------------------------------------------------
    // MININEC-exakt vector ψ (GOSUB 102)
    // -------------------------------------------------
    Complex vectorPotentialMININEC(
        const Geometry& geom,
        int obsSeg,
        int srcSeg,
        double obsOffset,
        int srcMode
        ) const;

private:
    double k_;     // vågtal
    double srm_;   // small-radius-modifier (SRM)

    static constexpr double C_[10] = {
        1.38629436112,
        0.09666344259,
        0.03590092383,
        0.03742563713,
        0.01451196212,
        0.5,
        0.12498593397,
        0.06880248576,
        0.03328353460,
        0.00441787012
    };

    // ===== Indexhjälpare (BASIC P1/P2) =====
    int baseP1(const Geometry& geom, int obsSeg) const;
    int baseP2(const Geometry& geom, int srcSeg) const;

    // ===== Kärnan: exakt GOSUB 87 / 102 =====
    std::complex<double>
    psiGaussCore(
        const Geometry& geom,
        int obsSeg,
        int srcSeg,
        int P1,
        int P2,
        int P3,
        bool vectorMode
        ) const
    ;

    // ===== Wrapper: standardfall =====
    Complex psiGauss(
        const Geometry& geom,
        int obsSeg,
        int srcSeg,
        bool vectorMode
        ) const;

    // ===== Wrapper: MININEC-indexstyrd =====
    Complex psiGaussMININEC(
        const Geometry& geom,
        int obsSeg,
        int srcSeg,
        double obsOffset,
        int srcMode,
        bool vectorMode
        ) const;

    void kernel28(
        double& re,
        double& im,
        double dx,
        double dy,
        double dz,
        double A,      // trådradie
        double A2,     // A*A
        double I6      // exact kernel-term (0 eller ≠0)
        ) const;
};

#endif
