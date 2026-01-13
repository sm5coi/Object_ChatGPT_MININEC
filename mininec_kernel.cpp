#include <iostream>
#include "mininec_kernel.hpp"
#include "geometry.hpp"
#include "types.hpp"
#include <cmath>

namespace {
constexpr double PI = 3.14159265358979323846;
constexpr double C0 = 299792458.0;
}

//using Complex = std::complex<double>;

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
MininecKernel::MininecKernel(double k, double srm)
    : k_(k),
    srm_(srm)
{

}

inline int wireOf(const Geometry& geom, int seg)
{
    return geom.segments[seg].wire;
}

int MininecKernel::baseP1(const Geometry& geom, int obsSeg) const
{
    //
    return 2 * geom.segments[obsSeg].wire + obsSeg;

}

int MininecKernel::baseP2(const Geometry& geom, int srcSeg) const
{
    // C++-anpassad MININEC-formel
    return 2 * geom.segments[srcSeg].wire + srcSeg;
}

struct SegmentGeom {
    double x1, y1, z1;
    double x2, y2, z2;
    double length;
};

static void segmentMidpoint(
    const Geometry& g, int s,
    double& x, double& y, double& z)
{
    const auto& seg = g.segments[s];
    const auto& n1  = g.nodes[seg.n1];
    const auto& n2  = g.nodes[seg.n2];

    x = 0.5 * (n1.x + n2.x);
    y = 0.5 * (n1.y + n2.y);
    z = 0.5 * (n1.z + n2.z);
}

static SegmentGeom segmentGeom(const Geometry& g, int s)
{
    const auto& seg = g.segments[s];
    const auto& a = g.nodes[seg.n1];
    const auto& b = g.nodes[seg.n2];

    SegmentGeom sg;
    sg.x1 = a.x; sg.y1 = a.y; sg.z1 = a.z;
    sg.x2 = b.x; sg.y2 = b.y; sg.z2 = b.z;

    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double dz = b.z - a.z;
    sg.length = std::sqrt(dx*dx + dy*dy + dz*dz);

    return sg;
}

static double segmentLength(const Geometry& g, int s)
{
    const auto& seg = g.segments[s];
    const auto& a = g.nodes[seg.n1];
    const auto& b = g.nodes[seg.n2];

    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double dz = b.z - a.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}


Complex MininecKernel::scalarPotential(
    const Geometry& geom,
    int obsSeg,
    int srcSeg
    ) const
{
    return psiGauss(geom, obsSeg, srcSeg, false);
}


Complex MininecKernel::scalarPotentialMININEC(
    const Geometry& geom,
    int obsSeg,
    int srcSeg,
    double obsOffset,
    int srcMode
    ) const
{
    return psiGaussMININEC(
        geom,
        obsSeg,
        srcSeg,
        obsOffset,
        srcMode,
        false
        );
}

Complex MininecKernel::vectorPotentialMININEC(
    const Geometry& geom,
    int obsSeg,
    int srcSeg,
    double obsOffset,
    int srcMode
    ) const
{
    return psiGaussMININEC(
        geom,
        obsSeg,
        srcSeg,
        obsOffset,
        srcMode,
        true
        );
}

void MininecKernel::kernel28(
    double& re,
    double& im,
    double dx,
    double dy,
    double dz,
    double A,
    double A2,
    double I6
    ) const
{
    double D3 = dx*dx + dy*dy + dz*dz;

    double D;
    if (A <= srm_) {
        D = std::sqrt(D3);
    } else {
        D = std::sqrt(D3 + A2);
    }

    if (D < 1e-12)
        D = 1e-12;

    // --- exact kernel correction (MININEC I6 logic) ---
    if (I6 != 0.0) {
        double B = D3 / (D3 + 4.0 * A2);

        double W0 =
            C_[0] + B*(C_[1] + B*(C_[2] + B*(C_[3] + B*C_[4])));
        double W1 =
            C_[5] + B*(C_[6] + B*(C_[7] + B*(C_[8] + B*C_[9])));

        double V0 = (W0 - W1 * std::log(B)) * std::sqrt(1.0 - B);

        re += (V0 + 0.5 * std::log(D3 / (64.0 * A2)))
                  / (M_PI * A) - 1.0 / D;
    }

    // --- exp(-j k R) / R ---
    double phase = k_ * D;
    re += std::cos(phase) / D;
    im -= std::sin(phase) / D;
}

Complex MininecKernel::kernel28Value(
    double dx,
    double dy,
    double dz,
    double A,
    double A2,
    double I6
    ) const
{
    double re = 0.0;
    double im = 0.0;
    kernel28(re, im, dx, dy, dz, A, A2, I6);
    return Complex(re, im);
}


std::complex<double>
MininecKernel::psiGauss(
    const Geometry& geom,
    int obsSeg,
    int srcSeg,
    bool vectorMode
    ) const
{
    // Standardfall: observation = segmentmitt (obsSeg)
    // Source: segment srcSeg n1->n2

    const auto& sSrc = geom.segments[srcSeg];
    int P2 = sSrc.n1;
    int P3 = sSrc.n2;

    // Core ska använda segmentmitt => obsIsNode=false
    bool obsIsNode = false;

    // P1 används inte i detta läge (men vi skickar obsSeg för formalia)
    int P1 = obsSeg;

    return psiGaussCore(geom, obsSeg, srcSeg, P1, P2, P3, obsIsNode);
}

std::complex<double>
MininecKernel::psiGaussMININEC(
    const Geometry& geom,
    int obsSeg,
    int srcSeg,
    double obsOffset,
    int srcMode,
    bool vectorMode
    ) const
{


    const int nSeg = (int)geom.segments.size();

    auto clampSeg = [&](int s) -> int {
        if (s < 0) return 0;
        if (s >= nSeg) return nSeg - 1;
        return s;
    };

    int useObsSeg = clampSeg(obsSeg);
    int useSrcSeg = clampSeg(srcSeg);

    // --- P1 = observationsnod (M±1/2) ---
    const auto& sObs = geom.segments[useObsSeg];
    int P1 = (obsOffset >= 0.0) ? sObs.n2 : sObs.n1;

    // --- MININEC: srcMode=-1 betyder segmentet (N-1) ---
    int useSrcSeg2 = useSrcSeg;

    if (srcMode < 0)
    {
        // MININEC: (N-1, N) kräver att N-1 existerar
        if (useSrcSeg == 0)
            return Complex(0.0, 0.0);

        useSrcSeg2 = useSrcSeg - 1; // INGEN clamp här
    }


    const auto& sSrc = geom.segments[useSrcSeg2];
    int P2 = sSrc.n1;
    int P3 = sSrc.n2;

    std::cout
        << "psiGaussMININEC obsSeg=" << obsSeg
        << " srcSeg=" << srcSeg
        << " obsOffset=" << obsOffset
        << " srcMode=" << srcMode
        << " => P1=" << P1
        << " P2=" << P2
        << " P3=" << P3
        << "\n";

    bool obsIsNode = true;
    return psiGaussCore(geom, useObsSeg, useSrcSeg2, P1, P2, P3, obsIsNode);
}


std::complex<double>
MininecKernel::psiGaussCore(
    const Geometry& geom,
    int obsSeg,
    int srcSeg,
    int P1,
    int P2,
    int P3,
    bool vectorMode
    ) const
{
    double re = 0.0;
    double im = 0.0;

    // --------------------------------------------------
    // 1) Observationspunkt X1,Y1,Z1
    // --------------------------------------------------
    double X1, Y1, Z1;

    // VIKTIGT:
    // - vectorMode==true  => P1 är NODE-index (GOSUB 102 / nod-observation)
    // - vectorMode==false => observation är SEGMENTMITT (standardfall)
    //
    // OBS: INGEN "gissning" på P1 längre.
    if (vectorMode) {
        // Observation i nod P1
        const auto& n = geom.nodes[P1];
        X1 = n.x; Y1 = n.y; Z1 = n.z;
    } else {
        // Observation i segmentmitt (obsSeg)
        const auto& s = geom.segments[obsSeg];
        const auto& a = geom.nodes[s.n1];
        const auto& b = geom.nodes[s.n2];
        X1 = 0.5*(a.x + b.x);
        Y1 = 0.5*(a.y + b.y);
        Z1 = 0.5*(a.z + b.z);
    }

    // --------------------------------------------------
    // 2) Källpunkter X2 och X3 (NODER)
    // --------------------------------------------------
    const auto& n2 = geom.nodes[P2];
    const auto& n3 = geom.nodes[P3];

    double X2 = n2.x, Y2 = n2.y, Z2 = n2.z;
    double X3 = n3.x, Y3 = n3.y, Z3 = n3.z;

    // --------------------------------------------------
    // 3) Gauss-integral längs källsegmentet P2->P3
    // --------------------------------------------------
    double sx = X3 - X2;
    double sy = Y3 - Y2;
    double sz = Z3 - Z2;

    double Ls = std::sqrt(sx*sx + sy*sy + sz*sz);
    if (Ls < 1e-12)
        return {0.0, 0.0};

    static const double t[4] = {
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526
    };

    static const double wgt[4] = {
        0.3478548451374538,
        0.6521451548625461,
        0.6521451548625461,
        0.3478548451374538
    };

    Complex sum(0.0, 0.0);

    for (int g = 0; g < 4; ++g)
    {
        double u = 0.5 * (1.0 + t[g]);

        double X = X2 + u * sx;
        double Y = Y2 + u * sy;
        double Z = Z2 + u * sz;

        double dx = X - X1;
        double dy = Y - Y1;
        double dz = Z - Z1;

        double A  = geom.segments[srcSeg].radius;
        double A2 = A * A;
        double I6 = 0.0;   // tills vidare

        sum += wgt[g] * kernel28Value(dx, dy, dz, A, A2, I6);
    }

    sum *= (Ls * 0.5);
    return sum;
}


Complex MininecKernel::vectorPotential(
    const Geometry& geom, int obsSeg, int srcSeg) const
{
    // Skalärpotential
    Complex psiS = scalarPotential(geom, obsSeg, srcSeg);

    // --- källsegmentets längd ---
    const auto& seg = geom.segments[srcSeg];
    const auto& p1  = geom.nodes[seg.n1];
    const auto& p2  = geom.nodes[seg.n2];

    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;
    double L  = std::sqrt(dx*dx + dy*dy + dz*dz);

    // --- vågtal ---
    double freq = geom.frequencyHz;
    double c0   = 299792458.0;
    double k    = 2.0 * M_PI * freq / c0;

    // --- MININEC första ordningens vektorterm ---
    double scale = std::cos(0.5 * k * L);

    return psiS * scale;
}

Complex MininecKernel::selfImpedanceAnalytic(
    const Geometry& geom, int segIdx) const
{
    const auto& seg = geom.segments[segIdx];

    const auto& p1 = geom.nodes[seg.n1];
    const auto& p2 = geom.nodes[seg.n2];

    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;

    double L = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (L <= 0.0)
        return Complex(0.0, 0.0);

    double a = 1e-3;

    const double gamma = 0.5772156649015328606;

    double realPart = std::log(2.0 * L / a) - 1.0 - gamma;

    double freq = geom.frequencyHz;
    double c0   = 299792458.0;
    double k    = 2.0 * M_PI * freq / c0;

    double imagPart = 0.5 * M_PI - k * L;

    return Complex(realPart, imagPart);
}

Complex MininecKernel::gradPhiContributionMININEC(
    const Geometry& geom,
    int obsSeg,
    int srcSeg
    ) const
{
    // std::cout << "gradPhi called obs=" << obsSeg << " src=" << srcSeg << "\n";

    // --------------------------------------------
    // MININEC BASIC 274–311 (princip):
    // gradPhi ≈ ( ψ(+,+) - ψ(+,-) - ψ(-,+) + ψ(-,-) ) / (Δz_obs * Δz_src)
    //
    // där:
    //   obsOffset = +0.5 / -0.5  (väljer "änden" på observationssegmentet)
    //   srcMode   = +1 / -1     (väljer (N,N+1) eller (N-1,N))
    //
    // Här använder vi scalarPotentialMININEC (GOSUB 87).
    // --------------------------------------------

    const double Lobs = geom.segments[obsSeg].length(geom);
    const double Lsrc = geom.segments[srcSeg].length(geom);

    if (Lobs <= 0.0 || Lsrc <= 0.0)
        return Complex(0.0, 0.0);

    // ψ( obsOffset, srcMode )
    Complex psi_pp = scalarPotentialMININEC(geom, obsSeg, srcSeg, +0.5, +1);
    Complex psi_pm = scalarPotentialMININEC(geom, obsSeg, srcSeg, +0.5, -1);
    Complex psi_mp = scalarPotentialMININEC(geom, obsSeg, srcSeg, -0.5, +1);
    Complex psi_mm = scalarPotentialMININEC(geom, obsSeg, srcSeg, -0.5, -1);

    // “Mixed finite difference”
    Complex mixed = (psi_pp - psi_pm - psi_mp + psi_mm);

    // Skala enligt segmentlängder (MININEC använder en diskret derivata)
    return mixed / (Lobs * Lsrc);
}

Complex MininecKernel::vectorPotentialContributionMININEC(
    const Geometry& geom,
    int obsSeg,      // M
    int srcSeg,      // N
    double obsOffset // +0.5 eller -0.5
    ) const
{
    // ψ(M±1/2, N,   N+1)
    Complex psi_plus  = psiGaussMININEC(geom, obsSeg, srcSeg, obsOffset, +1, true);

    // ψ(M±1/2, N-1, N)
    Complex psi_minus = psiGaussMININEC(geom, obsSeg, srcSeg, obsOffset, -1, true);

    // MININEC: vector potential term uses (psi_plus - psi_minus) / Lsrc
    const double Lsrc = geom.segmentLength(srcSeg);
    if (Lsrc < 1e-12) return Complex(0.0, 0.0);

    return (psi_plus - psi_minus) / Lsrc;
}



