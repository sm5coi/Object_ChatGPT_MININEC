#include "mininec_kernel.hpp"
#include "geometry.hpp"
#include <cmath>

namespace {
constexpr double PI = 3.14159265358979323846;
constexpr double C0 = 299792458.0;
}

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
MininecKernel::MininecKernel(double k, double srm)
    : k_(k),
    srm_(srm)
{
    // MININEC-konstanter för exact kernel (C0..C9)
    // C_[0] = 1.38629436112;
    // C_[1] = 0.09666344259;
    // C_[2] = 0.03590092383;
    // C_[3] = 0.03742563713;
    // C_[4] = 0.01451196212;
    // C_[5] = 0.5;
    // C_[6] = 0.12498593397;
    // C_[7] = 0.06880248576;
    // C_[8] = 0.0332835346;
    // C_[9] = 0.00441787012;
}




// MininecKernel::MininecKernel(double frequency)
// {
//     double lambda = C0 / frequency;
//     k_  = 2.0 * PI / lambda;
//     w2_ = 15.0;   // ≈ Z0 / (8π)
// }


inline int wireOf(const Geometry& geom, int seg)
{
    return geom.segments[seg].wire;
}

int MininecKernel::baseP1(const Geometry& geom, int obsSeg) const
{
    // BASIC: P1 = 2*W%(I) + I - 1
    return 2 * geom.wireOfSegment[obsSeg] + obsSeg - 1;
}

int MininecKernel::baseP2(const Geometry& geom, int srcSeg) const
{
    // BASIC: P2 = 2*W%(J) + J - 1
    return 2 * geom.wireOfSegment[srcSeg] + srcSeg - 1;
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

MininecKernel::Complex
MininecKernel::scalarPotential(
    const Geometry& geom,
    int obsSeg,
    int srcSeg
    ) const
{
    return psiGauss(geom, obsSeg, srcSeg, false);
}

MininecKernel::Complex
MininecKernel::scalarPotentialMININEC(
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

MininecKernel::Complex
MininecKernel::vectorPotentialMININEC(
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


std::complex<double>
MininecKernel::psiGauss(
    const Geometry& geom,
    int obsSeg,
    int srcSeg,
    bool vectorMode
    ) const
{
    int P1 = baseP1(geom, obsSeg);
    int P2 = baseP2(geom, srcSeg);
    int P3 = P2 + 1;

    return psiGaussCore(geom, obsSeg, srcSeg, P1, P2, P3, vectorMode);
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
    int P1 = baseP1(geom, obsSeg)
    + (obsOffset > 0 ? 1 : -1);

    int P2, P3;
    if (srcMode > 0) {
        P2 = baseP2(geom, srcSeg);
        P3 = P2 + 1;
    } else {
        P3 = baseP2(geom, srcSeg);
        P2 = P3 - 1;
    }

    return psiGaussCore(geom, obsSeg, srcSeg, P1, P2, P3, vectorMode);
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

    if (vectorMode) {
        // GOSUB 102: nod
        X1 = geom.nodes[P1].x;
        Y1 = geom.nodes[P1].y;
        Z1 = geom.nodes[P1].z;
    } else {
        // GOSUB 87: mittpunkt av segment
        const auto& s = geom.segments[P1];
        const auto& a = geom.nodes[s.n1];
        const auto& b = geom.nodes[s.n2];

        X1 = 0.5 * (a.x + b.x);
        Y1 = 0.5 * (a.y + b.y);
        Z1 = 0.5 * (a.z + b.z);
    }

    // --------------------------------------------------
    // 2) Källpunkter X2 och X3
    // --------------------------------------------------
    double X2, Y2, Z2;
    double X3, Y3, Z3;

    // P2
    {
        const auto& s = geom.segments[P2];
        const auto& a = geom.nodes[s.n1];
        const auto& b = geom.nodes[s.n2];

        X2 = 0.5 * (a.x + b.x);
        Y2 = 0.5 * (a.y + b.y);
        Z2 = 0.5 * (a.z + b.z);
    }

    // P3
    {
        const auto& s = geom.segments[P3];
        const auto& a = geom.nodes[s.n1];
        const auto& b = geom.nodes[s.n2];

        X3 = 0.5 * (a.x + b.x);
        Y3 = 0.5 * (a.y + b.y);
        Z3 = 0.5 * (a.z + b.z);
    }

    // --------------------------------------------------
    // 3) Relativ vektor (DET DU FRÅGADE OM)
    // --------------------------------------------------
    double dx = X3 - X1;
    double dy = Y3 - Y1;
    double dz = Z3 - Z1;

    // --------------------------------------------------
    // 4) Gauss + kernel28 (kommer härnäst)
    // --------------------------------------------------
    int w = geom.segments[srcSeg].wire;
    double A  = geom.wireRadius[w];
    double A2 = A * A;
    double I6 = 0.0;   // tills vidare

    kernel28(re, im, dx, dy, dz, A, A2, I6);

    return {re, im};
}



