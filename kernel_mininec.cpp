#include "kernel_mininec.hpp"
#include <cmath>

namespace {
constexpr double PI = 3.14159265358979323846;
constexpr double C0 = 299792458.0;
}

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
MininecKernel::MininecKernel(double frequency)
{
    double lambda = C0 / frequency;
    k_  = 2.0 * PI / lambda;
    w2_ = 15.0;   // ≈ Z0 / (8π)
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

std::complex<double>
MininecKernel::scalarPotential(
    const Geometry& geom,
    int obsSeg,
    int srcSeg
    ) const
{
    if (obsSeg == srcSeg) {
        return w2_ * psiScalarSelf(geom, srcSeg);
    }

    return w2_ * psiGauss(geom, obsSeg, srcSeg, false);
}

std::complex<double>
MininecKernel::vectorPotential(
    const Geometry& geom,
    int obsSeg,
    int srcSeg
    ) const
{
    if (obsSeg == srcSeg) {
        return w2_ * psiVectorSelf(geom, srcSeg);
    }

    return w2_ * psiGauss(geom, obsSeg, srcSeg, true);
}

void MininecKernel::kernel28(
    double& re,
    double& im,
    double dx,
    double dy,
    double dz
    ) const
{
    double R = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (R < 1e-12) R = 1e-12;

    double phase = k_ * R;

    re += std::cos(phase) / R;
    im -= std::sin(phase) / R;
}


std::complex<double>
MininecKernel::psiGauss(
    const Geometry& geom,
    int obsSeg,
    int srcSeg,
    bool vectorCase
    ) const
{
    auto obs = segmentGeom(geom, obsSeg);
    auto src = segmentGeom(geom, srcSeg);

    // Gauss 2-punkt (kan utökas senare)
    constexpr double w[2] = {1.0, 1.0};
    constexpr double x[2] = {-0.5773502692, 0.5773502692};

    double re = 0.0;
    double im = 0.0;

    for (int i = 0; i < 2; ++i) {
        double t = 0.5 * (x[i] + 1.0);

        double xs = src.x1 + t * (src.x2 - src.x1);
        double ys = src.y1 + t * (src.y2 - src.y1);
        double zs = src.z1 + t * (src.z2 - src.z1);

        double xo = 0.5 * (obs.x1 + obs.x2);
        double yo = 0.5 * (obs.y1 + obs.y2);
        double zo = 0.5 * (obs.z1 + obs.z2);

        kernel28(
            re, im,
            xo - xs,
            yo - ys,
            zo - zs
            );
    }

    double scale = 0.5 * src.length;
    return scale * std::complex<double>(re, im);
}

std::complex<double>
MininecKernel::psiScalarSelf(const Geometry& geom, int s) const
{
    double S = segmentLength(geom, s);
    double a = geom.segments[s].radius;
    return { 2.0 * std::log(S / a), -k_ * S };
}

std::complex<double>
MininecKernel::psiVectorSelf(const Geometry& geom, int s) const
{
    double S = segmentLength(geom, s);
    double a = geom.segments[s].radius;
    return { std::log(S / a), -0.5 * k_ * S };
}


