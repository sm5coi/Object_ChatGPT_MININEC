#include "geometry.hpp"
#include <cmath>

Geometry Geometry::standardDipoleMININEC(
    double freqHz,
    int    nSeg,
    double radius
    )
{
    Geometry g;

    const double c0 = 299792458.0;
    const double lambda = c0 / freqHz;
    const double L = 0.5 * lambda;          // halvvågsdipol
    const double dz = L / nSeg;

    // --- skapa noder ---
    for (int i = 0; i <= nSeg; ++i) {
        double z = -0.5 * L + i * dz;
        g.nodes.push_back({0.0, 0.0, z});
    }

    // --- skapa segment ---
    for (int i = 0; i < nSeg; ++i) {
        Geometry::Segment seg;
        seg.n1 = i;
        seg.n2 = i + 1;
        seg.wire = 0;           // endast en wire
        seg.radius = radius;

        g.segments.push_back(seg);
        //g.wireOfSegment.push_back(0);
    }

    g.frequencyHz = freqHz;

    return g;
}

double Geometry::segmentLength(int s) const
{
    const auto& seg = segments[s];
    const auto& a = nodes[seg.n1];
    const auto& b = nodes[seg.n2];

    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double dz = b.z - a.z;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

Vec3 Geometry::Segment::unitDir(const Geometry& geom) const
{
    const auto& n1p = geom.nodes[n1];
    const auto& n2p = geom.nodes[n2];

    Vec3 d(
        n2p.x - n1p.x,
        n2p.y - n1p.y,
        n2p.z - n1p.z
        );

    double L = d.norm();
    assert(L > 0.0);
    return d / L;
}


bool Geometry::areAdjacentOnSameWire(int i, int j) const
{
    const auto& si = segments[i];
    const auto& sj = segments[j];

    if (si.wire != sj.wire)
        return false;

    return std::abs(i - j) <= 1;
}

double Geometry::Segment::length(const Geometry& geom) const
{
    const auto& p1 = geom.nodes[n1];
    const auto& p2 = geom.nodes[n2];

    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

