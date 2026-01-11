#include "geometry.hpp"
#include <cmath>

double Geometry::segmentLength(int seg) const
{
    const auto& s = segments.at(seg);
    const auto& a = nodes.at(s.n1);
    const auto& b = nodes.at(s.n2);

    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double dz = b.z - a.z;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

Geometry Geometry::standardDipoleMININEC()
{
    Geometry g;

    const double a = 0.004;
    const double h = 0.191;
    const double L = 0.309;

    const int Nv = 4;
    const int Nh = 6;

    // --- Noder ---
    g.nodes.push_back({0,0,0});
    for (int i = 1; i <= Nv; ++i)
        g.nodes.push_back({0,0,h*i/Nv});

    int top = g.nodes.size() - 1;
    for (int i = 1; i <= Nh; ++i)
        g.nodes.push_back({0, L*i/Nh, h});

    // --- Trådar ---
    g.wireRadius = { a, a };
    g.wireSegLen = { h/Nv, L/Nh };

    // --- Segment ---
    for (int i = 0; i < Nv; ++i)
        g.segments.push_back({i, i+1, 0});

    for (int i = 0; i < Nh; ++i)
        g.segments.push_back({top+i, top+i+1, 1});

    return g;
}


