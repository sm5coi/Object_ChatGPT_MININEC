#include "geometry.hpp"
#include <cmath>

Geometry Geometry::dipole(
    double length,
    double radius,
    int segmentCount
    ) {
    Geometry g;

    // Build nodes along Z-axis, centered at origin
    double dz = length / segmentCount;
    int nodeCount = segmentCount + 1;

    g.nodes.reserve(nodeCount);

    double z0 = -0.5 * length;
    for (int i = 0; i < nodeCount; ++i) {
        g.nodes.push_back({0.0, 0.0, z0 + i * dz});
    }

    // Build segments
    g.segments.reserve(segmentCount);
    for (int i = 0; i < segmentCount; ++i) {
        g.segments.push_back({
            i,        // node i
            i + 1,    // node i+1
            radius
        });
    }

    // Feed at center segment
    g.feedSegment = segmentCount / 2;

    return g;
}


