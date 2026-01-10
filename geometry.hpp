#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>

struct Node {
    double x;
    double y;
    double z;
};

struct Segment {
    int n1;           // start node index
    int n2;           // end node index
    double radius;    // wire radius
};

class Geometry {
public:
    std::vector<Node> nodes;
    std::vector<Segment> segments;

    int feedSegment = -1;

    static Geometry dipole(
        double length,
        double radius,
        int segmentCount
        );
};

#endif
