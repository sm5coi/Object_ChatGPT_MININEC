#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "Vec3.hpp"
#include <vector>

class Geometry {
public:
    struct Node {
        double x, y, z;
    };

    struct Segment {
        int n1, n2;     // nodindex
        int wire;       // wire-id
        double radius;

        Vec3 unitDir(const Geometry& geom) const;
        double length(const Geometry& geom) const;
    };

    std::vector<Node>    nodes;
    std::vector<Segment> segments;
    std::vector<int>     wireOfSegment;

    double frequencyHz = 0.0;

    bool areAdjacentOnSameWire(int i, int j) const;

    static Geometry standardDipoleMININEC(
        double freqHz,
        int    nSeg,
        double radius
        );

    double segmentLength(int s) const;
};


#endif
