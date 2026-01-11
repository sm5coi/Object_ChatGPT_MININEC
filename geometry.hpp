#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>

struct Node {
    double x;
    double y;
    double z;
};

struct Segment {
    int n1;
    int n2;
    int wire;   // index W
};


class Geometry {
public:
    // Topologi
    std::vector<Node> nodes;
    std::vector<Segment> segments;

    // === MININEC ===
    std::vector<double> wireRadius;   // A(W)
    std::vector<double> wireSegLen;   // S(W)

    // MININEC: vilken tråd varje segment ligger på
    // motsvarar W%(I) i BASIC
    std::vector<int> wireOfSegment;

    // Hjälpfunktioner (ren geometri)
    double segmentLength(int seg) const;

    // === DEN DU FRÅGAR EFTER ===
    static Geometry standardDipoleMININEC();
};

#endif
