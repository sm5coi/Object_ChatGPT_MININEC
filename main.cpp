#include <iostream>

#include "geometry.hpp"
#include "mininec_kernel.hpp"
#include "impedance_assembler.hpp"
#include "impedance_matrix.hpp"

// #include "geometry.hpp"
// #include "kernel.hpp"
// #include "impedance_assembler.hpp"
// #include "kernel_mininec.hpp"
// #include "mom_solver.hpp"


// Forward declaration (tillfälligt)
class MininecKernel;

int main()
{
    double length = 15.0;     // meters
    int segments  = 11;
    double freq   = 10e6;     // 10 MHz
    double c0   = 299792458.0;
    double lambda = c0 / freq;
    double k = 2.0 * M_PI / lambda;

    double radius = 0.004;      // 4 mm
    double srm = 2.0 * radius;  // exakt som MININEC

    MininecKernel kernel(k, srm);

    Geometry geom; // = Geometry::dipole(length, radius, segments);

    std::cout << "Dipole geometry\n";
    std::cout << "Nodes:    " << geom.nodes.size() << "\n";
    std::cout << "Segments: " << geom.segments.size() << "\n";
    //std::cout << "Feed segment index: " << geom.feedSegment << "\n\n";

    for (size_t i = 0; i < geom.segments.size(); ++i) {
        const auto& s = geom.segments[i];
        const auto& n1 = geom.nodes[s.n1];
        const auto& n2 = geom.nodes[s.n2];

        std::cout << "Segment " << i
                  << " : (" << n1.x << "," << n1.y << "," << n1.z << ") -> "
                  << "(" << n2.x << "," << n2.y << "," << n2.z << ")\n";
    }

    ImpedanceAssembler assembler(kernel);

    ImpedanceMatrix Z = assembler.build(geom);

    // std::cout << "\nZ-matrix (real part):\n";
    // for (int i = 0; i < Z.size(); ++i) {
    //     for (int j = 0; j < Z.size(); ++j) {
    //         std::cout << Z(i,j).real() << " ";
    //     }
    //     std::cout << "\n";
    // }

    // MoMSolver solver;

    // auto Zin = solver.inputImpedance(Z, geom.feedSegment);

    // std::cout << "\nInput impedance:\n";
    // std::cout << "Zin = "
    //           << Zin.real()
    //           << " + j"
    //           << Zin.imag()
    //           << " ohm\n";



    return 0;
}
