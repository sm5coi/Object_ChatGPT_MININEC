#include "geometry.hpp"
#include "mininec_kernel.hpp"
#include "impedance_assembler.hpp"
#include "impedance_matrix.hpp"

int main()
{
    double freq = 10e6;
    double c0 = 299792458.0;
    double k = 2.0 * M_PI * freq / c0;

    double radius = 0.001;
    double srm = 2.0 * radius;

    Geometry geom = Geometry::standardDipoleMININEC(
        freq,
        11,
        radius
        );

    MininecKernel kernel(k, srm);
    ImpedanceAssembler assembler(kernel);

    auto Z = assembler.build(geom);

    // välj mittsegment
    int feedSeg = geom.segments.size() / 2;

    // här kan du senare lösa MoM-systemet
    // eller inspektera Z(feedSeg, feedSeg)

    return 0;
}
