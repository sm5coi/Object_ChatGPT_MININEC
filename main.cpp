#include "geometry.hpp"
#include "mininec_kernel.hpp"
#include "impedance_assembler.hpp"
#include "impedance_matrix.hpp"
#include <iostream>


std::vector<Complex> solveLinearSystem(
    std::vector<std::vector<Complex>> A,
    std::vector<Complex> b)
{
    int N = b.size();

    for (int i = 0; i < N; ++i) {
        // pivot
        int pivot = i;
        for (int r = i + 1; r < N; ++r)
            if (std::abs(A[r][i]) > std::abs(A[pivot][i]))
                pivot = r;

        std::swap(A[i], A[pivot]);
        std::swap(b[i], b[pivot]);

        Complex diag = A[i][i];
        if (std::abs(diag) == 0.0)
            throw std::runtime_error("Singular matrix");

        // normalize row
        for (int j = i; j < N; ++j)
            A[i][j] /= diag;
        b[i] /= diag;

        // eliminate
        for (int r = 0; r < N; ++r) {
            if (r == i) continue;
            Complex f = A[r][i];
            for (int j = i; j < N; ++j)
                A[r][j] -= f * A[i][j];
            b[r] -= f * b[i];
        }
    }

    return b; // solution vector I
}





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

    int N = Z.size();

    std::cout << "\n geom.segments[i].wire = ";
    for (int i = 0; i <= N; ++i)
    {
        std::cout << geom.segments[i].wire << "  ";
    }
    std::cout << std::endl;

    // Kopiera Z → A (lokal, lösbar matris)
    std::vector<std::vector<Complex>> A(
        N, std::vector<Complex>(N));

    std::cout << "\n Z = ";
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            A[i][j] = Z(i, j);
            std::cout << Z(i,j) << "  ";
        }
        std::cout << std::endl;
    }
    std::vector<Complex> V(N, Complex(0.0, 0.0));
    V[feedSeg] = Complex(1.0, 0.0);

    auto I = solveLinearSystem(A, V);

    std::cout << "\n I = ";
    for (int i = 0; i <= N; ++i)
    {
        std::cout << I[i] << "  ";
    }
    std::cout << std::endl;

    // inimpedans
    Complex Z_in = Complex(1.0, 0.0) / I[feedSeg];

    std::cout << "\n Z_in = " << Z_in << std::endl;


    // här kan du senare lösa MoM-systemet
    // eller inspektera Z(feedSeg, feedSeg)

    return 0;
}
