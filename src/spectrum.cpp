#include "mixing.hpp"
#include <Eigen/Eigenvalues>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <algorithm>

namespace mix {
namespace {
Spectrum low_spectrum(const Sparse &l) {
    if (l.rows() <= 256) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver{Eigen::MatrixXd(l)};
        require(solver.info() == Eigen::Success, "dense Laplacian eigensolve failed");
        return {solver.eigenvalues()[1], solver.eigenvalues()[2], solver.eigenvectors().col(1)};
    }
    // A small negative shift makes the singular Laplacian safely invertible.
    Spectra::SparseSymShiftSolve<double> op(l);
    Spectra::SymEigsShiftSolver<decltype(op)> solver(op, 3, 24, -1e-7);
    solver.init();
    require(solver.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10,
                           Spectra::SortRule::SmallestAlge) == 3,
            "sparse Laplacian eigensolve failed");
    auto values = solver.eigenvalues();
    Vector v = solver.eigenvectors().col(1);
    require((l * v - values[1] * v).norm() < 1e-6, "poor Fiedler residual");
    return {values[1], values[2], v};
}
} // namespace
Spectrum connectivity(const Graph &graph) {
    if (!graph.directed())
        return low_spectrum(graph.laplacian());
    auto edges = graph.edges();
    for (auto &[i, j] : edges)
        if (i > j)
            std::swap(i, j);
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return low_spectrum(Graph(graph.size(), std::move(edges)).laplacian());
}
TransitionSpectrum symmetric_spectrum(const Sparse &p) {
    auto extrema = [](double second, double smallest) -> TransitionSpectrum {
        return {std::abs(smallest) > std::abs(second) ? smallest : second, second, smallest};
    };
    if (p.rows() <= 256) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Eigen::MatrixXd(p),
                                                              Eigen::EigenvaluesOnly);
        require(solver.info() == Eigen::Success, "kernel eigensolve failed");
        auto &values = solver.eigenvalues();
        return extrema(values[p.rows() - 2], values[0]);
    }
    Sparse l(p.rows(), p.cols());
    l.setIdentity();
    l -= p;
    auto low = low_spectrum(l);
    Spectra::SparseSymMatProd<double> op(l);
    Spectra::SymEigsSolver<decltype(op)> solver(op, 1, 24);
    solver.init();
    require(solver.compute(Spectra::SortRule::LargestAlge) == 1, "upper Laplacian eigensolve failed");
    return extrema(1 - low.second, 1 - solver.eigenvalues()[0]);
}
double symmetric_slem(const Sparse &p) {
    return std::abs(symmetric_spectrum(p).dominant);
}
namespace {
std::complex<double> dominant_mode(const Eigen::VectorXcd &values) {
    // Conjugates describe the same decay. Report the one with nonnegative imaginary part.
    Eigen::Index i;
    values.cwiseAbs().maxCoeff(&i);
    return {values[i].real(), std::abs(values[i].imag())};
}
} // namespace
TransitionSpectrum general_spectrum(const Sparse &p) {
    Eigen::EigenSolver<Eigen::MatrixXd> solver(Eigen::MatrixXd(p), false);
    require(solver.info() == Eigen::Success, "directed kernel eigensolve failed");
    auto values = solver.eigenvalues();
    Eigen::Index stationary;
    (values.array() - 1.).abs().minCoeff(&stationary);
    values[stationary] = 0;
    return {dominant_mode(values), {}, {}};
}
double general_slem(const Sparse &p) {
    return std::abs(general_spectrum(p).dominant);
}
TransitionSpectrum dhn_spectrum(int n, double c) {
    // Reflection interchanges Fourier modes k and -k on the 2n-cycle.
    // Each pair has trace 2a cos(theta) and determinant a*a - b*b.
    double b = c / n, a = 1 - b;
    Eigen::VectorXcd modes(2 * n - 1);
    modes[0] = b - a; // alternating mode; the constant mode was excluded
    for (int k = 1; k < n; ++k) {
        double theta = std::acos(-1.) * k / n;
        double sine = a * std::sin(theta);
        auto root = std::sqrt(std::complex<double>((b - sine) * (b + sine), 0));
        double center = a * std::cos(theta);
        modes[2*k-1] = center + root;
        modes[2*k] = center - root;
    }
    return {dominant_mode(modes), {}, {}};
}
} // namespace mix
