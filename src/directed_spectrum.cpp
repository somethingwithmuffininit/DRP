#include "mixing.hpp"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

namespace mix {
std::optional<TransitionSpectrum> sparse_directed_spectrum(const Sparse &p) {
    // Only materialize the Arnoldi basis for large, nonsymmetric kernels.
    Spectra::SparseGenMatProd<double> op(p);
    Spectra::GenEigsSolver<decltype(op)> solver(op, 4, 64);
    solver.init();
    if (solver.compute(Spectra::SortRule::LargestMagn, 1000, 1e-11) != 4)
        return {};
    auto values = solver.eigenvalues();
    auto vectors = solver.eigenvectors();
    Eigen::Index stationary, dominant;
    if ((values.array() - 1.).abs().minCoeff(&stationary) > 1e-8)
        return {};
    values[stationary] = 0;
    values.cwiseAbs().maxCoeff(&dominant);
    auto v = vectors.col(dominant);
    if ((p * v - values[dominant] * v).norm() > 1e-8 * v.norm())
        return {};
    return TransitionSpectrum{{values[dominant].real(), std::abs(values[dominant].imag())}, {}, {}};
}
} // namespace mix
