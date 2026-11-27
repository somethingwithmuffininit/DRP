#include "mixing.hpp"
#include <scs.h>


// Variables are edge weights w and the epigraph scalar s. With B_e the
// Laplacian of one edge, P = I - sum(w_e B_e). SCS expects b - A*x
// in the cone and packs lower triangles columnwise with sqrt(2) off diagonal.
namespace mix {
Kernel optimize_weights(const Graph &g) {
    require(!g.directed(), "FMMC requires undirected support");
    int n = g.size(), m = int(g.edges().size());
    int linear = m + n + 1, triangle = n * (n + 1) / 2, rows = linear + 2 * triangle;
    Vector b = Vector::Zero(rows), c = Vector::Zero(m + 1);
    c[m] = 1;
    std::vector<Entry> terms;
    b.segment(m, n).setOnes();
    for (int e = 0; e < m; ++e) {
        auto [i, j] = g.edges()[e];
        terms.emplace_back(e, e, -1);
        terms.emplace_back(m + i, e, 1);
        terms.emplace_back(m + j, e, 1);
    }
    terms.emplace_back(m + n, m, -1);
    for (int block = 0; block < 2; ++block) {
        const int offset = linear + block * triangle;
        auto row = [&](int i, int j) { return offset + j * n - j * (j - 1) / 2 + i - j; }; // i >= j
        double sign = block == 0 ? 1. : -1.;
        for (int j = 0; j < n; ++j) {
            terms.emplace_back(row(j, j), m, -1);
            b[row(j, j)] = sign * (1. - 1. / n);
            for (int i = j + 1; i < n; ++i)
                b[row(i, j)] = -sign * std::sqrt(2.) / n;
        }
        for (int e = 0; e < m; ++e) {
            auto [i, j] = g.edges()[e];
            terms.emplace_back(row(i, i), e, sign);
            terms.emplace_back(row(j, j), e, sign);
            terms.emplace_back(row(j, i), e, -sign * std::sqrt(2.));
        }
    }
    Sparse a(rows, m + 1);
    a.setFromTriplets(terms.begin(), terms.end());
    a.makeCompressed();
    ScsMatrix matrix{a.valuePtr(), a.innerIndexPtr(), a.outerIndexPtr(), rows, m + 1};
    ScsData data{};
    data.m = rows;
    data.n = m + 1;
    data.A = &matrix;
    data.b = b.data();
    data.c = c.data();
    scs_int sizes[] = {n, n};
    ScsCone cone{};
    cone.l = linear;
    cone.s = sizes;
    cone.ssize = 2;
    ScsSettings settings;
    scs_set_default_settings(&settings);
    settings.eps_abs = settings.eps_rel = 1e-7;
    settings.verbose = 0;
    Vector x = Vector::Zero(m + 1), y = Vector::Zero(rows), slack = Vector::Zero(rows);
    ScsSolution solution{x.data(), y.data(), slack.data()};
    ScsInfo info{};
    int status = scs(&data, &cone, &settings, &solution, &info);
    auto verified = [&] {
        return (status == SCS_SOLVED || status == SCS_SOLVED_INACCURATE) && x.allFinite() &&
               std::abs(info.res_pri) < 1e-5 && std::abs(info.res_dual) < 1e-5 &&
               std::abs(info.gap) < 1e-5;
    };
    int initial_iterations = 0;
    if (!verified() && status == SCS_SOLVED_INACCURATE &&
        x.allFinite() && y.allFinite() && slack.allFinite()) {
        // Near-degenerate modes on chorded cycles can exhaust the initial budget.
        // Continue from the primal/dual iterate and retain the same acceptance checks.
        initial_iterations = info.iter;
        settings.warm_start = 1;
        settings.max_iters *= 4;
        status = scs(&data, &cone, &settings, &solution, &info);
    }
    require(verified(), std::string("unverified SCS solution: ") + info.status +
            " primal=" + std::to_string(info.res_pri) + " dual=" + std::to_string(info.res_dual) +
            " gap=" + std::to_string(info.gap));
    // Make a usable stochastic kernel from solver roundoff and report the repair.
    Vector w = x.head(m).cwiseMax(0.);
    Vector sums = Vector::Zero(n);
    for (int e = 0; e < m; ++e) {
        auto [i, j] = g.edges()[e];
        sums[i] += w[e];
        sums[j] += w[e];
    }
    w /= std::max(1., sums.maxCoeff());
    double correction = (w - x.head(m)).cwiseAbs().maxCoeff();
    require(correction <= 1e-5, "excessive FMMC feasibility correction");
    auto result = weighted_kernel(g, w);
    double achieved = symmetric_slem(result.p);
    require(std::abs(achieved - x[m]) < 1e-5, "SCS objective disagrees with attained SLEM");
    result.diagnostics = {{"solver", "SCS"},
                          {"status", info.status},
                          {"iterations", initial_iterations + info.iter},
                          {"warm_start_retry", initial_iterations > 0},
                          {"objective", x[m]},
                          {"achieved_slem", achieved},
                          {"weight_correction", correction},
                          {"primal_residual", info.res_pri},
                          {"dual_residual", info.res_dual},
                          {"gap", info.gap}};
    return result;
}
} // namespace mix
