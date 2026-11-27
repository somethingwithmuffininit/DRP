#include "mixing.hpp"
#include "sample_transition.hpp"
#include <algorithm>
#include <numeric>
#include <queue>

namespace mix {
std::vector<int> choose_starts(int n, int limit, int count, int diagnostic, std::mt19937 &rng) {
    require(limit >= 0 && count > 0 && diagnostic >= 0 && diagnostic < n, "invalid start options");
    std::vector<int> ids(n);
    std::iota(ids.begin(), ids.end(), 0);
    if (n > limit && count < n) {
        std::swap(ids.front(), ids[diagnostic]);
        std::shuffle(ids.begin() + 1, ids.end(), rng);
        ids.resize(count);
        std::sort(ids.begin(), ids.end());
    }
    return ids;
}
namespace {
// A positive stationary distribution plus reachability from every vertex to
// vertex zero implies irreducibility. Level differences give the period.
int period(const Sparse &p) {
    std::vector<int> depth(p.rows(), -1);
    std::queue<int> pending;
    pending.push(0);
    depth[0] = 0;
    int divisor = 0;
    while (!pending.empty()) {
        int v = pending.front();
        pending.pop();
        for (Sparse::InnerIterator e(p, v); e; ++e)
            if (e.value() > 0) {
                int u = e.row();
                if (depth[u] < 0) {
                    depth[u] = depth[v] + 1;
                    pending.push(u);
                }
                divisor = std::gcd(divisor, std::abs(depth[v] + 1 - depth[u]));
            }
    }
    require(std::find(depth.begin(), depth.end(), -1) == depth.end(), "reducible transition kernel");
    return divisor;
}
Vector project(const Kernel &k, const Vector &p) {
    Vector result = Vector::Zero(k.positions);
    for (int i = 0; i < p.size(); ++i)
        result[k.position[i]] += p[i];
    return result;
}
struct Walk {
    Vector p, next;
    std::vector<int> first;
    double mass_error = 0, minimum = 1;
};
} // namespace
Json analyze(const Kernel &k, const Options &o, const std::vector<int> &starts) {
    k.validate();
    require(!o.eps.empty() && o.max_steps >= 0 && o.stride > 0, "invalid analysis options");
    for (double e : o.eps)
        require(e > 0 && e < 1, "epsilon must be in (0,1)");
    for (int t : o.snapshots)
        require(t >= 0, "negative snapshot time");
    require(!starts.empty() && std::is_sorted(starts.begin(), starts.end()) &&
                std::adjacent_find(starts.begin(), starts.end()) == starts.end(),
            "starts must be sorted and unique");
    auto diagnostic = std::find(starts.begin(), starts.end(), o.diagnostic);
    require(diagnostic != starts.end(), "diagnostic start must be measured");
    int diag = int(diagnostic - starts.begin());
    bool lifted = k.positions != k.p.rows();
    for (int i = 0; !lifted && i < k.p.rows(); ++i)
        lifted = k.position[i] != i;
    int spaces = lifted ? 2 : 1, count = int(starts.size()), thresholds = int(o.eps.size());
    int d = period(k.p);
    auto impossible = [&](int space, double e) {
        return space == 0 && d > 1 && e < 1. - 1. / d - 1e-14;
    };
    std::vector<Walk> walkers;
    for (int i : starts) {
        require(i >= 0 && i < k.p.rows(), "invalid start");
        walkers.push_back({Vector::Unit(k.p.rows(), i), Vector(k.p.rows()),
                           std::vector<int>(spaces * thresholds, -1)});
    }
    Vector pos_target = project(k, k.target);
    Json curves = Json::array(), snaps = Json::array(), crossings = Json::array(),
         individual = Json::array();
    std::vector<int> first(spaces * thresholds, -1);
    int horizon = 0;
    for (int t : o.snapshots)
        horizon = std::max(horizon, std::min(t, o.max_steps));
    constexpr int block = 32;
    std::vector<double> values(count * block * spaces);
    std::vector<double> least_probability(count * block);
    auto value = [&](int s, int t, int space) -> double & {
        return values[(s * block + t) * spaces + space];
    };
    bool done = false;
    int terminal = 0;
    for (int t0 = 0; !done;) {
        int len = int(std::min<long long>(block, 1LL + o.max_steps - t0));
        std::vector<Vector> snapshots(len);
#pragma omp taskgroup
        {
#pragma omp taskloop grainsize(4) shared(walkers, values, least_probability, snapshots, value, k, o, starts, pos_target)
            for (int s = 0; s < count; ++s) {
                auto &w = walkers[s];
                for (int b = 0; b < len; ++b) {
                    int t = t0 + b;
                    if (t > 0) {
                        w.next.noalias() = k.p.transpose() * w.p;
                        w.p.swap(w.next);
                    }
                    value(s, b, 0) = 0.5 * (w.p - k.target).lpNorm<1>();
                    if (lifted)
                        value(s, b, 1) = 0.5 * (project(k, w.p) - pos_target).lpNorm<1>();
                    w.mass_error = std::max(w.mass_error, std::abs(w.p.sum() - 1));
                    least_probability[s * block + b] = w.p.minCoeff();
                    w.minimum = std::min(w.minimum, least_probability[s * block + b]);
                    if (s == diag)
                        snapshots[b] = w.p;
                }
            }
        }
        for (int b = 0; b < len && !done; ++b) {
            terminal = t0 + b;
            bool complete = true, new_crossing = false;
            std::vector<double> maxima(spaces, 0);
            std::vector<double> minima(spaces, 1);
            for (int a = 0; a < spaces; ++a) {
                for (int s = 0; s < count; ++s) {
                    double distance = value(s, b, a);
                    require(std::isfinite(distance), "nonfinite trajectory");
                    maxima[a] = std::max(maxima[a], distance);
                    minima[a] = std::min(minima[a], distance);
                    for (int e = 0; e < thresholds; ++e) {
                        int &hit = walkers[s].first[a * thresholds + e];
                        if (hit < 0 && distance <= o.eps[e]) {
                            hit = terminal;
                            if (s == diag)
                                new_crossing = true;
                        }
                    }
                }
                for (int e = 0; e < thresholds; ++e) {
                    int &hit = first[a * thresholds + e];
                    if (hit < 0 && maxima[a] <= o.eps[e]) {
                        hit = terminal;
                        new_crossing = true;
                    }
                    if (hit < 0 && !impossible(a, o.eps[e]))
                        complete = false;
                }
            }
            done = terminal == o.max_steps || (!lifted && complete && terminal >= horizon);
            if (done || new_crossing || terminal % o.stride == 0) {
                for (int a = 0; a < spaces; ++a)
                    curves.push_back({{"step", terminal},
                                      {"space", a ? "position" : "full"},
                                      {"maximum", maxima[a]},
                                      {"minimum", minima[a]},
                                      {"diagnostic", value(diag, b, a)}});
            }
            if (done ||
                std::find(o.snapshots.begin(), o.snapshots.end(), terminal) != o.snapshots.end()) {
                auto &p = snapshots[b];
                double smallest = 1;
                for (int s = 0; s < count; ++s)
                    smallest = std::min(smallest, least_probability[s * block + b]);
                snaps.push_back({{"step", terminal},
                                 {"min_probability_over_starts", smallest},
                                 {"distribution", std::vector<double>(p.data(), p.data() + p.size())}});
            }
        }
        if (!done)
            t0 += len;
    }
    double error = 0, minimum = 1;
    for (int s = 0; s < count; ++s) {
        individual.push_back({{"start", starts[s]}, {"first_crossings", walkers[s].first}});
        error = std::max(error, walkers[s].mass_error);
        minimum = std::min(minimum, walkers[s].minimum);
    }
    require(error < 1e-9 && minimum >= 0, "probability conservation failed");
    for (int a = 0; a < spaces; ++a)
        for (int e = 0; e < thresholds; ++e) {
            int hit = first[a * thresholds + e];
            crossings.push_back(
                {{"space", a ? "position" : "full"},
                 {"epsilon", o.eps[e]},
                 {"status",
                  hit >= 0 ? "resolved" : (impossible(a, o.eps[e]) ? "periodic" : "unresolved")},
                 {"steps", hit < 0 ? Json(nullptr) : Json(hit)},
                 {"checked_through", terminal}});
        }
    return {{"period", d},
            {"scope", count == k.p.rows() ? "all" : "sampled"},
            {"starts", starts},
            {"diagnostic_start", o.diagnostic},
            {"crossings", crossings},
            {"per_start", individual},
            {"curves", curves},
            {"snapshots", snaps},
            {"max_mass_error", error},
            {"minimum_probability", minimum}};
}

Json trajectories(const Kernel &k, int draws, int reps, unsigned seed) {
    require(draws >= 2 && reps >= 2, "MCMC needs at least two draws and replicates");
    auto begin = Clock::now();
    TransitionRows rows = k.p;
    int cutoff = k.positions / 2;
    double truth = project(k, k.target).head(cutoff).sum();
    std::vector<double> means(reps);
#pragma omp taskloop shared(means, rows, k) grainsize(1)
    for (int r = 0; r < reps; ++r) {
        std::seed_seq sequence{seed, 0x4d434d43u, unsigned(r)};
        std::mt19937 generator(sequence);
        std::discrete_distribution<int> initial(k.target.data(), k.target.data() + k.target.size());
        int state = initial(generator), hits = 0;
        std::uniform_real_distribution<double> uniform;
        for (int t = 0; t < draws; ++t) {
            hits += k.position[state] < cutoff;
            state = sample_transition(rows, state, uniform(generator));
        }
        means[r] = double(hits) / draws;
    }
    double elapsed = seconds(begin);
    double mean = std::accumulate(means.begin(), means.end(), 0.) / reps;
    double variance = 0, mse = 0;
    for (double x : means) {
        variance += (x - mean) * (x - mean);
        mse += (x - truth) * (x - truth);
    }
    variance /= reps - 1;
    Json ess = variance > 0 ? Json(reps * truth * (1 - truth) / variance) : Json(nullptr);
    return {{"observable", "position < floor(n/2)"},
            {"initialization", "stationary"},
            {"draws", draws},
            {"replicates", reps},
            {"true_mean", truth},
            {"estimate", mean},
            {"rmse", std::sqrt(mse / reps)},
            {"means", means},
            {"ess", ess},
            {"seconds", elapsed},
            {"ess_per_second", ess.is_null() ? Json(nullptr) : Json(ess.get<double>() / elapsed)}};
}
} // namespace mix
