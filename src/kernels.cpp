#include "mixing.hpp"
#include <algorithm>
#include <numeric>

namespace mix {
std::vector<Edit> select_edges(const Graph &g, const Json &spec, const Spectrum &base,
                               std::mt19937 &rng) {
    const std::string mode = spec.at("policy");
    const int n = g.size();
    if (mode == "none")
        return {{mode, {}}};
    require(!g.directed(), "edge policies require undirected support");
    std::vector<Edit> result;
    bool anchored = spec.value("candidates", std::string("all")) == "anchored";
    require(anchored || spec.value("candidates", std::string("all")) == "all", "unknown candidate set");
    require(!anchored || (g.order() == "cycle" && n >= 4),
            "anchored candidates require an ordered cycle");
    auto eligible = [&](int i, int j) {
        return i >= 0 && i < j && j < n && !g.contains(i, j) &&
               (!anchored || (j == n - 1 && i >= 1 && i <= n - 3));
    };
    if (mode == "fixed" || mode == "halves") {
        if (mode == "halves" && n % 2)
            return {};
        Edge e = mode == "halves" ? Edge{n / 2 - 1, n - 1} : spec.at("edge").get<Edge>();
        if (e.first > e.second)
            std::swap(e.first, e.second);
        require(eligible(e.first, e.second), "fixed edge is ineligible");
        return {{mode, e}};
    }
    const int repeats = spec.value("samples", 1);
    require(repeats > 0, "invalid policy sample count");
    if (anchored && mode == "random" && g.edges().size() == std::size_t(n)) {
        std::uniform_int_distribution<int> k(2, n - 2);
        for (int s = 0; s < repeats; ++s)
            result.push_back({mode, Edge{k(rng) - 1, n - 1}});
        return result;
    }
    auto pick = [&](const auto &a) {
        return a[std::uniform_int_distribution<std::size_t>(0, a.size() - 1)(rng)];
    };
    // Uniform rejection has constant expected cost on sparse unrestricted graphs.
    // This avoids materializing the quadratic complement of a large real graph.
    if (mode == "random" && !anchored && 4LL * g.edges().size() <= 1LL * n * (n - 1)) {
        std::uniform_int_distribution<int> vertex(0, n - 1);
        while (int(result.size()) < repeats) {
            int i = vertex(rng), j = vertex(rng);
            if (i > j)
                std::swap(i, j);
            if (eligible(i, j))
                result.push_back({mode, Edge{i, j}});
        }
        return result;
    }
    if (mode == "min_degree") {
        std::vector<int> endpoints;
        int lowest = n;
        for (int i = 0; i < n; ++i) {
            if (g.degrees()[i] > lowest || g.degrees()[i] == n - 1)
                continue;
            if (anchored && i != n - 1 && !eligible(i, n - 1))
                continue;
            if (g.degrees()[i] < lowest)
                endpoints.clear();
            lowest = g.degrees()[i];
            endpoints.push_back(i);
        }
        require(!endpoints.empty(), "no eligible missing edges");
        for (int s = 0; s < repeats; ++s) {
            int i = pick(endpoints);
            std::vector<int> partners;
            for (int j = 0; j < n; ++j)
                if (eligible(std::min(i, j), std::max(i, j)))
                    partners.push_back(j);
            require(!partners.empty(), "no eligible partner");
            int j = pick(partners);
            result.push_back({mode, Edge{std::min(i, j), std::max(i, j)}});
        }
        return result;
    }
    require(mode == "random" || mode == "all" || mode == "fiedler" || mode == "optimal_connectivity",
            "unknown edge policy: " + mode);
    double best = -1;
    std::optional<Edge> selected;
    for (int i = 0; i < n; ++i)
        for (int j = anchored ? n - 1 : i + 1; j < n; ++j) {
            if (!eligible(i, j))
                continue;
            if (mode == "all" || mode == "random")
                result.push_back({mode, Edge{i, j}});
            else {
                double score = mode == "fiedler" ? std::abs(base.fiedler[i] - base.fiedler[j])
                                                 : connectivity(g.edited(Edge{i, j})).second;
                if (score > best + 1e-12) {
                    best = score;
                    selected = {i, j};
                }
            }
        }
    if (selected)
        result.push_back({mode, selected});
    require(!result.empty(), "no eligible missing edges");
    if (mode == "random") {
        std::vector<Edit> samples;
        for (int s = 0; s < repeats; ++s)
            samples.push_back(pick(result));
        return samples;
    }
    return result;
}

void Kernel::validate() const {
    int n = p.rows();
    require(n > 0 && p.cols() == n && target.size() == n && position.size() == std::size_t(n) &&
                positions > 0,
            "invalid kernel dimensions");
    require(target.allFinite() && target.minCoeff() > 0 && std::abs(target.sum() - 1) < 1e-10,
            "invalid target");
    for (int j = 0; j < p.outerSize(); ++j)
        for (Sparse::InnerIterator it(p, j); it; ++it)
            require(std::isfinite(it.value()) && it.value() >= 0, "invalid transition probability");
    for (int v : position)
        require(v >= 0 && v < positions, "invalid projection");
    require((p * Vector::Ones(n) - Vector::Ones(n)).lpNorm<Eigen::Infinity>() < 1e-11,
            "invalid row sums");
    require((p.transpose() * target - target).lpNorm<Eigen::Infinity>() < 1e-11,
            "target is not stationary");
    if (reversible) {
        for (int j = 0; j < p.outerSize(); ++j)
            for (Sparse::InnerIterator it(p, j); it; ++it)
                require(std::abs(target[it.row()] * it.value() - target[j] * p.coeff(j, it.row())) <
                            1e-11,
                        "reversible kernel violates detailed balance");
    }
}
static Kernel identity_kernel(int n, const std::vector<Entry> &entries, bool symmetric) {
    Kernel k{Sparse(n, n), Vector::Constant(n, 1. / n), std::vector<int>(n), n, symmetric};
    std::iota(k.position.begin(), k.position.end(), 0);
    k.p.setFromTriplets(entries.begin(), entries.end());
    k.p.prune(0.0);
    k.validate();
    return k;
}
Kernel weighted_kernel(const Graph &g, const Vector &weights) {
    require(!g.directed() && weights.size() == int(g.edges().size()), "invalid symmetric weight input");
    Vector sums = Vector::Zero(g.size());
    std::vector<Entry> entries;
    for (int e = 0; e < weights.size(); ++e) {
        auto [i, j] = g.edges()[e];
        double w = weights[e];
        require(w >= 0 && std::isfinite(w), "invalid edge weight");
        sums[i] += w;
        sums[j] += w;
        entries.emplace_back(i, j, w);
        entries.emplace_back(j, i, w);
    }
    for (int i = 0; i < g.size(); ++i) {
        require(sums[i] <= 1. + 1e-14, "edge weights exceed row budget");
        double holding = 1. - sums[i];
        entries.emplace_back(i, i, holding > 1e-14 ? holding : 0.);
    }
    return identity_kernel(g.size(), entries, true);
}
Kernel build_chain(const std::string &method, const Graph &base, const Graph &edited,
                   std::optional<Edge> edge, double dhn_c) {
    int n = base.size();
    if (method == "dhn") {
        require(!edge && !base.directed() && base.order() == "path" && int(base.edges().size()) == n - 1,
                "DHN requires an unedited, ordered undirected path");
        require(dhn_c > 0 && dhn_c < n, "DHN requires 0 < c < n");
        int states = 2 * n;
        Kernel chain{Sparse(states, states), Vector::Constant(states, 1. / states), {}, n, false};
        double flip = dhn_c / n;
        std::vector<Entry> entries;
        for (int x = 0; x < states; ++x) {
            entries.emplace_back(x, (x + 1) % states, 1. - flip);
            entries.emplace_back(x, (states - x) % states, flip);
            // Dana stores signed labels modulo 2n. States 0 and +1 both observe position 1.
            chain.position.push_back(x == 0 ? 0 : (x <= n ? x - 1 : states - x));
        }
        chain.p.setFromTriplets(entries.begin(), entries.end());
        chain.validate();
        return chain;
    }
    if (method == "clockwise") {
        require(base.order() == "cycle" && n >= 4 && int(base.edges().size()) == n,
                "clockwise requires a simple ordered cycle");
        std::vector<Entry> entries;
        for (int i = 0; i < n; ++i) {
            bool hub = edge && (i == edge->first || i == edge->second);
            entries.emplace_back(i, (i + 1) % n, 0.5);
            entries.emplace_back(i, i, hub ? 0.25 : 0.5);
            if (hub)
                entries.emplace_back(i, i == edge->first ? edge->second : edge->first, 0.25);
        }
        return identity_kernel(n, entries, false);
    }
    require(!edited.directed(), "reversible construction requires undirected support");
    if (method == "fmmc")
        return optimize_weights(edited);
    require(method == "max_degree" || method == "metropolis", "unknown chain method: " + method);
    auto &degree = edited.degrees();
    int maximum = *std::max_element(degree.begin(), degree.end());
    Vector weights(edited.edges().size());
    for (int e = 0; e < weights.size(); ++e) {
        auto [i, j] = edited.edges()[e];
        weights[e] = 1. / (method == "max_degree" ? maximum : std::max(degree[i], degree[j]));
    }
    return weighted_kernel(edited, weights);
}
} // namespace mix
