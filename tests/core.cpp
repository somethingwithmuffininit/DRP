#include "mixing.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

using namespace mix;
void near(double x, double y, double tolerance = 1e-11) {
    require(std::abs(x - y) <= tolerance, "mismatch: " + std::to_string(x) + " vs " + std::to_string(y));
}
void test_projected_first_crossings() {
    // Observe the two halves of a lazy directed four-cycle. The two starts
    // reach uniform observed positions at different times, then move away.
    Kernel chain{Sparse(4, 4), Vector::Constant(4, .25), {0, 0, 1, 1}, 2, false};
    for (int i = 0; i < 4; ++i) {
        chain.p.insert(i, i) = .5;
        chain.p.insert(i, (i + 1) % 4) = .5;
    }
    Options opts;
    opts.eps = {.1};
    opts.diagnostic = 1;
    opts.stride = 1;
    opts.max_steps = 6;
    auto a = analyze(chain, opts, {0, 1});
    require(a["per_start"][0]["first_crossings"][1] == 3 &&
            a["per_start"][1]["first_crossings"][1] == 1, "individual projected crossing");
    near(a["curves"][3]["diagnostic"], 0);   // position TV at t=1
    near(a["curves"][5]["diagnostic"], .25); // position TV at t=2
    require(a["crossings"][1]["steps"] == 6, "simultaneous projected crossing");
}
int main() {
    try {
        test_projected_first_crossings();
        std::mt19937 rng(42);
        auto bipartite = generate("bipartite", 5, Json::object(), rng);
        auto basic = build_chain("max_degree", bipartite, bipartite, {});
        near(symmetric_slem(basic.p), 2. / 3);
        auto optimized = optimize_weights(bipartite);
        near(symmetric_slem(optimized.p), 3. / 7, 1e-5);
        for (auto [i, j] : bipartite.edges())
            near(optimized.p.coeff(i, j), 2. / 7, 1e-5);
        // Port Caroline's independent dense reference and all-start relaxation check.
        for (int n = 4; n <= 12; ++n) {
            auto g = generate("cycle", n, Json::object(), rng);
            near(connectivity(g).second, 2. - 2. * std::cos(2. * std::acos(-1.) / n));
            for (int k = 1; k <= n - 2; ++k) {
                std::optional<Edge> chord;
                if (k > 1)
                    chord = Edge{k - 1, n - 1};
                auto edited = g.edited(chord);
                near(connectivity(edited).second, connectivity(g).second);
                auto chain = build_chain("clockwise", g, edited, chord);
                Eigen::MatrixXd reference = Eigen::MatrixXd::Zero(n, n);
                for (int j = 0; j < n; ++j) {
                    reference((j + n - 1) % n, j) = 0.5;
                    reference(j, j) = k > 1 && (j == k - 1 || j == n - 1) ? 0.25 : 0.5;
                }
                if (chord)
                    reference(k - 1, n - 1) = reference(n - 1, k - 1) = 0.25;
                for (int start = 0; start < n; ++start) {
                    Vector p = Vector::Unit(n, start), q = p;
                    for (int t = 0; t <= 40; ++t) {
                        near((p - q).norm(), 0);
                        near(.5 * (p - chain.target).lpNorm<1>(), (q - chain.target).cwiseMax(0.).sum());
                        p = (chain.p.transpose() * p).eval();
                        q = (reference.transpose() * q).eval();
                    }
                }
            }
        }
        // DHN's uniform path, using the one-based cyclic labeling of the 2n states.
        int n = 5;
        Kernel lifted{Sparse(2 * n, 2 * n), Vector::Constant(2 * n, 1. / (2 * n)), {}, n, false};
        std::vector<Entry> rows;
        for (int j = 1; j <= 2 * n; ++j) {
            lifted.position.push_back(j <= n ? j - 1 : 2 * n - j);
            rows.emplace_back(j - 1, j % (2 * n), 1. - 1. / n);
            rows.emplace_back(j - 1, (2 * n - j + 2 * n - 1) % (2 * n), 1. / n);
        }
        lifted.p.setFromTriplets(rows.begin(), rows.end());
        lifted.validate();
        near(lifted.p.coeff(n - 1, n), 0.8);
        near(lifted.p.coeff(n - 1, n - 1), 0.2);
        // Dana stores label x at x mod 2n. The paper fixture stores it at (x-1) mod 2n.
        auto path = generate("path", n, {}, rng);
        auto dana = build_chain("dhn", path, path, {});
        for (int i = 0; i < 2*n; ++i) {
            int x = (i+1) % (2*n);
            require(dana.position[x] == lifted.position[i], "DHN position-label mismatch");
            for (int j = 0; j < 2*n; ++j)
                near(dana.p.coeff(x, (j+1) % (2*n)), lifted.p.coeff(i, j));
        }
        auto tuned = build_chain("dhn", path, path, {}, 0.5);
        near(tuned.p.coeff(0, 0), 0.5/n);
        near(tuned.p.coeff(0, 1), 1.-0.5/n);
        auto starts = choose_starts(2 * n, 256, 32, 0, rng);
        Options o;
        o.max_steps = 40;
        o.stride = 1;
        auto measurement = analyze(dana, o, starts);
        require(measurement["curves"].size() == 82, "projected horizon");
        auto &curves = measurement["curves"];
        for (int t = 0; t <= 40; ++t) {
            require(curves[2 * t + 1]["maximum"].get<double>() <=
                        curves[2 * t]["maximum"].get<double>() + 1e-12,
                    "projection bound");
            if (t)
                require(curves[2 * t + 1]["maximum"].get<double>() <=
                            curves[2 * t - 1]["maximum"].get<double>() + 1e-12,
                        "worst projected TV contraction");
        }
        near(curves[1]["maximum"], 1. - 1. / n);
        auto g = Graph(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}});
        for (auto mode : {"random", "min_degree", "fiedler", "all", "optimal_connectivity"}) {
            auto edit = select_edges(g, {{"policy", mode}}, connectivity(g), rng);
            require(edit.size() == 1 && edit[0].edge == Edge{2, 3}, "missing-edge oracle");
        }
        auto cycle = generate("cycle", 32, {}, rng);
        auto ids = choose_starts(32, 256, 32, 0, rng);
        o.max_steps = 195;
        o.snapshots = {};
        auto check = [&](std::optional<Edge> chord) {
            return analyze(build_chain("clockwise", cycle, cycle.edited(chord), chord), o,
                           ids)["crossings"][0];
        };
        near(check({})["steps"], 195);
        near(check(Edge{15, 31})["steps"], 61);
        o.max_steps = 194;
        require(check({})["status"] == "unresolved", "crossing horizon");
        o.max_steps = 0;
        o.eps = {0.99};
        near(check({})["steps"], 0);
        std::cout << "Core checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
