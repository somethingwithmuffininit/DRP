#include "mixing.hpp"
#include "sample_transition.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace mix {
Json budget_trial(const Graph &graph, const std::string &method, const Json &policy,
                  const std::vector<double> &budgets, unsigned seed, double c, int split) {
    require(!budgets.empty() &&
                std::all_of(budgets.begin(), budgets.end(), [](double t) { return t > 0 && std::isfinite(t); }) &&
                std::is_sorted(budgets.begin(), budgets.end()) &&
                std::adjacent_find(budgets.begin(), budgets.end()) == budgets.end(),
            "budgets must be finite, positive and strictly increasing");
    require(split > 0 && split < graph.size(), "invalid observable partition");
    require(policy.at("policy") != "all" && policy.value("samples", 1) == 1,
            "a budget trial selects one edit");
    // Both observables are defined on the shared unedited graph.
    auto &degrees = graph.degrees();
    int rare = int(std::min_element(degrees.begin(), degrees.end()) - degrees.begin());
    auto begin = Clock::now();
    std::mt19937 edit_rng(seed);
    Spectrum base{};
    if (policy.at("policy") == "fiedler")
        base = connectivity(graph);
    auto edits = select_edges(graph, policy, base, edit_rng);
    require(edits.size() == 1, "budget policy must produce exactly one edit");
    auto edge = edits.front().edge;
    auto edited = graph.edited(edge);
    auto k = build_chain(method, graph, edited, edge, c);
    TransitionRows rows = k.p;
    // Start at physical vertex zero, heading right on Dana's path.
    int state = method == "dhn" ? 1 : 0;
    require(k.position[state] == 0, "budget starts must share physical vertex zero");
    double partition_truth = 0, rare_truth = 0;
    for (int i = 0; i < k.target.size(); ++i) {
        if (k.position[i] < split)
            partition_truth += k.target[i];
        if (k.position[i] == rare)
            rare_truth += k.target[i];
    }
    std::seed_seq stream{seed, 0x4d434d43u};
    std::mt19937 sampler(stream);
    std::uniform_real_distribution<double> uniform;
    double setup = seconds(begin);
    Json records = Json::array();
    std::uint64_t draws = 0, partition_hits = 0, rare_hits = 0;
    auto record = [&](double budget) {
        records.push_back({{"budget_seconds", budget}, {"draws", draws},
                           {"status", draws ? "ok" : (setup >= budget ? "setup_exhausted" : "no_draws")},
                           {"partition_estimate", draws ? Json(double(partition_hits)/draws) : Json(nullptr)},
                           {"rare_estimate", draws ? Json(double(rare_hits)/draws) : Json(nullptr)}});
    };
    size_t pending = 0;
    while (pending < budgets.size() && budgets[pending] <= setup)
        record(budgets[pending++]);
    while (pending < budgets.size()) {
        state = sample_transition(rows, state, uniform(sampler));
        int pos = k.position[state];
        bool in_partition = pos < split, at_rare = pos == rare;
        double elapsed = seconds(begin);
        // A transition completed after a deadline counts only for later budgets.
        while (pending < budgets.size() && elapsed > budgets[pending])
            record(budgets[pending++]);
        ++draws;
        partition_hits += in_partition;
        rare_hits += at_rare;
    }
    return {{"method", method}, {"policy", policy.at("policy")}, {"seed", seed},
            {"edge", edge ? Json(*edge) : Json(nullptr)}, {"setup_seconds", setup},
            {"initial_state", method == "dhn" ? 1 : 0}, {"initial_position", 0},
            {"partition_end_exclusive", split}, {"rare_vertex", rare},
            {"partition_truth", partition_truth}, {"rare_truth", rare_truth}, {"checkpoints", records}};
}
} // namespace mix
