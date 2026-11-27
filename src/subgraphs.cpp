#include "mixing.hpp"
#include <algorithm>
#include <numeric>

namespace mix {
std::vector<std::vector<int>> nested_orders(const Graph &parent, const std::string &kind, int reps,
                                            int max_size, std::mt19937 &rng) {
    int n = parent.size();
    require(!parent.directed() && max_size >= 3 && max_size <= n, "invalid subgraph parent or size");
    require(kind == "section" || kind == "connected", "unknown subgraph scheme");
    std::vector<int> roots(n);
    std::iota(roots.begin(), roots.end(), 0);
    require(reps > 0 && reps <= n, "nested replicates require distinct roots");
    std::shuffle(roots.begin(), roots.end(), rng);
    std::vector<std::vector<int>> orders(reps);
    if (kind == "section") {
        require(parent.order() == "cycle" && int(parent.edges().size()) == n && max_size < n,
                "sections require a proper segment of a simple ordered ring");
        for (int r = 0; r < reps; ++r)
            for (int j = 0; j < max_size; ++j)
                orders[r].push_back((roots[r] + j) % n);
        return orders;
    }
    std::vector<std::vector<int>> neighbors(n);
    for (auto [i, j] : parent.edges()) {
        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
    }
    for (int r = 0; r < reps; ++r) {
        std::vector<int> frontier{roots[r]};
        std::vector<bool> seen(n, false);
        seen[roots[r]] = true;
        while (int(orders[r].size()) < max_size) {
            require(!frontier.empty(), "subgraph exceeds connected component");
            auto index = std::uniform_int_distribution<size_t>(0, frontier.size()-1)(rng);
            int vertex = frontier[index];
            frontier[index] = frontier.back();
            frontier.pop_back();
            orders[r].push_back(vertex);
            for (int v : neighbors[vertex])
                if (!seen[v]) {
                    seen[v] = true;
                    frontier.push_back(v);
                }
        }
    }
    return orders;
}

Graph induced_graph(const Graph &parent, const std::vector<int> &ids, bool ordered_path) {
    std::vector<int> local(parent.size(), -1);
    for (int i = 0; i < int(ids.size()); ++i) {
        int v = ids[i];
        require(v >= 0 && v < parent.size() && local[v] == -1, "invalid parent vertex");
        local[v] = i;
    }
    std::vector<Edge> edges;
    for (auto [i, j] : parent.edges())
        if (local[i] >= 0 && local[j] >= 0)
            edges.emplace_back(local[i], local[j]);
    return Graph(int(ids.size()), std::move(edges), parent.directed(), ordered_path ? "path" : "");
}
} // namespace mix
