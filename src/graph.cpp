#include "mixing.hpp"
#include <algorithm>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <openssl/evp.h>
#include <set>
#include <sstream>

namespace mix {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::invalid_argument(message);
}
Graph::Graph(int n, std::vector<Edge> edges, bool directed, std::string order)
    : directed_(directed), order_(std::move(order)), edges_(std::move(edges)) {
    require(n >= 3, "graphs need at least three vertices");
    require(order_.empty() || order_ == "cycle" || order_ == "path", "invalid vertex ordering");
    degree_.resize(n);
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](int i) {
        while (i != parent[i]) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    for (auto &[i, j] : edges_) {
        require(i >= 0 && i < n && j >= 0 && j < n && i != j, "invalid edge endpoint or loop");
        if (!directed && i > j)
            std::swap(i, j);
        ++degree_[i];
        if (!directed)
            ++degree_[j];
        parent[root(i)] = root(j);
    }
    std::sort(edges_.begin(), edges_.end());
    require(std::adjacent_find(edges_.begin(), edges_.end()) == edges_.end(), "duplicate edge");
    for (int i = 1; i < n; ++i)
        require(root(i) == root(0), "disconnected graph");
    if (!order_.empty()) {
        for (int i = 0; i < n - 1; ++i)
            require(contains(i, i + 1), "ordering requires consecutive links");
        if (order_ == "cycle")
            require(contains(n - 1, 0), "missing closing link");
    }
}
bool Graph::contains(int i, int j) const {
    if (!directed_ && i > j)
        std::swap(i, j);
    return std::binary_search(edges_.begin(), edges_.end(), Edge{i, j});
}
Graph Graph::edited(std::optional<Edge> edge) const {
    if (!edge)
        return *this;
    auto edges = edges_;
    edges.push_back(*edge);
    return Graph(size(), std::move(edges), directed_, order_);
}
Sparse Graph::laplacian() const {
    require(!directed_, "combinatorial connectivity requires an undirected support");
    std::vector<Entry> terms;
    for (int i = 0; i < size(); ++i)
        terms.emplace_back(i, i, degree_[i]);
    for (auto [i, j] : edges_) {
        terms.emplace_back(i, j, -1);
        terms.emplace_back(j, i, -1);
    }
    Sparse l(size(), size());
    l.setFromTriplets(terms.begin(), terms.end());
    return l;
}
Graph import_edges(const std::filesystem::path &file, int n, bool directed, const std::string &order,
                   std::string &sha256) {
    std::ifstream input(file, std::ios::binary);
    require(bool(input), "cannot open " + file.string());
    std::string bytes{std::istreambuf_iterator<char>(input), {}};
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned length = 0;
    require(EVP_Digest(bytes.data(), bytes.size(), digest, &length, EVP_sha256(), nullptr) == 1,
            "SHA-256 failed");
    std::ostringstream hex;
    for (unsigned i = 0; i < length; ++i)
        hex << std::hex << std::setw(2) << std::setfill('0') << int(digest[i]);
    sha256 = hex.str();
    std::istringstream in(bytes);
    std::vector<Edge> edges;
    for (std::string line; std::getline(in, line);) {
        line = line.substr(0, line.find('#'));
        std::istringstream row(line);
        if ((row >> std::ws).eof())
            continue;
        int i, j;
        require(bool(row >> i >> j) && (row >> std::ws).eof(),
                "expected two integer endpoints in " + file.string());
        edges.emplace_back(i, j);
    }
    return Graph(n, std::move(edges), directed, order);
}
Graph generate(const std::string &family, int n, const Json &opts, std::mt19937 &rng) {
    require(n >= 3, "invalid graph size");
    std::vector<Edge> edges;
    if (family == "cycle" || family == "path") {
        for (int i = 1; i < n; ++i)
            edges.emplace_back(i - 1, i);
        if (family == "cycle")
            edges.emplace_back(0, n - 1);
        return Graph(n, std::move(edges), false, family);
    }
    if (family == "bipartite" || family == "barbell") {
        int split = opts.value("split", n / 2);
        require(split > 0 && split < n, "invalid partition size");
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                bool crosses = i < split && j >= split;
                if (crosses == (family == "bipartite"))
                    edges.emplace_back(i, j);
            }
        if (family == "barbell") {
            int bridges = opts.value("bridges", 1);
            require(bridges > 0 && bridges <= std::min(split, n - split),
                    "bridge count exceeds the available disjoint endpoints");
            for (int i = 0; i < bridges; ++i)
                edges.emplace_back(split - 1 - i, split + i);
        }
    } else if (family == "tree") {
        int k = opts.value("branching", 2);
        require(k >= 1, "invalid branching");
        for (int i = 1; i < n; ++i)
            edges.emplace_back((i - 1) / k, i);
    } else if (family == "er") {
        double p = opts.value("p", 0.6);
        require(p > 0 && p <= 1, "invalid ER probability");
        std::bernoulli_distribution link(p);
        for (int attempt = 0; attempt < 100; ++attempt) {
            edges.clear();
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    if (link(rng))
                        edges.emplace_back(i, j);
            try {
                return Graph(n, edges);
            } catch (const std::invalid_argument &) {
            }
        }
        throw std::runtime_error("failed to draw a connected ER graph in 100 attempts");
    } else if (family == "ba") {
        int m = opts.value("m", 3);
        require(m >= 1 && m < n, "invalid BA attachment count");
        std::vector<int> urn;
        auto add = [&](int i, int j) {
            edges.emplace_back(i, j);
            urn.push_back(i);
            urn.push_back(j);
        };
        for (int i = 0; i <= m; ++i)
            for (int j = i + 1; j <= m; ++j)
                add(i, j);
        for (int j = m + 1; j < n; ++j) {
            std::set<int> neighbors;
            std::uniform_int_distribution<std::size_t> pick(0, urn.size() - 1);
            while (int(neighbors.size()) < m)
                neighbors.insert(urn[pick(rng)]);
            for (int i : neighbors)
                add(i, j);
        }
    } else
        throw std::invalid_argument("unknown graph family: " + family);
    return Graph(n, std::move(edges));
}
} // namespace mix
