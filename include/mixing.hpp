#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <nlohmann/json.hpp>
#include <chrono>
#include <complex>
#include <filesystem>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace mix {
using Json = nlohmann::json;
using Vector = Eigen::VectorXd;
using Sparse = Eigen::SparseMatrix<double>;
using Entry = Eigen::Triplet<double>;
using Edge = std::pair<int, int>;
using Clock = std::chrono::steady_clock;
inline double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
void require(bool condition, const std::string &message);

class Graph {
  public:
    Graph(int n, std::vector<Edge> edges, bool directed = false, std::string order = "");
    int size() const { return int(degree_.size()); }
    bool directed() const { return directed_; }
    const auto &edges() const { return edges_; }
    const auto &degrees() const { return degree_; }
    const auto &order() const { return order_; } // IDs follow the declared order.
    bool contains(int i, int j) const;
    Graph edited(std::optional<Edge> edge) const;
    Sparse laplacian() const;

  private:
    bool directed_;
    std::string order_;
    std::vector<Edge> edges_;
    std::vector<int> degree_;
};
Graph generate(const std::string &family, int n, const Json &options, std::mt19937 &rng);
Graph import_edges(const std::filesystem::path &file, int n, bool directed, const std::string &order,
                   std::string &sha256);
std::vector<std::vector<int>> nested_orders(const Graph &parent, const std::string &kind, int replicates,
                                            int max_size, std::mt19937 &rng);
Graph induced_graph(const Graph &parent, const std::vector<int> &ids, bool ordered_path);
struct Spectrum {
    double second, third;
    Vector fiedler;
};
Spectrum connectivity(const Graph &graph);
struct TransitionSpectrum {
    std::complex<double> dominant; // representative nonstationary mode of largest modulus
    std::optional<double> second, smallest; // algebraic ordering for reversible chains
};
TransitionSpectrum symmetric_spectrum(const Sparse &p);
std::optional<TransitionSpectrum> sparse_directed_spectrum(const Sparse &p);
TransitionSpectrum general_spectrum(const Sparse &p);
TransitionSpectrum dhn_spectrum(int positions, double c);
double symmetric_slem(const Sparse &p);
double general_slem(const Sparse &p);
struct Edit {
    std::string policy;
    std::optional<Edge> edge;
};
std::vector<Edit> select_edges(const Graph &graph, const Json &policy, const Spectrum &base,
                               std::mt19937 &rng);

struct Kernel {
    Sparse p; // row-stochastic, p(i,j) is the probability of i -> j
    Vector target;
    std::vector<int> position;
    int positions;
    bool reversible;
    Json diagnostics = Json::object();
    void validate() const;
};
Kernel build_chain(const std::string &method, const Graph &base, const Graph &edited,
                   std::optional<Edge> edge, double dhn_c = 1.0);
Kernel weighted_kernel(const Graph &graph, const Vector &weights);
Kernel optimize_weights(const Graph &graph);
struct Options {
    std::vector<double> eps{0.25};
    std::vector<int> snapshots{0, 10, 100, 500};
    int max_steps = 100000, stride = 20, diagnostic = 0;
};
std::vector<int> choose_starts(int n, int all_limit, int sample_count, int diagnostic,
                               std::mt19937 &rng);
Json analyze(const Kernel &chain, const Options &options, const std::vector<int> &starts);
Json trajectories(const Kernel &chain, int draws, int reps, unsigned seed);
Json budget_trial(const Graph &graph, const std::string &method, const Json &policy,
                  const std::vector<double> &budgets, unsigned seed, double c, int split);
void run(const std::filesystem::path &config, const std::filesystem::path &output, int jobs);
} // namespace mix
