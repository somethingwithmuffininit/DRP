#include "mixing.hpp"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <memory>
#include <omp.h>

namespace mix {
namespace fs = std::filesystem;
namespace {
Json read_json(const fs::path &path) {
    std::ifstream in(path);
    require(bool(in), "cannot read " + path.string());
    return Json::parse(in);
}
void write_json(const fs::path &path, const Json &object) {
    auto staging = path;
    staging += ".partial";
    std::ofstream out;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out.open(staging);
    out << object.dump() << '\n';
    out.close();
    fs::rename(staging, path);
}
struct Instance {
    Graph graph;
    Spectrum spectrum;
    std::string id;
};
struct Case {
    std::shared_ptr<const Instance> instance;
    std::string method;
    Edit edit;
    unsigned sampling_seed;
    double selection_seconds;
    double dhn_c;
    int diagnostic;
};
void record_spectrum(Json &result, const Kernel &k, const Case &item) {
    auto start = Clock::now();
    TransitionSpectrum spectrum;
    if (k.reversible) {
        Sparse symmetric = k.p;
        for (int j = 0; j < symmetric.outerSize(); ++j)
            for (Sparse::InnerIterator e(symmetric, j); e; ++e)
                e.valueRef() *= std::sqrt(k.target[e.row()] / k.target[j]);
        spectrum = symmetric_spectrum(symmetric);
    } else if (item.method == "dhn") {
        spectrum = dhn_spectrum(k.positions, item.dhn_c);
    } else if (k.p.rows() <= 256) {
        spectrum = general_spectrum(k.p);
    } else {
        auto mode = sparse_directed_spectrum(k.p);
        if (!mode) {
            result["spectral_status"] = "eigensolve_unresolved";
            result["spectrum_seconds"] = seconds(start);
            return;
        }
        spectrum = *mode;
    }
    result["spectrum_seconds"] = seconds(start);
    double rho = std::abs(spectrum.dominant);
    require(std::isfinite(rho) && rho <= 1 + 1e-9, "invalid transition spectral radius");
    const bool periodic = result["analysis"]["period"].get<int>() > 1;
    rho = periodic ? 1. : std::min(1., rho);
    result["slem"] = rho;
    result["dominant_eigenvalue_real"] = spectrum.dominant.real();
    result["dominant_eigenvalue_imag"] = spectrum.dominant.imag();
    if (spectrum.second) {
        result["transition_lambda2"] = *spectrum.second;
        result["transition_lambda_min"] = *spectrum.smallest;
        result["spectral_gap"] = 1 - *spectrum.second;
    }
    result["absolute_spectral_gap"] = 1 - rho;
    // Infinite or numerically indistinguishable limits have an explicit status.
    if (rho == 1) {
        result["spectral_status"] = periodic ? "periodic" : "unit_within_precision";
        result["asymptotic_rate"] = 0.;
    } else if (rho <= 1e-15) {
        result["spectral_status"] = "zero_within_precision";
        result["spectral_time"] = 0.;
    } else {
        result["spectral_status"] = "convergent";
        double rate = -std::log(rho);
        result["asymptotic_rate"] = rate;
        result["spectral_time"] = 1 / rate;
    }
}
Json evaluate(const Case &item, const Json &settings, Options options) {
    options.diagnostic = item.diagnostic;
    auto &base = item.instance->graph;
    auto &spec = item.instance->spectrum;
    auto start = Clock::now();
    auto g = base.edited(item.edit.edge);
    Json result = {{"graph", item.instance->id},
                   {"method", item.method},
                   {"policy", item.edit.policy},
                   {"edge", item.edit.edge ? Json(*item.edit.edge) : Json(nullptr)},
                   {"vertices", g.size()},
                   {"lambda2_before", spec.second},
                   {"lambda3_before", spec.third},
                   {"start_seed", item.sampling_seed},
                   {"trajectory_seed", item.sampling_seed},
                   {"edge_selection_seconds", item.selection_seconds},
                   {"connectivity_space", "undirected_support"}};
    auto edited_spec = item.edit.edge ? connectivity(g) : spec;
    double lambda = edited_spec.second;
    result["lambda3_after"] = edited_spec.third;
    double delta = lambda - spec.second, gap = spec.third - spec.second;
    result["lambda2_after"] = lambda;
    result["delta_lambda2"] = delta;
    result["normalized_gain"] =
        gap > 1e-9 * std::max(1., spec.third) ? Json(delta / gap) : Json(nullptr);
    result["connectivity_seconds"] = seconds(start);
    if (item.method == "fmmc" && g.size() > settings.value("sdp_limit", 64)) {
        result["status"] = "skipped_sdp_limit";
        return result;
    }
    start = Clock::now();
    auto k = build_chain(item.method, base, g, item.edit.edge, item.dhn_c);
    if (item.method == "dhn")
        result["dhn_c"] = item.dhn_c;
    result["construction_seconds"] = seconds(start);
    result["setup_seconds"] = result["construction_seconds"].get<double>() + item.selection_seconds;
    result["solver"] = k.diagnostics;
    result["states"] = k.p.rows();
    result["positions"] = k.positions;
    result["stationary"] = std::vector<double>(k.target.data(), k.target.data() + k.target.size());
    result["state_to_position"] = k.position;
    result["reversible"] = k.reversible;
    result["status"] = "ok";
    if (item.method == "fmmc" || settings.value("export_kernel", false)) {
        Json entries = Json::array();
        for (int j = 0; j < k.p.outerSize(); ++j)
            for (Sparse::InnerIterator e(k.p, j); e; ++e)
                entries.push_back({e.row(), e.col(), e.value()});
        result["transitions"] = std::move(entries);
    }
    std::mt19937 generator(item.sampling_seed);
    auto starts = choose_starts(k.p.rows(), settings.value("all_starts_limit", 256),
                                settings.value("sample_starts", 32), options.diagnostic, generator);
    start = Clock::now();
    result["analysis"] = analyze(k, options, starts);
    result["propagation_seconds"] = seconds(start);
    record_spectrum(result, k, item);
    if (settings.contains("mcmc")) {
        auto m = settings.at("mcmc");
        result["mcmc"] = trajectories(k, m.at("draws"), m.at("replicates"), item.sampling_seed);
    }
    return result;
}
void run_budgets(const Json &cfg, const fs::path &output, int requested_jobs) {
    auto settings = cfg.at("wall_budget");
    auto budgets = settings.at("seconds").get<std::vector<double>>();
    int repeats = settings.at("repeats");
    int repetition = settings.value("only_repetition", -1);
    require(repetition >= -1 && repetition < repeats, "invalid budget repetition");
    require(repeats >= 2, "budget error estimates require repeated trials");
    require(fs::create_directories(output), "choose a fresh budget output directory");
    Json manifest = {{"kind", "wall_budget"}, {"status", "running"}, {"config", cfg},
                     {"jobs", 1}, {"requested_jobs", requested_jobs}, {"graphs", Json::array()}, {"run_count", 0}};
    write_json(output / "manifest.json", manifest);
    unsigned seed = cfg.value("seed", 42u);
    std::mt19937 graph_rng(seed), pairing(seed), order_rng(seed + 1u);
    int id = 0;
    try {
        for (const auto &study : cfg.at("studies")) {
            const auto &def = study.at("graph");
            for (int n : def.at("sizes")) {
                auto g = generate(def.at("family"), n, def, graph_rng);
                auto graph_id = "g" + std::to_string(manifest["graphs"].size());
                manifest["graphs"].push_back({{"id", graph_id}, {"name", study.at("name")},
                                              {"vertices", n}, {"generator", def}, {"edges", g.edges()}});
                std::vector<std::pair<std::string, Json>> cases;
                for (const auto &policy : study.at("policies"))
                    for (const auto &method : study.at("methods"))
                        cases.emplace_back(method.get<std::string>(), policy);
                int split = def.at("family") == "barbell" ? def.value("split", n / 2) : n / 2;
                for (int trial = 0; trial < repeats; ++trial) {
                    unsigned trial_seed = pairing();
                    std::shuffle(cases.begin(), cases.end(), order_rng);
                    // Advance both RNG streams even for blocks executed by another process.
                    if (repetition >= 0 && repetition != trial)
                        continue;
                    for (const auto &[method, policy] : cases) {
                        Json r;
                        if (method == "fmmc" && n > cfg.value("sdp_limit", 64)) {
                            r = {{"status", "skipped_sdp_limit"}, {"method", method}, {"policy", policy.at("policy")}};
                        } else {
                            r = budget_trial(g, method, policy, budgets, trial_seed, study.value("dhn_c", 1.), split);
                            r["status"] = "ok";
                        }
                        r["graph"] = graph_id;
                        r["trial"] = trial;
                        r["id"] = id;
                        write_json(output / ("run-" + std::to_string(id++) + ".json"), r);
                    }
                }
            }
        }
        require(id > 0, "empty budget experiment");
        manifest["status"] = "complete";
    } catch (const std::exception &e) {
        manifest["run_count"] = id;
        manifest["status"] = "failed";
        manifest["error"] = e.what();
        write_json(output / "manifest.json", manifest);
        throw;
    }
    manifest["run_count"] = id;
    write_json(output / "manifest.json", manifest);
    std::cout << id << " budget trials in " << output << '\n';
}
} // namespace
void run(const fs::path &config, const fs::path &output, int jobs) {
    require(jobs > 0, "jobs must be positive");
    Json cfg = read_json(config);
    if (cfg.contains("wall_budget")) {
        run_budgets(cfg, output, jobs);
        return;
    }
    Options opts;
    opts.eps = cfg.value("epsilons", opts.eps);
    opts.stride = cfg.value("stride", opts.stride);
    opts.max_steps = cfg.value("max_steps", opts.max_steps);
    opts.snapshots = cfg.value("snapshots", opts.snapshots);
    opts.diagnostic = cfg.value("diagnostic_start", 0);
    unsigned seed = cfg.value("seed", 42u);
    // Keep the paper's mt19937(seed) chord sequence. Generate random topologies
    // from a distinct stream so a graph draw cannot also determine its edit.
    std::seed_seq topology_seed{seed, 0x47524150u};
    std::mt19937 graph_rng(topology_seed), placement_rng(seed), start_rng(seed);
    std::vector<Case> cases;
    Json graphs = Json::array();
    for (const auto &study : cfg.at("studies")) {
        auto methods = study.at("methods").get<std::vector<std::string>>();
        require(!methods.empty(), "empty method list");
        auto def = study.at("graph");
        Json dataset = nullptr;
        std::optional<Graph> parent;
        Json fixtures = Json::array();
        const bool fixed_trials = def.contains("fixtures");
        if (fixed_trials) {
            require(!def.contains("dataset") && !def.contains("sizes") && !study.contains("replicates"),
                    "graph fixtures supply their own sizes and repetitions");
            auto file = config.parent_path() / def.at("fixtures").get<std::string>();
            auto records = read_json(file);
            require(records.is_array(), "graph fixtures must be an array of trials");
            for (auto &record : records)
                if (record.at("family") == def.at("family"))
                    fixtures.push_back(std::move(record));
            require(!fixtures.empty(), "no fixtures for the requested graph family");
        }
        if (def.contains("dataset")) {
            auto meta = config.parent_path() / def.at("dataset").get<std::string>();
            dataset = read_json(meta);
            std::string digest;
            parent = import_edges(meta.parent_path() / dataset.at("edges").get<std::string>(),
                                  dataset.at("vertices"), dataset.value("directed", false),
                                  dataset.value("order", std::string{}), digest);
            require(dataset.value("sha256", digest) == digest, "dataset checksum mismatch");
            require(dataset.value("edge_count", parent->edges().size()) == parent->edges().size(),
                    "dataset edge count mismatch");
            dataset["sha256"] = digest;
        }
        bool nested = def.contains("subgraph");
        require(!nested || parent.has_value(), "subgraph requires a parent dataset");
        auto sizes = fixed_trials ? std::vector<int>{fixtures.front().at("n").get<int>()}
                   : !parent || nested ? def.at("sizes").get<std::vector<int>>()
                                       : std::vector<int>{parent->size()};
        int reps = fixed_trials ? static_cast<int>(fixtures.size()) : study.value("replicates", 1);
        require(reps > 0 && !sizes.empty(), "invalid graph sizes or replicates");
        std::vector<std::vector<int>> orders;
        if (nested)
            orders = nested_orders(*parent, def.at("subgraph"), reps,
                                   *std::max_element(sizes.begin(), sizes.end()), graph_rng);
        for (int n : sizes)
            for (int rep = 0; rep < reps; ++rep) {
                const Json *fixture = fixed_trials ? &fixtures[rep] : nullptr;
                require(!fixture || fixture->at("n") == n, "fixture family mixes graph sizes");
                std::vector<int> ids;
                if (nested) {
                    require(n >= 3 && n <= int(orders[rep].size()), "invalid nested subgraph size");
                    ids.assign(orders[rep].begin(), orders[rep].begin() + n);
                }
                auto g = fixture ? Graph(n, fixture->at("edges").get<std::vector<Edge>>(),
                                         false, def.value("order", std::string{}))
                         : nested ? induced_graph(*parent, ids, def.at("subgraph") == "section")
                         : parent ? *parent : generate(def.at("family"), n, def, graph_rng);
                auto eig_start = Clock::now();
                auto base = connectivity(g);
                double eig_seconds = seconds(eig_start);
                auto id = "g" + std::to_string(graphs.size());
                auto instance = std::make_shared<Instance>(Instance{std::move(g), std::move(base), id});
                graphs.push_back({{"id", id},
                                  {"name", study.at("name")},
                                  {"vertices", n},
                                  {"replicate", rep},
                                  {"generator", def},
                                  {"source", dataset},
                                  {"edges", instance->graph.edges()}});
                if (nested)
                    graphs.back()["parent_vertex_ids"] = ids;
                unsigned starts = fixture ? fixture->at("start_seed").get<unsigned>() : start_rng();
                std::optional<std::mt19937> fixture_placements;
                if (fixture)
                    fixture_placements.emplace(fixture->at("policy_seed").get<unsigned>());
                auto &placements = fixture_placements ? *fixture_placements : placement_rng;
                for (const auto &policy : study.at("policies")) {
                    auto selection_start = Clock::now();
                    auto edits =
                        select_edges(instance->graph, policy, instance->spectrum, placements);
                    double selection = seconds(selection_start);
                    if (policy.at("policy") == "fiedler")
                        selection += eig_seconds;
                    if (!edits.empty())
                        selection /= edits.size();
                    for (auto &edit : edits)
                        for (auto &method : methods)
                            cases.push_back({instance, method, edit, starts, selection,
                                             study.value("dhn_c", 1.0),
                                             study.value("diagnostic_starts", Json::object()).value(method, opts.diagnostic)});
                }
            }
    }
    require(!cases.empty(), "empty experiment plan");
    require(fs::create_directories(output), "choose a fresh output directory: " + output.string());
    Json manifest = {{"version", 1}, {"status", "running"}, {"config", cfg},
                     {"jobs", jobs}, {"graphs", graphs},    {"run_count", cases.size()}};
    write_json(output / "manifest.json", manifest);
    std::vector<std::string> errors(cases.size());
    // One OpenMP team schedules independent cases and their child start tasks.
    // Running cases own their buffers and write disjoint files indexed by plan order.
    omp_set_dynamic(0);
#pragma omp parallel num_threads(jobs)
    {
#pragma omp single
        {
            for (std::size_t i = 0; i < cases.size(); ++i) {
#pragma omp task firstprivate(i) shared(cases, cfg, opts, errors, output)
                {
                    try {
                        auto r = evaluate(cases[i], cfg, opts);
                        r["id"] = i;
                        write_json(output / ("run-" + std::to_string(i) + ".json"), r);
                    } catch (const std::exception &e) {
                        errors[i] = e.what();
                    }
                }
            }
        }
    }
    manifest["failures"] = Json::object();
    for (std::size_t i = 0; i < errors.size(); ++i)
        if (!errors[i].empty())
            manifest["failures"][std::to_string(i)] = errors[i];
    manifest["status"] = manifest["failures"].empty() ? "complete" : "failed";
    write_json(output / "manifest.json", manifest);
    require(manifest["failures"].empty(), "failed runs: " + manifest["failures"].dump());
    std::cout << cases.size() << " runs in " << output << '\n';
}
} // namespace mix
