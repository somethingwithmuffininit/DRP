# MixBench

MixBench compares how changes to a graph or its transition probabilities affect discrete-time Markov-chain mixing. The primary measurement is the first step at which the maximum total variation (TV) distance across starting states reaches a chosen tolerance. The C++ executable builds the chains and advances complete probability distributions. An optional experiment samples actual MCMC trajectories.

## Report and notebook workflow

Open [the MixBench notebook](notebooks/mixbench.ipynb) for the methods, tests, plots and tables in one place. The [experiment report PDF](output/pdf/mixbench-experiment-report.pdf), [experiment write-up](EXPERIMENTS.md) and [presentation](presentation/mixing.pdf) are available for sharing.

The [editable report source](report/report.py) rebuilds the PDF from the saved figures. From the MixBench directory, with Lato fonts installed at `/usr/share/fonts/truetype/lato`:

```sh
python -m pip install reportlab pdfrw pypdf
python report/report.py
```

Use `--no-student-credits` to omit the report's credit annotations or `--output FILE` to save a separate PDF.

Use a Jupyter kernel with the Python requirements installed and Poppler available for displaying the plots. Select **Run All**. Review mode configures a release build, runs the 16 CTest suites, calculates the exact shuffle curve and replays five quick graph presets. It uses the report's measured plots and tables for the larger sweeps. The same notebook can run in a terminal:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/run_notebook.py notebooks/mixbench.ipynb
```

Review mode saves the build and test logs, JUnit XML, CSV tables, PDF plots and PNG previews under `results/notebook/review/`. The five replays and the exact riffle calculation also save individual run records and the shuffle CSV. Set `MIXBENCH_RESULTS` to a fresh directory to repeat all measurements, or reuse the default directory to retain complete matching runs. Use `MIXBENCH_BUILD_DIR` to select a separate build directory.

`MIXBENCH_NOTEBOOK_MODE=full` replays all graph presets, checks the original implementations and solver output, and saves fresh plots, a combined case table and a run inventory under `results/notebook/full/`. It uses the original student reference files in the sibling `Caroline`, `Dana` and `Danielle` directories. Set `MIXBENCH_RUN_TIMED=1` in full mode to include the 105 timed trials. The serial timing schedule takes about 105 minutes. On Linux, `MIXBENCH_BUDGET_CPUS` can assign five distinct physical cores, one per repetition. `MIXBENCH_JOBS` controls the worker count for the graph sweeps.

The [saved inventory](results/notebook/full/inventory.csv), when present locally, lists every graph preset. The report and presentation also include the recorded timed summaries. Run the timed sweep if you need fresh timings or the individual timed trial records. Generated results are ignored by Git.

To remove student credits, delete the notebook's Markdown cells tagged `student-credit`, uncomment `\def\HideMixBenchCredits{}` near the top of [mixing.tex](presentation/mixing.tex) and rebuild, or delete the **Student credit** text annotations from the report in a PDF editor.

The [LaTeX source bundle](presentation/latex-source.zip) rebuilds the slides with pdfLaTeX, Beamer and TikZ. Run `python tests/solver_retry.py` to repeat the standalone chorded-cycle SDP check, or pass `--run-json` with a saved solution to audit that run.

## Build and reproduce the students' results

From this directory, use a C++20 compiler, CMake 3.24+, Eigen 3.4+, OpenMP, BLAS/LAPACK, OpenSSL development headers, and Python 3. CMake fetches checksum-pinned sources for SCS, Spectra, and nlohmann/json on the first configure. For an offline build, supply extracted source directories with `-DFETCHCONTENT_SOURCE_DIR_SCS=...`, `-DFETCHCONTENT_SOURCE_DIR_SPECTRA=...`, and `-DFETCHCONTENT_SOURCE_DIR_JSON=...`.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target mixbench mix_tests -j 3
ctest --test-dir build --output-on-failure
build/mixbench run experiments/reproduce.json --jobs 4 --output results/my-first-run
python3 scripts/plot_results.py results/my-first-run
```

The plot command needs Matplotlib. Install `requirements.txt` in your preferred Python environment. Output directories must be fresh. Replotting an existing result regenerates its summaries and figures. A completed run has `"status": "complete"` in its manifest.

At TV tolerance 0.25, Danielle's five-vertex example takes 3 steps with maximum-degree weights and 2 with optimized weights. Their respective SLEMs are 2/3 and 3/7. On Caroline's 32-vertex cycle the clockwise baseline takes 195 steps and the equal-halves shortcut takes 61. All starting vertices enter these measurements.

On the included real graphs at the same tolerance, the fixed Glasgow shortcut reduces Caroline's all-start mixing time from 43 to 20 steps. On the unedited karate club, the maximum-degree, Metropolis, and optimized kernels take 48, 35, and 26 steps respectively. On the power grid the default 32-start sample gives lower bounds of 28,515 steps for maximum degree and 7,721 for Metropolis. The preset records the full power-grid SDP as skipped by its size limit.

## Choose an experiment

| Configuration under `experiments/` | Measurements |
|---|---|
| `link-policies.json` | Paired one-link choices on 40 archived graph realizations with fixed seeds |
| `matched-ring-links.json` | Opposite-side, random, fewest-neighbor and connectivity-based choices on one clockwise ring |
| `bridge-count.json` | Four bridge counts between two 8-cliques at four TV tolerances |
| `reproduce.json` | Danielle's complete bipartite example and Caroline's complete placement sweep |
| `dana-theorem1.json` | Dana’s lifted path at 100 positions, all 200 internal starts, 1,200 steps, and the minorization snapshot at step 400 |
| `method-scaling.json` | Matched size sweeps for path walks, cycle shortcuts, tree heuristics and small-graph connectivity optimization. Includes MCMC cost and error |
| `path-comparison.json` | Exhaustive starts for Dana and reversible controls on 10, 20, 50 and 100 positions at three TV tolerances |
| `cycle-placement.json` | Every allowed shortcut on the ordered 32-cycle, plus the no-shortcut baseline |
| `scaling.json` | Cycle and path sizes up to 256, odd and even cycles, four tolerances, and seeded random chords sampled with replacement |
| `real-networks.json` | Glasgow's 15-station ring with fixed shortcut and unedited controls, plus fixed-weight and optimized walks on Zachary's karate-club graph |
| `power-grid.json` | Sparse comparison on 4,941 vertices, 32 shared starting states, 40,000-step limit |
| `edge-comparison.json` | Exhaustive placements and reference policies on a cycle, two joined cliques, and the karate club |
| `paper-families.json` | Paired single-edge comparisons on connected Erdős–Rényi, preferential-attachment, and rooted tree instances |
| `mcmc.json` | Independent stationary trajectories, error against the known observable mean, effective sample size, throughput, and setup break-even |

The three real graphs ship as edge lists. Their sidecar JSON records the source, SHA-256, preprocessing, and vertex-ID mapping. The importer verifies the digest before analysis. Glasgow supplies an observed station order. Transition probabilities on all three datasets are experimental choices. [SPT station map](https://www.spt.co.uk/travel-with-spt/subway/maps-stations/) · [Network data and original citations](https://public.websites.umich.edu/~mejn/netdata/)

## Configure a comparison

Copy a preset and run it with `mixbench run CONFIG --jobs N --output DIR`. A study crosses its methods with its selected edge edits for each graph realization. Each method gets the same realized edge and the same sampled starts on the same state space.

```json
{
  "seed": 42,
  "epsilons": [0.25, 0.1],
  "max_steps": 20000,
  "studies": [{
    "name": "Random shortcut placement",
    "graph": {"family": "cycle", "sizes": [31, 32, 64]},
    "methods": ["clockwise", "max_degree", "fmmc"],
    "policies": [
      {"policy": "none"},
      {"policy": "random", "candidates": "anchored", "samples": 12}
    ]
  }]
}
```

- **Graph.** Choose `cycle`, `path`, `bipartite` or `barbell` (optional `split`), `er` (edge probability `p`, conditional on connectivity), `ba` (attachment count `m`, complete seed of size `m+1`), or `tree` (branching factor `branching`). Set study-level `replicates` for independent graph realizations. A dataset study uses `{"dataset":"../data/karate.json"}` relative to the config. Its sidecar points to a whitespace-separated list of zero-based endpoint pairs. Use `order:"cycle"` or `order:"path"` if IDs follow that order.
- **Fixed graph trials.** A graph can specify a `fixtures` JSON file relative to its preset and a `family` to select records from it. Each record supplies `family`, `n`, `edges`, `policy_seed` and `start_seed`. Records in a family share a size. Every policy on a record gets the same sampled starts, and its edits come from that record’s policy RNG.
- **Kernel.** `max_degree` assigns each edge probability 1/maximum degree. `metropolis` assigns 1/max(endpoint degrees). `fmmc` optimizes symmetric weights and permits holding at each vertex. `clockwise` uses Caroline's directed transitions on an ordered cycle with an optional reciprocal chord. Dana’s `dhn` kernel uses two internal states per position on an unedited ordered path. The presets target uniform observed positions. [Boyd, Diaconis & Xiao](https://web.stanford.edu/~boyd/papers/pdf/fmmc.pdf) · [Gerencsér, Definition 1](https://arxiv.org/pdf/1905.03223)
- **Edge policy.** Choose `none`, `fixed` with `edge:[i,j]`, `halves` for even cycles, `random`, or `all`. Caroline's preset uses `candidates:"anchored"`. Its edges are `[k-1,n-1]` with paper parameter 2 ≤ k ≤ n−2. The default candidate set contains every missing unordered pair. `min_degree`, `fiedler`, and `optimal_connectivity` implement the fewest-neighbors, Fiedler-coordinate and exhaustive connectivity policies. Random ties for the minimum degree use the run RNG. Spectral ties within numerical tolerance use the first eligible pair in vertex order. Exhaustive optimization computes an eigensolve for each candidate and is intended for small graphs. [Wang & Van Mieghem](https://www.nas.ewi.tudelft.nl/people/Piet/papers/Bionetics2008_AlgebraicConnectivity.pdf)
- **Workload.** Defaults are `epsilons:[0.25]`, `max_steps:100000`, `stride:20`, `snapshots:[0,10,100,500]`, and zero-based `diagnostic_start:0`. Distances are checked at every step including zero. `all_starts_limit:256` switches larger chains to `sample_starts:32` distinct starts including the diagnostic vertex. A sampled maximum gives a lower bound on the all-start mixing time. `sdp_limit:64` caps the full semidefinite solve. Override any of these in the top-level JSON.
- **Parallelism.** One OpenMP team schedules runs, blocks of starting states, and MCMC replicates. Each run writes a separate numbered file. Seeds, edits, and starting states are independent of job scheduling. The embedded SCS direct library disables its optional process-global interrupt handler. Independent solves have separate workspaces. An inaccurate iterate that fails the residual or gap check gets one warm-start retry with a larger iteration budget. For timing studies set `OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1`, hold the job count fixed, and run one configuration at a time. Use `--jobs 1` for isolated timing of competing methods.

## Read the results

The manifest contains the configuration, exact base edge lists, and source metadata. Each `run-ID.json` contains the edit, graph and internal state counts, solver diagnostics, spectral measurements, seeds, and analysis. For optimized chains the file includes the transition matrix as triples `[from,to,probability]`. Set `export_kernel:true` to save triples for other kernels.

- **Mixing time.** `analysis.crossings` reports a TV threshold, measurement space, status, and step. `unresolved` means the measured maximum still exceeds the threshold at the horizon. `periodic` identifies tolerances below a proven periodicity floor. Missing times serialize as null. The per-start first-crossing array follows the same space/epsilon order and uses −1 for an unreached threshold.
- **Spectra.** `slem` measures the transition kernel. `spectral_time` is the full chain’s asymptotic e-folding time, measured in transitions. `lambda2_before` and `lambda2_after` measure the unweighted combinatorial Laplacian of the undirected support. Reweighting a fixed graph preserves this connectivity. A repeated Fiedler eigenvalue makes Wang's normalized gain undefined, recorded as null. On an unweighted cycle every single chord preserves algebraic connectivity within eigensolver roundoff, although the TV crossing can change substantially.
- **Plots and tables.** `summary.csv` has one row per tolerance and measurement space. `vs_unedited` compares the same kernel before the edit. `vs_max_degree` compares with the maximum-degree kernel on the same edited graph. Ratios use matched resolved starts. Figures show TV by step, the selected start's distribution, spectral connectivity, placement and size sweeps, and the relationship between connectivity gains and mixing steps saved.
- **Sampling efficiency.** The MCMC preset measures the indicator that the observed position lies in the first half of the ID range. Starts come from stationarity to isolate sampling variance. Across independent replicates, effective draws are estimated from the exact observable variance divided by the variance of replicate means, multiplied by the replicate count. `mcmc.csv` reports the RMSE, effective draws per second, and the effective-draw count needed to recover added setup cost when throughput improves. These estimates depend on the observable, replicate count, and machine load.

## Compare performance across graph sizes

From the MixBench directory, these commands build the updated executable, run the size sweep and generate the comparisons. Use a fresh run directory. One worker gives isolated elapsed times for timing comparisons. A larger `--jobs` value runs simulations in parallel.

```sh
cmake --build build --target mixbench -j 3
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 build/mixbench run experiments/method-scaling.json --jobs 1 --output results/my-size-sweep
python -B scripts/plot_results.py results/my-size-sweep
```

The preset measures all internal starting states at its configured sizes. Paths and trees have 8–128 vertices. The cycle study includes even and odd sizes. Exhaustive optimization on two joined cliques uses 8–32 vertices. The FMMC solver limit is 64 vertices. Edit the `sizes` arrays to change the sweep. For Dana, the horizontal axis counts physical vertices and the internal chain has twice as many states.

| Output | How to read it |
|---|---|
| `summary.csv` | One row per run, TV space and tolerance, plus a status row for each skipped run. Includes matched speedups for the maximum, median and 90th-percentile crossing |
| `size-statistics.csv` | For each size and metric, median and interquartile range across graph realizations and edge draws. The case and graph counts show what contributes to each point |
| `size-fits.csv` | Descriptive exponent from fitting log(metric) against log(vertices) with equal weight per size. At least three positive, available size medians are required. Includes the fitted sizes, omissions, solver skips and log-scale R². An exponent near 2 describes approximately quadratic growth over those sizes |
| `size-figures/index.csv` | Index of figures for each study and tolerance covering mixing steps, spectra, speedups, cost and MCMC. Points show medians across runs at each size and shading shows the middle half of those runs |
| `figures/g*-methods.png` | Within-graph comparison of the median, 90th percentile, maximum-TV crossing and asymptotic rate for the displayed methods |

For comparisons across completed experiments, supply their directories explicitly. The command reads saved runs and writes `cases.csv`, `policies.csv`, the size tables and, with `--figures`, the size plots. Historical runs supply the fields they recorded. New spectral measurements require output from the updated executable.

```sh
python -B scripts/compare_results.py results/my-comparisons results/my-size-sweep results/my-first-run --figures
```

**Convergence and connectivity.** `transition_lambda2` and `transition_lambda_min` are the second-largest and smallest algebraic eigenvalues of a reversible transition matrix. Its ordinary gap is `1 − transition_lambda2`. For all kernels, `slem` is the largest modulus after excluding the stationary eigenvalue. The real and imaginary fields give one mode attaining that modulus. With tied modes the representative can vary. The absolute gap is `1 − slem`, the asymptotic log decay per step is `−log(slem)`, and its reciprocal is the e-folding time. A periodic chain has zero asymptotic decay. The status field identifies periodicity, numerical limits and an unresolved sparse eigensolve. Dana’s calculation uses paired Fourier modes on the lifted cycle. All spectral rates refer to internal states, including when a row displays position TV. Graph `lambda2` is the algebraic connectivity of the unweighted undirected support.

**Measured TV.** The observed log decay is `log(initial maximum TV / end maximum TV) / end step`, measured up to the first common crossing or the checked horizon. The row includes both distances and the end step. A zero-length window or end distance at the numerical floor leaves the rate blank. Start quartiles, mean and population standard deviation describe individual first crossings. Unreached starts rank above all observed crossings. A quantile is blank when its rank requires an unobserved crossing. The mean and standard deviation require every measured start to cross. A method can lead on the median and trail on the maximum.

**Interpreting scaling.** Each curve holds the experiment, graph generator and options, method, policy, tolerance, TV space and start scope fixed. Cycle parity and sampled start counts have separate series. Undefined or censored values produce gaps and the size tables give observation counts. A per-size aggregate requires every run at that size to supply the metric. Each random edge draw has equal weight. Estimated exponents describe this finite range. Single-size real networks have case comparisons and cannot establish a size trend on their own.

**Timing.** Setup includes kernel construction and the assigned edge-selection cost. Spectrum time measures the transition eigensolve. Propagation time covers all measured starts to the actual stopping horizon, including snapshots and any early stopping. `updates_per_second` counts logical start-transitions to the reported horizon. Use the same worker count, BLAS settings and machine for runtime comparisons. MCMC throughput and RMSE refer to the observable and trajectory budget in the manifest.

## Four comparison experiments

Each preset uses the normal executable and report command. For example, from the MixBench directory:

```sh
build/mixbench run experiments/edge-comparison.json --jobs 4 --output results/my-edge-study
python -B scripts/plot_results.py results/my-edge-study
```

| Preset | Comparison | Outputs to open first |
|---|---|---|
| `edge-comparison.json` | Every missing edge on a cycle, two joined cliques and the karate-club graph, at TV 0.25 and 0.1 | `edge-agreement.csv`, `edge-ranks.csv`, `edge-starts.csv` |
| `path-oracle.json` | Dana and maximum-degree walks from the same physical endpoint, over seven sizes from 8 to 512 | `path-scaling.png`, `path-controls.csv`, `path-sensitivity.csv` |
| `wall-budget.json` | Estimator error after 10 and 60 wall-clock seconds, including method setup | `budget-comparison.csv`, `budget-trials.csv` |
| `nested-networks.json` | Size sweeps within the Glasgow subway and Western power grid | `nested-realizations.csv`, `nested-variation.csv`, `size-figures/index.csv` |

**Edge rankings.** Rank 1 includes all tied optima. Connectivity gains within an absolute tolerance of 1e−10 of the leading value of a tie group share its rank. TV rankings use exact crossing steps. The agreement table shows the overlap of optimal edge sets and the range of extra steps among connectivity-optimal edges. The ranks table lists every edge and each heuristic draw relative to the fastest edge for both the median and maximum crossing. The per-start table gives each heuristic’s steps saved over the unedited graph and extra steps over the fastest edge from that vertex. A censored candidate leaves the affected oracle ranking and regret blank. To analyze a previously completed exhaustive sweep at its recorded tolerances, use `scripts/study_reports.py DIRECTORY --report edge-ranking` with the project Python interpreter.

**Path controls.** For (n) positions, the report checks the support Laplacian’s second eigenvalue against `4 sin²(π/(2n))` and the reversible transition matrix’s second eigenvalue against `cos(π/n)`, with absolute tolerance 1e−10. Both start at physical position zero. Study-level `diagnostic_starts` assigns DHN internal start 1 (heading right) and reversible start 0. The report uses position TV, plots steps divided by (n) or (n²), and refits the exponent while dropping successive smallest sizes. Each fit needs at least three resolved, positive crossing times. An oracle failure leaves the numerical check table on disk.

**Time budgets.** Use one worker and single-threaded BLAS for comparable timings:

```sh
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 build/mixbench run experiments/wall-budget.json --jobs 1 --output results/my-budget
python -B scripts/plot_results.py results/my-budget
```

The preset schedules 21 method/policy combinations, each with five trials. Allow approximately 105 minutes for the default 60-second endpoint. The 10-second reading comes from the same trajectory. Change `wall_budget.seconds` and `repeats` for shorter trials. Trials execute serially in seeded, shuffled order within each repetition. Methods share the graph, edit seed, sampling seed and physical initial position. The clock starts before method-specific edge selection, editing, optimization and sampler construction. Shared input graph generation and file I/O are outside the budget. The sampler counts completed transitions through each deadline. A setup call that overruns the budget is recorded when it returns. Observables are membership in the left clique (first half of positions on a cycle or path) and visits to the smallest-ID minimum-degree vertex of the base graph. The comparison table shows bias, RMSE, number of estimates, median draws and setup exhaustion. Bias and RMSE require an estimate from every planned repeat.

**Nested networks.** Four distinct starting stations define contiguous Glasgow sections of 4–14 stations. Six distinct seed vertices define power-grid samples of 12, 24, 48 and 96 vertices. Each grid sample grows by choosing uniformly from its current frontier. Each size uses a prefix of the same replicate’s growth order and retains every parent edge between selected vertices. The manifest records parent vertex IDs in local-ID order and the parent checksum. Methods on each graph share candidate rules and realized edits. The preset measures every start. The reports summarize repeated edge draws within each realization first, then show the spread across realizations of the same size. These are overlapping samples from one parent network. With the bundled SDP limit, FMMC at 96 vertices is recorded as skipped.

## Add a student method

The public C++ API is in `include/mixing.hpp`. A topology method returns an `Edit` from `select_edges`. A transition method returns a `Kernel` from `build_chain`. Both use the same immutable graph. Transition matrices are row-stochastic, and the simulator updates column distributions with the transpose. The validator checks mass, stationarity, and detailed balance for kernels declared reversible.

## Paired single-link trials

The [link policy preset](experiments/link-policies.json) uses the same runner and reports as the other graph experiments. It replays 40 archived graph fixtures and seeds from the original study. The dense random and preferential-attachment families each have ten 100-vertex graphs. The 255-vertex binary tree has twenty sets of policy draws and sampled starts. Each edit shares 16 starting vertices with the unedited maximum-degree walk. The preset uses exhaustive connectivity search on the two smaller families.

```sh
build/mixbench run experiments/link-policies.json --jobs 4 --output results/my-link-policies
python -B scripts/study_reports.py results/my-link-policies
python -B tests/link_policy_trials.py --results results/my-link-policies
```

`link-policy-summary.csv` reports mean normalized connectivity gain, mean TV steps saved and the number of slower trials. Positive steps saved means the edited walk mixed faster. The unedited and edited runs are paired by graph, method, threshold and start set. The [original submission comparison](results/aditi/REPORT.md) preserves the fixtures, measurements and independent Python verification.

## Reproduce Dana’s Theorem 1 checks

Set `methods: ["dhn"]` on a path with policy `none`. Study-level `dhn_c` defaults to 1, the value used in Theorem 1. The kernel accepts other rates in the range 0 < c < n. It stores Dana’s state labels modulo 2n. For the converse bound use internal start 0 and full-state TV. For the position comparison use internal start 1 (position 1, heading right), project out direction, and compare with `max_degree` on the path started at index 0. [Dana’s reference script](../Dana/theorem1.py)

Run the following from this directory. The first command saves ordinary MixBench results. The second checks all six of Dana’s sizes (10 through 500) against her actual Python functions using NumPy and SciPy. It prints CSV after each successful comparison and exits with an error on a mismatch.

```sh
build/mixbench run experiments/dana-theorem1.json --jobs 4 --output results/my-dana-run
python -B scripts/verify_dana.py --jobs 4
```

The verifier checks the transition matrices and projection exactly, all full-state per-start threshold crossings, the minimum and maximum TV curves over starts, the state-0 converse, both upper bounds through 12n, and both single-start position curves through TV 0.01. At step 4n, the measured minorization constant is 2n times the snapshot’s `min_probability_over_starts` over **all** 2n starts. For n = 500, the full-state worst-case TV crosses 0.1 at step 1,253. Starting at position 1, the lifted and reversible position curves cross 0.1 at steps 1,198 and 93,773 respectively.

## Run the regression suite

Open [the MixBench notebook](notebooks/mixbench.ipynb) in Jupyter and run all cells. It builds the project and shows what each test checks, its result and elapsed time, then walks through the experimental comparisons. The notebook uses the Python interpreter from the current kernel so the optional tests for NumPy, SciPy and Matplotlib are available. For the same suite from a terminal:

Configure using the Python environment containing the project requirements to enable every suite. From the MixBench directory, run

```sh
cmake -S . -B build -DPython3_EXECUTABLE="$(command -v python)"
cmake --build build --target mixbench mix_tests -j 3
OPENBLAS_NUM_THREADS=1 ctest --test-dir build --output-on-failure
```

The independently selectable CTest suites cover these comparisons.

| Suite | What it checks |
|---|---|
| `core` | Known weights and spectra, Caroline’s dense cycle transitions, edge policies, and a projected walk whose two starts first cross at steps 1 and 3 but first cross together at step 6 |
| `cli` | Serial and parallel agreement across methods, periodic and size-limited results, and dataset checksums |
| `dana_parity` | Dana’s actual Python implementation at 10 and 20 positions. Use the full verifier above for all six sizes |
| `dhn` | A dense walk built from physical positions and headings, including boundary reflections. Checks both TV spaces, every start’s crossings, snapshots, invalid inputs, and seeded MCMC mean and variance against finite-run predictions |
| `studies_edge_ranking`, `studies_path_oracle`, `studies_nested_networks`, `studies_wall_budget` | Independently runnable fixtures for tied optima, matched path starts, nested induced edges and paired deadlines |
| `riffle` | Eulerian class counts, repeated explicit small-deck shuffles and the 52-card reference probabilities |
| `bridges` | Nested bridge supports, an analytic eigenvalue formula, stochastic weights and invalid bridge counts |
| `link_policy_trials` | Replays one archived graph per family and compares chosen edges, matched starts, connectivity scores and TV crossings |
| `budget_repetitions` | Agreement of the serial and partitioned schedules for seeds, graph realizations, realized edges and trial order |
| `power_grid` | Diagnostic TV and snapshots from independent SciPy sparse propagation on the 4,941-vertex grid |
| `real_networks` | All 17 Glasgow and karate-club runs. Recomputes connectivity from the input edges and TV from dense transitions. For FMMC, checks feasibility and propagates the exported optimized weights |
| `edge_policies` | All missing edges on the 12-cycle, two joined cliques, and the karate-club graph. Checks each connectivity gain, the exhaustive optimum, and the heuristic selection rules |
| `plots` | Rendered figures and summary ratios for exhaustive, sampled, matched, and mismatched starts, including a threshold reached at step zero |
