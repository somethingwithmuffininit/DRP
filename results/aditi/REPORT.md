# Aditi’s one-link experiment

The C++ benchmark reused MixBench’s four edge policies. It ran ten Erdős–Rényi graphs and ten Barabási–Albert graphs with 100 vertices each, plus twenty repetitions on one depth-seven binary tree. Every strategy used the same graph and sixteen sampled start vertices within a trial. Each policy adds a single missing edge. The unedited graph provides the TV baseline.

- **Erdős–Rényi.** The normalized connectivity gain averaged 0.510 for Fiedler and 0.537 for the exhaustive rule. Every baseline and edited walk crossed TV 0.25 at step 2.
- **Barabási–Albert.** Fiedler averaged 0.825 normalized connectivity gain and saved 0.5 TV steps. The exhaustive rule averaged 0.906 and saved 0.4 steps. These data compare descriptive averages from ten graphs.
- **Binary tree.** The Fiedler rule joined two leaves and reduced the sampled-start crossing from 517 to 501 steps in each trial. Random and minimum-degree selection slowed 15 and 12 of 20 trials respectively. In every slowdown the maximum degree rose from 3 to 4, reducing every edge transition probability from 1/3 to 1/4 under the selected weighting rule. The repeated tree trials measure selection and start variation on one topology.

The original standalone Python run and the C++ study agreed on strategy ordering by mean normalized gain in each family. Their absolute means differ because the graph and edge draws differ.

The independent verifier compared 180 C++ rows with the original Python scoring functions on exactly the same graphs and edges. The largest absolute errors were approximately 1.12e−12 for the edited second Laplacian eigenvalue and 2.0e−12 for the normalized gain. A separate SciPy exhaustive oracle checked every optimum. The original standalone Python run uses a shared Python RNG, so its graph stream differs from the C++ study’s independently seeded fixtures. The [graph fixtures](../../experiments/fixtures/link-policy-graphs.json) are the source of truth for the paired comparisons.

[Trial CSV](trials.csv) · [Aggregate CSV](summary.csv) · [Per-row parity](parity.csv) · [Original standalone run](original.txt) · [Reproduction script](../../scripts/reproduce_aditi.sh)
