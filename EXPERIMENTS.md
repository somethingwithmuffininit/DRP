# MixBench experiments and results

I tested how quickly different random walks forget where they started, and whether faster mixing helps when we use a walk to estimate a probability. The clearest speedup came from Dana’s direction memory on long paths. Danielle’s optimized weights improved the slowest starts on the karate club graph. Caroline’s and Aditi’s link experiments showed that the location of a new link can matter as much as adding one. In the timed ring experiment, a shortcut shortened mixing but increased the error of an estimate.

## What I measured

Total variation (TV) is the largest difference between the walk’s current probability and its target probability for any event. For example, at TV 0.25, the probability of being in the left half of a network differs from its target by at most 25 percentage points. Unless specified otherwise, the step counts below use **TV 0.25** and report the first step when every tested start is within that distance. The target is uniform over the vertices or deck orders. The full power grid and the repeated link-policy study use sampled starts, so their measured maximum can be smaller than the all-start maximum. For Dana’s walk, position TV combines the two headings at each position, while full-state TV tracks both position and heading.

| Method | Change to the walk or graph |
| --- | --- |
| Maximum-degree baseline | Give each link probability one divided by the largest number of neighbors anywhere in the graph. Use the remaining probability to stay put. |
| Metropolis baseline | Use the same rule locally, dividing by the larger neighbor count at the two ends of each link. |
| Danielle’s optimized weights (FMMC) | Keep the links and choose symmetric move probabilities to make the slowest probability imbalance fade as fast as possible. |
| Dana’s lift | Give each position a left- and right-heading state. The walker usually continues in its current direction and can occasionally reverse early. |
| Caroline’s shortcut | Start with a walker that moves clockwise half the time and stays half the time. At either end of a new two-way link, use one quarter of the total probability to cross and one quarter to stay. |
| Aditi’s link rules | Compare a random link, a link from a vertex with fewest neighbors, the Fiedler rule, and a search over every eligible link. |

Algebraic connectivity measures how strongly the graph is joined. It is the second-smallest eigenvalue of the graph Laplacian, which compares vertex values with their neighbors. The Fiedler vector gives one coordinate per vertex for that eigenvalue. Aditi’s Fiedler rule links the eligible pair whose coordinates differ most. The exhaustive rule computes the actual connectivity score after each possible edit. Danielle instead uses the transition matrix. The magnitude of each nonstationary eigenvalue tells us how much of a corresponding probability imbalance survives one step.

## Card shuffling and the original examples

- **Exact riffle shuffle calculation.** I computed TV under the standard random riffle model for decks of 8, 16, 32 and 52 cards. For 52 cards, TV fell from 0.6135 at six riffles to 0.3341 at seven and 0.1672 at eight. Eight was the first crossing of 0.25. The calculation counts deck orders exactly. [MixBench notebook](notebooks/mixbench.ipynb) and [riffle curves](presentation/figures/riffle.pdf).
- **Danielle’s five-vertex graph.** Each of two vertices connects to all three vertices in the other group. Optimizing the probability on each link from 1/3 to 2/7 lowered the slowest remaining eigenvalue magnitude from 2/3 to 3/7 and shortened the crossing from 3 to 2 steps. Maximum degree and Metropolis both took 3 steps. The largest difference between the optimized C++ and Python matrix entries was about 8.7 × 10⁻¹¹. [Reproduction preset](experiments/reproduce.json) and [MixBench notebook](notebooks/mixbench.ipynb).
- **Caroline’s 32-vertex ring.** I fixed one end of the shortcut, tried all 29 eligible placements for the other end, and checked all 32 starts. The original took 195 steps and the opposite-side link took 61. Four other placements tied it. Caroline’s implementation agreed on all 960 start-by-placement crossings, including the original ring. Probability snapshots at steps 10 and 100 show how a shortcut sends part of the distribution to a different part of the ring. [Placement preset](experiments/cycle-placement.json), [placement sweep](presentation/figures/placement.pdf) and [probability snapshots](presentation/figures/cycle-mass.pdf).

## Direction memory on a line

**Same path and every starting position.** I compared three ordinary walks with Dana’s lift for paths of 10, 20, 50 and 100 positions. The lift uses both initial headings. The preset also tested TV 0.1 and 0.01. [Path comparison preset](experiments/path-comparison.json).

| Positions | Maximum degree | Metropolis | Optimized weights | Lift, position TV | Lift, full-state TV |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 19 | 19 | 19 | 13 | 14 |
| 20 | 76 | 76 | 76 | 27 | 30 |
| 50 | 474 | 474 | 474 | 70 | 77 |
| 100 | 1,898 | 1,898 | Above solver limit | 141 | 155 |

**One matched endpoint and increasing path length.** I started both walks at the same end, with the lifted walker heading into the line, and tested seven lengths from 8 to 512. At TV 0.1 the 512-position lift took 1,227 steps to mix positions and maximum degree took 98,328. The fitted size exponents were 1.039 and 2.000, respectively. Over the measured range, doubling the length roughly doubled the lifted time and quadrupled the ordinary time. The analytic eigenvalue checks passed at all seven sizes. [Endpoint preset](experiments/path-oracle.json) and [scaling figure](presentation/figures/path.pdf).

**Dana’s theorem and reference code.** I checked the finite-size upper and lower bounds and reproduced Dana’s Python curves for 10, 20, 50, 100, 200 and 500 positions. The largest curve difference was about 4.2 × 10⁻¹⁴. At 100 positions the lift took 248 steps for full-state TV 0.1 and 238 for position TV 0.1 over every start. In separate endpoint checks at 200 and 500 positions, the lifted versus ordinary position crossings at TV 0.1 were 478 versus 15,004 and 1,198 versus 93,773. At 500 positions the worst-start full-state crossing was 1,253. [Theorem preset](experiments/dana-theorem1.json) and [MixBench notebook](notebooks/mixbench.ipynb).

## Shortcut placement and bottlenecks

**Caroline’s size sweep.** I compared the unedited ring, the opposite-side link where available, and six seeded placements per size on rings from 15 to 256 vertices at four TV thresholds. At 256 vertices the original took 12,435 steps, the opposite-side link took 3,109, and the six random links ranged from 1,086 to 9,393. Four of those six beat the opposite-side link’s time. The companion sweeps tested reweighting on paths and odd and even rings. On the 32-vertex even ring, maximum degree and Metropolis alternated between the two parity classes and never met 0.25. Optimized weights reached it in 56 steps. On the 64-position path, all three weight rules took 777 steps. [Scaling preset](experiments/scaling.json) and [ring scaling figure](presentation/figures/cycle-scaling.pdf).

**Caroline and Aditi on the same ring.** I ran 73 cases with all 32 starts and Caroline’s clockwise transition rule. The 29 allowed links tied on algebraic connectivity, but mixing times spanned 61 to 172 steps. The exhaustive connectivity rule picked the first tied link in the implementation. The ring also has a tie in vertex degree, so the two random-policy sample medians describe these draws. [Matched ring preset](experiments/matched-ring-links.json) and [presentation, slide 15](presentation/mixing.pdf#page=15).

| One-link rule | Steps to TV 0.25 | Range across draws |
| --- | ---: | --- |
| Original ring | 195 | |
| Caroline’s opposite-side link | 61 | |
| Aditi’s Fiedler rule | 63 | |
| Random, median of 20 | 71.5 | 61–151 |
| Fewest neighbors, median of 20 | 64 | 61–172 |
| Exhaustive connectivity rule | 172 | |

**Every missing link on three networks.** I ranked the edited graphs by connectivity and by mixing from the slowest start, keeping tied winners. The two objectives agreed on all winning links for the graph of two joined cliques. On the ring, connectivity tied every link, and a connectivity-optimal choice could cost 11 extra steps. On the karate club graph, the best-connectivity link tied for fastest at 32 steps, but its median start took three steps more than the median-start optimum. [Exhaustive link preset](experiments/edge-comparison.json) and [connectivity versus mixing plot](presentation/figures/edges.pdf).

| Network | Missing links tested | Best connectivity | Fastest mixing | Shared winners | Original → fastest |
| --- | ---: | ---: | ---: | ---: | --- |
| 12-vertex clockwise ring | 54 | 54 | 6 | 6 | 28 → 13 |
| Two groups of six, initially joined by one bridge | 35 | 25 | 25 | 25 | 17 → 9 |
| 34-member karate club | 483 | 1 | 16 | 1 | 48 → 32 |

**Adding more bridges.** I connected two groups of eight, each internally fully connected, with 1, 2, 4 and 8 links between distinct pairs. As the bottleneck opened, both the crossing time and the advantage of optimizing the weights shrank. [Bridge-count preset](experiments/bridge-count.json) and [bridge-count plot](presentation/figures/edges.pdf).

| Bridges | Maximum degree | Metropolis | Optimized weights |
| ---: | ---: | ---: | ---: |
| 1 | 28 | 28 | 15 |
| 2 | 14 | 14 | 7 |
| 4 | 7 | 7 | 4 |
| 8 | 3 | 3 | 2 |

## Which link rules work across graph families

**Small synthetic graphs.** I compared the original graph, random, fewest-neighbor and Fiedler links at 50 and 100 vertices for dense random graphs and graphs that grow around hubs, with four realizations per size. I also tested fixed binary trees at 31 and 63 vertices. At size 100, every policy had a median of 3 steps on the dense graphs. On the hub-growing graphs, the median fell from 41 to 38 with Fiedler. On the 63-vertex tree the original took 120, the fewest-neighbors draw took 148, and Fiedler took 111. [Graph-family preset](experiments/paper-families.json) and [family curves](presentation/figures/families.pdf).

**Repeated single-link trials.** I tested ten dense 100-vertex graphs, ten 100-vertex preferential-attachment graphs and twenty repetitions on one 255-vertex binary tree. The original and edited walks shared 16 sampled starts in each trial. All dense walks reached the threshold in two steps, even though the edits increased connectivity. In the preferential-attachment graphs, the exhaustive rule won on mean scaled connectivity gain (0.906 versus Fiedler’s 0.825), while the mean TV savings were 0.4 and 0.5 steps respectively. The tree results were larger.

| Rule on the 255-vertex tree | Mean steps saved | Trials that slowed |
| --- | ---: | ---: |
| Random | −117.8 | 15 of 20 |
| Fewest neighbors | −95.5 | 12 of 20 |
| Fiedler | +16.0 | 0 of 20 |

Fiedler joined two leaves and brought the crossing from 517 to 501 steps. In all 27 slowdowns, the new edge raised maximum degree from 3 to 4, which reduced every link’s move probability from 1/3 to 1/4. The exhaustive search was omitted on the large tree. The original study checked all 180 rows against Aditi’s scoring code. I replayed those cases through the shared runner and reproduced the links, connectivity scores and TV crossings. [Shared link-policy preset](experiments/link-policies.json) and [original fixtures and reference results](results/aditi/REPORT.md).

**Method and size sweep.** I ran 405 cases across paths (8–128 positions), odd and even rings (8–128 vertices), trees (8–128), and pairs of fully connected groups (8–32 total). The tree and two-group sweeps used eight draws per size for each randomized policy. The run measured TV at 0.25 and 0.1, crossing times by start, eigenvalue decay rates, setup cost and short stationary trajectories. On fixed paths, the maximum-degree and optimized walks track closely in decay and crossing steps at the measured sizes, while optimized weights take more setup time. In the tree and two-group plots, the effect of selecting a new link varies with graph size and policy. The campaign recorded 393 measured cases, 12 solver size skips and 20 periodic threshold records. The setup curves reflect the parallel run environment. [Method-scaling preset](experiments/method-scaling.json), [link-policy curves](presentation/figures/method-scaling.pdf) and [setup and decay curves](presentation/figures/work-scaling.pdf).

## Real networks

**Karate club friendships.** The network contains 34 members and 78 friendship ties from Zachary’s study. I kept the same graph and tested all member starts. Original worst-start times were 48 (maximum degree), 35 (Metropolis) and 26 (optimized). The corresponding median-start times were 21, 15 and 21.5. A single Fiedler link shortened the worst-start times to 32, 25 and 17. A random link gave 46, 34 and 29. That last edit made the optimized walk slower than its unedited 26-step run. [Real-network preset](experiments/real-networks.json), [data description](data/karate.json) and [start-by-start curves](presentation/figures/starts.pdf).

**Glasgow Subway.** I used the published Outer Circle station order, treating each of the 15 stations as a vertex and consecutive stations as links. Then I tried one hypothetical Kinning Park–Cowcaddens shortcut. Results are worst-start steps for the modeled network. [Station data description](data/glasgow.json) and [network schematic](presentation/mixing.pdf#page=18).

| Walk | Original | With shortcut |
| --- | ---: | ---: |
| Clockwise | 43 | 20 |
| Maximum degree | 43 | 16 |
| Metropolis | 43 | 21 |
| Optimized weights | 14 | 16 |

**Full Western States power grid.** I kept the original 4,941 vertices and 6,594 physical links and used the same 32 sampled starts. Maximum degree took 28,515 steps over that sample, while Metropolis took 7,721. Median starts took 18,797.5 and 5,154.5 respectively. The symmetric optimizer exceeded the 64-vertex solver cap. [Power-grid preset](experiments/power-grid.json) and [data description](data/power.json).

**Growing sections from real networks.** I built consecutive Glasgow sections of 4–14 stations from four starting stations. Each section is a line, and sections with the same station count have the same shape. At four stations, maximum degree took 2 steps and Dana’s lift took 4. A separate all-start check on a nine-station section gave 16 and 12 steps, respectively. I also grew six connected samples of the power grid to 12, 24, 48 and 96 vertices, retaining all original links between the sampled vertices. At 96 vertices the maximum-degree crossing ranged from 289 to 584 steps over the six samples. A Fiedler link moved the median from 395 to 270. Optimized weights ran through 48 vertices in this sweep. [Nested-network preset](experiments/nested-networks.json), [growth diagram](presentation/figures/nested-map.pdf) and [nested results](presentation/figures/nested.pdf).

## Sampling and estimation

**Same number of simulated steps.** I started 32 replicate trajectories at stationarity and collected 20,000 draws each. The quantity estimated was the fraction of visits in the first half of the vertex order, with true value 0.5. On the karate club network the root mean squared errors were about 0.0173 for maximum degree, 0.0201 for Metropolis and 0.0143 for optimized weights. On the ring, the shortcut reduced TV time from 195 to 61 steps but increased the estimate’s error from 0.0037 to 0.0102. Optimized weights also had the highest estimated sampling throughput on the karate graph in the saved run, about 3.49 million effective draws per second versus 2.21 million for maximum degree. An effective draw expresses precision in terms of independent samples. The estimated break-even point for the extra setup work was about 14 million effective draws in that run. [Stationary sampling preset](experiments/mcmc.json) and [error and throughput plots](presentation/figures/mcmc.pdf).

**Same wall-clock budget.** I ran 21 method-and-edit combinations with five paired trials each. The 10-second reading came from the same trajectory as the 60-second reading, and the clock included edge selection and setup. Each walk estimated the fraction of visits in one half of the graph (true value 1/2) or at a specified low-degree vertex (true value 1/32). None of the runs spent a full checkpoint on setup. The table gives root mean squared error in probability units for both estimates. [Timed preset](experiments/wall-budget.json), [MixBench notebook](notebooks/mixbench.ipynb) and [two-group benchmark plot](presentation/figures/budget.pdf).

| Walk on 32 positions | Half, 10 s | One vertex, 10 s | Half, 60 s | One vertex, 60 s |
| --- | ---: | ---: | ---: | ---: |
| Ordinary path | 4.07 × 10⁻⁴ | 4.95 × 10⁻⁵ | 1.51 × 10⁻⁴ | 2.84 × 10⁻⁵ |
| Path with direction memory | 6.88 × 10⁻⁵ | 7.06 × 10⁻⁶ | 2.91 × 10⁻⁵ | 4.35 × 10⁻⁶ |
| Clockwise ring | 4.47 × 10⁻⁵ | 1.17 × 10⁻⁵ | 1.75 × 10⁻⁵ | 3.24 × 10⁻⁶ |
| Clockwise ring with shortcut | 1.59 × 10⁻⁴ | 1.67 × 10⁻⁵ | 4.41 × 10⁻⁵ | 3.79 × 10⁻⁶ |

Direction memory also reduced the one-vertex error on the path at both deadlines. For the original two-group graph at 60 seconds, Metropolis had the smallest measured error for left-group visits (1.16 × 10⁻⁴) and optimized weights had the smallest for visits to the specified vertex (9.89 × 10⁻⁶). The ranking depends on the quantity being estimated. Small differences have limited precision with five trials on one machine.

## Run inventory and checks

The original notebook audit covered **15 presets and 2,248 planned records**, including 105 timed trials. Of those, 2,170 produced measurements and 78 exceeded the optimizer’s size limit. Across thresholds and measurement spaces, 52 records were classified as periodic and none remained unresolved at the checked horizon. The exact riffle calculation, the independent student-code checks, the 180 repeated link-policy cases, and the later 73-run matched ring comparison are additional to that historical inventory. Both link studies now run through the same MixBench workflow as the other presets. The later reruns reproduced the main crossing times for the shared paths, endpoint scaling, Glasgow and karate. [MixBench notebook](notebooks/mixbench.ipynb).

Two numerical checks helped interpret close results. On four-station Glasgow sections, rounding put an optimized walk’s step-3 TV at 0.2500000000000001, so the strict comparison returned step 4. The difficult optimized 64-ring case with shortcut (51, 63), using zero-based IDs, needed a warm-start retry. Its final solver residuals and objective gap passed the acceptance checks, and the spectrum was checked against the exported transition matrix. The current release suite passed all 16 tests, including the replay of archived link-policy fixtures. [MixBench notebook](notebooks/mixbench.ipynb) and [solver check](tests/solver_retry.py).

The [MixBench notebook](notebooks/mixbench.ipynb) explains and runs the regression checks, displays the experiment plots, and saves the results. The [project README](README.md) gives commands to rerun each preset and describes the saved outputs. The [20-minute presentation](presentation/mixing.pdf) gives the short version with the network diagrams and paper citations.
