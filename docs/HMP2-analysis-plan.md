# Analysis protocol

## Research question

In people with Crohn’s disease, how does the relative abundance of common gut genera change between consecutive stool collections when a participant’s whole-grain report moves toward or away from “yesterday”? Each direction has its own comparison with intervals in which the same person started in the same dietary category and stayed there.

The primary estimates describe within-person associations. The dietary question records when the participant last ate whole grains. It provides neither portion size nor the sequence of meals between kits. Illness, treatment and other dietary changes can influence both the report and the microbiome.

## Data and analytic cohorts

Link Crohn’s disease stool metagenomic records to taxonomic profile columns by exact external assay ID. Retain the 187 source rows whose terminal taxonomic rank is genus. Require finite, nonnegative abundances and a positive full-profile total. For each participant, order all observed metagenomic receipt days and pair only consecutive days with eligible endpoints, increasing nominal study week and a receipt gap of 7–28 days inclusive. A date with multiple assays, an unusable profile or a missing whole-grain answer interrupts the sequence. An interior kit can contribute to two intervals.

Code either “yesterday” response frequency as Y. Code “2–3 days ago”, “4–7 days ago” and “none in the past seven days” as O. A stable interval retains the same binary category, even if the exact response changes within O.

| Cohort | Observations | Participants |
| --- | ---: | ---: |
| Crohn’s disease stool metagenomics | 750 assays on 736 participant-days | 65 |
| Unambiguous days with a positive full profile | 722 days | — |
| Days also having a valid whole-grain report | 710 days | 64 |
| Consecutive eligible pairs, full profile | 497 intervals from 641 distinct endpoints | 52 |
| Pairs with positive common-panel abundance at both endpoints | 492 intervals | 52 |
| Primary pairs after requiring both comparison arms within person and starting category | 351 intervals | 35 |

Select common genera by detection on at least 361 of the 722 valid days, independently of dietary responses. The panel comprises 22 source taxa, including two unclassified genus bins. Three endpoint kits have zero abundance across the panel, removing five intervals. Before this exclusion 359 intervals in 36 people have same-start comparators. Requiring comparator support again removes three additional intervals. All restrictions preserve the original adjacency.

**Panel in reporting order.** Alistipes, Anaerostipes, Bacteroides, Bifidobacterium, Blautia, Clostridium, Collinsella, Dialister, Dorea, Escherichia, Eubacterium, Faecalibacterium, Flavonifractor, Fusicatenibacter, Lachnoclostridium, Lachnospiraceae_unclassified, Odoribacter, Parabacteroides, Roseburia, Ruminococcaceae_unclassified, Ruminococcus and Ruthenibacterium. Use the full taxonomic identifiers for joins.

## Primary analysis

For each kit, divide each genus abundance by the total across all 187 genera to obtain $p_g$. For a common genus $g$ in the fixed panel $G$, compute its log ratio relative to the panel geometric mean with offset $a=10^{-6}$.

$$
y_g=\log_2(p_g+a)-\frac{1}{22}\sum_{h\in G}\log_2(p_h+a).
$$

For each interval from earlier kit A to later kit B, calculate $\Delta y_g=y_{Bg}-y_{Ag}$. For participant $i$ and direction $d$, let $C_{igd}$ be the mean change over switched intervals minus the mean change over same-start stable intervals. Average $C_{igd}$ equally across people with at least one interval in each arm to obtain $\hat\theta_{gd}$. Keep the chronological sign in both directions.

| Direction | Switched arm | Stable arm | Participants |
| --- | ---: | ---: | ---: |
| Toward yesterday, starting O | O → Y, 57 intervals | O → O, 181 intervals | 29 |
| Away from yesterday, starting Y | Y → O, 55 intervals | Y → Y, 58 intervals | 21 |

Fifteen participants contribute in both directions. A positive contrast means a greater relative increase or a smaller relative decrease in the switched arm. On the multiplicative scale, $2^{\hat\theta_{gd}}$ compares the geometric relative fold changes in the two arms. Absolute microbial load is unmeasured. Contrasts of opposite sign alone do not establish biological reversibility.

### Sampling uncertainty and Bayesian magnitude bounds

Report both directions for every genus. Resample the 35 participants with replacement 5,000 times. Each sampled row carries all 22 contrasts and both directions it supports. Average eligible copies within direction. Use R’s Mersenne-Twister, Inversion and Rejection settings, seed `20261129`, and type-7 quantiles for pointwise 95% percentile intervals.

For a simultaneous band over the 44 primary coordinates, let $s_j$ be the bootstrap standard deviation for coordinate $j$. In each draw compute the maximum absolute deviation from the corresponding point estimate, standardized by $s_j$, across all 44 coordinates. Its 95th percentile $c$ gives the band $\hat\theta_j\pm cs_j$. An empty direction in a draw requires investigation before computing intervals. If any $s_j$ is zero or nonfinite, retain valid marginal intervals and report the joint band as unavailable.

The **participant-level Bayesian bootstrap** summarizes plausible magnitudes of the same mean contrast. Draw 20,000 shared vectors of Dirichlet$(1,\ldots,1)$ weights on the 35 people, using seed `20261122`. Renormalize weights over the contributors in each direction, then calculate the weighted means. Report $F=2^{Q_{0.95}(|\theta_{gd}|)}$ alongside each signed estimate. The range $[1/F,F]$ contains approximately 95% of posterior mass for that coordinate’s ratio under exchangeable participants, conditional on their observed intervals. These are pointwise magnitude bounds. The simultaneous band addresses multiplicity. [Rubin, The Bayesian Bootstrap](https://doi.org/10.1214/aos/1176345338)

### Sensitivity analyses

| Change | Analysis and support |
| --- | --- |
| Zero handling | Repeat all 44 primary contrasts with $a=10^{-4}$ and identical resampled rows. Report pointwise intervals and a separate simultaneous band on this scale. |
| Recall window | Move the boundary to within three days versus longer ago. Re-establish same-start support and report both directions with pointwise intervals. Toward comprises 52 switched and 99 stable intervals in 24 people. Away comprises 54 switched and 123 stable intervals in 23 people. The union has 328 intervals in 35 people. |
| Antibiotics | Require an answered medication screen and a negative antibiotic flag at both endpoints. Re-establish support at the yesterday boundary and report pointwise intervals. Toward comprises 39 switched and 115 stable intervals in 24 people. Away comprises 34 switched and 40 stable intervals in 16 people. The union has 228 intervals in 29 people. |

## Clinical and dietary co-change

For all four primary transition groups, report participants, intervals, receipt-gap distributions, exact starting recall categories and variable-specific missingness. Summarize both endpoint antibiotic reports, later-kit diarrhoea, procedure and hospital stay, paired diarrhoea transitions, Harvey–Bradshaw Index (HBI), liquid or very soft stools over 24 hours and transitions in the other 20 dietary items. The medication, diarrhoea, procedure and hospital questions cover the preceding two weeks. Depending on the gap, the two recall windows can overlap or leave days uncovered. A missing medication parent makes an exported negative antibiotic flag unknown. Flag a positive parent with all medication children negative. For HBI require a completed, finite, nonnegative integer score and distinguish repeated from distinct clinical visit numbers. Exclude blank, sentinel 999 and invalid liquid-stool counts from that measure.

Report each arm’s starting HBI and B−A change on the same complete pairs. For paired diarrhoea reports, summarize onset and resolution separately as percentages of all complete intervals in each arm. Calculate arm differences with the primary participant weighting, restricted to people with usable records in both arms for each measure. Report pointwise intervals from 5,000 shared whole-participant resamples of the primary union. For the other foods, show each item with its own denominator and the count of changes on intervals complete for all 20. These descriptions inform interpretation of concurrent symptoms and treatment. Missing clinical responses affect the corresponding summaries only.

## Forward prediction benchmark

Across all 497 intervals, ask whether adding the whole-grain report at B improves estimation of the concurrent 187-genus profile for a participant held out of training. Use root proportions $h_g=\sqrt{p_g}$ and predict $h_{Bg}-h_{Ag}$ separately for each genus.

| Model | Predictors |
| --- | --- |
| Persistence | Carry forward the observed profile at A. |
| Baseline ridge | Unpenalized intercept, the genus’s own $h_{Ag}$, elapsed weeks and Y status at A. |
| Dietary update ridge | Baseline inputs plus O → Y and Y → O indicators. |

Use four participant-disjoint outer folds. Sort the 65 source participant IDs under `LC_COLLATE=C.UTF-8`. Assign `sample(rep(1:4, length.out = length(ids)))` with the R random-number settings above and seed `20261121`, then join folds to eligible intervals. The four folds contain 14, 11, 14 and 13 eligible people and 146, 86, 123 and 142 intervals, respectively.

In each training set, give each participant equal total weight. Center and scale predictors using training weights normalized to sum to one, center the response, and assign zero slope to constant predictors. Within each outer training set, use the three remaining fold IDs in turn for inner validation. For each model, choose one ridge penalty across genera from 0.001, 0.01, 0.1 and 1 by minimizing the pooled equal-participant validation loss, preferring the larger penalty in an exact tie. Refit on all outer training participants.

For validation and scoring, add predicted changes to the earlier root profile, clip negative components to zero and normalize to unit Euclidean length. A zero norm defaults to the earlier profile. Score with squared Hellinger loss, $\tfrac12\sum_g(\hat h_{Bg}-h_{Bg})^2$. Report all three held-out losses and the paired gain (baseline loss minus updated loss), averaging intervals within person first. Include fold gains and penalties. Estimate a pointwise interval by resampling the stored participant gains within outer fold 5,000 times with seed `20261129`. It conditions on the fitted models and folds. The later dietary answer is known at scoring time, so gain measures contemporaneous predictive information.

## Extensions of the dietary question

**Whole-recall change and community turnover.** Among the 497 intervals, 412 in 52 participants have all 21 dietary items answered at both ends. Let $K$ count items whose yesterday status changes, including whole grain. Regress Bray–Curtis turnover $D=\tfrac12\sum_{g=1}^{187}|p_{Bg}-p_{Ag}|$ on $K/3$, elapsed weeks and participant fixed intercepts by weighted least squares, with each participant contributing total weight one. Report the change in turnover per three additional item switches. Fifty participants with 410 intervals show within-person variation in $K$. Repeat under the antibiotic restriction (286 intervals in 48 people before checking variation). Recheck model rank for each fit. Resample whole people, assigning separate intercepts to repeated copies. If the original slope is unidentified, report it as such. If a bootstrap fit cannot identify it, report the estimate and failed-draw count and withhold the interval.

**Stool butyrate.** Use the source-annotated HILIC-negative feature `HILn_QI36` at Crohn’s disease stool visits with unique collection and external assay IDs. All 265 source visits have finite positive processed intensity and 261 have valid whole-grain recall. For the 31 people with both states (63 Y visits, 101 O), subtract the mean log₂ intensity at O from the mean at Y within person, then average equally across people. Report the log₂ difference and its fold ratio. Intensity is in relative instrument units, and sample-level batch identifiers are unavailable in the BIOM metadata. This comparison describes a metabolite association and cannot establish production or mediation.

For both extensions use 5,000 whole-participant resamples, seed `20261129`, and pointwise 95% percentile intervals. Report the ridge gain, turnover slope and butyrate contrast regardless of the primary findings. The simultaneous band applies only to the prespecified primary family.
