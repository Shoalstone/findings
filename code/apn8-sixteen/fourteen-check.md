# Independent check, i = 3: is there a quadratic APN F on F_2^8 with distribution [0^226, 2^15, 4^14]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `fourteen_check.cpp`,
> `fourteen_check_levels.py` → `src/`; `runs/check14/spreads.txt` → `data/spreads.txt`; the other outputs in
> `runs/check14/` → `logs/*-check14-*.log` and `logs/check14-ctl-*.txt` (`.prog` files are resume state, not
> included); `runs/bin/fourteen_check*.exe` were builds of the one source, build it from `src/` (`README.md`);
> `runs/a/fourteen.txt` (Rowen's run) → `logs/fourteen-output.txt`. The known-function lists (`apn_8bit.txt` and the
> rest) aren't included; `README.md` says where they come from.

Subagent (Claude Opus 5.5), 2026-10-24, for Rowen. Written without reading `sixteen.cpp` or `runs/a/fourteen*.txt`
until the verdict below was in. Files: `fourteen_check.cpp` (all parts, controls), `fourteen_check_levels.py`
(Python recount of the first levels), scratch and outputs in `runs/check14/`. Utilities (form generators, radicals,
the LUT/Walsh verifier) adapted from the earlier checker's `sixteen_check.cpp`.
Status: done. Verdict in section 6: no such function.
Builds: `runs/bin/fourteen_check*.exe` are successive builds of the one source; the search, legality and generator
code didn't change between them (later edits added the control modes, the oracle candidates they use, and an APN
filter for list lines), so every run above goes through the same main path.

## 1. The reduction, checked

Setup as in `sixteen-check.md` §1: S8 an 8-dim space of alternating forms with every a ≠ 0 in the radical of exactly
one nonzero form. i = 3: 226 bent, 15 rank 6, 14 rank 4; 14·15 + 15·3 = 255. The 14 rank-4 radicals are pairwise
complementary (a partial spread of solids), D = the 45 points they miss, and the 15 rank-6 radicals (lines) partition D.

**Hole structure.** Rowen leaned on Kurz, *Divisible codes*, Table 9.3 (no maximal partial spread of 14 solids in
PG(7,2)). I checked it by computation rather than cite it, and it holds. The i = 2 by-hand route doesn't carry over:
- Character sums: a hyperplane contains at most one rank-4 radical, so |H ∩ D| = 127 − 98 − 8k ∈ {29, 21} and
  W(u) = Σ_{x∈D} (−1)^{u·x} ∈ {13, −3}, with 13 taken 45 times (Σ_{u≠0} W = −45). So Cay(F_2^8, D) is
  srg(256, 45, 16, 6) and D is a projective [45, 8, {16, 24}] two-weight set. Unlike i = 2, that doesn't pin D down:
  there are two such sets (Kurz §11, citing an enumeration [113]): three disjoint solids, and one other.
  (Side remark: Kurz's Example 9.14 calls the other one the Baer-type construction, but for q = 2, t = 4 that
  construction, F_4^* · F_2^4 inside F_4^4 = F_2^8, is λF_2^4 for the three λ ∈ F_4^*, i.e. exactly three disjoint
  solids. So, as I read it, the label in the text is a slip; the other code is the one in [99, Thm 4.1]. Nothing below
  depends on this.)
- Also cheap but not enough on its own: every nonzero coset of each rank-4 radical S holds exactly 3 holes (the other
  13 radicals meet the 5-space S + ⟨c⟩ in one point each).
- **Computation (`fourteen_check holes 11`).** Put three of the 14 solids at Y = {x = 0}, X = {y = 0}, Z = graph(I)
  (any ordered triple of pairwise complementary 4-spaces can be moved there). The other 11 are graph(N) with N in
  V5824 = {N : N, N + I invertible}, pairwise with N + N' invertible: an 11-clique of the graph on V5824. Up to
  conjugation by GL(4,2) (which fixes X, Y, Z): let c be the smallest conjugacy class present (V5824 is 5 classes,
  sizes 1680, 1344, 1344, 1344, 112), conjugate one class-c element to the class representative r_c; the other 10
  are a 10-clique in adj(r_c) ∩ {class ≥ c}. Every such 10-clique (863,186 in all) was checked to have a triangle in
  the common neighbourhood of all 11, i.e. three more solids completing a spread. **0 failures.** Count check: the
  weighted count Σ_c |class c| · Σ_K 1/m(K) (m = number of class-c elements) gives 410,906,496 11-cliques, which is
  exactly 364 · 1,128,864 = C(14,11) times the number of 14-cliques (spread sets {0, I, ...}, below): each 11-clique
  lies in exactly one 14-clique.
So **every partial spread of 14 solids in PG(7,2) extends to a spread**, the holes are three pairwise complementary
solids H1, H2, H3, and Kurz's Table 9.3 entry (0 maximal sets for r = 14) is confirmed independently.

**Rank-6 radicals.** With H1, H2, H3 = X, Y, Z, the lines inside D are the 105 lines inside one solid and 15
transversals {(x, 0), (0, x), (x, x)}, one through each point of D: for a ∈ X, b ∈ Y, a + b ∈ D only if a + b ∈ Z,
i.e. a = b in the identification. So a transversal is always possible; the 15 rank-6 radicals are t transversals
(over a point set T ⊂ F_2^4 \ 0, the same T in all three solids) plus lines of each solid partitioning the rest.
The search below allows all 120 lines; no case analysis on t.

**Normalisation.** H1, H2, H3 → X, Y, Z as above (constructively: basis a_i of H1, b_i the H2-part of a_i in
H2 ⊕ H3). The 14 rank-4 radicals become graph(N_j), {N_j} a 14-clique of V5824 ({0, I, N_j} a spread set with Y as
∞). The residual group is G = {diag(h, h)} × S3, order 120,960: diag(h, h) acts by conjugation, S3 permutes X, Y, Z
and is generated by the swap (x, y) ↦ (y, x) (graph(N) ↦ graph(N⁻¹)) and (x, y) ↦ (x, x + y) (graph(N) ↦
graph(N + I)). The S3 part matters: the holes are unordered.

## 2. Spread configurations (`fourteen_check spreads`)

14-cliques of V5824, reduced by the same minimal-class trick (13-cliques in adj(r_c) ∩ {class ≥ c}): 2,873
conjugation-reduced cliques; weighted count **1,128,864** (Rowen's spread-set count for i = 1, reached here
independently). Each reduced clique's canonical form (lexicographically least sorted image over all 120,960 elements
of G) gives **73 G-orbits**, with stabiliser orders from the same pass; the orbit sizes sum to 1,128,864, matching the
weighted count. One orbit is the Desarguesian (field) one. So up to GL(8,2) there are 73 configurations
(spread, unordered triple of hole solids), and the search runs once on each representative (`runs/check14/spreads.txt`).

## 3. Method: point-branching on the forms, radicals fixed

Per representative the radicals are known, so the search only chooses forms:
- State: the span S (all 2^d elements) and the set of points covered by radicals of its nonzero elements. Every
  element must be legal: bent, or rank 4 with radical exactly one of the 14 graph(N_j), or rank 6 with radical inside
  D; radicals pairwise disjoint. Radicals are computed from each form directly, never taken from how it was built.
- Branch on the first uncovered point. Off D, the form of S8 covering it has rank 4 and radical the spread element
  through it: 28 candidates (pullbacks of the 28 nondegenerate forms on V/R). In D it has rank 6 with radical one of
  the 8 lines of D through it: 8 · 13,888 candidates. Complete (each solution's covering form is a candidate) and
  each solution is reached along one path. Points off D come first.
- All points covered at d < 8: all-bent completion over the 2^(28−d) coset representatives. d = 8 with points
  uncovered: dead.
- Every node where all 210 points off D are covered is logged (the analogue of Rowen's SHORT spans); any legal span of
  dimension ≤ 8 containing a form on each of the 14 radicals would produce one (follow it down the tree).

## 4. Results

`fourteen_check search spreads.txt asc.prog asc` (job 20261024-195115-check14-asc, 23 s wall on E-cores):
**0 solutions.** Nodes by span dimension 0..8: 73 / 2,044 / 57,232 / 740,398 / 1,196,814 / 64,384 / 845 / 0 / 0;
5.8e7 candidate forms tried. Nothing reached dimension 7; **no node ever had every point off D covered** (so no
in-D branching and no bent completion ever ran).

Cross-runs, each **0 solutions**, nothing at dimension 7, no node with D's complement covered:
- descending point order (a different branching point at almost every node): 73 / 2,044 / 57,232 / 743,146 /
  1,182,510 / 47,172 / 265.
- random point order (seed 7, off-D points shuffled, then D shuffled): 73 / 2,044 / 57,232 / 741,508 / 1,189,451 /
  48,006 / 270.
- **no S3 reduction**: all 2,873 conjugation-reduced 14-cliques instead of the 73 G-orbit representatives (so the
  result doesn't lean on the S3 part or the canonical-form computation, only on the minimal-class trick, which the
  weighted count checks): 2,873 / 80,444 / 2,252,432 / 29,126,873 / 46,816,683 / 2,164,216 / 23,703; 2.25e9 tries,
  14 min.

## 5. Controls

- `selftest` over all 73 representatives (0 failures): 28 / 13,888 nondegenerate forms on F_2^4 / F_2^6,
  |GL(4,2)| = 20160, |V5824| = 5824; per representative |D| = 45, the 14 radicals and D partition V \ 0, the clique
  is pairwise adjacent, D holds exactly 120 lines (105 inside a solid, 15 transversals), 8 through each hole; each
  radical's 28 generated forms are distinct with radical exactly that radical; for every 8th representative, every
  line's 13,888 forms are distinct, rank 6, with radical exactly that line (16.7M forms checked). Brute force over
  all 2^28 alternating forms: exactly 28 have radical R_0 of representative 0.
- S3: the 6 vertex permutations are distinct and closed under composition; in the controls the same group is built
  from the two linear maps and must carry configurations onto representatives (it always did).
- First levels recounted independently (`fourteen_check_levels.py`: forms by the block formula P^T d P, P = [N | I],
  d over the 28 nondegenerate 4x4 alternating matrices found by their Pfaffian, radicals by brute force over all 256
  vectors, its own legality test): spread 0, ascending order, nodes at dimension 0..4 = 1 / 28 / 784 / 10,234 /
  15,917, identical to the C++.
- **Known functions through the same code path** (`ctl`, `ctlbatch`). Take a known quadratic APN function with
  k ≥ 2 rank-4 components. Its rank-4 radicals are moved to X, Y (and graph(I) if k ≥ 3), a spread through them is
  found (a 14-clique of V5824 containing the remaining ones), three other spread elements (choice varied per run)
  are made the holes and mapped to X, Y, Z, and then an element of G is searched for that carries the 14
  remaining spread elements onto one of the 73 representatives in `spreads.txt`. That configuration-matching step
  is itself a test of section 2: every configuration must land on a representative. Then the main search runs on
  that representative, from the empty span, with the main generators, legality test, point order and completion,
  relaxed only by: a rank-6 radical outside D is legal if it is one of the known function's own, and those forms are
  also offered as candidates at their points. (The known function has k ≤ 4 of the 14 rank-4 forms; the other 10
  radicals' points are covered by its rank-6 forms, which is why the relaxation is needed.) Pass = the known space is
  among the spaces found. Every space found is re-verified by the LUT/derivative/Walsh check.
  - First run, BL22 line 12631 (186/65/4; 3 of its rank-6 radicals inside D, 62 outside): from the empty span the
    search finds exactly one space, the known one (APN, 186/65/4), in 0.1 s, at dimension 8.
  - BL22 (`apn_8bit.txt`, every function with ≥ 2 rank-4 components: 480 functions, 3 hole choices each): **1,440 of
    1,440 re-found** (924 with k = 2, 438 with k = 3, 78 with k = 4), 0 configurations failing to match a
    representative, 0 other spaces found; 20 of the 73 representatives hit.
  - BL22 again with a randomised spread through the radicals and 24 hole choices per function, for the 172 with
    k ≥ 3: **4,128 of 4,128 re-found** (3,504 with k = 3, 624 with k = 4), 68 of the 73 representatives hit.
  - BLP22 (`blp22_8bit.txt`, 37 functions with k = 2, 6 hole choices): 222 of 222.
  - The 3.8M (`new_apns.txt`, all 7,199 with k ≥ 2: 7,136 with k = 2 and the 63 with k = 3, one hole choice each):
    **7,199 of 7,199**.
  - Yu–Perrin QAM (`yuperrin_qam.txt`, 4 with k = 2, 6 choices): 24 of 24. Yu–Wang–Li + Edel–Pott
    (`yuwangli_qam.txt`, 12 with k = 2, 6 choices): 72 of 72; its line 3 isn't an APN space as read (the
    list's known odd line) and was skipped.
  - Across all control runs: no configuration failed to land on a representative, no known space was missed, and no
    other space was ever found. Representatives hit at least once: 68 of 73 (not 18, 22, 48, 49, 68; the code path
    is the same for all of them).

## 6. Verdict

**No quadratic APN function F_2^8 → F_2^8 has amplitude distribution [0^226, 2^15, 4^14]** (i = 3, 14 rank-4
components). It rests on:
- the hole structure, proved here by computation rather than cited (section 1): every 11-clique checked, and the
  count identity 410,906,496 = C(14,11) · 1,128,864;
- the configuration count (section 2): 73 G-orbits whose sizes sum to the weighted count 1,128,864, and the
  no-S3 cross-run over all 2,873 conjugation-reduced cliques giving the same answer;
- complete point-branching (each solution's covering form is a candidate at every step; generators checked against
  brute force), run in four shapes (ascending, descending, random order; no S3), all 0;
- positive controls through the same code path: every known function with ≥ 2 rank-4 components in every list,
  re-found from the empty span after being carried onto a representative, with no stray spaces.
The failure is early and uniform: in no configuration does any span reach dimension 7, and no legal span of
dimension ≤ 8 contains forms on all 14 rank-4 radicals (no node ever had every point off D covered). The rank-6
part of the problem (the in-D branching with 13,888-way candidates, and the transversal-line question) never comes
into play. Confidence: high.

## 7. Comparison with Rowen's run (read after the above)

`runs/a/fourteen.txt` (job 20261024-165047-fourteen, still running when last read at ~22:00): 140,000 of the 1,128,864
spread sets done (each with all 14 choices of the third hole), 1.93e9 nodes, **no SHORT lines**. That matches
this search exactly where the two can be compared: SHORT (a span of the 14 rank-4 forms covering all their radicals
in < 8 dims) is the event logged here as "every point off D covered", and it never occurs in any configuration. His
node rate, ~980 per (spread set, hole choice) with the form on graph(I) fixed, is the same order as mine, ~1,000 per
(representative, first form). His run is ~12% through after 4.6 h on shared E-cores, so roughly 30 h more; the orbit
reduction here (73 configurations instead of 15.8M ordered ones) covers the same cases in 23 s. Whether to let it
finish as the second independent search is his call; if he'd rather have his own
engine finish quickly, `runs/check14/spreads.txt` lists the 73 representatives (holes X, Y, Z = graph(I); his
normalisation puts graph(I) among the rank-4 radicals instead, so they'd need moving).
