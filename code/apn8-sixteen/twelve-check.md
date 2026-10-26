# Independent check, i = 5: is there a quadratic APN F on F_2^8 with distribution [0^218, 2^25, 4^12]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `twelve_check.cpp`,
> `twelve_check_levels.py` → `src/`; `runs/check12/orbits.txt` → `data/orbits12.txt`; `runs/check12/configs.txt` →
> `data/configs12.txt`; `runs/check12/maxc.txt` (5.6 MB) isn't included (`twelve_check parta` writes it, `README.md`);
> the other outputs in `runs/check12/` → `logs/check12-<name>-output.txt` (`<name>.out`), `logs/check12-ctl-*.txt`,
> `logs/check12-levels-output.txt` (`levels.txt`), and `logs/check12-sp-output.summary.txt` (`sp.out`'s header and
> summary; its 6,643 per-configuration lines are left out); the share runner's job logs (`check12-asc` and so on) →
> `logs/2026…-check12-*.log` (`.prog` files are resume state, not included); `runs/check14/spreads.txt` →
> `data/spreads.txt`; `runs/check13/configs.txt` → `data/configs.txt`; `runs/bin/twelve_check*.exe` were builds of
> the one source, build it from `src/` (`README.md`). Rowen's side: `pspread.cpp`, `cover.cpp` → `src/`;
> `runs/i5/holes_iso.py` → `src/holes_iso12.py`; `runs/i5/p12.txt` → `data/p12.txt`; the other outputs in `runs/i5/` →
> `logs/pspread1*-output.txt` and `logs/twelve-output.txt`; README "Twelve" is Rowen's working notes (not published),
> summarised in `README.md` here and in the write-up. The supplementary runs it once described as still going have
> finished; their results are folded into the text (marked as Rowen's) and listed in `README.md`. The known-function lists aren't included (`README.md`).

Subagent (Claude Opus 5.5), 2026-10-25, for Rowen ("search B" for i = 5). Written without reading `pspread.cpp`,
`cover.cpp`, `cover_prep.py`, README's i = 5 / "Twelve" parts, `runs/i5/` or any of Rowen's i = 5 output until the
verdict below was in. Files: `twelve_check.cpp` (all parts, controls), `twelve_check_levels.py` (Python recount of the
first levels), scratch and outputs in `runs/check12/` (`orbits.txt` = Part A's orbit list, `maxc.txt` = Part A's
maximal 9-cliques, `configs.txt`, `*.out`, `*.prog`, `ctl-*.txt`). Built on `thirteen_check.cpp`: utilities, V5824,
clique enumeration, point-set isomorphism, exact covers, canonical forms over frames, the forms search and the
controls carried over, generalised to any number of members.
Status: verdict in (section 6: no such function). Four supplementary runs were queued behind the share runner's
holds when this was written; all have since finished and agree (sections 4 and 5, last one 2026-10-26 04:38; added by Rowen).
Builds: `runs/bin/twelve_check.exe` ran Part A, configs, every search and the ordinary controls;
`twelve_check2.exe` is the same source plus the in-D control mode (section 5) and per-batch stats in `ctlbatch`.
`twelve_check3.exe` adds `maprowen` (section 7, after the verdict). Search, legality, generator and canonical-form
code are identical in all three.

## 1. The reduction

Setup as in `sixteen-check.md` §1. i = 5: 218 bent, 25 rank 6, 12 rank 4; 12·15 + 25·3 = 255. The 12 rank-4
radicals are pairwise complementary (a partial spread of 12 solids), D = the 75 points they miss, and the 25 rank-6
radicals (lines) partition D. Every point of a radical R_j can only be covered by a rank-4 form with radical exactly
R_j (rank-6 radicals lie in D), so all 12 rank-4 forms are forced by coverage, and the search needs nothing beyond
legality and covering every point.

A hyperplane contains at most one R_j, so it meets D in 43 or 35 points: D is a projective [75, 8, {32, 40}]
two-weight set, and (as at i = 4) the number of lines inside D is the same for every candidate: 325 in all 376
configurations below (13 through each hole).

**Normalisation.** Three of the 12 solids go to Y, X, Z = graph(I); the other 9 are graph(N), N in V5824, pairwise
N + N' invertible: a 9-clique. GL(4,2) conjugation fixes X, Y, Z and acts on the cliques.

## 2. Configurations up to GL(8,2), by computation

**Part A (`twelve_check parta`).** Every 9-clique up to conjugation by the minimal-class trick (weights 1/m kept exact
as multiples of 2520): **5,269,153 reduced 9-cliques**, weighted **2,811,561,984**. For each, Q = the common
neighbourhood in V5824 (= every solid disjoint from all 12 members), and each clique is canonicalised.
- Extension counts, exact (each ties to an earlier, separately computed number):
  - Σ (number of 5-cliques in Q) = **2,259,985,728 = C(14,9) · 1,128,864** (completions to a spread; 1,128,864 =
    the spread sets {0, I, ...}, fourteen-check §2).
  - Σ #{q in Q : P + q is a maximal 13} = **625,584,960 = 10 · 62,558,496** (thirteen-check's non-extendable count).
  - Σ #{q in Q : P + q extends to a spread} = **11,299,928,640 = 10 · 1,129,992,864** (thirteen-check's extendable
    count). Σ #{q : P + q neither} = 0.
  So extension is indeed not unique here: weighted, 2,258,825,184 9-cliques extend to a spread, short of
  C(14,9) · 1,128,864 = 2,259,985,728 (thirteen-check §8's numbers, reproduced), and the counts above count every
  extension. By |Q|: 0: 49,896,000; 1: 427,996,800; 2: 74,844,000; 5: 2,209,764,480; 6: 47,900,160; 10: 1,160,544.
- Types: **spr** (extends to a spread), **ext13** (Q nonempty, but every P + q is a maximal 13), **max** (Q empty).
  Weighted 2,258,825,184 / 502,840,800 / 49,896,000 (reduced 4,128,284 / 1,040,086 / 100,783), thirteen-check's
  `holesk 9` sizing exactly.
- **376 orbits: 280 spr, 84 ext13, 12 max.** Canonical form as thirteen-check §2 (`canonP`: frames = ordered member
  triples → (X, Y, Z), then diag(h, h); key = least image under the least class vector; |Aut P| = frames giving the
  key). **Per orbit**, the weighted number of 9-cliques landing on it equals 1320 · 20160 / |Aut P| exactly, for all
  376 (a stronger form of the sum identity: a missed orbit, a wrong |Aut| or a mis-merged pair would break it), and
  every reduced clique on an orbit had the same invariants (|Aut|, type, |Q|, extension counts). |Aut P|: spr 1..3600,
  ext13 1..48, max 2..24; 28 orbits have trivial Aut.
- The stored keys are canonical: each orbit's form recomputed from a random GL(8,2) image with shuffled members gives
  the same key and |Aut| (376 of 376; `configs`, and again in `selftest` with a frame carrying it back).

**Hole sets (`twelve_check configs`).** The 376 hole sets fall into **42 classes up to GL(8,2), 5 of them containing no
solid** (the maximal partial spreads), by isomorphism tests on the 42 distinct values of an invariant (pair degrees,
plane counts, hyperplane sizes), which happens to separate them all. **This is exactly Kurz's Table 9.3 entry for
r = 12 (42 point sets, 5 maximal).** By type: 21 hole sets of spr configurations (five disjoint solids), 16 of ext13,
5 of max; the type of every orbit agrees with whether its hole set contains a solid.
For every hole-set representative H, all exact covers of its complement by 12 solids, each canonicalised: every cover
lands on an orbit with that hole set, every such orbit is reached, and Σ_orbits |Aut H| / |Aut P| = number of covers,
for all 42 (|Aut H| by backtracking, independent of the canonical form). The five maximal ones: |Aut H| 32, 96, 32,
16, 2 with 4, 84, 4, 4, 1 covers and 1, 8, 1, 1, 1 orbits. (Kurz doesn't print r = 12 sets, so the comparison is the
count only.)

**Independent sources (`configs`).** Fourteen-check's 73 (spread, holes X Y Z) orbit reps × C(14,2) choices of two
more holes (6,643 configurations, which contain every spread-extendable 12-partial spread up to GL(8,2)): all land on
spr orbits and reach **all 280**. Thirteen-check's 31 non-extendable 13-configurations × 13 members removed (403):
none lands on a max orbit (as it can't), 37 on spr, 366 on ext13, reaching **all 84** ext13 orbits.

## 3. Method: point-branching on the forms, radicals fixed

As `thirteen-check.md` §3, with 12 members: every span element must be bent, rank 4 with radical exactly a member, or
rank 6 with radical inside D; radicals pairwise disjoint and computed from each form directly. Branch on the first
uncovered point: off D, the 28 forms on its member; in D, the 13,888 forms on each line of D through it disjoint from
what's covered. Off-D points come first, so rank-6 branching only starts at a node where all 12 rank-4 forms are
present (logged as "p1", the near-misses). All points covered at d < 8: all-bent completion; d = 8 with points
uncovered: dead. Complete (the covering form of each solution is always a candidate), each solution on one path.

## 4. Results

`twelve_check search configs.txt asc.prog asc` (job 20261025-092923-check12-asc, **20 s**): **0 solutions.**
Nodes by span dimension 0..8: 376 / 10,528 / 294,784 / 4,002,613 / 7,761,928 / 724,602 / 8,217 / 0 / 0
(spr 280 / 7,840 / 219,520 / 2,981,155 / 5,768,135 / 499,380 / 6,210; ext13 84 / 2,352 / 65,856 / 893,778 /
1,761,554 / 209,552 / 1,942; max 12 / 336 / 9,408 / 127,680 / 232,239 / 15,670 / 65); 3.6e8 candidate forms tried.
Nothing reached dimension 7, and **no node anywhere had all 12 rank-4 forms** (no near-miss, so no rank-6 branching
and no bent completion ever ran). The expected cost of the 13,888-way branching never arises: the failure is early
and uniform, as at i = 3.

Cross-runs, each **0 solutions**, nothing at dimension 7, no node with every off-D point covered:
- descending point order: 376 / 10,528 / 294,784 / 3,970,567 / 7,369,230 / 380,832 / 3,848 (20 s).
- random point order (seed 7): 376 / 10,528 / 294,784 / 3,979,633 / 7,462,326 / 442,121 / 5,272 (19 s).
- **without this classification, spread-extendable side**: fourteen-check's 73 spread reps × 91 (6,643
  configurations): 6,643 / 186,004 / 5,208,112 / 70,266,988 / 133,848,444 / 9,697,414 / 117,642 (8 min).
- **without this classification, ext13 side**: thirteen-check's 31 maximal 13s × 13 (403): 403 / 11,284 / 315,952 /
  4,288,357 / 8,486,696 / 1,618,431 / 9,352 (32 s).
- **without the canonical form, maximal side**: all 100,783 maximal 9-cliques as Part A enumerated them (only
  GL(4,2)-conjugation reduced). Finished 2026-10-26 04:38 (job check12-maxc, cancelled after 9,593 configurations, resumed as check12-maxc2;
  `maxc.out` + `maxc2.out`, ~2.6 h of running time): 100,783 / 2,821,924 / 79,013,872 / 1,065,082,714 /
  1,985,285,822 / 174,884,697 / 929,739, 0 solutions, nothing at dimension 7, no node with every off-D point
  covered. The max side was already covered by the 12 classified orbits and by the exact covers of the five maximal
  hole sets (section 2); this run checks it again without the classification.

## 5. Controls

- `selftest` (376 configurations, 0 failures): 28 / 13,888 nondegenerate forms, |GL(4,2)| = 20160, |V5824| = 5824;
  |D| = 75 and 12 members each; pairwise complementarity of every clique; 325 lines in D; each member's 28 forms
  distinct with radical exactly that member; for every 16th configuration all 325 · 13,888 line forms distinct, rank
  6, radical exactly the line (108M generated forms checked in all); brute force over all 2^28 forms: exactly 28 with
  radical R_5 of the last configuration; canonical-form invariance under a random GL(8,2) image (section 2).
- First levels recounted independently (`twelve_check_levels.py`, adapted from thirteen-check's: forms
  B = π(x)ᵀ d π(y) in a basis adapted to the radical, d over the 28 Pfaffian-nonzero 4×4 matrices, radicals by brute
  force, its own legality): configurations 0 (spr), 280 (ext13), 364 (max) give 1 / 28 / 784 / 10,804 / 21,176,
  1 / 28 / 784 / 10,711 / 19,628 and 1 / 28 / 784 / 10,594 / 18,535 through dimension 4, identical to the C++
  (`runs/check12/levels.txt`).
- **Known functions through the same code path** (`ctlbatch`), both paths for every function with k ≥ 2 rank-4
  components, as thirteen-check §5. Extendable path: radicals → X, Y (and graph(I)), a random spread through them, five
  other spread elements as holes, the 12-set canonicalised and matched to a listed representative (must be spr), the
  function carried there by the canonical frame. Non-extendable path: the radicals embedded into the members of an
  ext13 or max representative (cycled through all 96; for k = 4 a member 4-tuple of the right class). Then the main
  search from the empty span, relaxed only by: a rank-6 radical outside D is legal if it's the known function's own,
  and those forms are offered at their points. Pass = the known space is among those found; every space found is
  re-verified by LUT / derivatives / Walsh.
  - BL22 (`apn_8bit.txt`, all 480 with k ≥ 2, 3 choices per path): **1,440 of 1,440** on each path (924 / 438 / 78
    with k = 2 / 3 / 4), 0 other spaces, 0 setup failures; representatives hit: 210 of 280 spr, **96 of 96** ext13 +
    max.
  - BLP22 (37 with k = 2, 6 choices): 222 of 222 and 222 of 222. Yu–Perrin QAM (4, 6 choices): 24 of 24 and 24 of 24.
    Yu–Wang–Li + Edel–Pott (12, 6 choices): 72 of 72 and 72 of 72; line 3 isn't an APN space as read (the list's known
    odd line, as for the predecessors) and was skipped (the 6 + 6 "FAIL" lines in `ctl-yuwangli.txt`).
  - The 3.8M (`new_apns.txt`, the 7,199 with k ≥ 2, one choice per path): **7,199 of 7,199 re-found on each path**
    (7,136 / 63 with k = 2 / 3), 0 other spaces, 0 failures; representatives hit 241 spr and all 96 ext13 + max
    (`runs/check12/ctl-new.out` / `.txt`). *(Rowen, 2026-10-26 03:35, from the finished job check12-ctl.)*
  - Representatives: the batches finished by then hit 214 of the 280 spr and all 96 ext13 + max. The other 66 got a
    targeted run (`twelve_check2 ctlreps`, the first 40 BL22 functions with k ≥ 3 embedded into the members of each,
    `runs/check12/ctlreps.out`): **40 of 40 re-found for every one of the 66**, 0 other spaces. So every one of the
    376 configurations has had known functions carried onto it and re-found. *(Rowen, 2026-10-25 13:05, from the
    finished job check12-ctl2.)*
- **The rank-6 branching, exercised** (`twelve_check2 ctlbatch ... nond`). The main search never branched in D (no
  near-miss anywhere), and neither do the ordinary controls (a known function's out-of-D forms span all of it by the
  time every off-D point is covered). So a further control aims at that code path: take a rank-6 form g of the known
  function whose radical lies inside D, start from a random hyperplane W of the function not containing g, and put
  g's points first. Those points are covered in the function by g alone, so the search must find g among the main
  generator's 13,888 forms on its line, and pass legality and disjointness against all 128 elements of W. A trial on
  10 BL22 functions per path: 20 of 20 re-found, nothing else found, 1.06M in-D candidates tried. The full batch (BL22, all 480 with k ≥ 2, 3 choices per path, `runs/check12/ctl-ind-bl22.out`): **1,415 of 1,415**
  on the extendable path and **1,419 of 1,419** on the other (46 trials skipped because the starting span was already
  the whole space), 0 other spaces, 0 setup failures; 78.8M and 78.3M in-D candidates tried, every control branching at
  dimension 7 as designed. *(Rowen, 2026-10-25 13:05.)*
- Across all finished controls: no configuration failed to land on a listed representative, no known space was
  missed, no other space was ever found.


## 6. Verdict

**No quadratic APN function F_2^8 → F_2^8 has amplitude distribution [0^218, 2^25, 4^12]** (i = 5, 12 rank-4
components). It rests on:
- the classification of 12-solid partial spreads in PG(7,2) by computation (section 2): all 5,269,153 reduced
  9-cliques canonicalised; the per-orbit identity holds for all 376 orbits; the extension counts reproduce exactly the
  spread count (C(14,9) · 1,128,864, counting multiple extensions) and thirteen-check's two 13-counts (× 10); 42 hole
  sets, 5 of them maximal, as in Kurz's Table 9.3, with every hole set's exact covers accounted for orbit by orbit;
  and two outside sources (fourteen-check's spread reps × 91, thirteen-check's maximal 13s × 13) landing on, and
  reaching, every spr and every ext13 orbit;
- complete point-branching in five finished shapes (ascending, descending, random; spread reps × 91; maximal 13s × 13,
  the last two independent of this classification), all 0 solutions, nothing at dimension 7, no span ever holding
  all 12 rank-4 forms;
- controls through the same code path on both paths, on every list (3,516 on the four smaller lists + 14,398 from the 3.8M, all re-found, nothing else
  found), every one of the 376 configurations reached (66 spr by a targeted batch, 2,640 of 2,640), and a dedicated
  control of the rank-6 branching, which the main search never needed (2,834 of 2,834 on BL22);
- a separate first-levels recount agreeing node for node through dimension 4.
Like i = 3, the failure is early and uniform: no span reaches dimension 7 and no near-miss appears, so the 13,888-way
rank-6 branching expected to dominate never runs. Confidence: high. The maximal-clique cross-run without the canonical form (section 4, finished 2026-10-26) agrees: 0.

## 7. Comparison with Rowen's search (read after the above)

README "Twelve: i = 5" and `runs/i5/`, read after the above:
- **Classification.** Rowen's `pspread.cpp` (K = 12, 72 min) finds the same **376 configurations, 12 / 84 / 280**
  by extension (his w = 0 / 1 / 5 = my max / ext13 / spr), the same reduced counts (100,783 / 1,040,086 / 4,128,284)
  and weighted counts (49,896,000 / 502,840,800 / 2,258,825,184), orbit sums matching. Mapped directly
  (`twelve_check3 maprowen`, his `p12.txt` canonicalised by my code, `runs/check12/maprowen.out`): **a bijection**
  onto my 376, with |Aut| and type agreeing on every one. His maximal hole sets: 5 classes, |Aut H| 2, 96, 32, 32,
  16, covers 1, 84, 4, 4, 4: the same as mine. (His note that plane counts alone give 4 classes there; the invariant
  here, with pair degrees and hyperplane sizes added, separates all 42.) The two classifications share the method's
  shape (minimal-class reduction, canonical form over frames with the least class vector), which came down from
  thirteen-check's write-up; the checks on top differ (here: per-orbit identities, extension counts, covers of all
  42 hole sets, outside sources; his: K = 17, 14, 13 reproduced).
- **Search.** His `cover3` on all 376: **0 solutions**, no near-miss, 12,303,628 nodes in 49 s. Mine: 0 solutions,
  no near-miss, 12,802,672 nodes (roots not counted, as in his) in 20 s. Nodes by dimension 1..8, his 10,528 / 294,784 / 4,002,712 / 7,479,748 /
  509,108 / 6,747 / 1 / 0, mine (asc) 10,528 / 294,784 / 4,002,613 / 7,761,928 / 724,602 / 8,217 / 0 / 0: equal
  through dimension 2, then apart (different branching points; at i = 4 the two trees matched exactly only where the
  orders did). His tree has one node at dimension 7 and none of my five shapes has any. That isn't a conflict: a
  point-branching tree visits only spans whose elements cover the first uncovered points in its own order, so which
  legal 7-dimensional spans appear depends on the order, while each solution is reached in either.
- Same answer: no function, by two searches and two classifications that agree configuration for configuration.
