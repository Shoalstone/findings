# Independent check, i = 4: is there a quadratic APN F on F_2^8 with distribution [0^222, 2^20, 4^13]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `thirteen_check.cpp`,
> `thirteen_check_levels.py`, `cover.cpp`, `cover_prep.py`, `i4_holes.py` → `src/`; `runs/check13/configs.txt` →
> `data/configs.txt`; `runs/check13/kurz.txt` → `data/kurz13.txt` (the same file as Rowen's `runs/i4/kurz13.txt`);
> `runs/check13/holes10.txt` and `covers.txt` aren't included (regenerated in seconds, `README.md`); the other
> outputs in `runs/check13/` → `logs/*-check13-*.log` and `logs/check13-ctl-*.txt` (`.prog` files are resume state,
> not included); `runs/check14/spreads.txt` → `data/spreads.txt`; `runs/bin/thirteen_check*.exe` were builds of the
> one source, build it from `src/` (`README.md`). Rowen's side: `runs/i4/` → `data/` and `logs/*-i4-*.log`,
> `logs/i4-*-output.txt`, `logs/i3-reps-output.txt`; README "Thirteen" is Rowen's working notes (not published),
> summarised in `README.md` here and in the write-up. The known-function lists aren't included (`README.md`).

Subagent (Claude Opus 5.5), 2026-10-24, for Rowen. Written without reading `cover.cpp`, `cover_prep.py`,
`i4_holes.py`, README "Thirteen" or `runs/i4/` until the verdict below was in. Files: `thirteen_check.cpp` (all
parts, controls), `thirteen_check_levels.py` (Python recount of the first levels), scratch and outputs in
`runs/check13/` (`holes10.txt`, `configs.txt`, `covers.txt`, `*.prog`, `ctl-*.txt`, `kurz.txt` = Kurz's three
generator matrices cut from the paper text by `sed`, not retyped). Utilities (forms, radicals, V5824, clique
enumeration, verifier) copied from `fourteen_check.cpp`.
Status: done. Verdict in section 6: no such function.
Builds: `runs/bin/thirteen_check{,2..6}.exe` are successive builds of the one source. 2 added the `covers` /
`searchmem` input, 3 fixed a control-only bug (below), 4 `mapsp`, 5 `holesk` (sizing for other i, section 8),
6 `ctlreps`; search, legality and generator code didn't change.

## 1. The reduction, checked

Setup as in `sixteen-check.md` §1. i = 4: 222 bent, 20 rank 6, 13 rank 4; 13·15 + 20·3 = 255. The 13 rank-4 radicals
are pairwise complementary (a partial spread of solids), D = the 60 points they miss, and the 20 rank-6 radicals
(lines) partition D. In a solution every point of a radical R_j can only be covered by a rank-4 form with radical
exactly R_j (rank-6 radicals lie in D), so all 13 forms are forced by coverage; the search needs nothing beyond
legality and covering every point.

**Character sums.** A hyperplane contains at most one R_j, so it holds 36 or 28 holes (60 and 195 hyperplanes,
confirmed on every hole set below): D is a projective [60, 8, {24, 32}] two-weight set. As at i = 3 that doesn't pin
D down, and here it pins down even less than it looks: Cay(F_2^8, D) is strongly regular, so the line count (200,
10 through each hole) is the same for every candidate. Plane counts do separate them (below).

**Hole classification by computation (`thirteen_check holes`, `configs`).** Put three of the 13 solids at Y, X,
Z = graph(I); the other 10 are graph(N), N in V5824, pairwise N + N' invertible: a 10-clique. Up to GL(4,2)
conjugation by the minimal-class trick (as in `fourteen-check.md` §1; weights 1/m kept exact as multiples of 2520):
**2,358,714 reduced 10-cliques** (per class 1,339,160 / 740,385 / 126,321 / 72,768 / 80,080), 4 s. For each, Q =
the common neighbourhood in V5824 (the solids disjoint from all 13):
- Q nonempty: always |Q| = 4 and Q a 4-clique. Weighted count of these 10-cliques **1,129,992,864 = C(14,10) ·
  1,128,864**, and also = Σ (number of 4-cliques in Q): every extendable 13-partial spread lies in exactly one spread.
- Q nonempty without a 4-clique: **0** (agrees with i = 3: a partial spread of 14 always extends).
- Q empty (maximal partial spread): 141,165 reduced, weighted **62,558,496**.

Kurz's three printed sets (§9, after Table 9.3) checked from scratch: 60 distinct nonzero columns each, rank 8,
two-weight, no solid inside (tested against all 200,787 solids of PG(7,2)), pairwise non-isomorphic, stabilisers in
GL(8,2) of order **576, 14,400, 96** (backtracking over images of a basis, pruned by membership and per-point
pair/plane counts): his orders. They contain 12, 0 and 16 planes. Exact covers of each complement by 13 solids
(from the 1,264 / 1,280 / 1,280 solids inside it): **552, 10,440, 12**. So each really is the hole set of a maximal
partial spread of 13, and
  1716 · 20160 · covers / |Aut H| = 33,153,120 + 25,081,056 + 4,324,320 = **62,558,496**,
exactly the weighted count of non-extendable 10-cliques above (1716 ordered triples of members, 20160 = the
pointwise stabiliser of (X, Y, Z)). Directly as well: the canonical form (section 2) of every one of the 141,165
non-extendable cliques equals that of a Kurz cover, 0 exceptions, all 31 cover orbits reached; on a sample of 142,
an explicit isomorphism to the matching Kurz set was found every time.
The extendable hole sets fall into three types (four disjoint solids), with stabilisers **120, 288, 4320**, Kurz's
three four-solid codes; each holds 60 planes.

So Kurz's statement holds as he gives it: exactly six hole sets occur for 13 solids, three unions of four solids
(extendable) and his three printed sets (maximal). Table 9.3's r = 13 entry (6 point sets, 3 maximal) agrees. I
found nothing wrong in the passage. Not checked because nothing here uses it: that there are 12
[60, 8, {24, 32}] codes in all, and the Aut-720 one containing a solid. Kurz prints only hole sets; the counts of
partial spreads per hole set (552 / 10,440 / 12 covers, 15 / 15 / 1 orbits) aren't in his text (Rowen's
`cover_prep.py` found the same cover counts, section 7).

**Normalisation.** Unlike i = 3, three of the radicals (not holes) sit at X, Y, Z, and the canonical form ranges over
all choices of which three.

## 2. Configurations up to GL(8,2) (`thirteen_check configs`)

Canonical form of a 13-partial spread P: over every frame (ordered triple of members (A, B, C) → (X, Y, Z) by an
explicit basis change, then diag(h, h), h ∈ GL(4,2); these are all g ∈ GL(8,2) doing that), the sorted list of the
other 10 members as V5824 vertices. Key = the least such list over the frames whose ordering has the least sorted
conjugacy-class vector (an intrinsic condition, so the key is an orbit invariant; it's an image of P, so equal keys
mean equivalent configurations); only h moving a least-class member onto its class rep need be tried. Every frame
giving the key lies in that set, and those frames are in bijection with Aut(P), which gives |Aut P|.
Sources: (a) extendable: my own 14-cliques (2,873 conjugation-reduced, weighted 1,128,864) × 14 choices of the
fourth hole (holes X, Y, Z, graph(M), radicals the other 13): 40,222; (b) the 141,165 non-extendable Part A cliques;
(c) the 11,004 Kurz covers.
**190 configurations: 159 extendable, 31 non-extendable.** Identities:
- Σ 1716 · 20160 / |Aut P| = 1,129,992,864 (extendable) and 62,558,496 (non-extendable): both equal Part A's weighted
  counts. A missed orbit would leave a shortfall.
- Per Kurz set, Σ_orbits |Aut H| / |Aut P| = 552, 10,440, 12 (= the cover counts; links the backtracking |Aut H| to
  the canonical-form |Aut P|).
- By hole type: Kurz 1 / 2 / 3: 15 / 15 / 1 configurations; four solids with |Aut H| 288 / 4320 / 120: 34 / 58 / 67.
  |Aut P| ranges over 1..720 (six extendable configurations have trivial Aut).
- `selftest`: every configuration's canonical form recomputed from a random GL(8,2) image with the members shuffled
  gives the same key and |Aut| and a frame carrying it back (190 of 190).
- A second source for the extendable side (`mapsp`): fourteen-check's 73 (spread, 3 holes) G-orbit reps × 14 choices
  of the fourth hole (1,022, an independent route to every extendable orbit) canonicalise onto listed extendable
  representatives every time and hit all 159.

## 3. Method: point-branching on the forms, radicals fixed

As `fourteen-check.md` §3, with the members generic (masks, not graphs): every span element must be bent, rank 4 with
radical exactly a member, or rank 6 with radical inside D; radicals pairwise disjoint; radicals computed from each form
directly. Branch on the first uncovered point: off D, the 28 forms on its member; in D, the 13,888 forms on each line
of D through it disjoint from what's covered (every configuration has 200 lines in D). All points covered at d < 8:
all-bent completion; d = 8 with points uncovered: dead. Complete, each solution on one path. Nodes where every point
off D is covered (all 13 rank-4 forms present) are logged.

## 4. Results

`thirteen_check search configs.txt asc.prog asc` (job 20261024-222133-check13-asc, 111 s on E-cores): **0 solutions.**
Nodes by span dimension 0..8: 190 / 5,320 / 148,960 / 1,972,468 / 3,468,491 / 783,269 / 4,006 / 0 / 0 (extendable
159 / 4,452 / 124,656 / 1,650,663 / 2,900,658 / 622,874 / 3,569; non-extendable 31 / 868 / 24,304 / 321,805 /
567,833 / 160,395 / 437); 1.8e8 candidate forms tried. Nothing reached dimension 7.

**Four near-misses**, all in the most symmetric four-solid hole type (|Aut H| = 4320): configuration 45 (one) and
154 (three) each have a 6-dim span holding all 13 rank-4 forms (rank distributions 13 rank-4 + 50 bent, three times,
and once 13 rank-4 + 8 rank-6 + 42 bent). At each, every rank-6 form on every hole line through the first uncovered
hole is illegal (486,080 forms on 35 lines; correction by Rowen 2026-10-25 after a review: the 585,900 first given here
was the whole dimension-6 tries counter, which also counts 28-form tries at the other 3,565 dim-6 nodes), so they're dead. These are the i = 4 analogue of the i = 2 SHORT spans. The
in-D branching runs only there.

Cross-runs, each **0 solutions**, nothing at dimension 7, and the same four near-miss spaces (compared as subspaces):
- descending point order: 190 / 5,320 / 148,960 / 1,957,240 / 3,277,832 / 130,389 / 938.
- random point order (seed 7): 190 / 5,320 / 148,960 / 1,961,530 / 3,324,270 / 170,351 / 1,770.
- **without my classification, extendable side**: fourteen-check's 73 (spread, 3 holes) G-orbit reps
  (`runs/check14/spreads.txt`) × 14 choices of the fourth hole, 1,022 configurations covering every extendable
  orbit (with repeats): 1,022 / 28,616 / 801,248 / 10,543,898 / 18,202,346 / 1,152,254 / 12,084; four near-miss
  nodes again (configurations 737: three, 937: one; `mapsp` puts them on my 154 and 45), 12 min.
- **without Part A or the canonical form, non-extendable side**: all 11,004 Kurz covers in Kurz's coordinates
  (`thirteen_check covers`, `searchmem`): 11,004 / 308,112 / 8,627,136 / 115,373,106 / 207,621,076 / 42,791,038 /
  177,828; no near-miss node (none of the 31 non-extendable configurations has one), 53 min.

## 5. Controls

- `selftest` (190 configurations, 0 failures): 28 / 13,888 nondegenerate forms, |GL(4,2)| = 20160, |V5824| = 5824;
  |D| = 60, 200 lines, pairwise complementarity of the clique; each member's 28 forms distinct with radical exactly
  that member; for every 16th configuration all 200 · 13,888 line forms distinct, rank 6, radical exactly the line
  (33.4M forms checked in all); brute force over all 2^28 forms: exactly 28 with radical R_5 of the last
  configuration; canonical-form invariance (section 2).
- First levels recounted independently (`thirteen_check_levels.py`: forms B = π(x)ᵀ d π(y) in a basis adapted to the
  radical, d over the 28 Pfaffian-nonzero 4×4 matrices, radicals by brute force, its own legality): configuration 0
  (extendable) 1 / 28 / 784 / 10,573 and 170 (non-extendable) 1 / 28 / 784 / 10,288, identical to the C++.
- **Known functions through the same code path** (`ctlbatch`), two ways for every function with k ≥ 2 rank-4
  components. Extendable: its radicals → X, Y (and graph(I)), a random spread through them, four other spread elements
  as holes, the 13-set canonicalised and matched to a representative in `configs.txt` (this tests section 2: every
  such configuration must land on a listed extendable one), the function carried there by the canonical frame.
  Non-extendable: its radicals embedded into the members of a non-extendable representative (cycled through all 31;
  for k = 4 via a member 4-tuple of the right conjugacy class). Then the main search from the empty span, relaxed only
  by: a rank-6 radical outside D is legal if it is the known function's own, and those forms are offered at their
  points. Pass = the known space is among those found; every space found is re-verified by LUT/derivatives/Walsh.
  - BL22 (`apn_8bit.txt`, all 480 with k ≥ 2, 3 choices per path): extendable **1,440 of 1,440**, non-extendable
    **1,440 of 1,440** (924 / 438 / 78 with k = 2 / 3 / 4 on each path), 0 other spaces, 0 setup failures;
    representatives hit: 124 of 159 extendable, 31 of 31 non-extendable.
  - BLP22 (`blp22_8bit.txt`, 37 with k = 2, 6 choices per path): 222 of 222 and 222 of 222 (reps hit 73 / 31).
    Yu–Perrin QAM (4 with k = 2, 6 choices): 24 of 24 and 24 of 24. Yu–Wang–Li + Edel–Pott (12 with k = 2, 6
    choices): 72 of 72 and 72 of 72; its line 3 isn't an APN space as read (the list's known odd line) and was
    skipped (the 12 "FAIL" lines in `ctl-yuwangli.txt`). 0 other spaces anywhere.
  - The 3.8M (`new_apns.txt`, all 7,199 with k ≥ 2: 7,136 with k = 2 and the 63 with k = 3, one choice per path):
    **7,199 of 7,199** and **7,199 of 7,199**, 0 other spaces.
  - Representatives hit by the batches: 141 of 159 extendable, 31 of 31 non-extendable. The other 18 (among them 154,
    a near-miss configuration) were then targeted directly (`ctlreps`: the first 40 BL22 functions with k ≥ 3
    embedded into the members of each, the non-extendable path's embedding, plus 45 again): **40 of 40 for each of
    the 19**, 0 other spaces. So every one of the 190 configurations has re-found a known function through the main
    search.
  - Across all controls: no configuration failed to land on a listed representative, no known space was missed, no
    other space was ever found.
  - Bug found and fixed on the way (build 3): in single-control verbose mode a worker thread could take the search
    task and record the solution in the global map instead of the control's, so the known space showed as "not
    found". Batch mode wasn't affected (no nested parallelism there). The same pattern is in `fourteen_check.cpp`'s
    verbose `ctl` mode; its batch controls in `fourteen-check.md` aren't affected.

## 6. Verdict

**No quadratic APN function F_2^8 → F_2^8 has amplitude distribution [0^222, 2^20, 4^13]** (i = 4, 13 rank-4
components). It rests on:
- the hole structure, by computation and not by citing Kurz (section 1): every 10-clique checked; the extendable count
  is exactly C(14,10) · 1,128,864 with a unique extension each, the non-extendable count is exactly what Kurz's three
  sets and their 552 / 10,440 / 12 covers give, and every non-extendable configuration lands on a Kurz cover orbit;
- the configuration count (section 2): 159 + 31 orbits whose sizes sum to Part A's weighted counts, and two cross-runs
  that don't use the classification at all (fourteen-check's spread reps × 14; all 11,004 Kurz covers);
- complete point-branching run in five shapes (ascending, descending, random; spread file × 14; unreduced covers),
  all 0 solutions, nothing at dimension 7;
- positive controls through the same code path on both sides (extendable and non-extendable), on every list
  (18,674 in all, every one re-found, nothing else found), reaching every one of the 190 configurations.
Like i = 3 the failure is early and uniform, with one difference: here, four 6-dim spans (in two configurations,
not reduced further by their automorphisms) do hold all 13 rank-4 forms, and they die at once when any rank-6 form is tried. Confidence: high.

## 7. Comparison with Rowen's search (read after the above)

README "Thirteen" and `runs/i4/`:
- Non-extendable: Rowen takes the three hole sets from Kurz (citing completeness), finds the same **552 / 10,440 / 12**
  covers, and searches all 11,004 unreduced with `cover.cpp`: **0 solutions**, 374,898,296 nodes, 59 min (job
  i4-max). This check now supplies the step he cited: the three printed sets are all there is (section 1).
- Extendable: holes Y, X, graph(I), graph(T), T one per conjugacy class (5), rank-4 radicals the other 13 of each
  spread set through {0, I, T}: 25,850 covers, weighted check 15,804,096 = 14 · 1,128,864 (the same number as my
  source (a), |Ω| = 14 · 1,128,864). That's complete for the extendable case (any four disjoint solids go to
  X, Y, Z, graph(T), T up to conjugation). Job i4-ext (finished 2026-10-24 ~23:45): **0 solutions**, 723,800 tasks,
  768,769,492 nodes, 6,995 s.
- Same answer everywhere. On the non-extendable side the two searches can be compared node for node: my unreduced
  run over the same 11,004 covers visits 374,909,300 nodes, exactly his 374,898,296 plus the 11,004 roots (his
  tasks, 308,112 = cover × first form, are my dimension-1 nodes). Two independently written searches, same tree.
- His runs don't log the "all rank-4 forms present" event, so the four near-miss spans (extendable side) can't be
  matched against his output; `cover.cpp` passes through them the same way (rank-6 branching at the first hole).
- Cost: his searches cover 11,004 + 25,850 unreduced configurations (× 28 first forms); the GL(8,2) reduction here
  gives 190, searched in 2 min. The 73 i = 3 reps re-run through `cover.cpp` (`runs/i4/i3-reps.out`: 0, 19 s) is a
  check of `cover.cpp`, as he notes.

## 8. For i = 5 and below

What carries over: nothing here needs Kurz. Part A (k-cliques, minimal-class trick, exact weights), canonical forms
over frames of the partial spread itself, and the identity Σ |frames| / |Aut P| = Part A's weighted count found the
31 non-extendable configurations on their own, from Part A's cliques; Kurz's sets only served as a cross-check. For
12 solids the frames are 1320 ordered triples × 20160.
Sizing, i = 5 (`thirteen_check holesk 9`, seconds): 5,269,153 reduced 9-cliques. By how far the 12-partial spread
extends (largest clique in the common neighbourhood), weighted: maximal 49,896,000 (100,783 reduced); by exactly one
solid, to a maximal 13, 502,840,800; to a spread 2,258,825,184. That last is **less than C(14,9) · 1,128,864 =
2,259,985,728**, so at i = 5 some 12-partial spreads extend to a spread in more than one way. The extension stops
being unique, and the count check there has to count extensions (5-cliques in Q), as `holes` does with 4-cliques
here. Likewise for the maximal-13 extensions: if each were unique the weighted count would be 10 · 62,558,496.
Expect roughly 10^3 configurations for i = 5 (extendable ones from (spread, 5 holes), or straight from Part A's
canonical forms, ~2 h CPU at this rate). The forms search will cost more: fewer rank-4 forms means more in-D
branching at 13,888 forms a line, and near-miss spans (here 4) will become common.
