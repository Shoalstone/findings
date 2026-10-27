# Independent check, i = 6: is there a quadratic APN F on F_2^8 with distribution [0^214, 2^30, 4^11]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `eleven_check.cpp`,
> `eleven_check_levels.py` → `src/`; `runs/check11/parta-c*.txt` → `data/parta11-c*.txt`; `runs/check11/orbits.txt` →
> `data/orbits11.txt`; `runs/check11/configs.txt` → `data/configs11.txt`; `runs/check11/maxc-o798.txt` →
> `data/maxc11-o798.txt`; Part A's `maxc-c*.txt` (36 MB) and the sample `maxc-s50.txt` aren't included
> (`eleven_check parta` writes the first, the second is every 50th line of them, `README.md`); the other outputs in
> `runs/check11/` → `logs/check11-<name>-output.txt` (`<name>.out`), `logs/check11-ctl-*.txt`,
> `logs/check11-levels-output.txt` (`levels.txt`), and `logs/check11-{s12,s13,ssp,maxs50}-output.summary.txt` (headers
> and summaries; their per-configuration lines are left out); the share runner's job logs (`check11-asc` and so on) →
> `logs/2026…-check11-*.log`, whose commands list the 64 representatives of `ctlreps` ("the job record"); `.prog` files
> are resume state, not included, apart from their `P1` lines (the near-miss spans) for the three main runs →
> `logs/check11-{asc,desc,rand7}-p1.txt`; `runs/check12/configs.txt` → `data/configs12.txt`; `runs/check14/spreads.txt`
> → `data/spreads.txt`; `runs/check13/configs.txt` → `data/configs.txt`; `runs/papers/2112.11763.txt` is the text of
> Kurz's paper (arXiv 2112.11763, not included); `runs/bin/eleven_check.exe` was rebuilt from the one source, build it
> from `src/` (`README.md`). Rowen's side: `pspread.cpp`, `cover.cpp` → `src/` (`cover3` is a build of `cover.cpp`);
> `runs/i6/p11.txt` → `data/p11.txt`; `runs/i6/p11.out` → `logs/pspread11-output.txt`; `runs/i6/cover11.out` →
> `logs/eleven-output.txt`; README "Eleven: i = 6" is Rowen's working notes (not published), summarised in `README.md`
> here and in the write-up; `pspread3.cpp`, whose header section 7 mentions, is a later version of `pspread.cpp` with
> checkpoints, not used for eleven and not included. The known-function lists aren't included (`README.md`).

Subagent (Claude Opus 5.5), 2026-10-26, for Rowen ("search B" for i = 6). Written without reading `pspread*.cpp`,
`cover*.cpp`, `cover_prep.py`, README's i = 6 / "Eleven" parts, `runs/i6/`, `runs/i7/` or the `pspread*` logs until
the verdict below is in. Files: `eleven_check.cpp` (all parts, controls), `eleven_check_levels.py` (Python recount of
the first levels), outputs in `runs/check11/`. Built on `twelve_check.cpp`: utilities, V5824, clique enumeration,
point-set isomorphism, exact covers, canonical forms over frames, the forms search and the controls carried over,
generalised to 11 members; Part A reorganised (below).

**Status (2026-10-27 00:05): final verdict, none exists** (0 solutions everywhere, every control clean, every
configuration reached by a control). Comparison with search A in section 7: same 817, same answer, same near-misses.

## 1. The reduction

Setup as in `sixteen-check.md` §1. i = 6: 214 bent, 30 rank 6, 11 rank 4; 11·15 + 30·3 = 255. The 11 rank-4
radicals are pairwise complementary (a partial spread of 11 solids), D = the 90 points they miss, and the 30 rank-6
radicals (lines) partition D. Every point of a radical R_j can only be covered by a rank-4 form with radical exactly
R_j (rank-6 radicals lie in D), so all 11 rank-4 forms are forced by coverage, and the search needs nothing beyond
legality and covering every point.

A hyperplane contains at most one R_j, so it meets D in 127 - 77 - 8t = 50 or 42 points: D is a projective
[90, 8, {40, 48}] two-weight set, Cay(F_2^8, D) is strongly regular with eigenvalues 90, 10, -6, so λ = 90 + 10 - 6
- 60 = 34 and every hole lies on 17 lines inside D: 510 lines in every candidate (the search code counts 510 for each
configuration it builds).

**Normalisation.** Three of the 11 solids go to Y, X, Z = graph(I); the other 8 are graph(N), N in V5824, pairwise
N + N' invertible: an 8-clique. GL(4,2) conjugation fixes X, Y, Z and acts on the cliques.

## 2. Classification plan (how it differs from i = 5)

- **Part A (`eleven_check parta <c>`, one job per minimal class c, so an eviction loses one class).** Every 8-clique
  up to conjugation by the minimal-class trick (weights 1/m exact as multiples of 2520), each canonicalised
  (`canonP`, frames = 990 ordered member triples × diag(h, h)). Per clique only cheap things: |Q| (Q = the common
  neighbourhood = every solid disjoint from all 11) and the key. Maximal ones (Q empty) are written out (`maxc-c*.txt`).
- **`merge`**: per-orbit identity (weighted cliques on an orbit = 990 · 20160 / |Aut P| exactly) and, once per orbit
  rather than per clique, the extension invariants: the largest clique in Q (must be 0, 1, 2 or 6: an 11 extends to
  nothing, a maximal 12, a maximal 13 or a spread), the number of spreads through P (6-cliques in Q), and for every
  q in Q which of twelve-check's 376 orbits P + q lands on. Identities:
  - Σ_K |Q(K)| over all 8-cliques = 9 × twelve-check's weighted 9-cliques (2,811,561,984), from Part A alone;
  - Σ_P (990 · 20160 / |Aut P|) · #spreads through P = C(14,8) · 1,128,864;
  - **per 12-orbit O** (all 376): Σ_P (990 · 20160 / |Aut P|) · #{q : P + q ≅ O} = 9 · 1320 · 20160 / |Aut O|, which
    ties this classification to twelve-check's orbit by orbit;
  - the type from the clique size agrees with the types of the 12-orbits reached.
- **`configs`**: hole sets up to GL(8,2) against Kurz's Table 9.3 for r = 11 (**303 point sets, 58 maximal**,
  `runs/papers/2112.11763.txt`), exact covers of every hole set's complement by 11 solids (Σ |Aut H| / |Aut P| =
  covers, all orbits reached), and three outside families: twelve-check's 376 × 12 deletions (every 11 that extends),
  fourteen-check's 73 spread reps × C(14,3) (every spread-extendable 11), thirteen-check's 31 maximal 13s × C(13,2).

## 3. Method: point-branching on the forms, radicals fixed

As `twelve-check.md` §3 with 11 members and |D| = 90 (17 lines of D through each hole, 13,888 forms each).

## 4. Classification results (merge, configs)

**Part A** (all five classes finished, `runs/check11/parta-c*.txt`, each ending in `END`): reduced 8-cliques by
minimal class 0..4: 5,917,664 / 3,003,006 / 736,536 / 364,396 / 361,740 (10,383,342 in all; 8,611 / 4,950 / 815 /
832 / 195 s), of them maximal (Q empty) 599,168 / 86,195 / 17,880 / 7,125 / 0 = 710,368 (`maxc-c*.txt`).

**`merge`** (`runs/check11/merge.out`, 6 s): **all identities hold.**
- weighted 8-cliques 5,902,275,456 (an exact multiple of 2520 before division);
- Σ_K |Q(K)| = 25,304,057,856 = 9 × 2,811,561,984 (twelve-check's weighted 9-cliques, recomputed from its 376
  configurations' |Aut| here);
- **817 orbits**; no clique-level inconsistency (|Aut| or |Q| differing within an orbit); per-orbit identity
  (weighted cliques = 990 · 20160 / |Aut P|) holds for all 817, summing to the weighted total;
- largest clique in Q: 0 for 67 orbits, 1 for 43, 2 for 254, 6 for 453 (nothing else), so **spr 453, ext13 254,
  ext12 43, max 67** (weighted 3,377,565,576 / 1,760,081,400 / 342,619,200 / 422,009,280); type from the clique size
  agrees with the types of the 12-orbits reached for every orbit; every P + q lands on one of twelve-check's 376;
- Σ_P (990 · 20160 / |Aut P|) · #spreads through P = 3,389,978,592 = C(14,8) · 1,128,864;
- per 12-orbit link identity: holds for all 376 of twelve-check's orbits.

**`configs`** (`runs/check11/configs.out`, 19 s; writes `configs.txt`, 817 configurations in type order: spr 0..452,
ext13 453..706, ext12 707..749, max 750..816):
- canonical form recomputed from a random GL(8,2) image of every orbit (members shuffled): 0 mismatches (key and |Aut|);
- **hole sets up to GL(8,2): 303, of them 58 containing no solid: exactly Kurz's Table 9.3 for r = 11** (checked in
  `runs/papers/2112.11763.txt`: "# point sets ... 303", "# maximal point sets ... 58"); 302 distinct hole invariants,
  one pair separated only by the isomorphism test; type "max" ⇔ hole set contains no solid, for every orbit;
- exact covers of every hole set's complement by 11 solids: 2,134 covers in all; for each of the 303, every cover
  canonicalises onto a listed orbit with that hole set, all its orbits are reached, and Σ |Aut H| / |Aut P| = number
  of covers: 0 mismatches. With Kurz's count this makes the 817 complete independently of Part A's enumeration;
- outside sources, every configuration on a listed orbit: twelve-check's 376 × 12 (4,512) reach all 453 spr, 254
  ext13, 43 ext12; fourteen-check's spread reps × 364 (26,572) reach all 453 spr; thirteen-check's maximal 13s × 78
  (2,418) reach all 254 ext13 (and 82 spr).

## 5. Search results

Main search over the 817 configurations, three point orders (jobs 20261026-210508-check11-asc, -desc,
20261026-210509-check11-rand7, ~1,200 s each in a shared room), **0 solutions in each**, no bent completion ever ran:

| order | nodes by dimension 0..8 | near-misses (dim 6 / 7) | tries in D |
|---|---|---|---|
| asc | 817 / 22,876 / 640,528 / 8,993,779 / 19,647,718 / 1,938,778 / 24,532 / 38 / 0 | 22 / 37 | 6,944,000 |
| desc | 817 / 22,876 / 640,528 / 8,935,159 / 19,095,956 / 1,335,017 / 18,128 / 40 / 0 | 22 / 36 | 7,318,976 |
| rand 7 | 817 / 22,876 / 640,528 / 8,944,444 / 19,234,261 / 1,513,135 / 18,258 / 38 / 0 | 22 / 37 | 6,860,672 |

By type (asc): spr 453 configurations, 37 nodes at dim 7, near-misses 14 / 36; ext13 254, 1 at dim 7, near-misses
8 / 1; ext12 (43) and max (67) never reach dim 7 and never hold all 11 rank-4 forms (deepest: dim 6, 654 and 248
nodes). Near-misses are spans containing all 11 rank-4 forms, which then fail in D (58-59 distinct such spans per
order, logged as P1 lines in the .prog files). `asc.out`, `desc.out`, `rand7.out`, `.prog`.

Cross-runs that don't use this classification (finished, 0 solutions):
- `search12` (twelve-check's 376 × 12 = 4,512: every 11 that extends), job 20261026-175950-check11-s12, 5,220 s of
  running: nodes by dimension 0..8 4,512 / 126,336 / 3,537,408 / 49,655,433 / 109,538,964 / 10,777,400 / 165,744 / 111
  / 0; **near-misses** (every off-D point covered, i.e. all 11 rank-4 forms present): 97 at dimension 6, 103 at
  dimension 7; 25,637,248 candidate forms tried at points of D (17.4M at dim 6, 8.2M at dim 7); no bent completion;
  **0 solutions** (`s12.out`).
- `search13` (thirteen-check's 31 maximal 13s × 78 = 2,418), job -s13, 3,010 s: 2,418 / 67,704 / 1,895,712 /
  26,629,650 / 60,146,644 / 8,913,624 / 96,436 / 72 / 0; near-misses 22 at dim 6, 61 at dim 7; 10.1M tries in D;
  **0 solutions** (`s13.out`).

Unlike i = 5, the rank-6 branching does run here: spans holding all 11 rank-4 forms appear at dimensions 6 and 7, and
die there in D (nothing reaches 8). So the in-D control (section 6) matters for this verdict.

First levels recounted independently (`eleven_check_levels.py`, asc order, through dimension 4, `levels.txt`):
config 0 (spr) 1 / 28 / 784 / 11,152 / 24,772; 453 (ext13) 1 / 28 / 784 / 11,206 / 23,367; 707 (ext12) 1 / 28 /
784 / 10,939 / 24,007; 750 (max) 1 / 28 / 784 / 11,074 / 22,957: **identical to the C++ asc run** node for node on
all four (`asc.out`).

`selftest` (`selftest.out`, 817 configurations, **0 failures**): 28 / 13,888 nondegenerate forms, |GL(4,2)| = 20160,
|V5824| = 5824; |D| = 90, 11 members, 510 lines in D for all 817; pairwise complementarity; each member's 28 forms
distinct with radical exactly it; every 16th configuration's 510 × 13,888 line forms distinct, rank 6, radical
exactly the line (368,561,396 generated forms checked); brute force over 2^28: exactly 28 forms with radical R_5 of the
last configuration; canonical form invariant under a random GL(8,2) image, with |Aut| and frame.

- `searchsp` (fourteen-check's spread reps × 364 = 26,572: every 11 that extends to a spread), job
  20261026-175950-check11-ssp, resumed once from `ssp.prog` (26,572 configurations in all, each once): nodes by
  dimension 0..8 26,572 / 744,016 / 20,832,448 / 291,080,164 / 635,317,280 / 54,887,313 / 860,221 / 495 / 0;
  near-misses 263 at dim 6, 450 at dim 7; 82.9M tries in D (43.4M / 39.6M); no bent completion; **0 solutions**
  (`ssp.out`).
- `searchmax` on the 1-in-50 sample of Part A's maximal cliques (14,208, below), job
  20261026-210924-check11-maxs50: 14,208 / 397,824 / 11,139,072 / 154,792,269 / 321,538,731 / 14,093,008 / 17,710 /
  0 / 0; no near-miss (no span ever holds all 11 rank-4 forms, as for the max type in the main runs), no tries in D;
  **0 solutions** (`maxs50.out`).

**searchmax: sampled, not the full list.** The plan said all 710,368 maximal 8-cliques (×7 twelve-check's 100,783).
Measured at 0.35-1.7 s per configuration in the shared room that is days of box time, for a cross-check of a side
that is already covered three ways (the 67 max orbits searched in three orders; exact covers of all 58 maximal hole
sets accounted orbit by orbit; Kurz's count). So `maxc-s50.txt` = every 50th line of `maxc-c0..c3.txt` concatenated
(14,208 configurations: 11,984 / 1,724 / 357 / 143 by class; the smallest max orbit has 85 reduced cliques, so every
orbit very likely appears). Rowen can run the full list later if wanted: `searchmax <concatenated maxc> <prog> asc`,
resumable.

Orbit coverage of the sample (new mode `orbitcov`, canonicalises each listed clique and looks it up among the 817,
`maxs50-orbits.out`): all 14,208 on a listed max orbit, **66 of the 67 max orbits hit**; the one missed is orbit /
configuration 798, one of the two smallest (85 reduced cliques, expected 1.7 in the sample; its twin 799 got some).
Rather than the full list (~50 × the sample's 2.2 h, around 100 h of box time, to re-search orbits already covered),
a targeted run (job 20261026-235204-check11-max798): `orbitcov` over all 710,368 maximal cliques
(`maxall-orbits.out`): every one on a listed max orbit, all 67 hit, and the per-orbit tallies equal Part A's reduced
clique counts (`n` in `orbits.txt`) for every one of the 67; the 85 on orbit 798 went to `maxc-o798.txt`, and
`searchmax` on those (`max798.out`): 85 / 2,380 / 66,640 / 934,228 / 2,128,630 / 407,037 / 5,956 / 0 / 0, no
near-miss, **0 solutions**. So every max orbit has had raw Part A cliques (not canonical forms) searched too.

## 6. Controls

- **The in-D control** (`ctlbatch ... nond` on BL22, all 480 functions with k ≥ 2, 3 choices per path, job
  20261026-210845-check11-ctlnond, `ctl-ind-bl22.out` / `.txt`): **1,438 of 1,438** on the extendable path and
  **1,440 of 1,440** on the other (2 trials skipped: start span already the whole space), 0 other spaces, 0 setup
  failures; every control started at dimension 7 and found g among the 13,888 forms on its line: 98,965,888 and
  99,785,280 in-D candidates tried. Representatives reached: 318 on the extendable path, 364 of 364 (all ext13 +
  ext12 + max) on the other. This is the code path the near-misses exercise.
- Ordinary controls from the empty span (jobs 20261026-210844-check11-ctl and -210845-check11-ctlnew; `ctl-*.out` /
  `.txt`), same set-up as `twelve-check.md` §6 (extendable path: carried onto an spr configuration; other path:
  embedded into an ext13 / ext12 / max one):
  - BL22 (`apn_8bit.txt`, all 480 with k ≥ 2, 3 choices per path): **1,440 of 1,440** on each path (924 / 438 / 78
    with k = 2 / 3 / 4), 0 other spaces, 0 setup failures; representatives 319 (ext) and 364 (nonext).
  - BLP22 (37 with k = 2, 6 choices): **222 of 222** on each path; Yu–Perrin QAM (4, 6 choices): **24 of 24** on each.
  - Yu–Wang–Li + Edel–Pott (13 lines with k ≥ 2, 6 choices): **72 of 72** re-found on each path, 0 not found, 0 other
    spaces; the 6 + 6 setup failures are all line 3, "forms not an APN space (line not quadratic APN)" (the FAIL lines
    in `ctl-yuwangli.txt`), the same odd line twelve-check skipped. Checked here independently in Python: line 3 as
    read has differential uniformity 2 (APN) but a nonzero third derivative, so it isn't quadratic and has no form
    space to carry; out of scope for the control, not a miss. (The job exits 3 because of these 12.)
  - The 3.8M (`new_apns.txt`, the 7,199 with k ≥ 2, one choice per path): **7,199 of 7,199 re-found on each path**
    (7,136 / 63 with k = 2 / 3), 0 other spaces, 0 setup failures; representatives 383 (ext) and 364 (nonext).
  - Representatives, both paths and all five files together: **753 of 817**: every ext13, ext12 and max configuration
    (364) and 389 of the 453 spr. The other 64 (spr 93, 112, ..., 452, listed in the job record) get the targeted
    run, `ctlreps` with the first 40 BL22 functions with k ≥ 3 embedded into each one's members (`ctlreps.out`, job
    20261026-235204-check11-ctlreps, exit 0): **40 of 40 re-found for every one of the 64**, 0 other spaces. So
    every one of the 817 configurations has had known functions carried onto it and re-found.

## Verdict (final, 2026-10-27 00:05, subagent Claude Opus 5.5 for Rowen, continuing the session that wrote 1-6)

**There is no quadratic APN function on F_2^8 with Walsh distribution [0^214, 2^30, 4^11].** 0 solutions in the
complete search over all 817 classes of 11-solid partial spreads (three point orders), and 0 in four cross-runs that
don't rest on this classification: twelve-check's 376 × 12, thirteen-check's maximal 13s × 78, fourteen-check's
spread reps × 364, and Part A's raw maximal cliques (a 1-in-50 sample plus every clique of the one max orbit the
sample missed, so all 67 max orbits). The classification is backed by every identity in `merge` (including the
orbit-by-orbit link to twelve-check's 376), by Kurz's 303 / 58 with all exact covers accounted for, and by three
outside sources. The search is backed by controls on every one of the 817 configurations (BL22, BLP22, Yu–Perrin,
Yu–Wang–Li, the 3.8M, plus the targeted ctlreps: every known quadratic APN function tried was re-found, 0 misses, the
only setup failures being Yu–Wang–Li's non-quadratic line 3), and the rank-6 branching in D, which does run here
(59 near-misses per order), by the 2,878 / 2,878 in-D control. Signed before reading anything of search A.

## 7. Comparison with search A (read after the verdict above was signed and committed)

README "Eleven: i = 6" and `runs/i6/` (`p11.txt`, `p11.out`, `cover11.out`), `pspread3.cpp`'s header:
- **Classification.** Search A's `pspread 11` finds the same **817 configurations, 67 / 43 / 254 / 453** by extension
  (its w = 0 / 1 / 2 / 6 = max / ext12 / ext13 / spr here), the same reduced counts (710,368 / 619,505 / 3,324,960 /
  5,728,509, summing to Part A's 10,383,342; its 710,368 maximal = the `maxc` files here) and the same weighted counts
  (422,009,280 / 342,619,200 / 1,760,081,400 / 3,377,565,576), orbit sums matching. Mapped directly (`eleven_check
  mapa`, its `p11.txt` read as 11 solids by their bases and canonicalised by this code, `mapa.out`): **a bijection**
  onto the 817 here, with type and |Aut| agreeing on every one. As at i = 5, the two classifications share the
  method's shape (minimal-class reduction, canonical form over frames with the least class vector), so the matching
  reduced counts say both enumerations are complete for the same reduction rather than that two different reductions
  agree; the checks on top are independent (here: per-orbit identities, the link to twelve-check's 376 orbit by orbit,
  Kurz's 303 / 58 and exact covers of every hole set, three outside families).
- **Search.** A's `cover3` on all 817: **0 solutions**, 30,505,655 nodes, by dimension 1..8 22,876 / 640,528 /
  8,993,626 / 19,289,387 / 1,533,912 / 25,288 / 38 / 0. Here (asc) 22,876 / 640,528 / 8,993,779 / 19,647,718 /
  1,938,778 / 24,532 / 38 / 0: equal through dimension 2, then apart, as at i = 5 (different representatives and
  branching points).
- **Near-misses agree configuration by configuration.** A logged 58 (22 at dim 6, 36 at dim 7) on 11 configurations;
  mapped onto ours they are exactly the 11 where the runs here have theirs (248, 299, 314, 316, 429, 435, 452, 541,
  558, 609, 649), with the same count on each (6, 6, 6, 27, 1, 1, 3, 3, 1, 1, 3). desc matches A's 58 exactly; asc
  and rand 7 have one more, a second distinct span at 558 (59). That isn't a conflict: a point-branching tree only
  visits partial spans whose forms cover the first uncovered points in its own order, so which near-misses appear
  depends on the order, while a solution would be reached in any.
- **What this check adds beyond A:** controls (A had none at i = 6; every one of the 817 has had known functions
  carried onto it and re-found here, and the in-D path the near-misses exercise has its own 2,878 / 2,878), and four
  cross-runs that don't go through the 817.
- Same answer: no function, by two searches and two classifications that agree configuration for configuration,
  near-misses included.

## State (2026-10-27 00:20)

Done. Verdict final (none exists), comparison with search A written (section 7). `eleven_check.exe` was rebuilt for
two comparison-only modes (`orbitcov`, `mapa`); the search code is unchanged. Nothing running.
