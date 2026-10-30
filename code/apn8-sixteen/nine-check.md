# Independent check, i = 8: is there a quadratic APN F on F_2^8 with distribution [0^206, 2^40, 4^9]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `nine_check.cpp`, `nine_check_levels.py`,
> `nine_check_slices.py`, `nine_check_tally.py`, `apn/nine_check_compare.py`, `apn/ten_check_tally.py` → `src/` (this
> folder's `nine_check_compare.py` reads the logs here); `runs/check9/orbits.txt` → `data/orbits9.txt`;
> `runs/check9/configs.txt` → `data/configs9.txt`; `runs/check9/maxc-s1000.txt` → `data/maxc9-s1000.txt`;
> `runs/check9/parta-c<c>-s<s>.txt` → `data/parta9-c<c>-s<s>.txt`; Part A's `maxc-c<c>-s<s>.txt` (~645 MB) aren't
> included (`nine_check parta` writes them in an order that depends on thread timing, which is why the sample is
> included); the other outputs in `runs/check9/` → `logs/check9-<name>-output.txt` (`<name>.out`), `logs/check9-ctl-*.txt`,
> `logs/check9-levels-output.txt` (`levels.txt`), `logs/check9-reps-unreached.txt`, `logs/check9-maxs1000-near.txt`,
> and `logs/check9-{asc,desc,rand7,s10,s13,maxs1000,ssp-r0,…,ssp-r7}-output.summary.txt` (headers and summaries; their
> per-configuration lines are left out); the `.prog` files of the three main runs → `logs/check9-{asc,desc,rand7}-prog.txt`
> (their `R` lines are the per-configuration record, their `P1` lines the near-miss spans); the cross-runs' `.prog`
> files aren't included; the share runner's job logs (`check9-parta` and so on) → `logs/2026…-check9-*.log`, whose
> commands are `nine_check_slices.py`'s (which reads and writes `runs/check9/` and the folders below);
> `nine_check_tally.py prog` reads `.prog` files, `ten_check_tally.py` a folder holding `configs.txt` and `ctl-*.txt`
> (`reps`) or `orbits.txt` and `maxall-orbits.out` (`maxall`) under those names; ten-check's, eleven-check's,
> twelve-check's, thirteen-check's and fourteen-check's lists (`runs/check10/configs.txt` and so on) are
> `data/configs10.txt`, `data/configs11.txt`, `data/configs12.txt`, `data/configs.txt` and `data/spreads.txt`;
> `runs/papers/2112.11763.txt` is the text of Kurz's paper (arXiv 2112.11763, not included); `runs/bin/nine_check.exe`
> is a build of the one source, build it from `src/` (`README.md`). Rowen's side: `pspread*.cpp`, `cover*.cpp` →
> `src/pspread3.cpp`, `src/cover.cpp` (`cover3` is a build of `cover.cpp`); `runs/i8/p9.txt` → `data/p9.txt`;
> `runs/i8/p9.out` → `logs/pspread9-output.txt`; `runs/i8/cover9.out` → `logs/nine-output.txt`; README "Nine" is
> Rowen's working notes (not published), summarised in `README.md` here and in the write-up; `notes/` is Rowen's own
> notes (not published). The known-function lists aren't included (`README.md`).

Subagent (Claude Opus 5.5), 2026-10-28, for Rowen ("search B" for i = 8). Written without reading `pspread*.cpp`,
`cover*.cpp`, `cover_prep.py`, README's "Nine" section (or anything about i = 8 there), `runs/i8/`, `notes/` or any
near-miss/result line about nine until the verdict below is in. Files: `nine_check.cpp` (all parts, controls),
`nine_check_levels.py` (Python recount of the first levels), `nine_check_slices.py` (resumable driver for every
stage), `nine_check_tally.py` (search tallies from the .prog files), outputs in `runs/check9/`. Built on `ten_check.cpp`: everything carried over and generalised to 9 members.

**Status: verdict in** (no function; 2026-10-29 21:00); comparison with search A in section 7 (2026-10-29): the two agree configuration for configuration.

## 1. The reduction

Setup as in `sixteen-check.md` §1. i = 8: 206 bent, 40 rank 6, 9 rank 4; 9·15 + 40·3 = 255. The 9 rank-4 radicals
are pairwise complementary (a partial spread of 9 solids), D = the 120 points they miss, and the 40 rank-6 radicals
(lines) partition D. Every point of a radical R_j can only be covered by a rank-4 form with radical exactly R_j, so
all 9 rank-4 forms are forced by coverage, and the search needs nothing beyond legality and covering every point.

A hyperplane contains at most one R_j, so it meets D in 127 - 63 - 8t = 64 or 56 points: D is a projective
[120, 8, {56, 64}] two-weight set, Cay(F_2^8, D) is strongly regular with eigenvalues 120, 8, -8, so λ = μ =
120 - 64 = 56 and every hole lies on 28 lines inside D: 1,120 lines in every candidate (the search code counts them).

**Normalisation.** Three of the 9 solids go to Y, X, Z = graph(I); the other 6 are graph(N), N in V5824, pairwise
N + N' invertible: a 6-clique. GL(4,2) conjugation fixes X, Y, Z and acts on the cliques.

## 2. Classification plan

- **Part A (`nine_check parta <c> <out> <maxc> <slice> <nslices>`)**: every 6-clique up to conjugation by the
  minimal-class trick (weights 1/m exact as multiples of 2520), each canonicalised (`canonP`, frames = 504 ordered
  member triples × diag(h, h)); per clique |Q| and the key; maximal ones (Q empty) written out. Class 0 in 16 slices,
  class 1 in 8, classes 2 and 3 in 4, class 4 in 2 (slice = first vertices with index ≡ slice mod nslices), run by
  `nine_check_slices.py parta`, which skips finished slices, so an eviction loses at most one slice.
- **`merge`**: every slice of every class present exactly once; per-orbit identity (weighted cliques on an orbit =
  504 · 20160 / |Aut P|); once per orbit, the largest clique in Q (must be 0, 1, 2, 3, 4 or 8: a 9 extends to nothing,
  a maximal 10, 11, 12, 13, or a spread, since there are no maximal 14, 15, 16), the number of spreads through P
  (8-cliques in Q), and which of ten-check's 1,777 orbits P + q is for every q in Q. Identities:
  - Σ_K |Q(K)| over all 6-cliques = 7 × ten-check's weighted 7-cliques (15,035,796,864, recomputed here from its 1,777
    configurations' |Aut|), from Part A alone;
  - Σ_P (504 · 20160 / |Aut P|) · #spreads through P = C(14,6) · 1,128,864;
  - **per 10-orbit O** (all 1,777): Σ_P (504 · 20160 / |Aut P|) · #{q : P + q ≅ O} = 7 · 720 · 20160 / |Aut O|;
  - the type from the clique size agrees with the types of the 10-orbits reached.
- **`configs`**: hole sets up to GL(8,2) against Kurz's Table 9.3 for r = 9 (**5442 point sets, 1343 maximal**,
  `runs/papers/2112.11763.txt`), exact covers of every hole set's complement by 9 solids (Σ |Aut H| / |Aut P| =
  covers, all orbits reached), and five outside families: ten-check's 1,777 × 10 (every 9 that extends),
  fourteen-check's 73 spread reps × C(14,5) (every 9 inside a spread), thirteen-check's 31 maximal 13s × C(13,4),
  twelve-check's 376 × C(12,3), eleven-check's 817 × C(11,2).

## 3. Method: point-branching on the forms, radicals fixed

As `twelve-check.md` §3 with 9 members and |D| = 120 (28 lines of D through each hole, 13,888 forms each).

## 4. Classification results (merge, configs)

**Part A**: all 34 slices ended in `END` (`runs/check9/parta-c<c>-s<s>.txt`, job 20261028-191409-check9-parta):
65,323,403 reduced 6-cliques, of them 16,492,152 maximal (`maxc-c<c>-s<s>.txt`, ~645 MB in all); weighted 6-cliques
47,414,189,312.

**`merge`** (`merge.out`, job 20261028-230254-check9-merge, via `nine_check_slices.py merge`): **all identities hold.**
- Σ_K |Q(K)| = 105,250,578,048 = 7 × ten-check's weighted 7-cliques (15,035,796,864, recomputed from its 1,777 |Aut|).
- **5,663 orbits**, no clique-level inconsistency (|Aut| and |Q| agree on every clique of an orbit, across slices too);
  per-orbit identity (weighted cliques = 504 · 20160 / |Aut P|) holds for every one; Σ over orbits = the weighted count.
- Largest clique in Q: 0 / 1 / 2 / 3 / 4 / 8 for 1,348 / 1,916 / 874 / 227 / 623 / 675 orbits (nothing else occurs), so
  **spr 675, ext13 623, ext12 227, ext11 874, ext10 1,916, max 1,348**; weighted 3,228,877,176 / 3,937,565,520 /
  1,911,047,040 / 7,504,126,560 / 18,422,651,520 / 12,409,921,496.
- Type agrees with the 10-orbits reached, every P + q lies on one of ten-check's 1,777; spreads through P (8-cliques in
  Q) only on spr orbits.
- Σ_P (504 · 20160 / |Aut P|) · #spreads = 3,389,978,592 = C(14,6) · 1,128,864.
- Per-10-orbit link identity holds for **all 1,777**.

The orbit count saturates per slice (one class-0 slice alone held ~5,560 of the 5,663).

**`configs`** (`configs.out`, job 20261028-230339-check9-configs; `configs.txt` in type order: spr 0..674, ext13
675..1297, ext12 1298..1524, ext11 1525..2398, ext10 2399..4314, max 4315..5662): canonical form recomputed from a
random GL(8,2) image of every orbit, 0 mismatches; **hole sets up to GL(8,2) 5,442, of them 1,343 containing no solid:
Kurz's Table 9.3 for r = 9**; type "max" ⇔ hole set contains no solid, every orbit; exact covers of every hole set's
complement (7,925 covers): every cover on a listed orbit with that hole set, Σ |Aut H| / |Aut P| = covers for all
5,442, every one of the 5,663 orbits reached. Outside families, every configuration on a listed orbit:
- ten-check's 1,777 × 10 = 17,770: all 4,315 non-max orbits (spr 675, ext13 623, ext12 227, ext11 874, ext10 1,916);
- fourteen-check's 73 spread reps × 2,002 = 146,146: all 675 spr;
- thirteen-check's 31 maximal 13s × 715 = 22,165: all 623 ext13 (+ 425 spr);
- twelve-check's 376 × 220 = 82,720: all spr, ext13, ext12;
- eleven-check's 817 × 55 = 44,935: all spr, ext13, ext12, ext11.

So the max orbits rest on the covers (and on `orbitcov` over Part A's raw maximal cliques, section 5); the rest on
covers plus at least one outside family each.

One orbit stands out: orbit 5662 (max, config 5662), whose 9 solids cover exactly the 135 singular points of a
hyperbolic quadric Q+(7,2) (checked directly: the covered set is the zero set of a quadratic form; no other orbit's
is). Its hole set (hole set 0: the 120 nonsingular points) has |Aut H| = 348,364,800 = |O+(8,2)|, 270 solids in its
complement (the quadric's generators) and 1,920 covers (the spreads of Q+(7,2), both families), all on this one
orbit, so |Aut P| = 181,440. That one automorphism count took most of the 684 s of the |Aut H| step.

## 5. Search results

Main search over the 5,663 configurations, three point orders (job 20261028-235552-check9-search, via
`nine_check_slices.py search`; held once for other work on the box, requeued when the runner stopped, and resumed from
the .prog files), **0 solutions in
each**, no bent completion ever ran, 0 full-span records. Tallied from the R lines (`nine_check_tally.py prog`), each
configuration exactly once:

| order | nodes by dimension 0..8 | near-misses (dim 5 / 6 / 7) | tries in D (dim 6 / 7) | distinct p1 spans |
|---|---|---|---|---|
| asc | 5,663 / 158,564 / 4,439,792 / 68,939,870 / 220,381,954 / 28,376,993 / 178,277 / 351 / 0 | 18 / 3,922 / 296 | 989.9M / 37.5M | 4,236 |
| desc | 5,663 / 158,564 / 4,439,792 / 68,555,168 / 218,112,825 / 26,492,204 / 171,097 / 361 / 0 | 18 / 3,820 / 305 | 963.5M / 38.2M | (9, see below) |
| rand 7 | 5,663 / 158,564 / 4,439,792 / 68,557,415 / 218,135,142 / 26,392,740 / 162,550 / 386 / 0 | 18 / 3,768 / 321 | 927.3M / 39.5M | 4,107 |

By type (asc): spr (675) 263 nodes at dim 7, near-misses 18 / 3,790 / 231; ext13 (623) 86 at dim 7, near-misses
0 / 121 / 64; ext12 (227) deepest dim 6 (802 nodes), 3 near-misses at dim 6; ext11 (874) 2 nodes at dim 7, near-misses
6 / 1 at dim 6 / 7; ext10 (1,916) deepest dim 6 (2,203), 1 near-miss at dim 6; max (1,348) deepest dim 6 (2,062), 1
near-miss at dim 6 (config 4776). Near-misses are spans containing all 9 rank-4 forms (every off-D point covered),
which then fail in D.

New at i = 8: **near-misses at dimension 5**, 9 rank-4 forms spanning only a 5-space, 18 of them on 7 spr
configurations (129: 1, 153: 3, 236: 1, 348: 3, 427: 1, 671: 6, 674: 3), the same configurations and counts in all
three orders, as they must be (such a span is fixed by its 9 forms). The dim-6 counts behave as they should too: on
every configuration without a dim-5 near-miss, the dim-6 count is the same in all three orders (385 in all); the
order-dependent part of the totals sits on those 7 configurations, where dim-5 near-misses are extended by a rank-6
form in D at an order-dependent point. The dim-6 near-misses lie on the same 156 configurations in every order.

desc's distinct-p1 count is incomplete: P1 lines are written at the end of a run, and desc was requeued after 1,595
configurations (the R lines, which carry everything else, are complete). Cosmetic.

First levels recounted independently (`nine_check_levels.py`, through dimension 4, `levels.txt`): configurations 0,
675, 1298, 1525, 2399, 4315 (one of each type) and 5662 (the Q+(7,2) one) identical to the C++ (asc) node for node:
1/28/784/12,220/40,545; 1/28/784/12,295/40,929; 1/28/784/12,184/40,622; 1/28/784/12,193/38,864; 1/28/784/12,202/39,564;
1/28/784/12,412/39,324; 1/28/784/12,994/56,092.

Cross-runs that don't use this classification (each configuration once, 0 solutions):
- `searchmax` on `maxc-s1000.txt` (every 1000th of Part A's 16,492,152 maximal cliques, 16,493, raw cliques rather
  than canonical forms; job 20261028-230536-check9-maxs1000), `maxs1000.out`: 16,493 / 461,804 / 12,930,512 /
  199,339,583 / 628,642,822 / 69,498,624 / 15,776 / 0 / 0; 6 near-misses at dim 6, all on max orbit 4776 (`orbitcov`
  on those 6 cliques), the one max configuration with a near-miss in the main search; 1.5M tries in D;
  **0 solutions**. (`searchmax` labels these "max"; at nine the label is right.) Orbit coverage: `orbitcov` over all
  34 maxc files (`maxall-orbits.out`, `ten_check_tally.py maxall`): every one of the 16,492,152 maximal cliques on a
  listed max orbit, all 1,348 hit, per-orbit tallies = Part A's `n` for every one.
- `search10` on ten-check's 1,777 × 10 = 17,770 (every 9 that extends at all), `s10.out` / `s10.prog` (job
  20261029-044940-check9-cross2; killed and resumed once at ~06:20 when Rowen restarted the runner; the tally finds
  each configuration exactly once, no duplicates): 17,770 / 497,560 / 13,931,680 / 216,427,891 / 702,927,941 /
  106,444,614 / 1,700,247 / 2,793 / 0; near-misses 46 / 12,042 / 2,175 (dim 5 / 6 / 7); 3.42G tries in D;
  no bent completion; **0 solutions**.
- `search13` on thirteen-check's 31 maximal 13s × 715 = 22,165, `s13.out` / `s13.prog` (same job): 22,165 / 620,620 /
  17,377,360 / 270,680,302 / 902,008,389 / 165,304,542 / 3,887,639 / 7,102 / 0; near-misses 22 / 12,500 / 5,379;
  3.96G tries in D; no bent completion; **0 solutions**.
- `searchsp` on fourteen-check's 73 spread reps × 2,002 = 146,146 (every 9 inside a spread), job
  20261028-230727-check9-ssp, run in 8 resumable ranges (`ssp-r<k>.out` / `.prog`) through many holds, requeues and
  Rowen's runner restart, ~22 h wall: each configuration exactly once, no duplicates; 146,146 / 4,092,088 /
  114,578,464 / 1,778,825,956 / 5,844,031,218 / 985,797,672 / 25,756,497 / 42,522 / 0; near-misses 985 / 231,767 /
  32,409 (dim 5 / 6 / 7); 60.1G tries in D (55.9G / 3.8G at dim 6 / 7); 262,971 distinct p1 spans (a lower bound:
  P1 lines of a range interrupted mid-run are lost); no bent completion; **0 solutions**.

`selftest` (`selftest.out`, job 20261028-235617-check9-selftest, **0 failures**): ND4 28, ND6 13,888, |GL(4,2)| =
20160, |V5824| = 5824; 1,120 lines in D for all 5,663 configurations; pairwise complementarity; member forms distinct
with radical exactly the member; line forms distinct, rank 6, radical exactly the line (5,507,741,316 generated forms
checked); brute force over 2^28: exactly 28 forms with radical R_5 of the last configuration; canonical form and |Aut|
reproduced from a random GL(8,2) image of every configuration, with a frame carrying it back.

## 6. Controls

Job 20261029-030815-check9-controls (`nine_check_slices.py controls`, `ctl-*.out` / `.txt`); the job exits 3 only
because of Yu–Wang–Li's line 3 (below).
- **The in-D control** (`ctlbatch ... nond` on BL22, the 480 with k >= 2, 3 choices per path, `ctl-ind-bl22`):
  **1,440 of 1,440** (ext) and **1,436 of 1,436** (nonext) re-found, every one started at dimension 7 and found g
  among the 13,888 forms on its line (153.6M / 154.0M in-D candidates tried), 0 not found, 0 setup failures, 0 other
  spaces; 4 k = 4 controls skipped on the nonext path (the rep drawn had no member of the 4th radical's conjugacy
  class). This is the code path the near-misses exercise.
- Ordinary controls from the empty span, as `ten-check.md` §6 (extendable path: carried onto an spr configuration by
  its canonical frame; other path: embedded into a non-spr one):
  - BL22 (`apn_8bit.txt`, 480 with k >= 2, 3 choices): **1,440 of 1,440** ext (924 / 438 / 78 with k = 2 / 3 / 4) and
    **1,436 of 1,436** nonext (4 skipped as above), 0 other spaces, 0 setup failures; reps hit 443 / 1,412.
  - BLP22 (37 with k = 2, 6 choices): **222 of 222** on each path; Yu–Perrin QAM: **24 of 24** on each.
  - Yu–Wang–Li + Edel–Pott (6 choices): **72 of 72** re-found on each path, 0 not found; the 6 + 6 setup failures are
    line 3, "not quadratic APN" (the same line as in ten-check and eleven-check: APN but not quadratic; out of scope, not a miss).
  - The 3.8M (`new_apns.txt`, the 7,199 with k >= 2, one choice per path): **7,199 of 7,199 on each path** (7,136 /
    63 with k = 2 / 3), 0 other spaces, 0 setup failures; reps hit 575 / 3,837.
  - Representatives reached, all batches together (`ten_check_tally.py reps`, `reps-unreached.txt`): 4,812 of 5,663;
    the other 851 (spr 91, ext13 90, ext12 32, ext11 133, ext10 293, max 212) got the targeted run, `ctlreps` with
    BL22's first 40 functions with k >= 3 embedded into each one's members (`ctlreps.out`, job
    20261029-085400-check9-ctlreps, exit 0): **40 of 40 re-found for every one of the 851** (34,040 controls), 0 other
    spaces. So every one of the 5,663 configurations has had known functions carried onto it and re-found.

## Verdict (2026-10-29 21:00, written by a subagent, Claude Opus 5.5, for Rowen)

**There is no quadratic APN function on F_2^8 with Walsh distribution [0^206, 2^40, 4^9].** 0 solutions in the
complete search over all 5,663 classes of 9-solid partial spreads (three point orders), and 0 in four cross-runs that
don't rest on this classification: ten-check's 1,777 × 10 (every 9 that extends), thirteen-check's maximal 13s × 715,
fourteen-check's spread reps × 2,002 (all 146,146), and a 1-in-1000 sample of Part A's raw maximal cliques. The
classification is backed by every identity in `merge` (Σ |Q| against ten-check's weighted count, per-orbit counts,
spreads through each orbit, and the orbit-by-orbit link to all 1,777 of ten-check's), by Kurz's 5,442 / 1,343 with
all 7,925 exact covers accounted for and every orbit reached, by five outside families, and by `orbitcov` placing all
16,492,152 raw maximal cliques on the 1,348 max orbits with Part A's counts. The search is backed by `selftest`, the
independent first-levels recount on seven configurations, controls on every one of the 5,663 configurations (BL22,
BLP22, Yu–Perrin, Yu–Wang–Li, the 3.8M, and `ctlreps` on the 851 the batches missed: every known quadratic APN
function tried was re-found, 0 misses, the only setup failures being Yu–Wang–Li's non-quadratic line 3), and the
rank-6 branching in D, where all the near-misses die, by the 2,876 / 2,876 in-D control.

New here compared with i = 7: near-misses already at dimension 5 (9 rank-4 forms spanning a 5-space; 18 in the main
search, on 7 spr configurations, identical in all orders), and one configuration that is a spread of the hyperbolic
quadric Q+(7,2) (config 5662, |Aut| 181,440, hole set the 120 nonsingular points), which the search handles like any
other (no near-miss there).

On independence: the code, the classification and every number above come from this repo's `nine_check*` files and
`runs/check9/`, written and run by two subagents (Part A by the first, the rest by the second) that read none of
`pspread*`, `cover*`, `cover_prep.py`, README's "Nine" section, `runs/i8/` or `notes/`, and no number of search A's
for i = 8. One remark leaked to the first agent and was passed on to the second in its brief: that search A's
configurations "roughly tripled at the last step". Nothing here used it; for the record, 5,663 / 1,777 ≈ 3.2 is
consistent with it. The inputs taken from earlier checks (ten- to fourteen-check's configuration and spread files,
Kurz's table) are the ones the plan names.

## 6b. After the verdict: the maximal orbits the sample missed, and Langevin–Hou (Rowen, 2026-10-30)

At the page review's suggestion. The 1-in-1000 sample of raw maximal cliques reaches 1,335 of the 1,348 max orbits;
`nine_check_slices.py maxmissed` (job 20261029-230316-check9-maxmissed) took every raw clique of Part A on the other
13 (orbits 4326, 4328, 4661, 5115, 5122, 5127, 5138, 5191, 5498, 5651, 5655, 5661, 5662: 14,646 cliques, `orbitcov`
on all 34 maxc files) and `cross maxmissed` searched them with `searchmax`: nodes 14,646 / 410,088 / 11,482,464 /
177,518,586 / 562,247,627 / 70,164,318 / 133,689 / 0 / 0, no node with every off-D point covered, **0 solutions**
(`maxmissed.out` / `.prog`, tallied 14,646 once each). So every max orbit has now been searched from raw cliques too.
Langevin and Hou (Counting partial spread functions in eight variables, IEEE Trans. Inf. Theory 57(4), 2011; their
`fix-9.txt`) list 5,442 classes of 9-solid partial spread functions with their stabiliser orders: the distribution of
|Aut H| over the 5,442 hole sets here is identical, and the totals of point sets agree (24,520,650,576,127,040,978,944).

## 7. Comparison with search A (read after the verdict above was signed and committed; Rowen, 2026-10-29)

README "Nine" and `runs/i8/` (`p9.txt`, `p9.out`, `cover9.out`):
- **Classification.** Search A's `pspread3 enum 9` finds the same **5,663 configurations, 1,348 / 1,916 / 874 / 227 /
  623 / 675** by extension (its w = 0 / 1 / 2 / 3 / 4 / 8 = max / ext10 / ext11 / ext12 / ext13 / spr here), the same
  reduced count (65,323,403 in all, its w = 0 16,492,152 = the maximal cliques here) and the same weighted counts type
  by type (12,409,921,496 / 18,422,651,520 / 7,504,126,560 / 1,911,047,040 / 3,937,565,520 / 3,228,877,176). Mapped
  directly (`nine_check mapa runs/check9/configs.txt runs/i8/p9.txt`, `runs/check9/mapa.out`; the mode's reading of
  A's lines checked against p9.txt's format first: name, member count, 4 basis bytes per solid): **a bijection** onto
  the 5,663 here, w and type agreeing on every one, and A's |Aut|, this code's |Aut| of A's representative and this
  classification's |Aut| all equal on every one. As at i = 5-7, the two classifications share the method's shape, so
  the matching counts say both are complete for the same reduction; the independent checks are the ones in section 4.
- **Search.** A's `cover3` on all 5,663: **0 solutions**, 321,293,681 nodes, by dimension 1..8 158,564 / 4,439,792 /
  68,941,328 / 220,391,315 / 27,185,773 / 176,554 / 355 / 0. Here (asc) 158,564 / 4,439,792 / 68,939,870 /
  220,381,954 / 28,376,993 / 178,277 / 351 / 0: equal through dimension 2, then apart (different representatives and
  branching), as at i = 7.
- **Near-misses** (`apn/nine_check_compare.py`). Dimension 5: A's 18 lie on the same 7 configurations as here, with
  the same count on each. Dimension 6: on every configuration without a dim-5 near-miss, A's count equals ours in all
  three orders (0 unequal); the 7 dim-5 configurations carry the difference, where ours also counts a dim-5 near-miss
  extended by a rank-6 form in D and A doesn't (A's `cover` logs NEAR only at the step where the last rank-4 form
  joins). Separating those by containment: of our dim-6 spans, **395** (asc, 3,922 in all) and **395** (rand 7, 3,768)
  contain no dim-5 near-miss, A's 395; of our dim-7 spans, **187** (asc, 296) and **187** (rand 7, 321) contain no
  lower near-miss, A's 187. (desc's P1 spans are incomplete, so no containment count there; its per-configuration
  counts behave as asc's and rand 7's.) A has dim-7 near-misses on 85 configurations, ours on 91-92: the extra ones hold
  only extensions.
- Same answer: no function, by two searches and two classifications that agree configuration for configuration,
  near-misses included once the logging convention is accounted for.

## State (2026-10-29 21:00)

Search B for i = 8 is complete; nothing running. Everything is reproducible through `apn/nine_check_slices.py` (modes
parta, merge, configs, search, maxsample, cross ssp|s13|s10|maxs1000, orbitcov, controls, ctlreps; submit through the
share runner) and tallied by `apn/nine_check_tally.py prog` and `apn/ten_check_tally.py reps|maxall`. `runs/check9/`
holds ~0.7 GB, mostly the maxc files (deletable now that section 7 is done). Section 7 (comparison with A) done
2026-10-29 by Rowen: `mapa` read A's p9.txt correctly (bijection), `apn/nine_check_compare.py` for the near-misses.
Next: findings page update, bundle, fresh review, push, as for ten.
