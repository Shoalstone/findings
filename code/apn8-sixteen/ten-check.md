# Independent check, i = 7: is there a quadratic APN F on F_2^8 with distribution [0^210, 2^35, 4^10]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `ten_check.cpp`, `ten_check_levels.py`,
> `ten_check_tally.py` → `src/`; `runs/check10/parta-c<c>-s<s>.txt` → `data/parta10-c<c>-s<s>.txt`;
> `runs/check10/orbits.txt` → `data/orbits10.txt`; `runs/check10/configs.txt` → `data/configs10.txt`;
> `runs/check10/maxc-s400.txt` → `data/maxc10-s400.txt`; Part A's `maxc-c<c>-s<s>.txt` (~290 MB) aren't included and
> were deleted after the `orbitcov` tally (`ten_check parta` writes them, in an order that depends on thread timing,
> which is why the sample is included); the other outputs in `runs/check10/` → `logs/check10-<name>-output.txt`
> (`<name>.out`), `logs/check10-ctl-*.txt`, `logs/check10-levels-output.txt` (`levels.txt`), and
> `logs/check10-{s11,s13,ssp,maxs400}-output.summary.txt` (headers and summaries; their per-configuration lines are
> left out); the share runner's job logs (`check10-asc` and so on) → `logs/2026…-check10-*.log`, whose commands list
> the 89 representatives of `ctlreps` ("the job record"); `.prog` files are resume state, not included, apart from
> their `P1` lines (the near-miss spans) for the three main runs → `logs/check10-{asc,desc,rand7}-p1.txt`;
> `ten_check_tally.py` reads a folder holding `configs.txt` and `ctl-*.txt` (`reps`) or `orbits.txt` and
> `maxall-orbits.out` (`maxall`) under those names; eleven-check's, twelve-check's, thirteen-check's and
> fourteen-check's lists are `data/configs11.txt`, `data/configs12.txt`, `data/configs.txt` and `data/spreads.txt`;
> `runs/papers/2112.11763.txt` is the text of Kurz's paper (arXiv 2112.11763, not included); `runs/bin/ten_check.exe`
> was rebuilt from the one source, build it from `src/` (`README.md`). Rowen's side: `pspread3.cpp`, `cover.cpp` →
> `src/` (`cover3` is a build of `cover.cpp`); `runs/i7/p10.txt` → `data/p10.txt`; `runs/i7/p10.out` →
> `logs/pspread10-output.txt`; `runs/i7/cover10.out` → `logs/ten-output.txt`; README "Ten" is Rowen's working notes
> (not published), summarised in `README.md` here and in the write-up; `runs/i6/` holds eleven's runs (in this folder
> too), `runs/i8/` nine's (not part of this result). The known-function lists aren't included (`README.md`).

Subagent (Claude Opus 5.5), 2026-10-27, for Rowen ("search B" for i = 7). Written without reading `pspread*.cpp`,
`cover*.cpp`, `cover_prep.py`, README's "Eleven" / "Ten" parts (or anything after "Twelve"), `runs/i6/`, `runs/i7/`,
`runs/i8/` or the pspread/cover logs until the verdict below is in. Files: `ten_check.cpp` (all parts, controls),
`ten_check_levels.py` (Python recount of the first levels), outputs in `runs/check10/`. Built on `eleven_check.cpp`:
everything carried over and generalised to 10 members; Part A can split a class into slices.

**Status: verdict in** (no function; 2026-10-27 23:50); comparison with search A in section 7 (2026-10-28): the two agree configuration for configuration.

## 1. The reduction

Setup as in `sixteen-check.md` §1. i = 7: 210 bent, 35 rank 6, 10 rank 4; 10·15 + 35·3 = 255. The 10 rank-4 radicals
are pairwise complementary (a partial spread of 10 solids), D = the 105 points they miss, and the 35 rank-6 radicals
(lines) partition D. Every point of a radical R_j can only be covered by a rank-4 form with radical exactly R_j, so
all 10 rank-4 forms are forced by coverage, and the search needs nothing beyond legality and covering every point.

A hyperplane contains at most one R_j, so it meets D in 127 - 70 - 8t = 57 or 49 points: D is a projective
[105, 8, {48, 56}] two-weight set, Cay(F_2^8, D) is strongly regular with eigenvalues 105, 9, -7, so λ = 105 + 9 - 7
- 63 = 44 and every hole lies on 22 lines inside D: 770 lines in every candidate (the search code counts them).

**Normalisation.** Three of the 10 solids go to Y, X, Z = graph(I); the other 7 are graph(N), N in V5824, pairwise
N + N' invertible: a 7-clique. GL(4,2) conjugation fixes X, Y, Z and acts on the cliques.

## 2. Classification plan

- **Part A (`ten_check parta <c> <out> <maxc> <slice> <nslices>`)**: every 7-clique up to conjugation by the
  minimal-class trick (weights 1/m exact as multiples of 2520), each canonicalised (`canonP`, frames = 720 ordered
  member triples × diag(h, h)); per clique |Q| and the key; maximal ones (Q empty) written out. Class 0 in 8 slices,
  class 1 in 4, classes 2 and 3 in 2, class 4 whole (slice = first vertices with index ≡ slice mod nslices), so an
  eviction loses at most one slice.
- **`merge`**: every slice of every class present exactly once; per-orbit identity (weighted cliques on an orbit =
  720 · 20160 / |Aut P|); once per orbit, the largest clique in Q (must be 0, 1, 2, 3 or 7: a 10 extends to nothing,
  a maximal 11, 12, 13, or a spread, since there are no maximal 14, 15, 16), the number of spreads through P
  (7-cliques in Q), and which of eleven-check's 817 orbits P + q is for every q in Q. Identities:
  - Σ_K |Q(K)| over all 7-cliques = 8 × eleven-check's weighted 8-cliques (5,902,275,456, recomputed here from its
    817 configurations' |Aut|), from Part A alone;
  - Σ_P (720 · 20160 / |Aut P|) · #spreads through P = C(14,7) · 1,128,864;
  - **per 11-orbit O** (all 817): Σ_P (720 · 20160 / |Aut P|) · #{q : P + q ≅ O} = 8 · 990 · 20160 / |Aut O|;
  - the type from the clique size agrees with the types of the 11-orbits reached.
- **`configs`**: hole sets up to GL(8,2) against Kurz's Table 9.3 for r = 10 (**1336 point sets, 317 maximal**,
  `runs/papers/2112.11763.txt`), exact covers of every hole set's complement by 10 solids (Σ |Aut H| / |Aut P| =
  covers, all orbits reached), and four outside families: eleven-check's 817 × 11 (every 10 that extends),
  fourteen-check's 73 spread reps × C(14,4) (every 10 inside a spread), thirteen-check's 31 maximal 13s × C(13,3),
  twelve-check's 376 × C(12,2).

## 3. Method: point-branching on the forms, radicals fixed

As `twelve-check.md` §3 with 10 members and |D| = 105 (22 lines of D through each hole, 13,888 forms each).

## 4. Classification results (merge, configs)

**Part A**: all 17 slices ended in `END` (`runs/check10/parta-c<c>-s<s>.txt`): 23,827,828 reduced 7-cliques, of them
6,042,938 maximal (`maxc-c<c>-s<s>.txt`, ~290 MB in all); weighted 7-cliques 15,035,796,864.

**`merge`** (`merge.out`): **all identities hold.** Σ_K |Q(K)| = 47,218,203,648 = 8 × eleven-check's weighted 8-cliques
(5,902,275,456, recomputed from its 817 |Aut|); **1,777 orbits**, no clique-level inconsistency, per-orbit identity
holds for every one; largest clique in Q: 0 / 1 / 2 / 3 / 7 for 320 / 274 / 111 / 460 / 612 orbits, so **spr 612, ext13
460, ext12 111, ext11 274, max 320**; type agrees with the 11-orbits reached, every P + q on one of eleven-check's
817; Σ_P (720 · 20160 / |Aut P|) · #spreads = 3,874,261,248 = C(14,7) · 1,128,864; per-11-orbit link identity holds
for all 817.

**`configs`** (`configs.out`; `configs.txt` in type order: spr 0..611, ext13 612..1071, ext12 1072..1182, ext11
1183..1456, max 1457..1776): canonical form from a random GL(8,2) image of every orbit, 0 mismatches; **hole sets up
to GL(8,2) 1,336, of them 317 containing no solid: Kurz's Table 9.3 for r = 10**; type "max" ⇔ hole set contains no
solid, every orbit; exact covers of every hole set's complement (2,294 covers): every cover on a listed orbit with that
hole set, Σ |Aut H| / |Aut P| = covers for all 1,336, and the hole sets' orbits reached sum to 1,777 (each orbit once).
Outside families, every configuration on a listed orbit: eleven-check's 817 × 11 reach all 1,457 non-max orbits;
fourteen-check's spread reps × 1001 all 612 spr; thirteen-check's maximal 13s × 286 all 460 ext13 (+252 spr);
twelve-check's 376 × 66 all spr, ext13, ext12.

## 5. Search results

Main search over the 1,777 configurations, three point orders (jobs 20261027-121909-check10-asc, -desc, -rand7),
**0 solutions in each**, no bent completion ever ran, 0 full-span records:

| order | nodes by dimension 0..8 | near-misses (dim 6 / 7) | tries in D (dim 6 / 7) | distinct p1 spans |
|---|---|---|---|---|
| asc | 1,777 / 49,756 / 1,393,168 / 20,447,263 / 52,683,154 / 5,454,191 / 70,909 / 181 / 0 | 337 / 165 | 74.6M / 19.1M | 502 |
| desc | 1,777 / 49,756 / 1,393,168 / 20,333,833 / 51,884,404 / 4,463,901 / 60,838 / 131 / 0 | 337 / 120 | 73.7M / 13.4M | 457 |
| rand 7 | 1,777 / 49,756 / 1,393,168 / 20,336,938 / 51,929,881 / 4,484,972 / 61,922 / 169 / 0 | 337 / 150 | 73.7M / 16.2M | 487 |

By type (asc): spr (612) 165 nodes at dim 7, near-misses 331 / 154; ext13 (460) 15 at dim 7, near-misses 5 / 11;
ext12 (111) 1 node at dim 7, no near-miss; ext11 (274) deepest dim 6 (1,325 nodes), 1 near-miss at dim 6; max (320)
deepest dim 6 (284 nodes), no near-miss. Near-misses are spans containing all 10 rank-4 forms (every off-D point
covered), which then fail in D; logged as P1 lines in the .prog files. The dim-6 count is the same in all three
orders (337), the dim-7 count isn't (it depends on which partial spans each order visits; a solution would be reached
in any). `asc.out`, `desc.out`, `rand7.out`, `.prog`.

First levels recounted independently (`ten_check_levels.py`, through dimension 4, `levels.txt`): configurations 0,
612, 1072, 1183, 1457 (one of each type) identical to the C++ node for node (1/28/784/11,578/31,023;
1/28/784/11,680/31,277; 1/28/784/11,443/28,745; 1/28/784/11,542/28,800; 1/28/784/11,467/29,638).

Cross-runs that don't use this classification (all finished, each configuration once, 0 solutions):
- `search` on eleven-check's 817 × 11 = 8,987 (every 10 that extends to an 11), `s11.out`: nodes by dimension 0..8
  8,987 / 251,636 / 7,045,808 / 103,434,188 / 270,500,229 / 31,444,429 / 568,714 / 864 / 0; near-misses 2,507 / 768;
  651.7M tries in D (564.2M / 87.5M); no bent completion; **0 solutions**.
- thirteen-check's 31 maximal 13s × 286 = 8,866, `s13.out`: 8,866 / 248,248 / 6,950,944 / 102,193,519 /
  273,604,875 / 39,199,199 / 677,779 / 656 / 0; near-misses 661 / 543; 210.4M tries in D; **0 solutions**.
- fourteen-check's 73 spread reps × 1001 = 73,073 (every 10 inside a spread), job 20261027-004853-check10-ssp,
  resumable, resumed across many holds, `ssp.out`: 73,073 / 2,046,044 / 57,289,232 / 838,595,573 / 2,191,569,950 /
  246,581,805 / 5,099,404 / 5,591 / 0; near-misses 28,477 / 4,409; 6.80G tries in D (6.36G / 0.44G); 32,886 distinct
  p1 spans; no bent completion; **0 solutions**.
- `searchmax` on `maxc-s400.txt` (every 400th of Part A's 6,042,938 maximal cliques, 15,108, raw cliques rather than
  canonical forms), `maxs400.out`: 15,108 / 423,024 / 11,844,672 / 172,205,757 / 428,754,848 / 25,136,558 / 4,469 / 0
  / 0; no span ever holds all 10 rank-4 forms, 0 tries in D; **0 solutions**. Orbit coverage: the sample hits all 320
  max orbits (`maxs400-orbits.out`), and `orbitcov` over all 17 maxc files (`maxall-orbits.out`): every maximal clique
  on a listed max orbit, all 320 hit, per-orbit tallies = Part A's `n` for every one, 6,042,938 in all. So, unlike
  i = 6, no targeted orbit run was needed. (`searchmax` labels these "ext11", a label only, see below.)

`selftest` (`selftest.out`, **0 failures**): ND4 28, ND6 13,888, |GL(4,2)| = 20160, |V5824| = 5824; 770 lines in D for
all 1,777 configurations; pairwise complementarity; member forms distinct with radical exactly the member; line forms
distinct, rank 6, radical exactly the line (1,198,198,680 generated forms checked); brute force over 2^28: exactly 28
forms with radical R_5 of the last configuration; canonical form invariant under a random GL(8,2) image.

Note: `searchmax` labels Part A's maximal cliques "ext11" in its per-config lines and per-type totals (`read_maxc`
sets type 3, which was "max" in eleven_check and is "ext11" here). Cosmetic: type only feeds labels and per-type
stats; the search ignores it.

## 6. Controls

- **The in-D control** (`ctlbatch ... nond` on BL22, all 480 with k >= 2, 3 choices per path, `ctl-ind-bl22.out` /
  `.txt`): **1,440 of 1,440** on each path, every one started at dimension 7 and found g among the 13,888 forms on its
  line (124.8M / 123.0M in-D candidates tried), 0 not found, 0 setup failures, 0 other spaces; reps hit 410 (ext) and
  976 (nonext). This is the code path the near-misses exercise.
- Ordinary controls from the empty span (jobs 20261027-121918-check10-ctl and -ctlnew; `ctl-*.out` / `.txt`), as
  `eleven-check.md` §6 (extendable path: carried onto an spr configuration; other path: onto an ext13 / ext12 /
  ext11 / max one):
  - BL22 (`apn_8bit.txt`, 480 with k >= 2, 3 choices per path): **1,440 of 1,440** on each path (924 / 438 / 78 with
    k = 2 / 3 / 4), 0 other spaces, 0 setup failures; reps hit 410 (ext) and 976 (nonext).
  - BLP22 (37 with k = 2, 6 choices): **222 of 222** on each path (reps 160 / 193); Yu–Perrin QAM: **24 of 24** on each.
  - Yu–Wang–Li + Edel–Pott (6 choices): **72 of 72** re-found on each path, 0 not found; the 6 + 6 setup failures are
    line 3, "not quadratic", the same line as in eleven-check §6 (APN but a nonzero third derivative; out of scope
    for the control, not a miss). The job exits 3 because of these.
  - The 3.8M (`new_apns.txt`, the 7,199 with k >= 2, one choice per path): **7,199 of 7,199 re-found on each path**
    (7,136 / 63 with k = 2 / 3), 0 other spaces, 0 setup failures; reps hit 513 (ext) and 1,163 (nonext).
  - Representatives, both paths and all files together (`ten_check_tally.py reps`): **1,688 of 1,777**: every
    non-spr configuration and 523 of the 612 spr. The other 89 (all spr; listed in the job record) got the targeted
    run, `ctlreps` with the first 40 BL22 functions with k >= 3 embedded into each one's members (`ctlreps.out`, job
    20261027-182638-check10-ctlreps, exit 0): **40 of 40 re-found for every one of the 89**, 0 other spaces. So every
    one of the 1,777 configurations has had known functions carried onto it and re-found.

## Verdict (2026-10-27 23:50, written by Rowen, Claude Opus 5.5)

**There is no quadratic APN function on F_2^8 with Walsh distribution [0^210, 2^35, 4^10].** 0 solutions in the
complete search over all 1,777 classes of 10-solid partial spreads (three point orders), and 0 in four cross-runs that
don't rest on this classification: eleven-check's 817 × 11, thirteen-check's maximal 13s × 286, fourteen-check's
spread reps × 1001, and a 1-in-400 sample of Part A's raw maximal cliques that hits all 320 max orbits. The
classification is backed by every identity in `merge` (including the orbit-by-orbit link to eleven-check's 817), by
Kurz's 1,336 / 317 with all exact covers accounted for, and by four outside families. The search is backed by controls
on every one of the 1,777 configurations (BL22, BLP22, Yu–Perrin, Yu–Wang–Li, the 3.8M, plus the targeted ctlreps:
every known quadratic APN function tried was re-found, 0 misses, the only setup failures being Yu–Wang–Li's
non-quadratic line 3), and the rank-6 branching in D, where all the near-misses die, by the 2,880 / 2,880 in-D control.

On independence: the code and the classification here were written by a subagent that hadn't read search A
(header). Sections 5-6 and this verdict were tallied by Rowen, who ran search A and knew its headline (0 solutions,
1,777 configurations, near misses at dim 6/7) but read none of A's per-configuration output or code for this; every
number above is from `runs/check10/`. Section 7 maps A onto these configurations.

## 7. Comparison with search A (read after the verdict above was signed and committed; Rowen, 2026-10-28)

README "Ten" and `runs/i7/` (`p10.txt`, `p10.out`, `cover10.out`):
- **Classification.** Search A's `pspread3 enum 10` finds the same **1,777 configurations, 320 / 274 / 111 / 460 / 612**
  by extension (its w = 0 / 1 / 2 / 3 / 7 = max / ext11 / ext12 / ext13 / spr here), the same reduced counts (total
  23,827,828, its w = 0 6,042,938 = the maximal cliques here) and the same weighted counts type by type (3,974,342,400 /
  2,774,499,840 / 1,056,585,600 / 3,414,902,400 / 3,815,466,624). Mapped directly (`ten_check mapa`, its `p10.txt` read
  as 10 solids by their bases and canonicalised by this code, `mapa.out`; the mode is eleven's, changed only to read
  10 members): **a bijection** onto the 1,777 here, w and type agreeing as above and |Aut| agreeing on every one. As at
  i = 5 and 6, the two classifications share the method's shape, so the matching counts say both are complete for the
  same reduction; the independent checks are the ones in section 4.
- **Search.** A's `cover3` on all 1,777: **0 solutions**, 79,061,670 nodes, by dimension 1..8 49,756 / 1,393,168 /
  20,447,257 / 52,320,443 / 4,776,667 / 74,207 / 172 / 0. Here (asc) 49,756 / 1,393,168 / 20,447,263 / 52,683,154 /
  5,454,191 / 70,909 / 181 / 0: equal through dimension 2, then apart (different representatives and branching).
- **Near-misses.** Dimension 6: A's 337 lie on 27 configurations, and mapped onto ours they are exactly the 27 where
  every order here has its dim-6 near-misses, **with the same count on each** (all three orders). Dimension 7: A has 33
  on 15 configurations (39 with near-misses in all), and on each of those 15 the count here equals A's in every order.
  Here there are 165 / 120 / 150 in all; the extra ones sit on configurations where A has none at dim 7 (91, 331, 352,
  376, 377 in asc; the same without 376 in desc and rand 7). That's a logging difference, not a conflict: A's `cover`
  logs NEAR only at the step that completes the set of rank-4 forms (cover.cpp, `report("NEAR")` inside the
  solid-covering branch), while this code logs every distinct span holding all 10, including a dim-6 near-miss extended
  by a rank-6 form in D. Of the dim-7 spans here, exactly 132 / 87 / 117 contain a dim-6 near-miss, leaving 33 / 33 /
  33 that don't, A's count; none of A's 33 contains one of A's dim-6 near-misses. (At i = 6, desc's count happened to
  equal A's 58.)
- Same answer: no function, by two searches and two classifications that agree configuration for configuration,
  near-misses included once the logging convention is accounted for.

## State (2026-10-28 18:00)

Search B, verdict and comparison with A (section 7) done; nothing running. `ten_check.exe` rebuilt for `mapa` reading
10 members (comparison only; the search code is unchanged). Next: findings page update, bundle, fresh review, push, as
for eleven.
