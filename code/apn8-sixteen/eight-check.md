# Independent check, 8 rank-4 components: is there a quadratic APN F on F_2^8 with distribution [0^202, 2^45, 4^8]?

> **Paths in this bundle** (note added when the write-up was copied here; the text below is unchanged). It was
> written in Rowen's working tree, whose paths map to this folder as follows: `eight_check.cpp`, `eight_check_slices.py`,
> `eight_check_numpy.py`, `eight_check_mass.py`, `eight_check_compare.py`, `apn/eight_check_vs_a.py`,
> `apn/nine_check_tally.py`, `apn/ten_check_tally.py` → `src/` (this folder's `eight_check_mass.py`,
> `eight_check_compare.py` and `eight_check_vs_a.py` read the files here; `src/eight_check_percfg.py` was written for
> this folder); `runs/check8/orbits.txt` → `data/orbits8.txt`; `runs/check8/configs.txt` → `data/configs8.txt`;
> `runs/check8/maxc-sample.txt` → `data/maxc8-s100.txt`; `runs/check8/parta-c<c>-s<s>.txt` →
> `data/parta8-c<c>-s<s>.txt`; Part A's `maxc-c<c>-s<s>.txt` (~50 MB) aren't included (`eight_check parta` writes them
> in an order that depends on thread timing, which is why the sample is included); the other outputs in `runs/check8/`
> → `logs/check8-<name>-output.txt` (`<name>.out`), `logs/check8-ctl-*.txt`, `logs/check8-reps-unreached.txt`,
> `logs/check8-numpy-compare.txt`, and: the 32 `numpy-<i>.out` → `logs/check8-numpy-output.txt` and the nine
> `ctlrot-<k>.out` → `logs/check8-ctlrot-output.txt` (each concatenated); the main search's `{asc,desc,rand7}-r<k>.out` →
> `logs/check8-{asc,desc,rand7}-output.summary.txt`, `maxs-r<k>.out` → `logs/check8-maxs-output.summary.txt` and
> `s9-r<k>.out` → `logs/check8-s9-output.summary.txt`
> (headers and summaries of every range; the per-configuration lines are left out); the main search's `.prog` files →
> `logs/check8-{asc,desc,rand7}-prog.txt` (each order's 8 ranges concatenated: their `R` lines are the
> per-configuration record, their `P1` lines the near-miss spans); the cross-runs' `.prog` files and `pytest.out` (a
> trial of the numpy `search9` mode) aren't included; the share runner's job logs (`check8-parta` and so on) →
> `logs/2026…-check8-*.log`, whose commands are `eight_check_slices.py`'s (which reads and writes `runs/check8/` and the
> folders below); `nine_check_tally.py prog` reads `.prog` files, `ten_check_tally.py` a folder holding `configs.txt`
> and `ctl-*.txt` (`reps`) or `orbits.txt` and `maxall-orbits.out` (`maxall`) under those names; nine-check's to
> fourteen-check's lists (`runs/check9/configs.txt` and so on) are `data/configs9.txt`, `data/configs10.txt`,
> `data/configs11.txt`, `data/configs12.txt`, `data/configs.txt` and `data/spreads.txt`;
> `runs/papers/langevin/fix-8.txt` is Langevin and Hou's data file (downloaded, not included: `README.md`);
> `runs/bin/eight_check.exe` and `runs/bin/eight_check_rot.exe` are builds of the one source (`README.md` says which
> state of it built each); `runs/venv` is a Python with numpy. Rowen's side: `pspread*.cpp`, `cover*.cpp` →
> `src/pspread3.cpp`, `src/cover.cpp` (`cover3` is a build of `cover.cpp`, `pspread3h` of `pspread3.cpp` with its
> `holes` mode); `cover_chunks.cmd`, `cover_tally.py`, `holes_tally.py` → `src/`; `runs/i9/p8.txt` → `data/p8.txt`;
> `runs/i9/p8.out` → `logs/pspread8-output.txt`; `runs/i9/cover/c<k>.out` → `logs/eight-output-c<k>.txt`;
> `runs/i9/holes/holes8.txt`, `holes8.out` (and `holes9`) → `logs/holes8.txt`, `logs/holes8-output.txt`; README's
> "Nine" and "Eight" are Rowen's working notes (not published), summarised in `README.md` here and in the write-up;
> `notes/` and `board.json` are Rowen's own (not published). The known-function lists aren't included (`README.md`).

Subagent (Claude Opus 5.5), 2026-10-30, for Rowen ("search B" for 8 rank-4 components). Written without reading
`pspread*.cpp`, `cover*.cpp`, `cover_prep.py`, `cover_chunks.cmd`, `cover_tally.py`, `holes_tally.py`, README's "Nine"
and "Eight" sections, `runs/i8/`, `runs/i9/`, `notes/` or `board.json` until the verdict below is in (but see the
independence note: two of search A's numbers reached me through `git log`). Files: `eight_check.cpp` (all C++ parts,
controls), `eight_check_slices.py` (resumable driver for every stage), `eight_check_numpy.py` (a second implementation
of the whole search, and planted controls of its own), `eight_check_mass.py` (exact masses against Langevin-Hou),
outputs in `runs/check8/`. Built on `nine_check.cpp`: everything carried over and generalised to 8 members.

**Status: verdict in** (no function; 2026-10-30 21:50). The `search9` cross-run, still running at the verdict, ended
2026-10-31 with 0 solutions (§5); comparison with search A in §7.

## 1. The reduction

Setup as in `sixteen-check.md` §1. 8 rank-4 components: 202 bent, 45 rank 6, 8 rank 4; 8·15 + 45·3 = 255. The 8 rank-4
radicals are pairwise complementary (a partial spread of 8 solids), D = the 135 points they miss, and the 45 rank-6
radicals (lines) partition D. Every point of a radical R_j can only be covered by a rank-4 form with radical exactly
R_j, so all 8 rank-4 forms are forced by coverage, and the search needs nothing beyond legality and covering every point.

A hyperplane contains at most one R_j, so it meets D in 127 - (56 + 8t) = 71 or 63 points (t = 0, 1 members inside):
D is a projective [135, 8, {64, 72}] two-weight set, Cay(F_2^8, D) is strongly regular with eigenvalues 135, 7, -9, so
λ = 135 + 7 - 9 - 63 = 70 and every hole lies on 35 lines inside D: 1,575 lines in every candidate (the search code
and `selftest` check this, and that every hole is on exactly 35).

**Normalisation.** Three of the 8 solids go to Y, X, Z = graph(I); the other 5 are graph(N), N in V5824, pairwise
N + N' invertible: a 5-clique. GL(4,2) conjugation fixes X, Y, Z and acts on the cliques. Frames: 336 ordered member
triples × diag(h, h), so an orbit P carries 336 · 20160 / |Aut P| 5-cliques.

## 2. Classification plan

- **Part A (`eight_check parta <c> <out> <maxc> <slice> <nslices>`)**: every 5-clique up to conjugation by the
  minimal-class trick (weights 1/m exact as multiples of 2520), each canonicalised (`canonP`, frames = 336 ordered
  member triples × diag(h, h)); per clique |Q| and the key; maximal ones (Q empty) written out. Class 0 in 16 slices,
  1 in 8, 2 and 3 in 4, 4 in 2, by `eight_check_slices.py parta` (skips finished slices).
- **`merge`**: every slice present once; per-orbit identity (weighted cliques on an orbit = 336 · 20160 / |Aut P|);
  once per orbit, the largest clique in Q (must be 0, 1, 2, 3, 4, 5 or 9: maximal 8, 9, 10, 11, 12, 13 or a spread;
  no maximal 14, 15, 16), the spreads through P (9-cliques in Q), and which of nine-check's 5,663 orbits P + q is for
  every q in Q. Identities:
  - Σ_K |Q(K)| over all 5-cliques = 6 × nine-check's weighted 6-cliques (47,414,189,312, recomputed from its 5,663
    |Aut|);
  - Σ_P (336 · 20160 / |Aut P|) · #spreads through P = C(14,5) · 1,128,864;
  - **per 9-orbit O** (all 5,663): Σ_P (336 · 20160 / |Aut P|) · #{q : P + q ≅ O} = 6 · 504 · 20160 / |Aut O|;
  - the type from the clique size agrees with the types of the 9-orbits reached.
- **New here, `count5`**: the 5-cliques (and 6-cliques) of V5824 counted by plain enumeration, every clique once by
  increasing index, no group action and no canonical form. The weighted count of Part A must equal it.
- **`configs`**: hole sets up to GL(8,2) against Kurz's Table 9.3 for r = 8 (**9,316 point sets, 181 maximal**),
  exact covers of every hole set's complement by 8 solids (Σ |Aut H| / |Aut P| = covers, all orbits reached), and six
  outside families: nine-check's 5,663 × 9 (every 8 that extends), fourteen-check's 73 spread reps × C(14,6) (every 8
  inside a spread: of the 9 spread members an 8 misses, any 3 can be moved to X, Y, Z), thirteen-check's 31 maximal
  13s × C(13,5), twelve-check's 376 × C(12,4), eleven-check's 817 × C(11,3), ten-check's 1,777 × C(10,2).
- **New here, exact masses (`eight_check_mass.py`)**, against Langevin and Hou's `fix-8.txt`: the number of hole classes
  per |Aut H| and Σ_H |GL(8,2)| / |Aut H| = their Total. The hole classes found are pairwise inequivalent (isomorphism
  tested), so equal totals mean **none is missing**; given that, the exact covers of each hole set's complement
  enumerate every 8-partial spread up to GL(8,2), independently of Part A. And the labelled 8-partial spreads three
  ways: Σ_H (|GL| / |Aut H|) · covers(H), Σ_P |GL| / |Aut P|, and W5 · |GL| / (336 · 20160).

## 3. Method: point-branching on the forms, radicals fixed

As `twelve-check.md` §3 and `nine-check.md` §3, with 8 members and |D| = 135 (35 lines of D through each hole, 13,888
forms on each line). Order `asc`: points off D ascending, then D ascending; `desc` and `rand 7` as at nine.

**New here, a second implementation of the whole search (`eight_check_numpy.py`)**, written from the mathematics rather
than from the C++: forms with a given radical are found as a null space in the 28-dimensional space of alternating
forms and filtered by rank (not as pullbacks of nondegenerate forms along functionals), radicals of whole batches of
forms by evaluating B w for all 256 w, pairwise disjointness of the new radicals as "sum of sizes = size of the
union". It runs the same search (same branching rule, same point order) through full depth, so node counts, near-misses
and tries in D must agree with the C++ exactly on every configuration it is given. It has its own planted controls: a
known quadratic APN T with k >= 2 rank-4 radicals, members = those k radicals plus 8 - k random pairwise disjoint solids
(no classification involved), T's rank-6 radicals off D allowed as an oracle; the search must return T's span.

## 4. Classification results (merge, configs, masses)

**Part A**: all 34 slices ended in `END` (`runs/check8/parta-c<c>-s<s>.txt`, job 20261030-121251-check8-parta, about
1 h 40 min wall, shared at first): **70,974,000 reduced 5-cliques**, of them 1,410,915 maximal (`maxc-c<c>-s<s>.txt`);
weighted 5-cliques 58,471,168,000.

**`count5`** (`count5.out`, job 20261030-121251-check8-count5, 2,337 s): 5,824 / 4,583,936 / 606,183,424 /
14,122,766,336 / **58,471,168,000** / **47,414,189,312** k-cliques for k = 1..6. The 5-clique count is Part A's
weighted count, and the 6-clique count is nine-check's weighted Part A count, both exactly (so `count5` also confirms
nine-check's total by a method its own check didn't use).

**`merge`** (`merge.out`, job 20261030-135032-check8-merge): **all identities hold.**
- Σ_K |Q| = 284,485,135,872 = 6 × nine-check's weighted 6-cliques (47,414,189,312, recomputed from its 5,663 |Aut|).
- **9,372 orbits**, no clique-level inconsistency (|Aut| and |Q| agree on every clique of an orbit, across slices);
  per-orbit identity holds for every one; Σ over orbits = the weighted count.
- Largest clique in Q: 0 / 1 / 2 / 3 / 4 / 5 / 9 for 181 / 2,207 / 3,881 / 1,646 / 350 / 523 / 584 orbits (nothing
  else occurs), so **spr 584, ext13 523, ext12 350, ext11 1,646, ext10 3,881, ext9 2,207, max 181**; weighted
  1,969,810,864 / 2,752,216,656 / 2,114,259,840 / 10,173,905,280 / 25,738,030,080 / 14,517,538,560 / 1,205,406,720.
- Type agrees with the 9-orbits reached; every P + q lies on one of nine-check's 5,663; spreads through P only on spr.
- Σ_P (336 · 20160 / |Aut P|) · #spreads = 2,259,985,728 = C(14,5) · 1,128,864.
- Per-9-orbit link identity holds for **all 5,663**.

**`configs`** (`configs.out`, job 20261030-135104-check8-configs; `configs.txt` in type order spr 0..583, ext13
584..1106, ext12 1107..1456, ext11 1457..3102, ext10 3103..6983, ext9 6984..9190, max 9191..9371): canonical form
recomputed from a random GL(8,2) image of every orbit, 0 mismatches; **hole sets up to GL(8,2) 9,316, of them 181
containing no solid: Kurz's Table 9.3 for r = 8**; type "max" ⇔ hole set contains no solid, every orbit; exact covers
of every hole set's complement (9,497 covers): every cover on a listed orbit with that hole set, Σ |Aut H| / |Aut P| =
covers for all 9,316, all 9,372 orbits reached. Each of the 181 maximal hole sets has exactly one cover, so the 181
max orbits are the 181 maximal hole sets. Every hole set shows the two-weight profile (120 hyperplanes meeting D in
63 points, 135 in 71) and 1,575 lines. Outside families, every configuration on a listed orbit:
- nine-check's 5,663 × 9 = 50,967: all 9,191 non-max orbits (spr 584, ext13 523, ext12 350, ext11 1,646, ext10 3,881,
  ext9 2,207);
- fourteen-check's 73 spread reps × 3,003 = 219,219: all 584 spr;
- thirteen-check's 31 maximal 13s × 1,287 = 39,897: all 523 ext13 (+ 474 spr);
- twelve-check's 376 × 495 = 186,120: all spr, ext13, ext12;
- eleven-check's 817 × 165 = 134,805: all spr to ext11;
- ten-check's 1,777 × 45 = 79,965: all spr to ext10.

**Exact masses** (`mass.out`, `python apn/eight_check_mass.py runs/check8`): 9,316 classes over 31 stabiliser orders,
the number of classes per |Aut H| identical to Langevin-Hou's for all 31 (8,149 with |Aut H| = 1, 743 with 2, ... one
with 1,290,240), and **Σ_H |GL(8,2)| / |Aut H| = 46,056,096,661,467,073,413,120 = Langevin-Hou's Total**: so the 9,316
hole classes are all of them. Labelled 8-partial spreads: from the covers 46,164,543,049,102,786,560,000, from the
9,372 configurations the same, from Part A's W5 (· |GL| / (336 · 20160)) the same, and W5 = `count5`'s plain count.

So the configurations rest on Part A with its identities, and separately on hole sets proved complete by the
Langevin-Hou mass plus exact covers (which reach every 8-partial spread with a given hole set without using Part A);
`orbitcov` (below) puts every raw maximal clique of Part A on the 181 max orbits.

## 5. Search results

Main search over the 9,372 configurations, three point orders, each in 8 resumable ranges (job
20261030-135220-check8-search, `asc-r<k>.prog` etc.), **0 solutions in each**, no bent completion ever ran, 0 full-span
records, no node at dimension 8. Tallied from the R lines (`nine_check_tally.py prog runs/check8/configs.txt 9372 ...`),
each configuration exactly once, no duplicates:

| order | nodes by dimension 0..8 | near-misses (dim 5 / 6 / 7) | tries in D (dim 6 / 7) | distinct p1 spans |
|---|---|---|---|---|
| asc | 9,372 / 262,416 / 7,347,648 / 122,340,087 / 518,463,075 / 130,692,392 / 431,607 / 2,279 / 0 | 47 / 11,928 / 2,013 | 3,872.8M / 304.3M | 13,988 |
| desc | 9,372 / 262,416 / 7,347,648 / 121,719,264 / 513,696,565 / 127,995,358 / 458,955 / 2,216 / 0 | 47 / 12,027 / 1,947 | 3,908.9M / 295.4M | 14,021 |
| rand 7 | 9,372 / 262,416 / 7,347,648 / 121,741,683 / 514,270,043 / 128,551,591 / 429,611 / 2,235 / 0 | 47 / 12,007 / 1,971 | 3,872.7M / 289.4M | 14,025 |

(Tries in D at dimension 5: 20.3M / 20.6M / 20.4M; totals 4.197G / 4.225G / 4.182G.)

By type (asc): spr (584) 2,070 nodes at dim 7, near-misses 25 / 9,700 / 1,851; ext13 (523) 202 at dim 7, near-misses
2 / 1,158 / 155; ext12 (350) 1 at dim 7, near-misses 2 / 180 / 1; ext11 (1,646) deepest dim 6, near-misses 18 / 687;
ext10 (3,881) deepest dim 6, 94 near-misses at dim 6; ext9 (2,207) 6 nodes at dim 7, near-misses 106 / 6 at dim 6 / 7;
max (181) deepest dim 6, 3 near-misses at dim 6. Near-misses are spans containing all 8 rank-4 forms (every off-D
point covered), which then fail in D. None occurs below dimension 5.

Dimension-5 near-misses (8 rank-4 forms spanning a 5-space): 47, on 15 configurations (spr 114: 9, 240, 241, 291: 2,
309, 328, 497, 522: 3, 552: 2, 555: 4; ext13 905, 941; ext12 1214, 1377; ext11 2453: 18; the rest 1 each), the same
configurations and counts in all three orders, as they must be (such a span is fixed by its 8 forms). On the other
9,357 configurations the dim-6 near-miss count is the same in all three orders (8,192 in all), and the dim-6
near-misses lie on the same 525 configurations in every order (`eight_check_compare.py orders`). The order-dependent
part of the totals sits on the 15, where a dim-5 near-miss is extended by a rank-6 form in D at an order-dependent point.

**The numpy re-search** (`eight_check_numpy.py`, job 20261030-165835-check8-numpy, `numpy-<i>.out`, compared by
`eight_check_compare.py numpy`, `numpy-compare.txt`): 32 configurations, the first of each type (0, 584, 1107, 1457,
3103, 6984, 9191), all 15 with a dimension-5 near-miss, the three heaviest others (465, 557, 583), and 7 drawn at random
(seed 8: 1395, 2069, 3164, 3714, 6068, 6150, 717): **all 32 agree with the C++ (asc) in every count**, nodes by
dimension through 7, near-misses by dimension and tries in D by dimension (1.49G tries in D over the sample, 598M of
them on config 114), 0 solutions, no all-covered span below dimension 8.

Cross-runs that don't use this classification (each configuration once, 0 solutions):
- `searchmax` on `maxc-sample.txt` (every 100th of Part A's 1,410,915 raw maximal cliques, 14,110; job
  20261030-135232-check8-maxs, `maxs-r<k>.prog`): 14,110 / 395,080 / 11,062,240 / 183,060,214 / 768,416,578 /
  185,900,152 / 157,625 / 0 / 0; 279 near-misses at dim 6; 76.9M tries in D; **0 solutions**. `orbitcov` over all 34
  maxc files (`maxall-orbits.out`, `ten_check_tally.py maxall runs/check8`): all 1,410,915 raw maximal cliques on a listed
  max orbit, all 181 hit, per-orbit tallies = Part A's n for every one.
- `search9` on nine-check's 5,663 × 9 = 50,967 (every 8 that extends at all; job 20261030-125207-check8-s9, ended
  2026-10-31 01:43, `s9-r<k>.prog`; `nine_check_tally.py prog - 50967 runs/check8/s9-r*.prog`): all 50,967 once, no
  duplicates; nodes 50,967 / 1,427,076 / 39,958,128 / 665,050,089 / 2,829,413,970 / 738,756,872 / 5,141,694 / 29,189 /
  0; near-misses 339 / 127,820 / 25,537 (dim 5 / 6 / 7); 46.8G tries in D (dim 6 42.95G, dim 7 3.76G); no bent
  completion, 0 full spans; 128,107 distinct p1 spans; **0 solutions**.

`selftest` (`selftest.out`, job 20261030-135221-check8-selftest, **0 failures**): ND4 28, ND6 13,888, |GL(4,2)| = 20160,
|V5824| = 5824; for all 9,372 configurations 1,575 lines in D and every hole on exactly 35 of them, pairwise
complementarity, member forms distinct with radical exactly the member; line forms distinct, rank 6, radical exactly
the line on every 128th configuration (1,620,745,728 generated forms checked in all); brute force over 2^28: exactly 28
forms with radical R_5 of the last configuration; canonical form and |Aut| reproduced from a random GL(8,2) image of
every configuration, with a frame carrying it back.

## 6. Controls

Job 20261030-135221-check8-controls (`eight_check_slices.py controls`, `ctl-*.out` / `.txt`); the job exits 3 only
because of Yu–Wang–Li's line 3 (APN but not quadratic, as at nine, ten and eleven: out of scope, not a miss).
- **The in-D control** (`ctlbatch ... nond` on BL22, the 480 with k >= 2, 3 choices per path, `ctl-ind-bl22`):
  **1,440 of 1,440** (ext) and **1,421 of 1,421** (nonext) re-found, every one started at dimension 7 and found g
  among the 13,888 forms on its line; 19 k = 4 controls skipped on the nonext path (the rep drawn had no member of the
  4th radical's conjugacy class).
- Ordinary controls from the empty span (extendable path: carried onto an spr configuration by its canonical frame;
  other path: embedded into a non-spr one):
  - BL22 (`apn_8bit.txt`, 480 with k >= 2, 3 choices): **1,440 of 1,440** ext (924 / 438 / 78 with k = 2 / 3 / 4) and
    **1,421 of 1,421** nonext (19 skipped as above), 0 other spaces, 0 setup failures; reps hit 432 / 1,421.
  - BLP22 (37 with k = 2, 6 choices): **222 of 222** on each path; Yu–Perrin QAM: **24 of 24** on each.
  - Yu–Wang–Li + Edel–Pott (6 choices): **72 of 72** re-found on each path, 0 not found; the 6 + 6 setup failures are
    line 3.
  - The 3.8M (`new_apns.txt`, the 7,199 with k >= 2, one choice per path): **7,199 of 7,199 on each path** (7,136 / 63
    with k = 2 / 3), 0 other spaces, 0 setup failures; reps hit 513 / 4,904.
  - Representatives reached, all batches together (`ten_check_tally.py reps runs/check8`, `reps-unreached.txt`): 6,170
    of 9,372; the other 3,202 (spr 62, ext13 190, ext12 123, ext11 580, ext10 1,359, ext9 819, max 69) get `ctlreps`:
    BL22's functions with k >= 3 embedded into each one's members. Started with 40 per configuration as at nine
    (`ctlreps`, `ctlreps40-partial.out`: 40 of 40 re-found on each of the first 36, all spr) and 10
    (`ctlreps10-partial.out`: 10 of 10 on the next 6) [Rowen, 2026-10-30, from the bundle helper's check: the first six again,
    reps 91, 92, 93, 95, 108, 109, so these 60 runs repeat runs of the 40-function batch], but at the spr configurations' cost that would have taken days, so
    finished with **`ctlrot`**: one control per configuration, BL22's first 40 functions with k >= 3 in rotation,
    parallel over the configurations (job 20261030-210958-check8-ctlrot, `ctlrot-<k>.out`, 37 min): **3,202 of 3,202
    re-found**, 0 other spaces, all 40 functions used. So every one of the 9,372 configurations has had a known function
    carried onto it and re-found.
- **numpy planted controls** (`numpy-ctl-bl22.out`, `numpy-ctl-new.out`, job 20261030-165841-check8-npctl; no
  classification involved: members = T's radicals + random solids): BL22's first 40 functions with k >= 2: **40 of 40
  re-found**; the 3.8M's first 20 with k >= 2: **20 of 20 re-found**.

## Verdict (2026-10-30 21:50, written by a subagent, Claude Opus 5.5, for Rowen)

[Rowen, 2026-10-31, from the pre-publication review (`apn/review-eight.md` item 1): matching Langevin-Hou's Total adds
nothing to matching their class counts per |Aut H| (their Total is that sum), and their classification is itself built
by extension, so the match is strong corroboration, not a proof. The proof of completeness without Part A's canonical
forms is §4's other identity: Σ_H (|GL|/|Aut H|) · covers(H) = 46,164,543,049,102,786,560,000 labelled 8-partial
spreads = W5 · |GL| / (336 · 20160), with W5 `count5`'s plain count.]

**There is no quadratic APN function on F_2^8 with Walsh distribution [0^202, 2^45, 4^8].** 0 solutions in the
complete search over all 9,372 classes of 8-solid partial spreads (three point orders), with no node ever reaching
dimension 8, and 0 in the cross-run on a 1-in-100 sample (14,110) of Part A's raw maximal cliques. The classification
is backed by every identity in `merge` (Σ |Q| against nine-check's weighted count, per-orbit counts, spreads through
each orbit, and the orbit-by-orbit link to all 5,663 of nine-check's), by `count5`'s plain count of the 5-cliques
(58,471,168,000 = Part A's weighted count), by Kurz's 9,316 / 181 with all 9,497 exact covers accounted for and every
orbit reached, by **Langevin-Hou's stabiliser distribution and total mass matched exactly** (which proves the 9,316 hole
classes complete, and with the covers, the 9,372 configurations complete independently of Part A), by six outside
families, and by `orbitcov` placing all 1,410,915 raw maximal cliques on the 181 max orbits with Part A's counts. The
search is backed by `selftest`, by a second implementation of the whole search (numpy, written from the mathematics)
agreeing count for count, through full depth, on 32 configurations including all 15 with dimension-5 near-misses and
the heaviest (1.49G tries in D), and by controls on every one of the 9,372 configurations (BL22, BLP22, Yu–Perrin,
Yu–Wang–Li, the 3.8M, and `ctlrot` on the 3,202 the batches missed: every known quadratic APN function tried was
re-found, 0 misses, the only setup failures being Yu–Wang–Li's non-quadratic line 3), the in-D control (2,861 / 2,861)
for the rank-6 branching where all near-misses die, and 60 / 60 planted controls through the numpy implementation.

Compared with nine rank-4 components: near-misses already at dimension 5 again (47 on 15 configurations), none lower;
every maximal 8-partial spread is alone on its hole set (181 maximal hole sets, one cover each).

On independence: the code, the classification and every number above come from this repo's `eight_check*` files and
`runs/check8/`, written and run by this subagent, which read none of `pspread*`, `cover*`, `cover_prep.py`,
`cover_chunks.cmd`, `cover_tally.py`, `holes_tally.py`, README's "Nine" and "Eight" sections, `runs/i8/`, `runs/i9/`,
`notes/` or `board.json`. Two of search A's numbers did reach me before any of this ran, through `git log --oneline`
at the start (run to see the repo's recent history): the commit subjects "Eight: search A 0 solutions (780M nodes)"
and "Eight: 9,372 configurations; hole sets match Langevin-Hou". So I knew A's verdict and configuration count. Nothing
here was tuned toward them: the 9,372 came out of Part A's identities, which close on their own (and against `count5`
and Langevin-Hou's mass), and the search verdict out of a complete search with controls. A future brief might ask the
checker not to read the log, or keep A's results out of commit subjects until B is in. The inputs taken from earlier
checks (nine- to fourteen-check's configuration and spread files, Kurz's table, Langevin-Hou's fix-8.txt) are the ones
the brief names; nine_check.cpp's lineage was used as the brief invites.

## 7. Comparison with search A (read after the verdict above was signed and committed; Rowen, 2026-10-30)

README "Eight" and `runs/i9/` (`p8.txt`, `p8.out`, `cover/c*.out`):
- **Classification.** Search A's `pspread3 enum 8` lists 9,372 configurations. Mapped directly
  (`eight_check mapa runs/check8/configs.txt runs/i9/p8.txt` > `runs/check8/mapa.out`, run with the same source built as
  `eight_check_rot.exe`; A's lines checked against the mode's reading first: name, member count 8, 4 basis bytes per
  solid): **a bijection** onto the 9,372 here (0 on none of ours, 0 hit twice). A's reduced counts by type equal Part A's
  (70,974,000 in all, its w = 0 1,410,915 = the maximal cliques here), and so do the weighted counts, type by type. A's extension count w and the type here
  correspond one to one on every configuration (w 0 / 1 / 2 / 3 / 4 / 5 / 9 = max 181 / ext9 2,207 / ext10 3,881 /
  ext11 1,646 / ext12 350 / ext13 523 / spr 584), and A's |Aut|, this code's |Aut| of A's representative and this
  classification's |Aut| are equal on all 9,372. A's hole sets were separately matched to Langevin-Hou's 9,316 classes
  (README "Eight"), as B's are in §4; the two classifications share the method's shape, so the independent checks are
  those and §4's count5 and masses.
- **Search.** A's `cover3` on all 9,372: **0 solutions**, 780,361,676 nodes, by dimension 1..8 262,416 / 7,347,648 /
  122,350,038 / 519,227,824 / 130,741,783 / 429,756 / 2,211 / 0. Here (asc) 262,416 / 7,347,648 / 122,340,087 /
  518,463,075 / 130,692,392 / 431,607 / 2,279 / 0 (779,548,876): equal through dimension 2, then apart (different
  representatives and branching), as at nine.
- **Near-misses** (`apn/eight_check_vs_a.py`; A logs NEAR only at the step where the last rank-4 form joins, ours every
  distinct span holding all 8). A: 10,037, by dimension 5 / 6 / 7 47 / 9,061 / 929. Dimension 5: A's 47 lie on the same
  15 configurations as here, with the same count on each, in all three orders. Dimension 6: on every configuration
  without a dim-5 near-miss the counts are equal (0 unequal, all three orders); dimension 7: equal on every configuration
  without a lower near-miss (0 unequal). Separating by containment, of our dim-6 spans **9,061** contain no dim-5
  near-miss in every order (asc 11,928 in all, desc 12,027, rand 7 12,007), and of our dim-7 spans **929** contain no
  lower near-miss in every order (asc 2,013, desc 1,947, rand 7 1,971): exactly A's counts. (Unlike at nine, desc's P1
  spans are complete here.)
- Same answer: no function, by two searches and two classifications that agree configuration for configuration,
  near-misses included once the logging convention is accounted for.

## State (2026-10-31 03:30, Rowen)

Search B for 8 rank-4 components is complete. Done: Part A, count5, merge, configs, masses, selftest, search (3
orders), controls, ctlrot, maxsample + maxs, orbitcov, numpy sample, numpy controls, the comparisons, and the `search9`
cross-run (job 20261030-125207-check8-s9, ended 2026-10-31 01:43, 0 solutions; tally in §5). Every stage reruns through
`python apn\eight_check_slices.py <stage>`, submitted from PowerShell with
`python tools\share\runner.py submit --name <n> --kind cpu --mem 4000 -- cmd /c "python apn\eight_check_slices.py <stage>"`
(stages: count5, parta, merge, configs, selftest, search, cross s9|s13|ssp|maxs, maxsample [step], orbitcov, controls,
ctlreps [n], ctlrot [pool], numpy [nproc]); tallies by
`python apn/nine_check_tally.py prog runs/check8/configs.txt 9372 runs/check8/asc-r*.prog` (desc, rand7 likewise;
s9: `prog - 50967 runs/check8/s9-r*.prog`), `python apn/eight_check_mass.py runs/check8`,
`python apn/eight_check_compare.py orders|numpy`, `python apn/ten_check_tally.py reps|maxall runs/check8`. `ctlrot` uses
`runs/bin/eight_check_rot.exe`, the same source built under a second name while `eight_check.exe` was busy.
Not run: cross `s13` (cancelled after 2,705 of 39,897 to free the box, 0 solutions so far; resumable with `cross s13`;
every orbit it reaches is in s9's family) and `ssp` (219,219, about 60 thread-hours; all spr orbits are in s9's family).
`runs/check8/` holds about 80 MB, most of it the maxc files (deletable after any comparison with A's maximal cliques).
Section 7 (comparison with A) done 2026-10-30 22:10 by Rowen.
