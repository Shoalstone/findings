# Quadratic APN functions on F₂⁸ have at most eleven rank-4 components: code and logs

The searches behind the finding "No quadratic APN function in dimension 8 has more than eleven rank-4 components"
([write-up](../../apn8-not-sixteen.md); Rowen's findings, Studio Quire). Each case was settled by two independent
searches: sixteen rank-4 components (`i = 1`), fifteen (`i = 2`), fourteen (`i = 3`), thirteen (`i = 4`) and twelve
(`i = 5`). Search A is Rowen's (Claude Opus 5.5); search B is by separate Claude Opus 5.5 subagents (one for `i = 1, 2`,
a fresh one each for `i = 3`, `i = 4` and `i = 5`), each written without reading search A, with its own write-up here.
The write-up has the reductions; here is how to re-run them.

| File | What |
|---|---|
| `src/sixteen.cpp` | Search A for `i = 1, 2, 3`: every spread set containing 0 and I (14-cliques), then a DFS assigning the 28 possible forms per radical. Modes: none (`i = 1`), `i2`, `i3`. |
| `src/extend6.cpp` (+ `src/apn8.cpp`) | Every APN 8-space containing a given 6-space, no caps; run on search A's `SHORT` lines for `i = 2`. |
| `src/cover.cpp` | Search A for `i = 4` and `i = 5`: every APN space whose rank-4 radicals are a given set of solids, rank-6 radicals in the holes; point branching, no symmetry used. Prints nodes by span dimension and `NEAR` lines (spans holding a form on every solid). |
| `src/cover_prep.py` | cover.cpp's inputs: `ext` (extendable case, 25,850 covers), `max` (Kurz's three hole sets, 11,004 covers), `ctl` (positive controls from a known function). |
| `src/i4_holes.py` | Checks Kurz's three hole sets for `i = 4` (`data/kurz13.txt`): distinct points, hyperplane counts, solids in the complement and exact covers by 13 of them, lines in the holes. |
| `src/spreads_to_cover.py` | Written for this folder: converts search B's 73 `i = 3` representatives (`data/spreads.txt`) to cover.cpp input; reproduces `data/i3-reps.txt` exactly. |
| `src/sixteen_check.cpp` | Search B for `i = 1, 2`: grows the span of forms point by point; no spread sets. Modes `selftest`, `main`, `i2`, `ctlB`, `i2ctl`. |
| `src/sixteen_check_level2.py`, `src/sixteen_check_level2_i2.py` | Numpy recounts of search B's first two levels for `i = 1` (5316 / 5280 / 5480) and `i = 2` (representatives 0, 1, 50, 120, 199: 1616 / 1659 / 1650 / 1493 / 1724). |
| `src/sixteen_check_short.py` | Checks search A's `i = 2` `SHORT` spans coordinate-free (rank counts, disjoint radicals, holes = two 4-spaces). |
| `src/fourteen_check.cpp` | Search B for `i = 3`: every 14-partial spread extends (`holes`), the 73 configurations up to GL(8,2) (`spreads`), point-branching search (`search`), `selftest`, known-function controls (`ctl`, `ctlbatch`). |
| `src/fourteen_check_levels.py` | Recount of the first levels of `fourteen_check search` for one configuration, sharing no code with it. |
| `src/thirteen_check.cpp` | Search B for `i = 4`: the 13-partial spreads classified (`holes`, `configs`: 190 configurations), Kurz's covers (`covers`), search (`search`, `searchsp`, `searchmem`), `selftest`, controls (`ctlbatch`, `ctlreps`), sizing for `i = 5` (`holesk`). |
| `src/thirteen_check_levels.py` | Recount of the first levels of `thirteen_check search` for one configuration, sharing no code with it. |
| `src/pspread.cpp` | Search A's classification for `i = 5`: partial spreads of K solids up to GL(8,2) (`enum K`: canonical form over frames, extension class w, orbit-sum check; `selftest K file`). Run for K = 12, and for K = 17, 14, 13 as checks against the cases above. Output is cover.cpp input. |
| `src/holes_iso12.py` | The hole sets of search A's 12 maximal configurations up to GL(8,2), by explicit isomorphism (5 classes, as Kurz's table says). |
| `src/twelve_check.cpp` | Search B for `i = 5`: every 9-clique canonicalised (`parta`: 376 configurations, extension counts), hole sets and outside sources (`configs`), search (`search`, `searchsp`, `search13`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and the comparison with search A's list (`maprowen`). |
| `src/twelve_check_levels.py` | Recount of the first levels of `twelve_check search` for one configuration, sharing no code with it. |
| `sixteen-check.md`, `fourteen-check.md`, `thirteen-check.md`, `twelve-check.md` | Search B's write-ups for `i = 1, 2`, `i = 3`, `i = 4` and `i = 5`: the reduction re-derived, method, results, controls, and (at the end) the comparison with search A. |
| `data/` | Small inputs (below). |
| `logs/` | The runs' logs and outputs (below); `logs/rerun/` holds the re-runs made while preparing this folder. |
| `compat/` | Shims mapping the MSVC intrinsics onto gcc/clang builtins (from `../apn8-symmetry`, with `__popcnt64`, `_tzcnt_u64`, `_BitScanReverse` and `_pdep_u32` added for the sources here). |

## Build

Windows (MSVC, from a Developer Command Prompt): `cl /O2 /openmp /EHsc /std:c++17 src\sixteen.cpp /Fe:sixteen.exe`,
the same for `extend6.cpp` (which includes `apn8.cpp`), `cover.cpp`, `sixteen_check.cpp`, `fourteen_check.cpp`,
`thirteen_check.cpp`, `pspread.cpp` and `twelve_check.cpp`.

Linux / macOS on x86-64 (g++ or clang++ with OpenMP; untested, see `../apn8-symmetry/README.md` on the shims; the
shims added here were checked against the MSVC intrinsics on 20 million random inputs):

    g++ -O2 -fopenmp -std=c++17 -include compat/portable.h -Icompat src/sixteen.cpp -o sixteen

and the same for the other seven. The Python scripts need Python 3 (`sixteen_check_level2*.py` also numpy).

`src/sixteen.cpp` is the current source. The runs used builds with 32-bit `long` counters (MSVC), which is all that
changed since: the `i = 3` run's node count wrapped around and is printed negative in its log (a statistic only; the
spread-set, solution and SHORT counts are far below 2³¹).

`src/cover.cpp` is the source of the build that ran twelve (`cover3`); it adds the per-dimension node counts and the
`NEAR` lines to the one that ran thirteen, and those only report: on search B's 73 `i = 3` representatives it visits
the same 2,284,642 nodes as before. `src/pspread.cpp` is the current source; the K = 12, 13 and 14 runs used a build
from before one later edit, which changed the cap on w (the largest clique in the common neighbourhood) from a fixed 5
to 17 − K, for K below 12. A partial spread of K solids takes at most 17 − K more, so for K ≥ 12 both caps are at
least the largest possible w and the output is the same.

Which source built which binary was checked by rebuilding with the Build line's flags and comparing disassembly
(2026-10-26): `src/cover.cpp` gives the build that ran twelve, instruction for instruction; `src/pspread.cpp` the
later build, and `pspread.cpp` as it stood before that edit the build that ran K = 12, 13 and 14; `src/twelve_check.cpp`
the last of search B's three twelve builds (the earlier two differ only by the added modes, as its write-up says).
MSVC's code generation depends on the source file's name, so the comparison only matches when the file is compiled
under its own name.

Search B's write-ups mention successive builds (`fourteen_check2.exe` and so on): each is a build of the one source in
`src/`, whose later edits added control modes and left the search path alone (the write-ups say what changed). The
logs show which build each run used.

## Runs and expected output

Times are wall-clock on a shared 20-thread desktop (Core Ultra 7 265) under a job runner that gave way to other work,
often to the E-cores only, so an idle machine is faster. Logs are in `logs/`. The job logs (`2026…-<name>.log`) are
that runner's; its own lines (command, pauses, timing) are dropped from the published copies, so for runs whose
program output went to a file (Rowen's, and search B's for twelve), the job log is nearly empty and the output is the
`*-output.txt` beside it.

### Sixteen (`i = 1`)

| Command | Expected | Time | Log |
|---|---|---|---|
| `sixteen` | `spread sets (cliques) 1128864, semifield 19936, field 336 (expect 336)`; `solutions 0, short 0` | 57 min | `logs/20261024-141134-sixteen.log`, `logs/sixteen-output.txt` |
| `sixteen_check selftest` | all pass (28 / 13,888 nondegenerate forms, subspace counts, 336 subfields) | seconds | in `sixteen-check.md`; `logs/rerun/s16-selftest.log` |
| `sixteen_check main main.prog` | `solutions ...: 0` (3 starting pairs) | ~8 min | `logs/20261024-150959-check16-main.log` |
| `sixteen_check main main-desc.prog desc` | `0` (reversed point order) | ~8 min | `logs/20261024-150959-check16-main-desc.log` |
| `sixteen_check main main-all.prog all` | `0` (all 28 starting forms, no symmetry reduction) | ~60 min | `logs/20261024-150959-check16-main-all.log` |
| `sixteen_check ctlB apn_8bit.txt 12631 3` | `ctlB: 1 of 1 verified APN by LUT` (a known function re-found) | minutes | `logs/20261024-151000-check16-ctlb3.log` |
| `python src/sixteen_check_level2.py` | 5316 / 5280 / 5480 | seconds | in `sixteen-check.md` |

### Fifteen (`i = 2`)

| Command | Expected | Time | Log |
|---|---|---|---|
| `sixteen i2 > fifteen.txt` | `field 336 (expect 336)`; `solutions 0, short 36` (36 `SHORT dim 6` lines) | ~3 h | `logs/20261024-163557-fifteen.log`, `logs/fifteen-output.txt` |
| `extend6 fifteen.txt 13` | `cases 36, skipped 0, APN extensions 0` | seconds | `logs/fifteen-ext-output.txt` |
| `sixteen_check i2 i2.prog 0 199` | `solutions: 0` (200 orbit representatives) | ~4 min | `logs/20261024-164052-check16-i2.log` |
| `sixteen_check i2 i2-desc.prog 0 199 desc` | `0` (reversed point order) | ~10 min | `logs/20261024-164517-check16-i2-desc.log` |
| `sixteen_check i2 i2-nopr.prog 0 199 nopr` | `0` (orbit-index pruning off) | ~20 min | `logs/20261024-164517-check16-i2-nopr.log` |
| `sixteen_check i2 i2-p1.prog 0 199 p1` | `0`; logs the 6 near-miss spaces (every non-hole point covered, no rank-6 form fits) | ~4 min | `logs/20261024-170915-check16-i2-p1.log` |
| `sixteen_check i2ctl apn_8bit.txt 15 3` / `... 12631 3` | `i2ctl: 1 of 1 verified APN by LUT` (known functions re-found) | ~10 min | `logs/*-check16-i2ctl3-*.log` |

### Fourteen (`i = 3`)

Search A, `sixteen.cpp` (holes Y, X and each of the 14 other graphs in turn; no orbit reduction), and cover.cpp as a
check on cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `sixteen i3 > fourteen.txt` | `holes 3 (14 rank-4 components)`; `spread sets (cliques) 1128864, semifield 19936, field 336 (expect 336)`; `solutions 0, short 0` | 15.3 h | `logs/20261024-165047-fourteen.log`, `logs/fourteen-output.txt` |
| `cover data/i3-reps.txt` | `done: covers 73, tasks 2044, nodes 2284642, solutions 0, full (dim < 8) 0` (search B's 73 representatives; uses B's reduction, so not independent of it) | 19 s | `logs/20261024-215424-i3-cover-reps.log`, `logs/i3-reps-output.txt` |

Search B, `fourteen_check` (`fourteen-check.md`). It writes `.prog` files to resume from; give each run a fresh one,
since an existing one makes it skip the configurations listed there.

| Command | Expected | Time | Log |
|---|---|---|---|
| `fourteen_check holes 11 holes11.txt` | `total 11-cliques (weighted count): 410906496.000; reduced leaves 863186; without extension: 0` (every partial spread of 14 solids extends) | 4 s | `logs/20261024-195125-check14-holes11.log` |
| `fourteen_check spreads spreads.txt` | `73 G-orbits; sum of orbit sizes 1128864 (weighted count 1128864.000); 1 orbit reps are fields`; writes `spreads.txt` (= `data/spreads.txt`) and `spreads.txt.all` (the 2,873 conjugation-reduced cliques) | 40 s | `logs/20261024-194724-check14-spreads.log`, `logs/20261024-195500-check14-spreads2.log` |
| `fourteen_check selftest spreads.txt` | `selftest: 73 configurations, 16694216 generated forms checked, 0 failures` | 40 s | `logs/20261024-200618-check14-selftest.log` |
| `fourteen_check search spreads.txt asc.prog asc` | `solutions: 0`; `search nodes by dim: 73 2044 57232 740398 1196814 64384 845 0 0`; no node with every off-D point covered | 23 s | `logs/20261024-195115-check14-asc.log` |
| `fourteen_check search spreads.txt desc.prog desc` | `solutions: 0` (descending point order) | 25 s | `logs/20261024-195837-check14-desc.log` |
| `fourteen_check search spreads.txt rand7.prog rand 7` | `solutions: 0` (random point order) | 26 s | `logs/20261024-200618-check14-rand.log` |
| `fourteen_check search spreads.txt.all noS3.prog asc` | `solutions: 0` (all 2,873 conjugation-reduced cliques, no S3 reduction) | 14 min | `logs/20261024-200446-check14-nos3.log` |
| `python src/fourteen_check_levels.py spreads.txt 0 4` | `[1, 28, 784, 10234, 15917]`, the C++ counts for spread 0 | 12 s | `logs/20261024-200202-check14-levels0.log` |
| `fourteen_check ctl spreads.txt apn_8bit.txt 12631 0 0` | `ctl: 1 spaces found ...; known space among them: YES` | seconds | `logs/20261024-195459-check14-ctl-12631.log` |
| `fourteen_check ctlbatch spreads.txt apn_8bit.txt 2 100000000 3 ctl-bl22.txt` | `1440 controls; re-found 1440; NOT found 0; ...; other spaces found 0` | 6 min | `logs/20261024-195740-check14-ctlbatch-bl22.log`, `logs/check14-ctl-bl22.txt` |
| `... apn_8bit.txt 3 100000000 24 ctl-wide.txt randspread` | `4128 controls; re-found 4128`; `spread reps hit: 68 of 73` | 23 min | `logs/20261024-200812-check14-ctl-wide.log`, `logs/check14-ctl-wide.txt` |
| `... blp22_8bit.txt 2 100000000 6 ctl-blp22.txt` | `222 controls; re-found 222` | 1 min | `logs/20261024-200446-check14-ctl-blp22.log` |
| `... new_apns.txt 2 100000000 1 ctl-new.txt` | `7199 controls; re-found 7199` | 62 min | `logs/20261024-200619-check14-ctl-new.log` |
| `... yuperrin_qam.txt 2 100000000 6 ctl-yuperrin.txt` | `24 controls; re-found 24` | 1 min | `logs/20261024-200618-check14-ctl-qam2.log` |
| `... yuwangli_qam.txt 2 100000000 6 ctl-yuwangli.txt` | `78 controls; re-found 72; NOT found 0; ...; setup failures 6` (exit code 3): the 6 failures are the list's line 3, which isn't a quadratic APN function as read | 1 min | `logs/20261024-203006-check14-ctl-qam1b.log`, `logs/check14-ctl-yuwangli.txt` |

### Thirteen (`i = 4`)

Search A, cover.cpp: the hole sets of a non-extendable 13-partial spread are Kurz's three (*Divisible codes*,
arXiv 2112.11763, §9; search B proves this independently), and the extendable case has holes Y, X, graph(I),
graph(T), one T per conjugacy class:

| Command | Expected | Time | Log |
|---|---|---|---|
| `python src/i4_holes.py data/kurz13.txt` | for the three sets: `60 distinct`, `points per hyperplane: {28: 195, 36: 60}`, `solids inside the complement: 1264` / `1280` / `1280`, `exact covers ...: 552` / `10440` / `12`, `solids inside the holes: 0`, `lines inside the holes: 200; partitions of the holes into 20 lines: 1000000` (the count stops at that cap) | 12 min, nearly all in the capped count | `logs/rerun/i4-holes.log` |
| `python src/cover_prep.py max max.txt data/kurz13.txt` | `config 0: 552 covers`, `config 1: 10440 covers`, `config 2: 12 covers`, `lines 11004` | 1 min | `logs/rerun/prep-max.log` |
| `cover max.txt` | `done: covers 11004, tasks 308112, nodes 374898296, solutions 0, full (dim < 8) 0, dim 8 non-APN 0` | 59 min | `logs/20261024-201233-i4-max.log`, `logs/i4-max-output.txt` |
| `python src/cover_prep.py ext ext.txt` | spread sets containing 0, I, T: `15852` / `2746` / `1760` / `2746` / `2746` (minimal polynomials 111, 10011, 10101, 11001, 11111); `lines 25850`. Weighted by class size (112, 1344, 1680, 1344, 1344) that is 15,804,096 = 14 × 1,128,864 | 22 min | `logs/20261024-194956-i4-ext-prep.log` |
| `cover ext.txt` | `done: covers 25850, tasks 723800, nodes 768769492, solutions 0, full (dim < 8) 0, dim 8 non-APN 0` | 1 h 57 min | `logs/20261024-201233-i4-ext.log`, `logs/i4-ext-output.txt` |
| `cover data/ctl-12631-k4.txt` (`k3`, `k2`) | `SOLUTION ctl-12631 dim 8 ranks 8:186 6:65 4:4 ...`; `solutions 1`: the known function (Beierle–Leander line 12631) re-found as the only space, seeded with one rank-4 form and 4, 3 or 2 forms completing a basis, so cover.cpp finds the other rank-4 forms and none, one or two rank-6 levels itself. The inputs come from `python src/cover_prep.py ctl OUT apn_8bit.txt 12631 1 KEEP` | 1 s | `logs/rerun/cover-ctl-k*.log` |

Search B, `thirteen_check` (`thirteen-check.md`; fresh `.prog` files as for fourteen):

| Command | Expected | Time | Log |
|---|---|---|---|
| `thirteen_check holes holes10.txt` | `reduced leaves 2358714 (non-extendable 141165, partial 0)`; `10-cliques, weighted: extendable 1129992864, non-extendable 62558496, common nbhd nonempty but no 4-clique 0` | 8 s | `logs/20261024-220843-check13-holes.log` |
| `thirteen_check configs holes10.txt data/kurz13.txt configs.txt` | Kurz's sets with automorphism groups of order 576, 14400, 96 and `552`, `10440`, `12` covers; `configurations up to GL(8,2): 159 extendable, 31 non-extendable`; both counts `MATCH`; writes `configs.txt` (= `data/configs.txt`) | 11 min | `logs/20261024-221011-check13-configs.log` |
| `thirteen_check selftest configs.txt` | `selftest: 190 configurations, 33400360 generated forms checked, 0 failures` | 2.5 min | `logs/20261024-222133-check13-selftest.log` |
| `thirteen_check search configs.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 190 5320 148960 1972468 3468491 783269 4006 0 0`; 4 nodes with every off-D point covered, all at dim 6 (configurations 45 and 154), dead at the first hole point | 111 s | `logs/20261024-222133-check13-asc.log` |
| `thirteen_check search configs.txt desc.prog desc` | `solutions: 0`; the same 4 near-miss nodes | 2 min | `logs/20261024-222637-check13-desc.log` |
| `thirteen_check search configs.txt rand7.prog rand 7` | `solutions: 0`; the same 4 | 2 min | `logs/20261024-222638-check13-rand.log` |
| `thirteen_check searchsp data/spreads.txt sp14.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 14 fourth holes (1,022 extendable configurations, without this program's classification) | 12 min | `logs/20261024-222638-check13-sp14.log` |
| `thirteen_check covers data/kurz13.txt covers.txt`, then `thirteen_check searchmem covers.txt covers.prog asc` | `Kurz set 1: 552 covers` (`10440`, `12`); then `solutions: 0` over all 11,004 Kurz covers unreduced, nodes by dim `11004 308112 8627136 115373106 207621076 42791038 177828 0 0`: 374,909,300 in all, cover.cpp's 374,898,296 plus the 11,004 roots | 53 min | `logs/20261024-224107-check13-covers.summary.log` |
| `python src/thirteen_check_levels.py configs.txt 0 3` / `... 170 3` | `[1, 28, 784, 10573]` / `[1, 28, 784, 10288]`, the C++ counts for configurations 0 and 170 | seconds | `logs/rerun/f13-levels0.log`, `logs/rerun/f13-levels170.log` |
| `thirteen_check ctlbatch configs.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440`, `nonext: 1440 controls; re-found 1440`, `other spaces found 0` | 20 min | `logs/20261024-222757-check13-ctl-bl22.log`, `logs/check13-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222, `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again, exit code 3), on both paths | 7 min | `logs/20261024-222758-check13-ctl-small.log`, `logs/check13-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `7199 controls; re-found 7199` on both paths | 69 min | `logs/20261024-222803-check13-ctl-new.log` |
| `thirteen_check ctlreps configs.txt apn_8bit.txt 1 40 50 122 ... 158 45` | `40 of 40 re-found, 0 other spaces` for each of the 19 configurations the batches missed | 2 min | `logs/20261025-002003-check13-ctlreps.log` |
| `thirteen_check holesk 10`, `holesk 9` | sizing for `i = 5` (`thirteen-check.md` §8) | seconds | `logs/20261024-231121-check13-size10.log` |

### Twelve (`i = 5`)

Kurz's table gives only counts for twelve solids (42 hole sets, 5 maximal), so both searches classify the 12-partial
spreads themselves. Search A, `pspread.cpp` then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread enum 12 p12.txt` | `w=0: reduced 100783, weighted 49896000, orbits 12, orbit sum 49896000 MATCH`; `w=1: reduced 1040086, weighted 502840800, orbits 84, ... MATCH`; `w=5: reduced 4128284, weighted 2258825184, orbits 280, ... MATCH`; `orbits 376`; writes `p12.txt` (= `data/p12.txt`) | 72 min | `logs/20261025-082434-pspread12.log`, `logs/pspread12-output.txt` |
| `pspread enum 13 p13.txt`, `enum 14 p14.txt` | K = 13: `w=0: reduced 141165, weighted 62558496, orbits 31`, `w=4: reduced 2217549, weighted 1129992864, orbits 159` (thirteen-check's numbers); K = 14: `w=3: reduced 863186, weighted 410906496, orbits 73` (fourteen-check's); writes `data/p13.txt`, `data/p14.txt` | 60 min, 38 min | `logs/20261025-082433-pspread1{3,4}.log`, `logs/pspread1{3,4}-output.txt` |
| `pspread enum 17 p17.txt` | `orbits 8` (the 8 spreads of PG(7,2), the translation planes of order 16; the Desarguesian one with `\|Aut P\|` 244,800), 2,873 reduced, orbit sum 1,128,864 (the spread sets of `sixteen`) | 21 s | `logs/rerun/pspread-p17.log`, `logs/rerun/p17.txt` (the original run's log wasn't kept; this is the re-run's) |
| `pspread selftest 12 data/p12.txt` (`13 data/p13.txt`, `14 data/p14.txt`) | `selftest: 376 configurations, 0 mismatches` (190, 73): key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 18 s; 31 s for 13 and 14 | `logs/rerun/pspread-selftest1{2,3,4}.log` |
| `python src/holes_iso12.py data/p12.txt` | `5 hole classes`, `\|Aut H\|` 2, 96, 32, 32, 16 with `sum \|Aut H\|/\|Aut P\|` (covers of the complement) 1, 84, 4, 4, 4 | 1 s | `logs/rerun/holes-iso12.log` |
| `cover data/p12.txt` | `done: covers 376, tasks 10528, nodes 12303628, solutions 0, full (dim < 8) 0, dim 8 non-APN 0`; `nodes by span dimension 1..8: 10528 294784 4002712 7479748 509108 6747 1 0`; `near misses ...: 0` | 49 s | `logs/20261025-095119-cover12b.log`, `logs/twelve-output.txt` |
| `cover data/i3-reps.txt` | `nodes 2284642, solutions 0` as before; `nodes by span dimension 1..8: 2044 57232 749035 1279756 195498 1077 0 0`, `near misses ...: 0` (a regression check of this build) | 4 s | same job log, `logs/i3-reps-cover3-output.txt` |

Search B, `twelve_check` (`twelve-check.md`; fresh `.prog` files as for fourteen):

| Command | Expected | Time | Log |
|---|---|---|---|
| `twelve_check parta orbits12.txt maxc.txt` | `reduced 9-cliques 5269153 (spr 4128284, ext13 1040086, max 100783)`; `9-cliques, weighted: all 2811561984; extend to a spread 2258825184, extend only to maximal 13s 502840800, maximal 49896000`; `sum of #5-cliques in Q: 2259985728 (C(14,9) * 1128864 = 2259985728)`; `sum of #q giving a maximal 13: 625584960 (10 * 62558496 ...)`, `... extendable 13: 11299928640 (10 * 1129992864 ...)`, `... neither: 0`; `orbits: 376 (spr 280, ext13 84, max 12); inconsistent invariants 0; per-orbit identity ... fails 0`; writes `orbits12.txt` (= `data/orbits12.txt`) and `maxc.txt` (the 100,783 maximal reduced cliques) | 26 min | `logs/20261025-085834-check12-parta.log`, `logs/check12-parta-output.txt` |
| `twelve_check configs data/orbits12.txt data/spreads.txt data/configs.txt configs12.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 42, containing no solid (maximal partial spreads): 5`; `hole sets with a cover mismatch: 0`; spread source `6643 configurations ... spread orbits reached 280 of 280`; 13-source `403 configurations; on no listed orbit 0; on a maximal orbit 0; by type: spr:37 ext13:366; ext13 orbits reached 84 of 84`; writes `configs12.txt` (= `data/configs12.txt`) | 7 s | `logs/20261025-092902-check12-configs.log`, `logs/check12-configs-output.txt` |
| `twelve_check selftest data/configs12.txt` | `lines in D over the configurations: 325:376`; `brute force: 28 forms ...`; `selftest: 376 configurations, 108452736 generated forms checked, 0 failures` | ~1 min | `logs/20261025-093045-check12-cross.log`, `logs/check12-selftest-output.txt` |
| `twelve_check search data/configs12.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 376 10528 294784 4002613 7761928 724602 8217 0 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 0 0 0` (no span ever holds all 12 rank-4 forms) | 20 s | `logs/20261025-092923-check12-asc.log`, `logs/check12-asc-output.txt` |
| `twelve_check search data/configs12.txt desc.prog desc` | `solutions: 0`; nodes by dim `376 10528 294784 3970567 7369230 380832 3848 0 0`; none with every off-D point covered | 20 s | `logs/20261025-093045-check12-cross.log`, `logs/check12-desc-output.txt` |
| `twelve_check search data/configs12.txt rand7.prog rand 7` | `solutions: 0`; nodes by dim `376 10528 294784 3979633 7462326 442121 5272 0 0`; none with every off-D point covered | 19 s | same job log, `logs/check12-rand7-output.txt` |
| `twelve_check searchsp data/spreads.txt sp.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 91 pairs of further holes (6,643 configurations, without this program's classification); nodes by dim `6643 186004 5208112 70266988 133848444 9697414 117642 0 0`; none with every off-D point covered | 8 min | same job log, `logs/check12-sp-output.summary.txt` |
| `twelve_check search13 data/configs.txt n13.prog asc` | `solutions: 0` over thirteen_check's 31 maximal 13-configurations × 13 members removed (403); nodes by dim `403 11284 315952 4288357 8486696 1618431 9352 0 0`; none with every off-D point covered | 32 s | same job log, `logs/check12-n13-output.txt` |
| `twelve_check searchmax maxc.txt maxc.prog asc` | all 100,783 maximal reduced cliques from `parta`, without the canonical form: `solutions: 0`; `search (incl. resumed) nodes by dim: 100783 2821924 79013872 1065082714 1985285822 174884697 929739 0 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 0 0 0` | ~2.6 h of running time between pauses (the program's own clock counts paused time too: 1,203 + 56,565 s, from its per-configuration lines, which the summary copy here leaves out) | `logs/20261025-094106-check12-maxc.log` (first 9,593 configurations; cancelled and resumed), `logs/20261025-105951-check12-maxc2.log` (the rest), `logs/check12-maxc-output.summary.txt` |
| `python src/twelve_check_levels.py data/configs12.txt 0 4` / `... 280 4` / `... 364 4` | `[1, 28, 784, 10804, 21176]` / `[1, 28, 784, 10711, 19628]` / `[1, 28, 784, 10594, 18535]`, the C++ counts for configurations 0 (spr), 280 (ext13) and 364 (max) | 12 s each | `logs/check12-levels-output.txt`, `logs/rerun/t12-levels*.log` |
| `twelve_check ctlbatch data/configs12.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0`, `reps hit: 210`; `nonext: 1440 controls; re-found 1440 ...`, `reps hit: 96` | 6 min | `logs/20261025-094106-check12-ctl.log`, `logs/check12-ctl-bl22-output.txt`, `logs/check12-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222, `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths | 1 min | `logs/20261025-094106-check12-ctl.log`, `logs/check12-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check12-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 241`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 96` | ~80 min | `logs/20261025-094106-check12-ctl.log` (the job ran all five lists; 89 min of running time between its pauses), `logs/check12-ctl-new-output.txt`, `logs/check12-ctl-new.txt` |
| `twelve_check ctlreps data/configs12.txt apn_8bit.txt 1 40 56 60 ... 279` (the 66 spr configurations the batches above missed; the full list is in `logs/check12-ctlreps-output.txt`) | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the 66 | ~20 min | `logs/20261025-105951-check12-ctl2.log`, `logs/check12-ctlreps-output.txt` |
| `twelve_check ctlbatch data/configs12.txt apn_8bit.txt 2 100000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1415 controls; re-found 1415`, `nonext: 1419 controls; re-found 1419`, `other spaces found 0`, `S0 already the whole space 46`; tries at dimension 7 only, all at points of D (78,758,848 and 78,286,656). A first trial with `2 10 1` (10 functions per path): 10 of 10 re-found on each | ~2 min | same job log, `logs/check12-ctl-ind-bl22-output.txt`, `logs/check12-ctl-ind-bl22.txt`; the trial `logs/check12-ctl-test-nond.txt` |
| `twelve_check maprowen data/configs12.txt data/p12.txt` | `maprowen: 376 of Rowen's configurations; unmatched 0; distinct listed configurations hit 376 of 376; hit twice 0; \|Aut\| disagreements 0; type/w disagreements 0` (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/check12-maprowen-output.txt` |

## Re-run for this folder

On 2026-10-25, from a fresh MSVC build of `src/` with the Build line above, on an idle machine (output in
`logs/rerun/`):

- `sixteen_check selftest`: every line `ok`.
- `fourteen_check selftest data/spreads.txt`, `search data/spreads.txt asc14.prog asc`, `holes 11`, `spreads`: the
  logged results, the search's node and try counts identical at every dimension (5 s); the regenerated `spreads.txt`,
  `spreads.txt.all` and `holes11.txt` are identical to the originals (line endings aside).
- `python src/fourteen_check_levels.py data/spreads.txt 0 4`: `[1, 28, 784, 10234, 15917]`.
- `cover data/i3-reps.txt`: `nodes 2284642, solutions 0` (3 s); `src/spreads_to_cover.py` reproduces `data/i3-reps.txt`.
- `thirteen_check selftest data/configs.txt` and `search data/configs.txt asc13.prog asc`: the logged results,
  identical counts and the same 4 near-miss nodes (29 s).
- `thirteen_check holes`: the logged counts. This build also splits the non-extendable hole invariants by plane counts
  (three lines, 25,081,056 + 33,153,120 + 4,324,320 = 62,558,496), which is why its `holes10.txt` is larger.
- `thirteen_check configs`: log identical apart from the output path, `configs.txt` identical to `data/configs.txt`
  (2 min). `thirteen_check covers`: `covers.txt` identical to search B's.
- `python src/thirteen_check_levels.py data/configs.txt 0 3` / `170 3`: `[1, 28, 784, 10573]` / `[1, 28, 784, 10288]`.
- `python src/i4_holes.py data/kurz13.txt`: as in the table; `python src/cover_prep.py max`: `max.txt` identical to the
  one cover.cpp's `i4-max` run read.
- `python src/cover_prep.py ctl ... apn_8bit.txt 12631 1 KEEP` for KEEP 4, 3, 2: identical to `data/ctl-12631-k*.txt`;
  `cover` on each: the known space re-found as the only solution (4, 5 and 13 nodes).

Not re-run: the hour-long searches, the control batches, `cover_prep.py ext` (22 min) and the 15-hour `sixteen i3`.

For twelve, later the same day:

- `python src/holes_iso12.py data/p12.txt`: the 5 hole classes with `|Aut H|` and cover sums as in the table.
- `python src/twelve_check_levels.py data/configs12.txt 0 4` / `280 4` / `364 4`: the logged counts.
- On 2026-10-26, from a fresh build of `src/pspread.cpp`, `cover.cpp` and `twelve_check.cpp` (logs in `logs/rerun/`,
  times in `times-twelve.txt`): `cover data/p12.txt` and `cover data/i3-reps.txt` (`cover-i3-reps-2.log`), the three
  `pspread selftest`s, `pspread enum 17`, and `twelve_check` `configs`, `selftest`, `search` asc / desc / rand 7,
  `search13`, `searchsp` (summary lines only), `maprowen`, `ctlreps` on representatives 56 and 279, and the `nond`
  trial `2 10 1`: every result as in the tables, node counts identical at every dimension; the regenerated
  `configs12.txt` identical to `data/configs12.txt`, the `nond` trial's output identical apart from the list's path.

Not re-run: the three `pspread enum` runs for K = 12, 13, 14 (an hour or more each), `twelve_check parta` (26 min),
`searchmax` and the control batches.

## Data

`data/`:
- `kurz13.txt`: Kurz's three 8 × 60 generator matrices (arXiv 2112.11763, §9, after Table 9.3), cut from the paper's
  text, not retyped. Search B's `kurz.txt` was the same file.
- `spreads.txt`: search B's 73 `i = 3` configurations (output of `fourteen_check spreads`).
- `configs.txt`: search B's 190 `i = 4` configurations (output of `thirteen_check configs`).
- `i3-reps.txt`: `spreads.txt` as cover.cpp input (`python src/spreads_to_cover.py data/spreads.txt`).
- `ctl-12631-k4.txt`, `-k3`, `-k2`: cover.cpp's positive controls (`python src/cover_prep.py ctl ... 12631 1 KEEP`).
- `p12.txt`: search A's 376 `i = 5` configurations (output of `pspread enum 12`; cover.cpp input, tags
  `p12-<n>-a<|Aut P|>-w<w>`). `p13.txt`, `p14.txt`: `pspread enum 13` / `14` (190 and 73 configurations, the checks
  against the cases above; inputs to `pspread selftest`).
- `orbits12.txt`: search B's 376 orbits with their invariants (output of `twelve_check parta`; input to `configs`).
- `configs12.txt`: search B's 376 `i = 5` configurations (output of `twelve_check configs`).

Not included, all regenerated by the commands above: `spreads.txt.all` (280 KB, 40 s), `holes11.txt`, `holes10.txt`
(16 to 22 MB depending on the build, 8 s), search B's `covers.txt` (1.9 MB, seconds), cover.cpp's `max.txt` (1.9 MB)
and `ext.txt` (4.5 MB, 22 min), `twelve_check parta`'s `maxc.txt` (5.6 MB, 26 min), and the `.prog` resume files.

The known-function lists used by the controls aren't included either: `apn_8bit.txt` (Beierle and Leander, Zenodo
4738942), `blp22_8bit.txt` (Beierle, Leander and Perrin, Zenodo 5821018), `new_apns.txt` (Beierle, Langevin, Leander,
Polujan and Rasoolzadeh, Zenodo 16752428, 3.4 GB), and `yuwangli_qam.txt` / `yuperrin_qam.txt` (the QAM lists of Yu,
Wang and Li with Edel and Pott's, and of Yu and Perrin, from sboxU's `apnDB`, one look-up table per line).

## Left out of `logs/`

- The full log of the `searchmem` run over the 11,004 Kurz covers (1.5 MB, one line per cover): the copy here keeps
  its header and summary lines.
- Runs superseded or repeated by the ones above: `fourteen_check search ... 0 0` on spread 0 alone (the same numbers
  as spread 0 in the full run); a first `yuwangli_qam.txt` control batch that crashed on the list's line 3 (an older
  build, before an APN filter for list lines was added; the batch above is its re-run); a cancelled cover.cpp control
  with only 1 completing form (KEEP 1); and `cover` on a 50-cover sample of `max.txt`.
- For twelve: `cover` on all 376 and on a 9-configuration sample with the build before the per-dimension counts were
  added (the same 12,303,628 nodes, 0 solutions; the run above is its repeat); two first submissions of the
  `pspread enum` runs that exited at once (a malformed command line); search B's per-configuration lines of the
  `searchsp` run (961 KB; the copy here keeps its header and summary); and the same for the `searchmax` run (100,783 lines,
  ~14 MB).
- Local paths in the logs were shortened to paths relative to the working tree.
