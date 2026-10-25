# Quadratic APN functions on F₂⁸ have at most twelve rank-4 components: code and logs

The searches behind the finding "No quadratic APN function in dimension 8 has more than twelve rank-4 components"
([write-up](../../apn8-not-sixteen.md); Rowen's findings, Studio Quire). Each case was settled by two independent
searches: sixteen rank-4 components (`i = 1`), fifteen (`i = 2`), fourteen (`i = 3`) and thirteen (`i = 4`). Search A
is Rowen's (Claude Opus 5.5); search B is by separate Claude Opus 5.5 subagents (one for `i = 1, 2`, a fresh one each
for `i = 3` and `i = 4`), each written without reading search A, with its own write-up here. The write-up has the reductions; here is how to re-run them.

| File | What |
|---|---|
| `src/sixteen.cpp` | Search A for `i = 1, 2, 3`: every spread set containing 0 and I (14-cliques), then a DFS assigning the 28 possible forms per radical. Modes: none (`i = 1`), `i2`, `i3`. |
| `src/extend6.cpp` (+ `src/apn8.cpp`) | Every APN 8-space containing a given 6-space, no caps; run on search A's `SHORT` lines for `i = 2`. |
| `src/cover.cpp` | Search A for `i = 4`: every APN space whose rank-4 radicals are a given set of solids, rank-6 radicals in the holes; point branching, no symmetry used. |
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
| `sixteen-check.md`, `fourteen-check.md`, `thirteen-check.md` | Search B's write-ups for `i = 1, 2`, `i = 3` and `i = 4`: the reduction re-derived, method, results, controls, and (at the end) the comparison with search A. |
| `data/` | Small inputs (below). |
| `logs/` | The runs' logs and outputs (below); `logs/rerun/` holds the re-runs made while preparing this folder. |
| `compat/` | Shims mapping the MSVC intrinsics onto gcc/clang builtins (from `../apn8-symmetry`, with `__popcnt64`, `_tzcnt_u64`, `_BitScanReverse` and `_pdep_u32` added for the sources here). |

## Build

Windows (MSVC, from a Developer Command Prompt): `cl /O2 /openmp /EHsc /std:c++17 src\sixteen.cpp /Fe:sixteen.exe`,
the same for `extend6.cpp` (which includes `apn8.cpp`), `cover.cpp`, `sixteen_check.cpp`, `fourteen_check.cpp` and
`thirteen_check.cpp`.

Linux / macOS on x86-64 (g++ or clang++ with OpenMP; untested, see `../apn8-symmetry/README.md` on the shims; the
shims added here were checked against the MSVC intrinsics on 20 million random inputs):

    g++ -O2 -fopenmp -std=c++17 -include compat/portable.h -Icompat src/sixteen.cpp -o sixteen

and the same for the other five. The Python scripts need Python 3 (`sixteen_check_level2*.py` also numpy).

`src/sixteen.cpp` is the current source. The runs used builds with 32-bit `long` counters (MSVC), which is all that
changed since: the `i = 3` run's node count wrapped around and is printed negative in its log (a statistic only; the
spread-set, solution and SHORT counts are far below 2³¹).

Search B's write-ups mention successive builds (`fourteen_check2.exe` and so on): each is a build of the one source in
`src/`, whose later edits added control modes and left the search path alone (the write-ups say what changed). The
logs show which build each run used.

## Runs and expected output

Times are wall-clock on a shared 20-thread desktop (Core Ultra 7 265) under a job runner that gave way to other work,
often to the E-cores only, so an idle machine is faster. Logs are in `logs/`. The job logs (`2026…-<name>.log`) are
that runner's; its own lines (command, pauses, timing) are dropped from the published copies, so for Rowen's runs,
whose program output went to a file, the job log is nearly empty and the output is the `*-output.txt` beside it.

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

## Data

`data/`:
- `kurz13.txt`: Kurz's three 8 × 60 generator matrices (arXiv 2112.11763, §9, after Table 9.3), cut from the paper's
  text, not retyped. Search B's `kurz.txt` was the same file.
- `spreads.txt`: search B's 73 `i = 3` configurations (output of `fourteen_check spreads`).
- `configs.txt`: search B's 190 `i = 4` configurations (output of `thirteen_check configs`).
- `i3-reps.txt`: `spreads.txt` as cover.cpp input (`python src/spreads_to_cover.py data/spreads.txt`).
- `ctl-12631-k4.txt`, `-k3`, `-k2`: cover.cpp's positive controls (`python src/cover_prep.py ctl ... 12631 1 KEEP`).

Not included, all regenerated by the commands above: `spreads.txt.all` (280 KB, 40 s), `holes11.txt`, `holes10.txt`
(16 to 22 MB depending on the build, 8 s), search B's `covers.txt` (1.9 MB, seconds), cover.cpp's `max.txt` (1.9 MB)
and `ext.txt` (4.5 MB, 22 min), and the `.prog` resume files.

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
- Local paths in the logs were shortened to paths relative to the working tree.
