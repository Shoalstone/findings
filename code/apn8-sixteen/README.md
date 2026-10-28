# Quadratic APN functions on F₂⁸ have at most nine rank-4 components: code and logs

The searches behind the finding "No quadratic APN function in dimension 8 has more than nine rank-4 components"
([write-up](../../apn8-not-sixteen.md); Rowen's findings, Studio Quire). Each case was settled by two independent
searches: sixteen rank-4 components (`i = 1`), fifteen (`i = 2`), fourteen (`i = 3`), thirteen (`i = 4`), twelve
(`i = 5`), eleven (`i = 6`) and ten (`i = 7`). Search A is Rowen's (Claude Opus 5.5); search B is by separate Claude
Opus 5.5 subagents (one for `i = 1, 2`, a fresh one each for `i = 3`, `i = 4`, `i = 5`, `i = 6` and `i = 7`), each
written without reading search A, with its own write-up here (for `i = 7`, Rowen tallied search B's results and wrote
its verdict, as that write-up says). The write-up has the reductions; here is how to re-run them.

| File | What |
|---|---|
| `src/sixteen.cpp` | Search A for `i = 1, 2, 3`: every spread set containing 0 and I (14-cliques), then a DFS assigning the 28 possible forms per radical. Modes: none (`i = 1`), `i2`, `i3`. |
| `src/extend6.cpp` (+ `src/apn8.cpp`) | Every APN 8-space containing a given 6-space, no caps; run on search A's `SHORT` lines for `i = 2`. |
| `src/cover.cpp` | Search A for `i = 4` to `i = 7`: every APN space whose rank-4 radicals are a given set of solids, rank-6 radicals in the holes; point branching, no symmetry used. Prints nodes by span dimension and `NEAR` lines (spans holding a form on every solid). |
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
| `src/pspread.cpp` | Search A's classification for `i = 5` and `i = 6`: partial spreads of K solids up to GL(8,2) (`enum K`: canonical form over frames, extension class w, orbit-sum check; `selftest K file`). Run for K = 12 and 11, and for K = 17, 14, 13 as checks against the cases above. Output is cover.cpp input. |
| `src/holes_iso12.py` | The hole sets of search A's 12 maximal configurations up to GL(8,2), by explicit isomorphism (5 classes, as Kurz's table says). |
| `src/twelve_check.cpp` | Search B for `i = 5`: every 9-clique canonicalised (`parta`: 376 configurations, extension counts), hole sets and outside sources (`configs`), search (`search`, `searchsp`, `search13`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and the comparison with search A's list (`maprowen`). |
| `src/twelve_check_levels.py` | Recount of the first levels of `twelve_check search` for one configuration, sharing no code with it. |
| `src/levels_all12.py` | Runs `twelve_check_levels.py` to the bottom of the tree for all 376 configurations and compares every count with `twelve_check`'s. |
| `src/eleven_check.cpp` | Search B for `i = 6`: every 8-clique canonicalised, one minimal class per run (`parta`), the classes merged with per-orbit identities and the link to `twelve_check`'s 376 (`merge`: 817 configurations), hole sets, exact covers and outside sources (`configs`), search (`search`, `search12`, `search13`, `searchsp`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and, for comparison only, `orbitcov` (which configuration each listed clique lies on) and `mapa` (search A's list mapped onto search B's). |
| `src/eleven_check_levels.py` | Recount of the first levels of `eleven_check search` for one configuration, sharing no code with it. |
| `src/pspread3.cpp` | Search A's classification for `i = 7`: `pspread.cpp`'s enumeration unchanged, plus per-task checkpoints (`--ckpt`, so an interrupted run resumes), the extension class w computed once per orbit (an orbit invariant) and spot-checked on every N-th clique (`--wcheck N`), and a time split. Run for K = 10, and for K = 12 and 14 as checks against `pspread.cpp`'s outputs. |
| `src/ten_check.cpp` | Search B for `i = 7`, built on `eleven_check.cpp`: every 7-clique canonicalised, a minimal class per run, split into slices (`parta`), the slices merged with per-orbit identities and the link to `eleven_check`'s 817 (`merge`: 1,777 configurations), hole sets, exact covers and outside sources (`configs`), search (`search`, `search11`, `search12`, `search13`, `searchsp`, `searchmax`), `selftest`, controls (`ctlbatch`, with `nond` for the rank-6 branching; `ctlreps`), and, for comparison only, `orbitcov` and `mapa`. |
| `src/ten_check_levels.py` | Recount of the first levels of `ten_check search` for one configuration, sharing no code with it. |
| `src/ten_check_tally.py` | Tallies over `ten_check`'s outputs: configurations no control reached (`reps`), and the `orbitcov` counts over all of Part A's maximal cliques against Part A's own (`maxall`). |
| `src/near_misses10.py` | Written for this folder: search A's and search B's near misses for `i = 7` compared configuration by configuration, from the logs here. |
| `sixteen-check.md`, `fourteen-check.md`, `thirteen-check.md`, `twelve-check.md`, `eleven-check.md`, `ten-check.md` | Search B's write-ups for `i = 1, 2`, `i = 3`, `i = 4`, `i = 5`, `i = 6` and `i = 7`: the reduction re-derived, method, results, controls, and (at the end) the comparison with search A. |
| `data/` | Small inputs (below). |
| `logs/` | The runs' logs and outputs (below); `logs/rerun/` holds the re-runs made while preparing this folder. |
| `compat/` | Shims mapping the MSVC intrinsics onto gcc/clang builtins (from `../apn8-symmetry`, with `__popcnt64`, `_tzcnt_u64`, `_BitScanReverse` and `_pdep_u32` added for the sources here). |

## Build

Windows (MSVC, from a Developer Command Prompt): `cl /O2 /openmp /EHsc /std:c++17 /arch:AVX2 src\sixteen.cpp /Fe:sixteen.exe`,
the same for `extend6.cpp` (which includes `apn8.cpp`), `cover.cpp`, `sixteen_check.cpp`, `fourteen_check.cpp`,
`thirteen_check.cpp`, `pspread.cpp`, `twelve_check.cpp`, `eleven_check.cpp`, `pspread3.cpp` and `ten_check.cpp`. The
runs' binaries were built with `/arch:AVX2`; the sources build without it too.

Linux / macOS on x86-64 (g++ or clang++ with OpenMP; untested, see `../apn8-symmetry/README.md` on the shims; the
shims added here were checked against the MSVC intrinsics on 20 million random inputs):

    g++ -O2 -fopenmp -std=c++17 -include compat/portable.h -Icompat src/sixteen.cpp -o sixteen

and the same for the other ten. The Python scripts need Python 3 (`sixteen_check_level2*.py` also numpy).

`src/sixteen.cpp` is the current source. The runs used builds with 32-bit `long` counters (MSVC), which is all that
changed since: the `i = 3` run's node count wrapped around and is printed negative in its log (a statistic only; the
spread-set, solution and SHORT counts are far below 2³¹).

`src/cover.cpp` is the source of the build that ran twelve, eleven and ten (`cover3`); it adds the per-dimension node
counts and the `NEAR` lines to the one that ran thirteen, and those only report: on search B's 73 `i = 3`
representatives it visits the same 2,284,642 nodes as before. `src/pspread.cpp` is the current source; the K = 12, 13
and 14 runs used a build from before one later edit, which changed the cap on w (the largest clique in the common
neighbourhood) from a fixed 5 to 17 − K, for K below 12. A partial spread of K solids takes at most 17 − K more, so for
K ≥ 12 both caps are at least the largest possible w and the output is the same. The K = 11 run used the later build:
there w goes up to 6 (an 11 that extends to a spread). `src/pspread3.cpp` ran K = 10; on K = 12 and 14 its output is
byte for byte `data/p12.txt` and `data/p14.txt` (Ten, below).

Which source built which binary was checked by rebuilding with the Build line's flags (`/arch:AVX2` included) and
comparing disassembly (2026-10-26, and for eleven 2026-10-27): `src/cover.cpp` gives the build that ran twelve and
eleven, instruction for instruction; `src/pspread.cpp` the later build (the one that ran K = 11), and `pspread.cpp` as
it stood before that edit the build that ran K = 12, 13 and 14; `src/twelve_check.cpp` the last of search B's three
twelve builds (the earlier two differ only by the added modes, as its write-up says); `src/eleven_check.cpp` the last
build of search B for eleven, which added the comparison-only modes `orbitcov` and `mapa` after the classification,
searches and controls had run (the source's history shows no other change since its first commit, made while Part A
was running). MSVC's code generation depends on the source file's name, so the comparison only matches when the file
is compiled under its own name. For `eleven_check.cpp` it also needed the original's compile line, from the source's
own folder with the object file named `eleven_check_` (`/Fo<dir>\eleven_check_`); compiled as `src\eleven_check.cpp`,
or to `eleven_check.obj`, the functions come out in a different order (the same instructions, counted by opcode).

For ten (2026-10-28) the binaries' object files were compared, function by function with the order ignored and the
compiler's hashed lambda names set aside (they vary from build to build): `src/pspread3.cpp` gives the build that ran
K = 10, 12 and 14, every function identical. `ten_check.cpp` changed once after its first commit (made as Part A was
submitted): the comparison-only `mapa` mode, which came from `eleven_check` reading 11 members, was changed to read
10 for the comparison with search A. The first commit's source gives the build that ran everything else (Part A to the
controls), and `src/ten_check.cpp` the rebuild that ran `mapa`, every function identical in both cases; the
executable that ran the rest was overwritten by that rebuild, its object file kept. In each comparison the only
difference is the size of a debug section, which records the build folder's path.

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
| `python src/levels_all12.py data/configs12.txt logs/check12-asc-output.txt OUT 8` | the same recount to dimension 8, which for twelve is search B's whole tree (it never branches at points of D here), for all 376 configurations: `levels recount: 376 of 376 configurations, 0 mismatches, 0 missing`, `total nodes (python): 12803048`, the sum of the C++ counts by dimension | ~1 h with 8 workers | `logs/check12-levels-all-output.txt` |
| `twelve_check ctlbatch data/configs12.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0`, `reps hit: 210`; `nonext: 1440 controls; re-found 1440 ...`, `reps hit: 96` | 6 min | `logs/20261025-094106-check12-ctl.log`, `logs/check12-ctl-bl22-output.txt`, `logs/check12-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222, `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths | 1 min | `logs/20261025-094106-check12-ctl.log`, `logs/check12-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check12-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 241`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 96` | ~80 min | `logs/20261025-094106-check12-ctl.log` (the job ran all five lists; 89 min of running time between its pauses), `logs/check12-ctl-new-output.txt`, `logs/check12-ctl-new.txt` |
| `twelve_check ctlreps data/configs12.txt apn_8bit.txt 1 40 56 60 ... 279` (the 66 spr configurations the batches above missed; the full list is in `logs/check12-ctlreps-output.txt`) | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the 66 | ~20 min | `logs/20261025-105951-check12-ctl2.log`, `logs/check12-ctlreps-output.txt` |
| `twelve_check ctlbatch data/configs12.txt apn_8bit.txt 2 100000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1415 controls; re-found 1415`, `nonext: 1419 controls; re-found 1419`, `other spaces found 0`, `S0 already the whole space 46`; tries at dimension 7 only, all at points of D (78,758,848 and 78,286,656). A first trial with `2 10 1` (10 functions per path): 10 of 10 re-found on each | ~2 min | same job log, `logs/check12-ctl-ind-bl22-output.txt`, `logs/check12-ctl-ind-bl22.txt`; the trial `logs/check12-ctl-test-nond.txt` |
| `twelve_check maprowen data/configs12.txt data/p12.txt` | `maprowen: 376 of Rowen's configurations; unmatched 0; distinct listed configurations hit 376 of 376; hit twice 0; \|Aut\| disagreements 0; type/w disagreements 0` (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/check12-maprowen-output.txt` |

### Eleven (`i = 6`)

Kurz's table again gives only counts for eleven solids (303 hole sets, 58 maximal), so both searches classify the
11-partial spreads themselves. Search A, `pspread.cpp` (the same program as for twelve) then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread enum 11 p11.txt` | `w=0: reduced 710368, weighted 422009280, orbits 67, orbit sum 422009280 MATCH`; `w=1: reduced 619505, weighted 342619200, orbits 43, ... MATCH`; `w=2: reduced 3324960, weighted 1760081400, orbits 254, ... MATCH`; `w=6: reduced 5728509, weighted 3377565576, orbits 453, ... MATCH`; `orbits 817`; writes `p11.txt` (= `data/p11.txt`) | ~2.1 h of running time between pauses, mostly on the E-cores (the program's clock, 65,072 s, counts paused time too) | `logs/20261025-095302-pspread11.log`, `logs/pspread11-output.txt` |
| `pspread selftest 11 data/p11.txt` | `selftest: 817 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 25 s | `logs/rerun/pspread-selftest11.log` (run for this folder) |
| `cover data/p11.txt` | `done: covers 817, tasks 22876, nodes 30505655, solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`; `nodes by span dimension 1..8: 22876 640528 8993626 19289387 1533912 25288 38 0`; `near misses ...: 58`: the `NEAR` lines, spans holding all 11 rank-4 forms with holes uncovered (22 at dim 6, 36 at dim 7, in 11 configurations), none of which reaches dim 8 | 89 s | `logs/20261026-092941-cover11.log`, `logs/eleven-output.txt` |

Search B, `eleven_check` (`eleven-check.md`; fresh `.prog` files as for fourteen). `parta` runs one minimal class per
job; the five outputs are `data/parta11-c*.txt`, its lists of maximal cliques aren't included (`Data` below):

| Command | Expected | Time | Log |
|---|---|---|---|
| `eleven_check parta C parta-cC.txt maxc-cC.txt` for C = 0 … 4 | `class C ...: N reduced 8-cliques, M maximal`, with N = `5917664`, `3003006`, `736536`, `364396`, `361740` (10,383,342 in all) and M = `599168`, `86195`, `17880`, `7125`, `0` (710,368); each file ends in `END`; `parta-cC.txt` = `data/parta11-cC.txt` | 8,611 / 4,950 / 815 / 832 / 195 s by the program's clock (about 2 h, 83 min, 14 min, 14 min and 3 min of running time, on the E-cores) | `logs/2026…-check11-pa*.log`, `logs/check11-parta-c*-output.txt` |
| `eleven_check merge orbits.txt data/configs12.txt data/parta11-c0.txt ... data/parta11-c4.txt` | `sum over 8-cliques of \|Q\| = 25304057856; 9 * (weighted 9-cliques, twelve-check) = 25304057856  OK`; `orbits 817; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 0:67 1:43 2:254 6:453`; `orbits by type: spr 453, ext13 254, ext12 43, max 67`; `... C(14,8) * 1128864 = 3389978592  OK`; `per 12-orbit O: ... fails for 0 of 376`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits11.txt`) | 6 s | `logs/check11-merge-output.txt` |
| `eleven_check configs orbits.txt data/configs12.txt data/spreads.txt data/configs.txt configs11.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 303, containing no solid (maximal partial spreads): 58`; `orbits whose type disagrees with 'hole set contains a solid': 0`; `hole sets with a cover mismatch: 0; covers in all 2134`; outside sources `twelve-check's 376 x 12: 4512 configurations; on no listed orbit 0; ... spr 453/453 ext13 254/254 ext12 43/43`, `fourteen-check's spread reps x 364: 26572 ...; on no listed orbit 0; ... spr 453/453`, `thirteen-check's maximal 13s x 78: 2418 ...; on no listed orbit 0; ... ext13 254/254`; writes `configs11.txt` (= `data/configs11.txt`; spr 0..452, ext13 453..706, ext12 707..749, max 750..816) | 19 s | `logs/check11-configs-output.txt` |
| `eleven_check selftest data/configs11.txt` | `lines in D over the configurations: 510:817`; `brute force: 28 forms with radical exactly R_5 of the last config (expect 28)`; `selftest: 817 configurations, 368561396 generated forms checked, 0 failures` | 4 min | `logs/check11-selftest-output.txt` |
| `eleven_check search data/configs11.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 817 22876 640528 8993779 19647718 1938778 24532 38 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 22 37 0`; `of which at points of D: 0 0 0 0 0 0 3930304 3013696 0` (the rank-6 forms tried in the holes); `distinct p1 spans this run: 59` | 18 min on the E-cores | `logs/20261026-210508-check11-asc.log`, `logs/check11-asc-output.txt`, near-miss spans `logs/check11-asc-p1.txt` |
| `eleven_check search data/configs11.txt desc.prog desc` | `solutions: 0`; nodes by dim `817 22876 640528 8935159 19095956 1335017 18128 40 0`; every off-D point covered `22` / `36` at dims 6 / 7; `distinct p1 spans this run: 58` | 18 min | `logs/20261026-210508-check11-desc.log`, `logs/check11-desc-output.txt`, `logs/check11-desc-p1.txt` |
| `eleven_check search data/configs11.txt rand7.prog rand 7` | `solutions: 0`; nodes by dim `817 22876 640528 8944444 19234261 1513135 18258 38 0`; `22` / `37`; `distinct p1 spans this run: 59` | 18 min | `logs/20261026-210509-check11-rand7.log`, `logs/check11-rand7-output.txt`, `logs/check11-rand7-p1.txt` |
| `eleven_check search12 data/configs12.txt s12.prog asc` | `solutions: 0` over twelve_check's 376 configurations × 12 members removed (4,512: every 11 that extends, without this program's classification); nodes by dim `4512 126336 3537408 49655433 109538964 10777400 165744 111 0`; every off-D point covered `97` / `103` at dims 6 / 7 | 76 min | `logs/20261026-175950-check11-s12.log`, `logs/check11-s12-output.summary.txt` |
| `eleven_check search13 data/configs.txt s13.prog asc` | `solutions: 0` over thirteen_check's 31 maximal 13-configurations × 78 pairs removed (2,418); nodes by dim `2418 67704 1895712 26629650 60146644 8913624 96436 72 0`; `22` / `61` | 50 min | `logs/20261026-175950-check11-s13.log`, `logs/check11-s13-output.summary.txt` |
| `eleven_check searchsp data/spreads.txt ssp.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 364 triples removed (26,572: every 11 that extends to a spread); nodes by dim `26572 744016 20832448 291080164 635317280 54887313 860221 495 0`; `263` / `450` | ~4.5 h of running time, mostly on the E-cores (evicted once and resumed from `ssp.prog`) | `logs/20261026-175950-check11-ssp.log`, `logs/check11-ssp-output.summary.txt` |
| `eleven_check searchmax maxc-s50.txt maxs50.prog asc`, where `maxc-s50.txt` is lines 1, 51, 101, … of `maxc-c0.txt` … `maxc-c3.txt` concatenated (14,208 maximal cliques) | `solutions: 0`; nodes by dim `14208 397824 11139072 154792269 321538731 14093008 17710 0 0`; no span holds all 11 rank-4 forms | ~2.1 h | `logs/20261026-210924-check11-maxs50.log`, `logs/check11-maxs50-output.summary.txt` |
| `eleven_check orbitcov orbits.txt configs11.txt maxc-s50.txt` | `orbitcov: 14208 cliques, 0 on no listed orbit, 0 on a non-max orbit; 66 orbits hit of 67 max orbits` (orbit 798 missed) | seconds | `logs/check11-maxs50-orbits-output.txt` |
| `eleven_check orbitcov orbits.txt configs11.txt maxc-all.txt maxc-o798.txt 798` (all four `maxc` files concatenated), then `eleven_check searchmax maxc-o798.txt max798.prog asc` | `orbitcov: 710368 cliques, 0 on no listed orbit, 0 on a non-max orbit; 67 orbits hit of 67 max orbits`, the per-orbit tallies equal to `n` in `orbits.txt`; `wrote 85 cliques` (= `data/maxc11-o798.txt`); then `solutions: 0`, nodes by dim `85 2380 66640 934228 2128630 407037 5956 0 0` | 3.5 min | `logs/20261026-235204-check11-max798.log`, `logs/check11-maxall-orbits-output.txt`, `logs/check11-max798-output.txt` |
| `python src/eleven_check_levels.py data/configs11.txt 0 4` / `453 4` / `707 4` / `750 4` | `[1, 28, 784, 11152, 24772]` / `[1, 28, 784, 11206, 23367]` / `[1, 28, 784, 10939, 24007]` / `[1, 28, 784, 11074, 22957]`, the C++ counts for one configuration of each type (spr, ext13, ext12, max) | 1 min for the four | `logs/20261026-210759-check11-levels.log`, `logs/check11-levels-output.txt` |
| `eleven_check ctlbatch data/configs11.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 319`; `nonext: 1440 controls; re-found 1440 ...`, `reps hit: 364` | ~1 h for the four lists | `logs/20261026-210844-check11-ctl.log`, `logs/check11-ctl-bl22-output.txt`, `logs/check11-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222, `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths (the job exits 3 because of the last) | (same job) | same job log, `logs/check11-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check11-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 383`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 364` | ~2.3 h | `logs/20261026-210845-check11-ctlnew.log`, `logs/check11-ctl-new-output.txt`, `logs/check11-ctl-new.txt` |
| `eleven_check ctlbatch data/configs11.txt apn_8bit.txt 2 100000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1438 controls; re-found 1438`, `nonext: 1440 controls; re-found 1440`, `other spaces found 0`, `S0 already the whole space 2`; every control starts at dimension 7 and finds its form among the 13,888 on its line (98,965,888 and 99,785,280 tries, all at points of D) | 5 min | `logs/20261026-210845-check11-ctlnond.log`, `logs/check11-ctl-ind-bl22-output.txt`, `logs/check11-ctl-ind-bl22.txt` |
| `eleven_check ctlreps data/configs11.txt apn_8bit.txt 1 40 93 112 ... 452` (the 64 spr configurations the batches above missed; the full list is in the job log's command) | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the 64 | 8 min | `logs/20261026-235204-check11-ctlreps.log`, `logs/check11-ctlreps-output.txt` |
| `eleven_check mapa data/configs11.txt data/p11.txt` | `mapa: 817 of A's configurations, 0 on none of ours, 817 of our 817 configurations hit, 0 hit more than once`, with type and `canon aut` agreeing with search A's tags on every line (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/check11-mapa-output.txt` |

The near misses agree between the two searches configuration by configuration: mapped through `mapa`, search A's
`NEAR` lines fall in search B's configurations 248, 299, 314, 316, 429, 435, 452, 541, 558, 609 and 649, 6, 6, 6, 27,
1, 1, 3, 3, 1, 1 and 3 of them, the counts search B's `desc` run has there (its per-configuration `p1` counts); `asc`
and `rand 7` have a second span at 558. Which partial spans a point-branching tree visits depends on its point order.

### Ten (`i = 7`)

Kurz's table gives counts for ten solids too (1,336 hole sets, 317 maximal), so both searches classify the 10-partial
spreads themselves once more. Many of these jobs ran side by side, so their times overlap; times are running time
between the runner's pauses. Search A, `pspread3.cpp` (`pspread.cpp`'s enumeration with checkpoints) then cover.cpp:

| Command | Expected | Time | Log |
|---|---|---|---|
| `pspread3 enum 12 p12.txt --ckpt p12.ckpt --wcheck 200`, `pspread3 enum 14 p14.txt --ckpt p14.ckpt --wcheck 50` | the checks of the new build: `w spot checks: 12411, mismatches 0`, `orbits 376`; `w spot checks: 17221, mismatches 0`, `orbits 73`; the outputs are byte for byte `data/p12.txt` and `data/p14.txt`. Resumed from a complete checkpoint the K = 14 run writes the same file at once, and from one holding part of the run (428 tasks and a torn last line, which it drops) the same file again | ~4.6 h for K = 12, beside the K = 10 run; 18 min for K = 14 | `logs/20261026-180053-p3-k12.log`, `logs/20261026-170942-p3-k14.log`, `logs/20261026-172744-p3-k14-resume.log`, `logs/20261026-172831-p3-k14-half.log`, `logs/pspread3-p12-output.txt`, `logs/pspread3-p14*-output.txt` |
| `pspread3 enum 10 p10.txt --ckpt p10.ckpt --wcheck 1000` | `w=0: reduced 6042938, weighted 3974342400, orbits 320, orbit sum 3974342400 MATCH`; `w=1: reduced 4255438, weighted 2774499840, orbits 274, ... MATCH`; `w=2: reduced 1707940, weighted 1056585600, orbits 111, ... MATCH`; `w=3: reduced 5835506, weighted 3414902400, orbits 460, ... MATCH`; `w=7: reduced 5986006, weighted 3815466624, orbits 612, ... MATCH`; `orbits 1777`; `w spot checks: 17563, mismatches 0`; writes `p10.txt` (= `data/p10.txt`, tags `p10-<n>-a<\|Aut P\|>-w<w>`) | ~5 h: 2.2 h on the E-cores before the runner restarted, then 2.8 h resumed from the checkpoint (114 tasks), mostly on an open box (the program's clock for that run, 10,610 s, counts paused time too) | `logs/20261026-180053-p3-k10.log`, `logs/pspread10-output.txt` |
| `pspread3 selftest 10 data/p10.txt` | `selftest: 1777 configurations, 0 mismatches`: key and `\|Aut P\|` unchanged after a random GL(8,2) map and shuffle | 21 s | `logs/rerun/pspread3-selftest10.log` (run for this folder) |
| `cover data/p10.txt` | `done: covers 1777, tasks 49756, nodes 79061670, solutions 0, full (dim < 8) 0, dim 8 non-APN 0, bad seeds 0`; `nodes by span dimension 1..8: 49756 1393168 20447257 52320443 4776667 74207 172 0`; `near misses ...: 370`: the `NEAR` lines, 337 at dim 6 (in 27 configurations) and 33 at dim 7 (in 15), 39 configurations in all, none reaching dim 8 | 202 s | `logs/20261026-235746-i7-cover10.log`, `logs/ten-output.txt` |

Search B, `ten_check` (`ten-check.md`; fresh `.prog` files as for fourteen). `parta` runs one slice of a minimal class
per job; the 17 outputs are `data/parta10-c*-s*.txt`, its lists of maximal cliques aren't included (`Data` below):

| Command | Expected | Time | Log |
|---|---|---|---|
| `ten_check parta C parta-cC-sS.txt maxc-cC-sS.txt S N` for class C = 0 … 4 in N = 8, 4, 2, 2, 1 slices | per slice `class C slice S/N ...: R reduced 7-cliques, M maximal`; by class R = 16,069,200, 5,277,651, 1,423,134, 621,223, 436,620 (23,827,828 in all) and M = 5,049,808, 785,215, 160,785, 46,290, 840 (6,042,938); each file ends in `END`; `parta-cC-sS.txt` = `data/parta10-cC-sS.txt` | 44 min from the first start to the last end on an open box, eleven jobs at a time; by the program's clock 1,407 to 2,548 s per slice of class 0, 596 to 1,747 s of class 1, at most 311 s for the rest | `logs/2026…-check10-pa*-s*.log`, `logs/check10-parta-c*-s*-output.txt` |
| `ten_check merge orbits.txt data/configs11.txt data/parta10-c0-s0.txt ... data/parta10-c4-s0.txt` (the 17 in class and slice order) | `eleven-check's configurations: 817; sum 990*20160/\|Aut\| = 5902275456 (its Part A: 5902275456)`; `reduced 7-cliques 23827828, of them maximal (Q empty) 6042938; weighted 7-cliques 15035796864`; `sum over 7-cliques of \|Q\| = 47218203648; 8 * (weighted 8-cliques, eleven-check) = 47218203648  OK`; `orbits 1777; clique-level inconsistencies (\|Aut\| or \|Q\|) 0 + 0; per-orbit identity ... fails 0; ... OK`; `largest clique in Q ...: 0:320 1:274 2:111 3:460 7:612`; `orbits by type: spr 612, ext13 460, ext12 111, ext11 274, max 320`; `type / Q / n7 inconsistencies 0; P + q on no eleven-check orbit 0`; `... C(14,7) * 1128864 = 3874261248  OK`; `per 11-orbit O: ... fails for 0 of 817`; `merge: all identities hold`; writes `orbits.txt` (= `data/orbits10.txt`) | under a minute with `configs` | `logs/20261027-045312-check10-merge-configs.log`, `logs/check10-merge-output.txt` |
| `ten_check configs orbits.txt data/configs11.txt data/configs12.txt data/spreads.txt data/configs.txt configs10.txt` | `canonical form recomputed from a random GL(8,2) image of every orbit: 0 mismatches`; `hole sets up to GL(8,2): 1336, containing no solid (maximal partial spreads): 317`; `orbits whose type disagrees with 'hole set contains a solid': 0`; `hole sets with a cover mismatch: 0; covers in all 2294`; outside sources `eleven-check's 817 x 11: 8987 configurations; on no listed orbit 0; ... spr 612/612 ext13 460/460 ext12 111/111 ext11 274/274 max 0/320`, `fourteen-check's spread reps x 1001: 73073 ...; on no listed orbit 0; ... spr 612/612`, `thirteen-check's maximal 13s x 286: 8866 ...; on no listed orbit 0; ... spr 252/612 ext13 460/460`, `twelve-check's 376 x 66: 24816 ...; on no listed orbit 0; ... spr 612/612 ext13 460/460 ext12 111/111`; writes `configs10.txt` (= `data/configs10.txt`; spr 0..611, ext13 612..1071, ext12 1072..1182, ext11 1183..1456, max 1457..1776) | (same job) | same job log, `logs/check10-configs-output.txt` |
| `ten_check selftest data/configs10.txt` | `selftest: ND4 28 (28), ND6 13888 (13888), GL4 20160 (20160), V5824 5824`; `lines in D over the configurations: 770:1777`; `brute force: 28 forms with radical exactly R_5 of the last config (expect 28)`; `selftest: 1777 configurations, 1198198680 generated forms checked, 0 failures` | 41 min, beside the searches below | `logs/20261027-121908-check10-selftest.log`, `logs/check10-selftest-output.txt` |
| `ten_check search data/configs10.txt asc.prog asc` | `solutions: 0`; `search (incl. resumed) nodes by dim: 1777 49756 1393168 20447263 52683154 5454191 70909 181 0`; `nodes with every off-D point covered, by dim: 0 0 0 0 0 0 337 165 0`; `of which at points of D: 0 0 0 0 0 0 74550784 19082112 0` (the rank-6 forms tried in the holes); `distinct p1 spans this run: 502` | 24 min, the three orders at once | `logs/20261027-121909-check10-asc.log`, `logs/check10-asc-output.txt`, near-miss spans `logs/check10-asc-p1.txt` |
| `ten_check search data/configs10.txt desc.prog desc` | `solutions: 0`; nodes by dim `1777 49756 1393168 20333833 51884404 4463901 60838 131 0`; every off-D point covered `337` / `120` at dims 6 / 7; at points of D `73661952` / `13360256`; `distinct p1 spans this run: 457` | 25 min | `logs/20261027-121909-check10-desc.log`, `logs/check10-desc-output.txt`, `logs/check10-desc-p1.txt` |
| `ten_check search data/configs10.txt rand7.prog rand 7` | `solutions: 0`; nodes by dim `1777 49756 1393168 20336938 51929881 4484972 61922 169 0`; `337` / `150`; at points of D `73717504` / `16221184`; `distinct p1 spans this run: 487` | 24 min | `logs/20261027-121909-check10-rand7.log`, `logs/check10-rand7-output.txt`, `logs/check10-rand7-p1.txt` |
| `ten_check search11 data/configs11.txt s11.prog asc` | `solutions: 0` over eleven_check's 817 configurations × 11 members removed (8,987: every 10 that extends, without this program's classification); nodes by dim `8987 251636 7045808 103434188 270500229 31444429 568714 864 0`; every off-D point covered `2507` / `768` at dims 6 / 7 | 53 min | `logs/20261027-004853-check10-s11.log`, `logs/check10-s11-output.summary.txt` |
| `ten_check search13 data/configs.txt s13.prog asc` | `solutions: 0` over thirteen_check's 31 maximal 13-configurations × 286 triples removed (8,866); nodes by dim `8866 248248 6950944 102193519 273604875 39199199 677779 656 0`; `661` / `543` | 47 min | `logs/20261027-004853-check10-s13.log`, `logs/check10-s13-output.summary.txt` |
| `ten_check searchsp data/spreads.txt ssp.prog asc` | `solutions: 0` over fourteen_check's 73 representatives × 1001 quadruples removed (73,073: every 10 inside a spread); nodes by dim `73073 2046044 57289232 838595573 2191569950 246581805 5099404 5591 0`; `28477` / `4409` | ~6.5 h of running time (resumed from `ssp.prog` across many holds) | `logs/20261027-004853-check10-ssp.log`, `logs/check10-ssp-output.summary.txt` |
| `ten_check orbitcov orbits.txt configs10.txt data/maxc10-s400.txt`, then `ten_check searchmax data/maxc10-s400.txt maxs400.prog asc`, where `maxc10-s400.txt` is lines 1, 401, 801, … of `maxc-c0-s0.txt` … `maxc-c4-s0.txt` concatenated in class and slice order (15,108 maximal cliques) | `orbitcov: 15108 cliques, 0 on no listed orbit, 0 on a non-max orbit; 320 orbits hit of 320 max orbits`; then `solutions: 0`, nodes by dim `15108 423024 11844672 172205757 428754848 25136558 4469 0 0`, no span holds all 10 rank-4 forms (`searchmax` labels these cliques `ext11` in its per-type lines, a label only, `ten-check.md` §5) | 99 min for the two | `logs/20261027-121919-check10-maxs400.log`, `logs/check10-maxs400-orbits-output.txt`, `logs/check10-maxs400-output.summary.txt` |
| `ten_check orbitcov orbits.txt configs10.txt maxc-cC-sS.txt` for each of the 17 files, appended to `maxall-orbits.out`; then `python src/ten_check_tally.py maxall DIR` on a folder holding `orbits.txt` and `maxall-orbits.out` | each `orbitcov: N cliques, 0 on no listed orbit, 0 on a non-max orbit`; then `orbitcov runs 17 (expect 17); cliques off the max orbits 0; max orbits 320, cliques 6042938 (Part A 6042938); orbits whose tally differs from Part A's n: 0` (so every maximal clique lies on a listed maximal configuration, and all 320 have cliques) | 46 min | `logs/20261027-121924-check10-maxcov.log`, `logs/check10-maxall-orbits-output.txt` |
| `python src/ten_check_levels.py data/configs10.txt 0 4` / `612 4` / `1072 4` / `1183 4` / `1457 4` | `[1, 28, 784, 11578, 31023]` / `[1, 28, 784, 11680, 31277]` / `[1, 28, 784, 11443, 28745]` / `[1, 28, 784, 11542, 28800]` / `[1, 28, 784, 11467, 29638]`, the C++ counts for one configuration of each type (spr, ext13, ext12, ext11, max) | 6 min for the five | `logs/20261027-121909-check10-levels.log`, `logs/check10-levels-output.txt` |
| `ten_check ctlbatch data/configs10.txt apn_8bit.txt 2 100000 3 ctl-bl22.txt` | `ext: 1440 controls; re-found 1440; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 410`; `nonext: 1440 controls; re-found 1440 ...`, `reps hit: 976` | 35 min for the four lists | `logs/20261027-121918-check10-ctl.log`, `logs/check10-ctl-bl22-output.txt`, `logs/check10-ctl-bl22.txt` |
| the same with `blp22_8bit.txt`, `yuperrin_qam.txt`, `yuwangli_qam.txt` (6 choices) | `222` of 222 (reps hit 160 / 193), `24` of 24, `72` of 78 with 6 setup failures (the list's line 3 again), on both paths (the job exits 3 because of the last) | (same job) | same job log, `logs/check10-ctl-{blp22,yuperrin,yuwangli}-output.txt`, `logs/check10-ctl-*.txt` |
| the same with `new_apns.txt` (1 choice) | `ext: 7199 controls; re-found 7199; NOT found 0; setup failures 0; other spaces found 0 (0 not APN)`, `reps hit: 513`; `nonext: 7199 controls; re-found 7199 ...`, `reps hit: 1163` | 1.7 h | `logs/20261027-121918-check10-ctlnew.log`, `logs/check10-ctl-new-output.txt`, `logs/check10-ctl-new.txt` |
| `ten_check ctlbatch data/configs10.txt apn_8bit.txt 2 100000 3 ctl-ind-bl22.txt nond` | the rank-6 branching exercised: `ext: 1440 controls; re-found 1440`, `nonext: 1440 controls; re-found 1440`, `other spaces found 0`; `started from S0 of dimension: 7:1440` on each path, every control finding its form among the 13,888 on its line (124,811,456 and 123,047,680 tries, all at dimension 7) | 3 min | `logs/20261027-121918-check10-ctlnond.log`, `logs/check10-ctl-ind-bl22-output.txt`, `logs/check10-ctl-ind-bl22.txt` |
| `python src/ten_check_tally.py reps DIR` on a folder holding `configs.txt` (= `data/configs10.txt`) and the `ctl-*.txt` files above | `reps hit 1688 of 1777; not reached 89 {'spr': 89}` and the 89 | seconds | `logs/rerun/t10-tally.log` (run for this folder) |
| `ten_check ctlreps data/configs10.txt apn_8bit.txt 1 40 131 132 ... 611` (those 89; the full list is in the job log's command) | `rep N (ext): 40 of 40 re-found, 0 other spaces` for each of the 89 | 33 min | `logs/20261027-182638-check10-ctlreps.log`, `logs/check10-ctlreps-output.txt` |
| `ten_check mapa data/configs10.txt data/p10.txt` | `mapa: 1777 of A's configurations, 0 on none of ours, 1777 of our 1777 configurations hit, 0 hit more than once`, with type (w = 0, 1, 2, 3, 7 for max, ext11, ext12, ext13, spr) and `canon aut` agreeing with search A's tags on every line (search A's list mapped onto search B's, after both verdicts) | seconds | `logs/20261028-174550-check10-mapa.log`, `logs/check10-mapa-output.txt` |

The near misses, compared by `python src/near_misses10.py` from the logs above: search A's 337 at dimension 6 lie in 27
configurations (mapped through `mapa`), and each of search B's three orders has its dimension-6 near misses in exactly
those 27, the same number in each. Search A's 33 at dimension 7 lie in 15 configurations, where every order of search B
has the same number. Search B's other dimension-7 near misses (132 / 87 / 117 in `asc` / `desc` / `rand 7`) all contain
one of its dimension-6 near misses: a rank-6 form fitted at a hole, which cover.cpp doesn't log, since it prints `NEAR`
only at the step where the last rank-4 form joins (none of search A's 33 contains one of its dimension-6 near misses).
Those sit in configurations 91, 331, 352 and 377 in every order, and 376 in `asc` too; each already has
dimension-6 near misses on both sides. The script prints these counts; the spans themselves are in different
coordinates on the two sides, so it compares counts, not spans.

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

For eleven, on 2026-10-27, from a fresh build of `src/cover.cpp`, `pspread.cpp` and `eleven_check.cpp` with the Build
line (the first two identical to the binaries that ran, by disassembly), as one job beside a long CPU job of
Rowen's (logs in `logs/rerun/`, times in `times-eleven.txt`): `cover data/p11.txt` (`cover-p11.log`,
100 s), `pspread selftest 11 data/p11.txt`, `eleven_check parta 4` (2.3 min), `merge`, `configs`, `selftest`, `search`
asc / desc / rand 7 (2 min each), the four `eleven_check_levels.py` recounts, `mapa`, `searchmax data/maxc11-o798.txt`,
and `ctlreps` on representatives 93 and 452: every result as in the tables. Node, try and near-miss counts are
identical at every dimension and per configuration; the `NEAR` lines of `cover` and the `P1` spans of the three
searches are the same sets as the originals'; the regenerated `parta-c4.txt` (apart from its timing), `orbits.txt` and
`configs11.txt` are identical to `data/parta11-c4.txt`, `data/orbits11.txt` and `data/configs11.txt`, and the other
logs to the originals apart from timings and output paths.

Not re-run: `pspread enum 11` (about 2 h), `eleven_check parta` for classes 0 to 3 (up to 2 h), the cross-runs
`search12`, `search13`, `searchsp` and `searchmax` on the sample (50 min to 4.5 h), and the control batches.

For ten, on 2026-10-28, from a fresh build of `src/cover.cpp`, `pspread3.cpp` and `ten_check.cpp` with the Build
line, as two jobs on the E-cores beside other work (logs in `logs/rerun/`, times in `times-ten.txt`): `cover
data/p10.txt` (`cover-p10.log`, 5.5 min), `pspread3 selftest 10 data/p10.txt`, `ten_check parta` for class 4 and for
class 3 slice 1 (`parta 3 ... 1 2`), `merge`, `configs`, the five `ten_check_levels.py` recounts, `mapa`, the `nond`
control batch (`ctl-ind-bl22`, 41 s), `ctlreps` on representatives 131 and 611, the first 500 configurations of
`search11` and `searchmax data/maxc10-s400.txt` and the first 800 of `searchsp` (configurations 0 to 499 or 799, given
as the last two arguments), `search` asc and desc (7 and 10 min) and `selftest` (13 min); and `ten_check_tally.py
reps` and `maxall` and `near_misses10.py` on the files here (`t10-tally.log`): every result as in the tables. Node,
try and near-miss counts are identical at every dimension and per configuration (in the slices, configuration by
configuration); the `NEAR` lines of `cover` and the `P1` spans of the two searches are the same sets as the
originals'; the regenerated `parta-c4-s0.txt`, `parta-c3-s1.txt`, `orbits.txt` and `configs10.txt` are identical to
`data/parta10-c4-s0.txt`, `data/parta10-c3-s1.txt` (apart from the timing on their `CLASS` lines),
`data/orbits10.txt` and `data/configs10.txt`, and the other logs to the originals apart from timings and
output paths. The two regenerated lists of maximal cliques hold all 57 lines of `data/maxc10-s400.txt` that fall in
their stretch of the concatenation (lines 1, 401, 801, … counted over the 17 files in class and slice order).

Not re-run: `pspread3 enum 10` (about 5 h), `ten_check parta` for the other slices (44 min on an open box for all 17),
the rest of `search11`, `search13`, `searchsp` and `searchmax` (47 min to 6.5 h), `search` rand 7, and the control
batches other than `nond`.

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
- `p11.txt`: search A's 817 `i = 6` configurations (output of `pspread enum 11`; cover.cpp input, tags
  `p11-<n>-a<|Aut P|>-w<w>`).
- `parta11-c0.txt` … `parta11-c4.txt`: search B's Part A for `i = 6`, one file per minimal class (outputs of
  `eleven_check parta`; inputs to `merge`).
- `orbits11.txt`: search B's 817 orbits with their invariants (output of `eleven_check merge`; input to `configs`).
- `configs11.txt`: search B's 817 `i = 6` configurations (output of `eleven_check configs`).
- `maxc11-o798.txt`: the 85 maximal 8-cliques on orbit 798 (output of `eleven_check orbitcov`; input to the last
  `searchmax`).
- `p10.txt`: search A's 1,777 `i = 7` configurations (output of `pspread3 enum 10`; cover.cpp input, tags
  `p10-<n>-a<|Aut P|>-w<w>`).
- `parta10-c0-s0.txt` … `parta10-c4-s0.txt`: search B's Part A for `i = 7`, one file per slice of a minimal class (17;
  outputs of `ten_check parta`; inputs to `merge`).
- `orbits10.txt`: search B's 1,777 orbits with their invariants (output of `ten_check merge`; input to `configs`).
- `configs10.txt`: search B's 1,777 `i = 7` configurations (output of `ten_check configs`).
- `maxc10-s400.txt`: every 400th of Part A's 6,042,938 maximal 7-cliques (input to `searchmax` and `orbitcov`).
  Included because it can't be cut again exactly: `parta` writes its maximal cliques in the order its threads finish.

Not included, all regenerated by the commands above: `spreads.txt.all` (280 KB, 40 s), `holes11.txt`, `holes10.txt`
(16 to 22 MB depending on the build, 8 s), search B's `covers.txt` (1.9 MB, seconds), cover.cpp's `max.txt` (1.9 MB)
and `ext.txt` (4.5 MB, 22 min), `twelve_check parta`'s `maxc.txt` (5.6 MB, 26 min), `eleven_check parta`'s
`maxc-c0.txt` … `maxc-c3.txt` (36 MB together, `maxc-c4.txt` is empty; 3.9 h of running time on the E-cores for the
five classes, 2 h of it class 0) and the 1-in-50 sample `maxc-s50.txt` taken from them (710 KB; how in the table),
`ten_check parta`'s `maxc-c*-s*.txt` (about 290 MB together; 44 min on an open box for all 17 slices; deleted after
the `orbitcov` tally, so `maxc10-s400.txt` is in `data/`), and the `.prog` resume files (search B's near-miss spans,
its `P1` lines, are in `logs/check11-*-p1.txt` and `logs/check10-*-p1.txt`).

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
- For eleven: search B's per-configuration lines of the `search12`, `search13`, `searchsp` and `searchmax` (sample)
  runs (4,512, 2,418, 26,572 and 14,208 lines, 7 MB together; the `*.summary.txt` copies keep the headers and summary
  lines); a test build's trial run before the real ones: Part A cut to the first few branches of each class (445
orbits), with `configs` and a control batch on that partial list (whose `FAIL` lines are functions landing on
configurations the partial list doesn't hold).
- For ten: search B's per-configuration lines of the `search11`, `search13`, `searchsp` and `searchmax` (sample) runs
  (8,987, 8,866, 73,073 and 15,108 lines, 15 MB together; the `*.summary.txt` copies keep the headers and summary
  lines). `logs/pspread10-output.txt` has both of the K = 10 run's starts (the program appends to it), and its job log
  both commands.
- Local paths in the logs were shortened to paths relative to the working tree.
