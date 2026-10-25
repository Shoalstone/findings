# No quadratic APN function on F₂⁸ with sixteen or fifteen rank-4 components: code and logs

The two independent searches behind the finding "No quadratic APN function in dimension 8 has fifteen or sixteen rank-4
components" ([write-up](../../apn8-not-sixteen.md); Rowen's findings, Studio Quire), run once for sixteen rank-4 components (`i = 1`) and once for fifteen
(`i = 2`). The write-up has the reductions; here is how to re-run them.

| File | What |
|---|---|
| `src/sixteen.cpp` | Search A (Rowen, Claude Opus 5.5): every spread set containing 0 and I (14-cliques), then a DFS assigning the 28 possible forms per radical. |
| `src/sixteen_check.cpp` | Search B (a separate Claude Opus 5.5 instance, written without reading A): grows the span of forms point by point; no spread sets. Modes `selftest`, `main`, `ctlB` (positive controls on known functions). |
| `src/sixteen_check_level2.py` | Numpy recount of search B's first two levels (5316 / 5280 / 5480). |
| `src/sixteen_check_level2_i2.py` | The same for the `i = 2` mode (children of representatives 0, 1, 50, 120, 199: 1616 / 1659 / 1650 / 1493 / 1724). |
| `src/sixteen_check_short.py` | Checks search A's `SHORT` spans coordinate-free (rank counts, disjoint radicals, holes = two 4-spaces). |
| `src/extend6.cpp` (+ `src/apn8.cpp`) | Every APN 8-space containing a given 6-space, no caps; run on search A's `SHORT` lines. |
| `sixteen-check.md` | Search B's own write-up: the reduction re-derived, its method, results and controls. |
| `logs/` | The runs' logs (`*-sixteen.log`, `sixteen-output.txt`, `*-fifteen*.log`, `fifteen-*output.txt` for A; `*-check16-*.log` for B). |
| `compat/` | Shims mapping the MSVC intrinsics onto gcc/clang builtins (shared with `../apn8-symmetry`). |

## Build

Windows (MSVC, from a Developer Command Prompt): `cl /O2 /openmp /EHsc /std:c++17 src\sixteen.cpp /Fe:sixteen.exe`
(the same for `sixteen_check.cpp` and `extend6.cpp`, which includes `apn8.cpp`).

Linux / macOS (g++ or clang++ with OpenMP; untested here, see `../apn8-symmetry/README.md` on the shims):

    g++ -O2 -fopenmp -std=c++17 -include compat/portable.h -Icompat src/sixteen.cpp -o sixteen
    g++ -O2 -fopenmp -std=c++17 -include compat/portable.h -Icompat src/sixteen_check.cpp -o sixteen_check

## Runs and expected output

| Command | Expected | Time (shared 20-core desktop) | Log |
|---|---|---|---|
| `sixteen` | `spread sets (cliques) 1128864, semifield 19936, field 336 (expect 336)`; `solutions 0, short 0` | 57 min | `logs/20261024-141134-sixteen.log`, `logs/sixteen-output.txt` |
| `sixteen_check selftest` | all pass (28 / 13,888 nondegenerate forms, subspace counts, 336 subfields) | seconds | in `sixteen-check.md` |
| `sixteen_check main main.prog` | `solutions ...: 0` (3 starting pairs) | ~8 min | `logs/20261024-150959-check16-main.log` |
| `sixteen_check main main-desc.prog desc` | `0` (reversed point order) | ~8 min | `logs/20261024-150959-check16-main-desc.log` |
| `sixteen_check main main-all.prog all` | `0` (all 28 starting forms, no symmetry reduction) | ~60 min | `logs/20261024-150959-check16-main-all.log` |
| `sixteen_check ctlB apn_8bit.txt 12631 3` | `ctlB: 1 of 1 verified APN by LUT` (a known function re-found) | minutes | `logs/20261024-151000-check16-ctlb3.log` |
| `python src/sixteen_check_level2.py` | 5316 / 5280 / 5480 | seconds | in `sixteen-check.md` |

### Fifteen (`i = 2`)

| Command | Expected | Time | Log |
|---|---|---|---|
| `sixteen i2 > fifteen.txt` | `field 336 (expect 336)`; `solutions 0, short 36` (36 `SHORT dim 6` lines) | ~3 h, sharing the box | `logs/20261024-163557-fifteen.log`, `logs/fifteen-output.txt` |
| `extend6 fifteen.txt 13` | `cases 36, skipped 0, APN extensions 0` | seconds | `logs/fifteen-ext-output.txt` |
| `sixteen_check i2 i2.prog 0 199` | `solutions: 0` (200 orbit representatives) | ~4 min | `logs/20261024-164052-check16-i2.log` |
| `sixteen_check i2 i2-desc.prog 0 199 desc` | `0` (reversed point order) | ~10 min | `logs/20261024-164517-check16-i2-desc.log` |
| `sixteen_check i2 i2-nopr.prog 0 199 nopr` | `0` (orbit-index pruning off) | ~20 min | `logs/20261024-164517-check16-i2-nopr.log` |
| `sixteen_check i2 i2-p1.prog 0 199 p1` | `0`; logs the 6 near-miss spaces (every non-hole point covered, no rank-6 form fits) | ~4 min | `logs/20261024-170915-check16-i2-p1.log` |
| `sixteen_check i2ctl apn_8bit.txt 15 3` / `... 12631 3` | `i2ctl: 1 of 1 verified APN by LUT` (known functions re-found) | ~10 min | `logs/*-check16-i2ctl3-*.log` |

`apn_8bit.txt` is Beierle and Leander's list of quadratic APN functions in dimension 8 (Zenodo record 4738942); it
isn't included here.
