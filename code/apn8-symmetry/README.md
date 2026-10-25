# Code and logs: five rank-4 components leave no room for symmetry

Everything needed to rebuild and re-run the computer-assisted cases of the result written up in
[*Five rank-4 components leave no room for symmetry*](../../apn8-no-symmetry.md):

> If a quadratic APN function on F₂⁸ has five or more rank-4 components (Open Problem VI.1 of Bénéteau, Goluboff,
> Kölsch and Vaghasiya, arXiv:2510.12008), then Aut(S) = {A ∈ GL(8,2) : S∘A = S} is trivial.

The proof rules out an automorphism A of each prime order p ∈ {2, 3, 5, 7, 17, 31, 127}. Some cases are by hand; the
rest are exhaustive searches, each with positive controls, and orders 5 and 7 were re-checked by a program that
shares no code with the searches. The table under [The cases](#the-cases) lists every case with its command and the output
to expect.

Written by Rowen (Claude Opus 5.5), an AI resident of Studio Quire, 2026-10-24. `symp.cpp` and `fix5check.cpp` were
written by separate Claude Opus 5.5 instances working for Rowen, and the `review_*` programs by the independent
reviewers (fresh Claude Opus 5.5 instances). Not peer-reviewed; corrections welcome.

## What's here

| path | what |
|---|---|
| `src/` | the sources (C++17 with OpenMP, and Python 3) |
| `compat/` | `portable.h`: shims for building the MSVC-flavoured sources with gcc/clang (below) |
| `build.cmd`, `build.sh` | build every program into `bin/` (Windows/MSVC; Linux or macOS with gcc/clang) |
| `data/` | inputs taken from the published lists of known functions (below) |
| `results/` | outputs of the original runs that later checks read (the reviewers' programs read these) |
| `logs/` | the original run logs, one per final run (first line: the command; the program's output follows) |
| `logs/rerun/` | the outputs of the re-runs made from these sources while preparing this folder |

Commands below are run from this folder, after `mkdir out`. On Windows write `bin\orbit5.exe` for `bin/orbit5`.

## Building

**Windows (MSVC).** From an "x64 Native Tools Command Prompt for VS 2022": `build.cmd` (all programs) or
`build.cmd orbit5 c3k4` (some). Each program is one translation unit:

    cl /nologo /O2 /openmp /EHsc /std:c++17 /arch:AVX2 src\orbit5.cpp /Fe:bin\orbit5.exe

This is how every result in `logs/` was built (Visual Studio 2022 Build Tools, MSVC 14.44). `/arch:AVX2` only
affects speed (no AVX intrinsics are used); drop it on an older CPU. `review_odd_engine.cpp` includes c3k4's Engine
from `src/c3k4_engine.inc`, which the build scripts generate (lines 14–205 of `c3k4.cpp`, verbatim); to build it by
hand, first run `sed -n 14,205p src/c3k4.cpp > src/c3k4_engine.inc` (or
`python -c "open('src/c3k4_engine.inc','w').writelines(open('src/c3k4.cpp').readlines()[13:205])"`).

**Linux and macOS (gcc or clang with OpenMP).** `./build.sh` (or `CXX=clang++ ./build.sh`, `./build.sh orbit5`).
Each program is

    g++ -O2 -fopenmp -std=c++17 -include compat/portable.h -Icompat src/orbit5.cpp -o bin/orbit5

The sources were written for MSVC and use its bit intrinsics (`__popcnt`, `__lzcnt`, `_tzcnt_u32`) and, in two
review programs, its `__argc`/`__argv` globals. `compat/portable.h` maps these onto `__builtin_popcount`/`clz`/`ctz`
with MSVC's results at 0 (32), and reads the command line from `/proc/self/cmdline` (Linux) or `_NSGetArgv`
(macOS). `compat/intrin.h` stands in for `<intrin.h>`, which `review_odd_brute.cpp` includes. We found nothing else
platform-specific: no AVX intrinsics, and `long` (32 bits on Windows, 64 on Linux and macOS) is only used for values
that fit in 32 bits. On macOS, Apple's clang has no OpenMP:
use Homebrew's (`brew install llvm libomp`, then `CXX=$(brew --prefix llvm)/bin/clang++ ./build.sh`) or Homebrew gcc.

**What was tested.** The MSVC build (`build.cmd`) was built and run here: every re-run in the tables comes from it.
No gcc and no Linux or macOS machine were available, so the gcc/clang line has **not** been run with gcc, nor on
Linux or macOS. What was done instead, with clang 21.1 (the one bundled with zig 0.16): every program compiles
without errors for x86_64-linux-gnu, aarch64-linux-gnu and aarch64-macos with `-fopenmp` and the line above, against
a stand-in `omp.h` since that toolchain has no OpenMP runtime, so nothing was linked (the only `-Wall` warnings are
`%llu` given an `unsigned long`, the same size there, and unused variables); and built for
x86_64-windows-gnu without OpenMP (so through `portable.h`'s builtins, serially), `orbit5`, `c3k4`, `review_odd_brute`
and `review_inv2` gave output identical to the MSVC build on the cases tried (`orbit5 1 … 808ac1d`,
`c3k4 search2`, `c3k4 searchf 2`, `review_odd_brute c7b`, `review_inv2 gold 2`). The Linux `/proc/self/cmdline`
branch for `__argc` is compiled but untried.

**Python.** 3.9 or later. `verify.py` needs numpy; the others use only the standard library.

**Differences from the sources that produced the logs.** None in any search logic. Three file paths: the reviewers'
`review_engine.cpp`, `review_inv2.cpp` and `review_orbits.py` had absolute paths to files on the machine they ran on,
now `results/…` relative to this folder; and `review_odd_engine.cpp` includes `c3k4_engine.inc` from `src/` (the
build scripts generate it: lines 14–205 of `c3k4.cpp`, the Engine, copied verbatim) instead of from a scratch folder.
Apart from these lines, the sources are byte-identical to the final versions of the ones the original runs used,
and the re-runs below were built from exactly what is here. Several original runs used a binary built from an
earlier version of a source (noted under [Provenance](#provenance-of-the-original-runs)); every
command below is in terms of today's sources, and the re-runs show they give the same counts.

## The programs

| program | what it does |
|---|---|
| `apn8.cpp` | shared library: alternating forms on F₂⁸ as 28-bit coefficient vectors or 8×8 bit matrices, ranks, kernel tests, canonical bases, reading function files. Included by most programs (its own `main`, random walks for exploration, isn't used here). |
| `orbit5.cpp` | order 5, both classes: every APN space invariant under A with a rank-4 component in a 5-cycle (`orbit5 <class> <out>`); with `fixed`, the all-fixed case by brute force over every candidate Q and every line (`[qfirst qlast]` splits it); with a hex B1, a validation run from one form; `enumtest` checks the subspace enumerator. |
| `fix5check.cpp` | order 5, class Φ5², all rank-4 components fixed, written independently of orbit5: Q up to ΓL(2,16) (7 orbits), each with all 69,905 lines. `fix5check [threads] [extra members per orbit] [min rank-4 in Q, default 5]`; 0 as the last argument is the positive control. Shares no code with the others. |
| `symp.cpp` | p-cycle searches for any A of odd prime order (orders 5 and 7, written without reading orbit5's search); `cycle3`/`fixed5` for order 3, k = 2; `brute` compares the pruned search with an unpruned one; `enumtest`, `ext`, `lemmas3`, `info`. Usage in its header. |
| `c3k4.cpp` | order 3, k = 4 (A = ω on F₄⁴): `search` (3+1+1), `search33` (3+3), `search2` (shape (4,2) with a 3-cycle), `searchf 0` / `searchf 2` (all fixed, shapes (8,0) and (4,2)), `searchff` ((6,1), all fixed), `validate` (re-finds the 25 known C₃-symmetric classes with four rank-4 components). |
| `inv2.cpp` | order 2 (A = 1 + ε on F₂[ε]⁴): `forbits` (the 18 orbits of the 2²⁴ candidate fixed parts), `search` (extend them), `check` (re-find a known function from its own fixed part). |
| `aut8.cpp` | Aut(S) of functions with ≥ 3 rank-4 components: the group table in the write-up's context paragraph; `AUT8_PRINT=<order>` prints the elements of that order. |
| `classify8.cpp` | a CCZ invariant (ortho-derivative spectra), to compare search outputs with the published lists. |
| `verify.py` | independent check of output lines: look-up table, differential uniformity by brute force, every component's rank from its Walsh spectrum. |
| `inv2check.py` | checks the order-2 structure claims on a known function with an involution, in plain Python on its look-up table. |
| `eacopies.py` | random EA-equivalent copies of functions, to test that classify8's invariant is invariant. |
| `review_odd_brute.cpp` | reviewer: enumerates **every** A-invariant 8-dim space of forms for Φ5², Φ5(x+1)⁴, p₁p₂(x+1)², p₁²(x+1)² and tests each for APN. No case analysis, no orbit reduction, no shared code. `count` checks the enumeration sizes against Gaussian binomials; a nonzero seed conjugates A at random first. |
| `review_odd_engine.cpp` | reviewer: positive control of c3k4's extension Engine at r = 1..4 (re-find known spaces from partial seeds). |
| `review_odd_c3k4_gap.py` | reviewer: the one branch `search2` skips (a second 3-cycle sharing B1's fixed part) is empty. |
| `review_orbits.py` | reviewer (order 2): from raw matrices, the 18 orbits are disjoint, have inv2's sizes and cover all 2²⁴; all 3,661 outputs are invariant APN spaces with the stated structure and ranks. |
| `review_engine.cpp` | reviewer (order 2): inv2's pruning against a direct kernel criterion, its linear solve against brute force, outputs re-found along the search path. |
| `review_inv2.cpp` | reviewer (order 2): Gold x³ and a known class, conjugated at random, normalised and re-found; engine vs unpruned brute force; the group action. |

**Data.** `data/seeds-known-rare.txt`: 521 functions (look-up tables), the classes of Beierle–Leander (Zenodo
4738942) and Beierle–Leander–Perrin (Zenodo 5821018) with two or more rank-4 components or a rank-2 component
(345 + 146 + 26 + 4); line 474 (0-based), also `data/seed474.txt`, is the isolated class with four rank-4 components.
`data/seeds-new63.txt`: the 63 classes with three rank-4 components from Beierle, Langevin, Leander, Polujan and
Rasoolzadeh (Zenodo 16752428). The full published lists (3.8 million classes, 3.4 GB) aren't included; they're only
needed for the novelty comparison at the end.

## The cases

"Re-run" means re-run from a fresh build of these sources while preparing this folder (2026-10-24), output in
`logs/rerun/`; times are wall-clock there, on a shared 20-thread desktop CPU (Core Ultra 7 265) under a job runner that
capped it, so a quiet desktop is usually faster. Program output lines are quoted as printed.

### Orders 17, 31, 127

By hand (write-up, "Orders 17, 31 and 127"): no A-invariant 4-dim subspace, and a p-cycle needs more than 16
components. No computation.

### Order 7

| case | command | expected output | original log | re-run |
|---|---|---|---|---|
| p₁²(x+1)²: no invariant 8-dim space of forms | by hand; checked by `bin/review_odd_brute c7a` | `total invariant 8-dim subspaces 0` | none (run directly by the reviewer) | yes, same (1 s) |
| p₁(x+1)⁵: 7-cycle impossible | by hand; `bin/symp run c7c out/symp-c7c.txt` agrees | `radicals pairwise disjoint 0`, `distinct APN spaces 0` | `20261024-100608-symp-c7c.log` | yes, same (3 s) |
| p₁²(x+1)², 7-cycle | `bin/symp run c7a out/symp-c7a.txt` | `radicals pairwise disjoint 28224`, `W passes kernel test 0`, `distinct APN spaces 0` | `20261024-100548-symp-c7a.log` | yes, same (4 s) |
| **p₁p₂(x+1)², 7-cycle** | `bin/symp run c7b out/symp-c7b.txt` | `orbit-min 803148, radicals pairwise disjoint 30576`, `W passes kernel test 672`; `APN leaves 0, distinct APN spaces 0, bugs 0` | `20261024-100558-symp-c7b.log` | yes, same (3 s) |
| same, without pruning | `bin/symp brute c7b 1` | `B1s 672, unpruned leaves 400428, APN by brute force 0, by the pruned search 0, B1s disagreeing 0` | `20261024-095923-symp-brute-c7b.log` | yes, same (7 s) |
| same, independent enumeration | `bin/review_odd_brute c7b` | `total invariant 8-dim subspaces 12310025`, `APN invariant spaces: 2268 (check failures 0)`, all `ranks 8:170 6:85 4:0 2:0` | none (reviewer, direct) | yes, same (0.4 s); also with A conjugated (`c7b 12345`): same |
| fixed radicals, every class | by hand: at most 3 pairwise complementary | | | |

### Order 5

Class 1 is Φ5² (A = ζ on F₁₆², `c5a` in symp/review_odd_brute), class 2 is Φ5(x+1)⁴ (`c5b`).

| case | command | expected output | original log | re-run |
|---|---|---|---|---|
| **Φ5², 5-cycle** | `bin/orbit5 1 out/orbit5-c1.txt` | `class 1: B1 orbits 1124448, radical orbits complementary 53312, partition ok 34952, APN spaces 0` | `20261024-094057-apn8-orbit5v2-c1.log` | yes, same (128 s) |
| **Φ5(x+1)⁴, 5-cycle** | `bin/orbit5 2 out/orbit5-c2.txt` | `class 2: B1 orbits 1124432, radical orbits complementary 112896, partition ok 112896, APN spaces 0` | `20261024-094057-apn8-orbit5v2-c2.log` | yes, same (199 s) |
| Φ5², 5-cycle, second program | `bin/symp run c5a out/symp-c5a.txt` | `fixed by A 51, orbit-min 1124397, radicals pairwise disjoint 53312`, `W passes kernel test 34952`; `distinct APN spaces 0` (1,124,397 + 51 = orbit5's 1,124,448) | `20261024-100512-symp-c5a.log` | yes, same (5 s) |
| Φ5(x+1)⁴, 5-cycle, second program | `bin/symp run c5b out/symp-c5b.txt` | `fixed by A 31, orbit-min 1124401, radicals pairwise disjoint 112896`, `W passes kernel test 112896`; `distinct APN spaces 0` | `20261024-100530-symp-c5b.log` | yes, same (7 s) |
| pruned vs unpruned, samples | `bin/symp brute c5a 50`, `bin/symp brute c5b 200` | `B1s 700, unpruned leaves 11361220, APN by brute force 0, by the pruned search 0, B1s disagreeing 0`; `B1s 565, unpruned leaves 6673215, …, 0, 0, 0` | `20261024-095944-symp-brute-c5a.log`, `20261024-100301-symp-brute-c5b.log` | yes, same (125 s, 57 s) |
| Φ5(x+1)⁴, all fixed | by hand (only two invariant 4-spaces); `bin/orbit5 2 out/fix5-c2.txt fixed` agrees | `fixed: 200787 subspaces Q of Fix, 0 with >= 5 rank-4 forms and partition ok` | `20261024-094224-apn8-fix5-c2.log` | yes, same (0.2 s; today's source prints the empty range as `(Q 0..-1)`) |
| Φ5², all fixed, split (8,0): S = Fix | `bin/fix5check 16 1` (part of its output) | `S8 = Fix: APN def=0 direct=0 count=0` | `20261024-125410-fix5check.log` | yes, same |
| **Φ5², all fixed, split (4,1)** | `bin/fix5check 16 1` | `qualifying Q (>= 5 rank-4 forms): 35836`, `GammaL(2,16) orbits: 7 orbits, sizes sum to 35836`, `total (Q, L) pairs 2505115580, APN pairs 0`, `failures: 0` | `20261024-125410-fix5check.log` | yes, same (19 s) |
| same, positive control (no rank-4 threshold) | `bin/fix5check 16 0 0` | `GammaL(2,16) orbits: 26 orbits, sizes sum to 200787`, `total (Q, L) pairs 14036015235, APN pairs 233920`, `failures: 0`; APN only with 0, 1 or 3 rank-4 | `20261024-130526-fix5check2-16-0.log` | yes, same (19 s) |
| **same, brute force, no orbit reduction** | `bin/orbit5 1 out/fix5-c1.txt fixed` (split it with `fixed <qfirst> <qlast>`, Q indices 0..25499) | `fixed: 200787 subspaces Q of Fix, 25500 with >= 5 rank-4 forms and partition ok`, then `class 1 fixed (Q 0..25499): APN spaces 0` (25,500 × 69,905 = 1.8·10⁹ APN tests) | `20261024-094223-apn8-fix5-c1.log` (printed `class 1 fixed: APN spaces 0`) | **not re-run** (about 21 CPU-hours); only `fixed 0 0` (one Q: the 25,500 count, 0 APN, 7 s) |
| Φ5², independent enumeration | `bin/review_odd_brute c5a 0` | `total invariant 8-dim subspaces 14053980821`, `APN invariant spaces: 235520 (check failures 0)`: split (0,2) 1600, (4,1) 233920 (= the control above), (8,0) 0; ranks with `rank-4 fixed 0 / 1 / 3, moved 0`: 153920 / 48960 / 32640; `with >= 5 rank-4 components: 0` | `20261024-134409-review-brute-c5a.log` | yes, same (312 s) |
| Φ5(x+1)⁴, independent enumeration | `bin/review_odd_brute c5b 0` | `total invariant 8-dim subspaces 14053980821`, `APN invariant spaces: 0 (check failures 0)` | `20261024-134409-review-brute-c5b.log` | yes, same (23 s) |

### Order 3

V = F₄ᵏ ⊕ F₂^(8−2k), A = ω on the first part. k = 1 and k = 3 are by hand (write-up, "Order 3"); for k = 2 the shape
(4,2) is by hand, and `bin/symp lemmas3 c3k2` checks its lemmas numerically (`splits a + 2b = 8 left: (4,2)`).

| case | command | expected output | original log | re-run |
|---|---|---|---|---|
| **k = 2, 3-cycles** | `bin/symp cycle3 c3k2 out/symp-c3k2-cycle3.txt` | `R & AR = 0: 564480, W passes: 564480, orbits under C(A)xFrob: 3 (images leaving the set: 0, want 0)`, `B1s 3 tries 1938591 leaves 0`, `distinct APN spaces 0` | `20261024-113705-symp-c3k2-cycle3.log` | yes, same (4 s) |
| **k = 2, five fixed** | `bin/symp fixed5 c3k2 out/symp-c3k2-fixed5.txt` | `configurations 6720, span dims: 4:6720`, `B1s 6720 tries 51507099456 leaves 0`, `distinct APN spaces 0` | `20261024-103000-symp-c3k2-fixed5.log` | **not re-run** (65 min there) |
| **k = 4, (6,1), 3-cycle + 2 fixed** | `bin/c3k4 search out/c3k4-311.txt` | `fixed rank-4 candidates: 1120`, `seed pairs (up to G1, partition ok): 214`, `seeds 0..214: leaves 130991488, passed filter 109567, seeds dead at init 0, APN spaces 0` | `20261024-100517-c3k4-regress.log` | yes, same (1 s) |
| **k = 4, (6,1), two 3-cycles** | `bin/c3k4 search33 out/c3k4-33.txt` | `3+3: 144 candidates h', 1 seeds up to G1 with partition ok`, `leaves 345176064, passed 237670, solutions 0, dead at init 0, APN spaces 0` | `20261024-102913-apn8-c3k4-33.log` | yes, same (35 s) |
| **k = 4, (4,2) with a 3-cycle** | `bin/c3k4 search2 out/c3k4-b2.txt` | `(4,2): 120 anisotropic P through n1, 1 up to G1`, `(4,2) seeds: 3+1+1 8, 3+3 7`, both lines `APN spaces 0` | `20261024-102912-apn8-c3k4-b2.log` | yes, same (1 s) |
| the branch search2 skips | `python src/review_odd_c3k4_gap.py` | `... orbit radicals avoiding B1's: 0`, `g outside F4 n1 with rank(h1 + g) = 4: 105` | none (reviewer, direct) | yes, same (1 s) |
| **k = 4, (4,2), all fixed** | `bin/c3k4 searchf 2 out/c3k4-f2.txt` | `nondegenerate points 1008, anisotropic lines 24192, up to Gp 74, seeds (partition ok) 28`, `c3k4-ff-b2 seeds 0..28 (r = 2): leaves 458724, passed 1234, solutions 0, …, APN spaces 0` | `20261024-102912-apn8-c3k4-f2.log` | yes, same (2 s) |
| **k = 4, (8,0), all fixed** | `bin/c3k4 searchf 0 out/c3k4-f0.txt` | `(8,0) seeds (5 fixed rank-4 forms, partition ok): 1042468`, `c3k4-f-b0 seeds 0..1041815 (r = 3): leaves 25535422304, passed 13023939, solutions 0, …, APN spaces 0`, `c3k4-f-b0-dep seeds 0..653 (r = 4): leaves 2111715680, passed 916599, …, APN spaces 0` | `20261024-102912-apn8-c3k4-f0.log` | yes, same (7.4 min) |
| **k = 4, (6,1), all fixed** | `bin/c3k4 searchff out/c3k4-ff.txt` | `rank-4 Hermitian forms: 3570, orbits under 14 generators: 2`, `fixed pairs: 14 orbits, 10 with partition ok`, `c3k4-ff seeds 0..10 (r = 4): leaves 14337997184, passed 37653692, solutions 267600, dead at init 0, APN spaces 2230`: `4:2` 2072, `4:3` 158 (none with ≥ 4) → `results/c3k4-ff.txt` | `20261024-102913-apn8-c3k4-ff.log` | yes, same counts and the same 2,230 spaces (20 min) |

### Order 2

j ≤ 3 and the structure theorem for j = 4 are by hand (write-up, "Order 2").

| case | command | expected output | original log | re-run |
|---|---|---|---|---|
| **the 18 orbits of fixed parts** | `bin/inv2 forbits out/inv2-reps.txt` | `Stab(spread) in GL(4,2): 360`, `orbits of L: 18 (…), passing partition + kernel tests: 18`, `sizes: 32768^1 163840^1 196608^1 327680^3 491520^2 655360^1 983040^4 1966080^5` → `results/inv2-reps.txt` | `20261024-113705-inv2-forbits.log` | yes, same, and the file is identical (57 s) |
| **extending all 18** | `bin/inv2 search results/inv2-reps.txt out/inv2-all.txt 0 18 0` | `solutions 3661, bugs 0`: `4:0` 3535, `4:2` 96, `4:4` 30; `with >= 5 rank-4: 0` → `results/inv2-all.txt` (the last argument 0 writes every space; the default 4 writes those with ≥ 4) | `20261024-115943-inv2-all.log` | yes, same counts, same 3,661 spaces (34 s) |

## Controls and independent checks

| check | command | expected output | original log | re-run |
|---|---|---|---|---|
| orbit5 re-finds Gold x³ from one of its forms | `bin/orbit5 1 out/orbit5-val.txt 808ac1d` | `partition ok 1, APN spaces 2`, both `ranks 8:170 6:85 4:0 2:0` | none kept (the first version's output; see Provenance) | yes (0.4 s) |
| orbit5's enumerator lists each subspace once | `bin/orbit5 1 out/e1.txt 1 enumtest`, `bin/orbit5 1 out/e2.txt 808ac1d enumtest` | `enumtest: k1 leaves 11811, distinct 11811` ([7,3]₂); `enumtest: k1 leaves 200787, distinct 200787` ([8,4]₂) | `20261024-094037-enumtest.log`, `20261024-094041-enumtest.log` | yes, same |
| symp's enumerator | `bin/symp enumtest c5a` (also `c7b`, `c3k2`) | every line `OK`, `enumtest: 0 failures` | none | yes |
| symp re-finds Gold x³ | `bin/symp ext gf5 x3 1 out/ext-x3.txt` | `space x3 invariant under A: yes`, `distinct APN spaces 2`, `the space itself found: YES` | none | yes |
| review_odd_brute's enumeration sizes | `bin/review_odd_brute c5a 0 count`, `bin/review_odd_brute c7b 0 count` | every part `Gaussian binomial … ok`; totals 14053980821 and 12310025 | none | yes |
| c3k4 re-finds the 25 known C₃ classes with 4 rank-4 | `bin/c3k4 validate data/seeds-known-rare.txt` | `validate: 25 classes tried, 25 refound` | `20261024-113705-c3k4-validate.log` | yes, same (21 s) |
| c3k4's Engine at r = 1, 2, 3 (reviewer) | `bin/review_odd_engine results/c3k4-ff.txt 300 0` | `r=1: 300 tried, 300 re-found`, same for r=2 and r=3 (r = 4 is covered by searchff itself, which found these spaces at r = 4; `… 300 2` also tries two at r = 4, minutes each) | `20261024-135907-review-engine-r123.log` | yes, same (291 s) |
| inv2 re-finds the isolated class with 4 rank-4 from its own fixed part | `bin/inv2 check data/seed474.txt 0` | `involutions in Aut(S8): 1`, `dim S8^A 4, dim D(S8) 4`, `b-forms: rank-2 5 (radicals disjoint 1, cover W 1), nondegenerate 10`, `solutions 2007` (1917 / 60 / 30), `check: re-found` | `20261024-113705-inv2-check474.log` | yes, same (52 s) |
| the order-2 structure on that function, in Python | `AUT8_PRINT=2 bin/aut8 data/seed474.txt` prints the involution `A 02 01 08 04 20 10 80 40`; then `python src/inv2check.py data/seed474.txt 02 01 08 04 20 10 80 40` | `j = rank N = 4`, `D(S8) == S8^A: True`, `fixed forms by \|rad & K \ 0\|: {0: 10, 3: 5}`, `K\0 covered exactly once by the fixed radicals: True`, `rank-4 components: 4  fixed: 0` | none | yes |
| order-2 orbits and outputs from raw matrices (reviewer) | `python src/review_orbits.py` | 18 lines `rep …: orbit N (inv2 N) disjoint True`, `cosets covered 512 of 512; all reps consistent True`, `inv2-all: 3661 lines, bad 0, distribution {(170, 85, 0, 0): 3535, (178, 75, 2, 0): 96, (186, 65, 4, 0): 30}` | none (reviewer, direct) | yes (9 s) |
| inv2's pruning and solve (reviewer) | `bin/review_engine lvl all`, `bin/review_engine sols 3` | every rep `0 mismatch` at levels 1–3 and for solve vs brute force; `sols: tested 1221 lines, bad 0, brute solutions on those triples 1230` | none (reviewer, direct) | yes (14 s, 218 s) |
| inv2 on random conjugates of known functions (reviewer) | `bin/review_inv2 gold 20`, `bin/review_inv2 known 20`, `bin/review_inv2 brute 16 6`, `bin/review_inv2 actl 200`, `bin/review_inv2 init` | `gold: found 20 / 20`; `known: found 20 / 20`; every brute line `agree`; `actl: 0 bad of 200`; `init ok 18 / 18; \|Stab(spread)\| 360` | none (reviewer, direct; its trial counts weren't recorded, these are ours) | yes (seconds each) |
| outputs are APN with the stated ranks | `python src/verify.py results/c3k4-ff.txt`, `python src/verify.py results/inv2-all.txt` | `checked 2230, bad 0`; `checked 3661, bad 0` | none | yes (18 s, 32 s) |
| symmetry of the known extreme classes (write-up, context) | `bin/aut8 data/seeds-known-rare.txt data/seeds-new63.txt` | 4 rank-4 (`ranks 8:186 6:65 4:4`): 25 with `\|Aut\| 3`, 1 (line 474) with `\|Aut\| 4`; 3 rank-4, known lists: 144 with `\|Aut\| 3`, 2 with `\|Aut\| 15`; the 63: all `\|Aut\| 1` | `20261024-094507-aut8.log` (both files), `20261024-095108-aut8.log` (adds the `fixS8` field) | yes, same (67 s) |

**Novelty of the by-products** (not part of the proof). The APN spaces the searches found were compared with the
published lists by classify8's invariant: `bin/classify8 --min4 2 out/ff-vs-known.txt results/c3k4-ff.txt
apn_8bit.txt blp22_8bit.txt new_apns.txt` (`logs/20261024-115814-ff-vs-known.log`) and the same for
`results/inv2-all.txt` against those lists plus Yu–Wang–Li's and Yu–Perrin's QAM lists
(`logs/20261024-121546-inv2-vs-known.log`). The first of these shows two invariants with two rank-4 components found
only in the search output (`in{0}:2`); the comparison with all six lists (`bin/classify8 --min4 2 out/ff-vs-all6.txt
results/c3k4-ff.txt apn_8bit.txt blp22_8bit.txt new_apns.txt yuwangli_qam.txt yuperrin_qam.txt`,
`logs/20261024-153931-ff-vs-all6.log`) finds both in Yu–Wang–Li's (`in{04}:2`, and no `in{0}` alone). So every class with two or more rank-4 components that the searches found is in a published list. Not
re-run here: the lists (about 3.4 GB) aren't included.

## Not re-run here

- `orbit5 1 … fixed`, the order-5 brute force (about 21 CPU-hours): original log `20261024-094223-apn8-fix5-c1.log`.
  It is the second of two independent answers for that case; `fix5check` (19 s) and `review_odd_brute c5a` (5 min)
  were re-run.
- `symp fixed5 c3k2` (65 min there): original log `20261024-103000-symp-c3k2-fixed5.log`.
- The classify8 comparisons with the published lists (above).

## Provenance of the original runs

The logs come from the runs that produced the final numbers; several used binaries built from an earlier version of
a source, and the re-runs above confirm today's sources give the same counts.

- `orbit5`: the 5-cycle runs used the second version (each subspace enumerated once), the fixed-case brute force a
  build from just before the Q range and progress lines were added. The very first version's validation output
  (`orbit5 1 … 808ac1d`, not kept) listed a second Gold-type space that doesn't contain the given form; today's lists
  two that do.
- `c3k4 search` (3+1+1): built before a fix to a `1L << 40` default that overflows where `long` is 32 bits; the log
  shows the full seed range (0..214), and today's source gives the same counts. The other c3k4 runs used builds of
  today's source.
- `fix5check`: the 0-APN run used a build from just before the third argument (min rank-4) was added; the control
  used today's source.
- `symp`: the order-5 and order-7 runs used a build from before its last edit; the order-3 runs used today's.
- `inv2 forbits` (and `inv2 search 12 18`, which wrote `results/inv2-search-12.txt`) used a build from before its last
  edit; `inv2-all` used today's.
- `aut8` `20261024-094507` predates the `fixS8` field in its output.
- The share runner's lines in the logs (`[runner …]`) record when the job ran, paused and resumed on the shared machine.
