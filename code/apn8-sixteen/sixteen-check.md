# Independent check: is there a quadratic APN F on F_2^8 with amplitude distribution [0^234, 2^5, 4^16]?

Subagent (Claude Opus 5.5), 2026-10-24, for Rowen. Written without reading `sixteen.cpp` or `runs/a/sixteen.txt`.
Files: `sixteen_check.cpp` (search + controls), `sixteen_check_level2.py` (numpy recount of the first two levels),
scratch in `runs/check16/`.
Status: done (2026-10-24). Verdict in section 5: no such function.

## 1. The reduction, checked

Setup. S8 = {B_b} is an 8-dim space of alternating forms (b -> B_b injective), and F is APN iff every a != 0 lies in
the radical of exactly one nonzero B_b. (For a != 0, {b : a in rad B_b} is the annihilator of Im(L_a),
L_a(x) = F(x+a)+F(x)+F(a)+F(0); APN <=> Im L_a has dim 7 <=> that set is one line {0, b}.) Conversely any 8-dim
space of alternating forms is the S8 of a quadratic F (F_i(x) = sum_{j<k} B_i(e_j,e_k) x_j x_k), so existence of F is
existence of S8. [0^234, 2^5, 4^16] = 234 nondegenerate forms, 5 of rank 6 (2-dim radical), 16 of rank 4 (4-dim);
5*3 + 16*15 = 255, so the radicals partition V \ 0.

Step 1 holds. Two 4-dim radicals meet in 0, so they're complementary. For a hyperplane H, each rank-4 radical lies in
H (15 points) or meets it in a 3-space (7 points); if t of them lie in H, H\0 has 127 - (112 + 8t) = 15 - 8t of the
15 uncovered points, so t in {0, 1} and |H ∩ holes| in {15, 7}. With T = holes ∪ {0} (16 points),
sum_{t in T} (-1)^{u.t} = 1 + (2h - 15) in {0, 16} for every u != 0. Parseval: 256 |T| = 16^2 (1 + #{u != 0 :
value 16}) gives exactly 15 nonzero u with u ⊥ T, so T^⊥ has 16 elements, span(T) has dim 4 = log2 |T|, and
T = span(T) is a 4-space R17. It meets every rank-4 radical in 0, so with them it's a spread (17 = 2^4 + 1), and the
five 2-dim radicals (inside T by construction) are a line spread of R17. (This is the partial-spread-of-deficiency-1
fact; the paper's Thm II.21 lists [2^5, 4^16] among the partition types, consistent.)

Step 2 holds: GL(8,2) is transitive on ordered triples of pairwise complementary 4-spaces (V = A ⊕ B, the third is the
graph of an isomorphism A -> B, change basis on B), the stabiliser of (X, Y, graph(I)) is {diag(g, g)}, and every
4-space complementary to Y is a graph {(x, Nx)}; complementary to X <=> N invertible, to graph(N') <=> N + N'
invertible.

Step 3 holds: a form with radical exactly W is the pullback of a nondegenerate form on V/W, and (x, y) -> Nx + y is a
surjection V -> F_2^4 with kernel graph(N); there are 28 nondegenerate alternating forms on F_2^4 (2^2 * 7).
The form with radical X is d(y, y'), and diag(g, g) acts on d transitively, so one fixed d is fine.
**One caution** (I can't see whether it applies, since I haven't read sixteen.cpp): the d-normalisation uses up the
same group diag(g, g) that also acts on spread sets by conjugation N -> g N g^{-1}. Fixing d is valid if the spread sets
are then enumerated in full (or up to conjugation by Sp(d) only, order 720). Fixing d *and* taking spread sets only up
to GL(4,2)-conjugacy would use the group twice and could miss solutions.

My normalisation (same facts, used in a different order): R17 = Y = {x = 0}; the radical X = {y = 0} carries
B0 = J(y, y') (J the standard symplectic form; diag(I, h) moves any d to J and fixes X, Y); the rank-4 radical
through the point (e0, e0) is graph(N1) for some invertible N1, and diag(N1, I) moves it to graph(I) fixing X, Y, B0;
finally diag(h, h), h in Sp(J), fixes X, Y, graph(I) and B0 and acts on the form M1 of graph(I) by congruence, with
orbits of sizes 1 (M1 = J), 15 and 12 on the 28 forms. So three starting pairs (B0, B1) cover everything, and the
search runs once from each (and once from all 28 as a check).

## 2. Method (how it differs from a spread-set clique search)

No spread sets are enumerated. The search grows the space of forms itself, one form at a time, and lets the spread
emerge as the radicals of what gets added:

- State: the span S (all 2^d elements), and the set of points already covered by radicals of its nonzero elements.
  Every element of S must be legal (nondegenerate; or rank 4 with radical meeting Y trivially; or rank 6 with
  radical inside Y) and radicals must be pairwise disjoint. Radicals are computed directly from each form (Gray code
  over all 256 vectors), never taken from how the form was built.
- Branching: take the first point not covered (points outside Y in ascending order, then points of Y) and branch
  over every legal form whose radical contains it. Outside Y that's a rank-4 form with radical graph(N), N x = y,
  N invertible, times the 28 nondegenerate forms on the quotient; on Y it's a rank-6 form with radical a line of Y
  through the point, times the 13,888 nondegenerate forms on V/L. In any solution the form covering that point is in
  S8, so this is complete, and each solution is reached by exactly one path. Each step raises dim S by one, so the
  depth is at most 6 after B0, B1.
- If every point is covered while dim S < 8, the rest must be all-bent: every coset representative of the
  complement (2^(28-d) of them) is tried. If dim S = 8 with points uncovered, the branch dies.
- Starts: B0 = J(y,y') and B1 = M1(x+y, x'+y') for the three Sp(J)-orbit representatives M1 (and, separately, all 28).

So the pruning is by the whole span (all 2^d elements legal and disjoint) rather than pairwise spread compatibility
first and forms second. Generators, legality test and verifier share nothing with sixteen.cpp, and any solution is
re-verified from scratch (LUT from the forms, every derivative 2-to-1, Walsh spectrum by fast transform).

## 3. Result

`sixteen_check main` (job 20261024-150959-check16-main, ~470 s wall on a shared box): **0 solutions**.
16,076 dim-3 nodes (5316 + 5280 + 5480 over the three M1), then by span dimension 4 / 5 / 6: 2,786,225 / 1,798,966 /
14,042 nodes; ~2.9e9 candidate forms tried. No branch ever reached dimension 7, no node had every point covered (so
the bent-extension step never ran). Progress file `runs/check16/main.prog`.

Cross-runs of the same search with a different tree:
- descending point order (`main ... desc`, a different branching point at every node): 16,080 dim-3 nodes,
  2,816,830 / 2,215,586 / 13,146 at dim 4 / 5 / 6, again nothing at dim 7, **0 solutions**.
- all 28 forms M1 instead of the 3 orbit representatives (`main ... all`, no symmetry reduction after the
  normalisation): 148,808 dim-3 nodes, 25,552,430 / 15,895,972 / 39,852 at dim 4 / 5 / 6, nothing at dim 7,
  **0 solutions** (~2.7e10 tries, ~60 min wall).

Comparison with Rowen's run (read only after the above): `runs/a/sixteen.txt` enumerates 1,128,864 spread sets
containing 0 and I (19,936 semifield, 336 field) and finds 0 solutions in the form search. Same verdict. The 336 is
the full count of F_16-subfields, so his spread sets are not reduced up to conjugacy, and the caution in section 1
about fixing d doesn't apply to it.

## 4. Controls

- `selftest` (all pass): 28 / 13,888 nondegenerate forms on F_2^4 / F_2^6; |GL(4,2)| = 20160, |Sp(4,2)| = 720;
  10,795 lines and 200,787 4-spaces of F_2^8, 127 / 11,811 through each point; every generated form has exactly the
  intended radical, 28 (13,888) distinct per radical; a brute force over all 2^28 alternating forms finds exactly 28
  with radical graph(I) and 112,881,664 = 2^12*7*31*127 nondegenerate ones; the F_16-subfields of M_4(F_2) number
  **336** and are spread sets.
- First two levels recounted independently (`sixteen_check_level2.py`: numpy, forms built from the block formula
  [[N^T M N, N^T M], [M N, M]], radicals by brute force): 5316 / 5280 / 5480 dim-3 children for M1 = 12 / 13 / 33,
  the same as the C++.
- Known functions through the main code path (`ctlB`): bring a known function with ≥ 2 rank-4 components into the
  normalised coordinates (X = a rank-4 radical, B0 = J-form, a second rank-4 radical = graph(I), Y a 4-space missing
  every rank-4 radical), relax the rules only by allowing rank-6 radicals outside Y, seed B0, B1 and some rank-6
  forms, and search with the same generators, legality, branching and completion. With seed dimension 6 it re-finds
  the known space every time (BL22 lines 12631, 12633, 12640, 12650, 12655, 12733 of type 186/65/4 and 15, 12534,
  12624 of type 182/70/3; each verified APN with the right spectrum by the LUT check), and finds nothing else.
  A side observation: every one of these has a 4-space Y missing all its rank-4 radicals that contains 4 or 5 of its
  rank-6 radicals; for 12631, 12650, 12655, 12733 and the three type-4 ones, five, i.e. Y is partitioned exactly as
  R17 would be in a [0^234, 2^5, 4^16] function.
- Seed dimension 5 (line 12631, job 20261024-151000-check16-ctlb3): 3596 dim-6 nodes, 2.4e9 tries, 110 dim-7
  nodes, exactly one space found, the known one (verified APN, 186/65/4), 377 s.
- "Control C" (line 12631, job ...-check16-ctlc): seed B0, B1 only, as in the main search; rank-4 forms (graph(N)
  generator) and in-Y rank-6 forms generated in full, only the outside-Y rank-6 candidates restricted to the known
  function's own (otherwise the relaxed tree is far too big). It found the known space, so the main code path does
  carry a real solution from (B0, B1) to dimension 8. Cancelled after ~70 min on the full box without finishing the
  tree (it has no progress counter), so "nothing else" isn't established there; it only matters as a positive control.

## 5. Verdict

**No quadratic APN function F_2^8 -> F_2^8 has amplitude distribution [0^234, 2^5, 4^16]**, by an exhaustive search
organised differently from Rowen's (forms grown point by point, no spread-set enumeration), run three ways (orbit
representatives; all 28 M1; reversed point order), each 0, and agreeing with Rowen's spread-set search. Confidence:
high, resting on (a) the reduction (section 1, checked by hand), (b) completeness of point-branching (each solution's
covering form is among the candidates; generators checked to produce every form with the given radical, against a
2^28 brute force for one radical), (c) the positive controls through the same code path.
Not covered here: the other open cases i = 2..12 of VI.1 (the R17 argument is specific to 16 rank-4 components).


---

# Case i = 2: [0^230, 2^10, 4^15] (15 rank-4, 10 rank-6)

Same terms: written without reading sixteen.cpp's i2 mode or `runs/a/fifteen.txt` until the verdict was in. Code:
`sixteen_check.cpp` modes `i2`, `i2count`, `i2ctl`, `i2probe` (exes `runs/bin/sixteen_check3.exe`, later `5`);
`sixteen_check_level2_i2.py`, `sixteen_check_short.py`.
Status: done (2026-10-24). Verdict in section 8: no such function.

## 6. Rowen's reduction for i = 2, checked

Step 1 holds, and is slightly simpler than stated: k = 2 can't happen at all (two complementary 4-spaces span V, so
no hyperplane contains both), so a hyperplane holds 22 or 14 of the 30 holes and W(u) = 2h - 30 is 14 or -2. Then
Σ_{u≠0} W = -30 gives 16a = 480, a = 30, b = 225, and Parseval checks out (30·196 + 225·4 = 6780 = 7680 - 900).

Step 2 holds. The eigenvalues of Cay(F_2^8, D) (D = holes) are the character sums: 30 once, 14 thirty times, -2
225 times. A connected regular graph with three eigenvalues is strongly regular, with μ = k + θτ = 30 - 28 = 2 and
λ = μ + θ + τ = 14: srg(256, 30, 14, 2), the parameters of L2(16). Shrikhande's theorem (L2(n) is determined by its
parameters for n ≠ 4) is a statement about all graphs with those parameters, so it applies to this Cayley graph as
stated. But the conclusion doesn't need it; only λ = 14 and μ = 2 at the vertex 0 are used:
- Γ(0) is the graph on D with x ~ y iff x + y ∈ D; it's 14-regular (λ).
- Let y ≠ z be neighbours of x in Γ(0) with y ≁ z. Then 0, x and x + y + z are all common neighbours of y and z
  in the Cayley graph (y + x, z + x ∈ D, and (x+y+z) + y = x + z, (x+y+z) + z = x + y ∈ D). μ = 2 forces
  x + y + z ∈ {0, x}: x + y + z = x gives y = z, and x + y + z = 0 gives z = x + y, but then y + z = x ∈ D, so y ~ z.
  Either way a contradiction, so x's 14 neighbours form a clique, x plus them is a 15-clique, and since every degree
  is 14 it's a whole component. So Γ(0) is two disjoint 15-cliques.
- In a clique C, x ≠ y gives x + y ∈ D adjacent to x (x + (x + y) = y ∈ D), so x + y ∈ C: C ∪ 0 is a 4-space.
  The two are disjoint, so complementary: D = (R ∪ R') \ 0.

Step 3 holds: each rank-4 radical misses D, so the 15 of them plus R, R' are 17 pairwise complementary 4-spaces (a
spread); a line {a, b, a+b} ⊂ D with a ∈ R, b ∈ R' would have a + b ∉ R ∪ R', so each rank-6 radical lies in R or
in R' (5 in each, since each splits its 15 points). Normalising R = Y, R' = X and one rank-4 radical to graph(I) is
triple transitivity; the form on graph(I) can then be fixed (diag(g, g) is transitive on its 28 forms). One addition:
the swap (x, y) ↦ (y, x) also preserves everything (X ∪ Y, graph(I), the form J(x+y, x'+y')), since R and R' play
symmetric roles here; it sends graph(N) to graph(N^-1). I use it.

## 7. Method for i = 2

The same point-branching engine with new rules: rank-4 radicals must miss X and Y; rank-6 radicals must lie inside
X or inside Y; X points get lines of X. Point order: points off X ∪ Y ascending, then Y, then X.
Start: B1 = J(x+y, x'+y') on graph(I). The residual group G = {diag(h, h) : h ∈ Sp(J)} × <swap> (order 1440) fixes
X ∪ Y, graph(I) and B1. The candidates for a second rank-4 form are the 28 forms on each graph(N) with N, N + I
invertible, N ≠ I: 5824 N (Rowen's count too), 163,072 forms, which split into **200 G-orbits** (sizes
1x12 1x24 1x40 3x60 8x120 3x144 4x180 16x240 3x288 10x360 37x480 39x720 74x1440, summing to 163,072). For each
orbit representative j, search every S8 containing B1 and rep_j in which every rank-4 element has orbit index ≥ j
(valid: in any solution take the rank-4 form other than B1 with the smallest orbit index, and apply the element of G
that moves it to its representative).

Checks so far:
- the 1440 group elements are distinct, fix B1 and preserve X ∪ Y; orbit sizes sum to the candidate count and divide
  1440;
- numpy recount (`sixteen_check_level2_i2.py`, block-formula forms, brute-force radicals) of the first branching for
  reps 0, 1, 50, 120, 199 (no pruning): 1616, 1659, 1650, 1493, 1724 children, identical to `i2count`;
- known functions through the i = 2 code path (`i2ctl`): many known functions have two complementary 4-spaces, each
  split by 5 of their rank-6 radicals and missing every rank-4 radical (most of the type-4 ones in BL22, and 11 of
  the 26 type-5). Normalised into these coordinates (that pair = X, Y; a rank-4 radical = graph(I) with form J; the
  minimal-orbit-index rank-4 form moved onto its representative, pruning on), relaxed only by allowing rank-6
  radicals outside X ∪ Y, seeded to dimension 6, the search re-finds the known space and nothing else for lines 15,
  12494, 12505, 12514 (182/70/3) and 12631 (186/65/4), each verified APN with the right spectrum.

## 8. Result for i = 2

`sixteen_check i2` over all 200 orbit representatives (job 20261024-164052-check16-i2, a few minutes):
**0 solutions**. Nodes by span dimension 2 / 3 / 4 / 5 / 6: 200 / 155,307 / 5,014,809 / 644,278 / 837; nothing at
dimension 7; the all-bent extension step never ran. Cross-runs, each **0 solutions**:
- descending order of the points off X ∪ Y: 200 / 159,536 / 5,138,508 / 600,417 / 516, nothing at 7;
- orbit-index pruning off (so the result doesn't lean on the minimal-index argument, only on the orbit reduction):
  200 / 312,842 / 15,083,757 / 2,748,440 / 4,752, nothing at 7.

Positive controls at seed dimension 5 (B1, the representative, 3 rank-6 forms; the rest found by the search): line 15
(182/70/3) and line 12631 (186/65/4) each re-found as the only space, 4408 / 4817 dim-6 nodes, ~2.9e9 / 3.1e9 tries,
500 / 540 s.

Comparison with Rowen (read after the above): `runs/a/fifteen.txt` (his run still going) lists "SHORT" spans,
6-dimensional with 15 rank-4 forms. `sixteen_check_short.py` checks them coordinate-free: 15 rank-4, 48 bent, no
rank 6, radicals disjoint, holes = two complementary 4-spaces, i.e. fully legal partial spaces with every non-hole
point covered. `sixteen_check i2probe` normalises each into my coordinates (all land on orbit representative 51) and
reruns that representative with logging: **all 20 are visited** by my search at the moment every point off X ∪ Y
is covered. A full run with logging (job 20261024-170915-check16-i2-p1) shows these are the only such nodes in the
whole search: **6 spaces, all from representative 51**, each 6-dimensional with 15 rank-4 and 48 bent forms, and
none takes even one rank-6 form (7 lines of Y through the first hole point x 13,888 forms each, all rejected). So
the two searches meet the same near-misses; mine says they're dead ends. (Rowen's run was still going when I read
it, so I can't compare final verdicts.)

Verdict: **no quadratic APN F_2^8 -> F_2^8 has amplitude distribution [0^230, 2^10, 4^15]**, on the same footing as
i = 1: the reduction (section 6, with an elementary replacement for Shrikhande), complete point-branching, three
differently-shaped runs at 0, positive controls through the same path, and agreement with Rowen's independent
search on the near-misses.
