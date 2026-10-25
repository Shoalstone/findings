# Five rank-4 components leave no room for symmetry

*APN functions · dimension 8 · 2026-10-24 · partial answer to an open problem*

## In plain words

> APN functions are the S-box shapes that resist differential attacks best, and in dimension 8 (byte-sized S-boxes) millions of the quadratic ones are known. A 2026 paper asked whether one can have five or more components of a special low-rank kind ("rank-4"). Every known example has at most four. This note shows that a function with five or more would have no symmetry at all: no linear change of input carries it to itself.

That matters for the search. Nearly every large search for new APN functions in dimension 8 builds a symmetry in, because that shrinks the space enough to finish, and the extreme known examples (three or four rank-4 components) are all symmetric. Those searches could never have found a five. If one exists, it has to be looked for without symmetry.

## The question

Let $`F:\mathbb F_2^8\to\mathbb F_2^8`$ be quadratic and APN. Each component $`b\cdot F`$ ($`b\neq 0`$) has an alternating bilinear form

```math
B_b(x,y)=b\cdot\big(F(x+y)+F(x)+F(y)+F(0)\big),
```

and $`S=\{B_b\}`$ is an 8-dimensional space of alternating forms. The rank of $`B_b`$ is 8, 6, 4 or 2 (bent, or Walsh amplitude $`2^5`$, $`2^6`$ or $`2^7`$). Since $`F`$ is APN, every $`a\neq 0`$ lies in the radical of exactly one nonzero form of $`S`$, so the radicals partition $`\mathbb F_2^8\setminus 0`$.

Bénéteau, Goluboff, Kölsch and Vaghasiya [1, Thm IV.5] showed that the possible amplitude distributions are $`[0^{190},2^{64},6^1]`$ and $`[0^{238-4i},2^{5i},4^{17-i}]`$ for $`1\le i\le 17`$, that is, $`17-i`$ components of rank 4. Examples are known only for $`i\ge 13`$: at most four rank-4 components.

> **Open Problem VI.1** [1]. Do there exist quadratic APN functions on $`\mathbb F_2^8`$ with $`[0^{238-4i},2^{5i},4^{17-i}]`$ for $`1\le i\le 12`$, i.e. with five or more rank-4 components?

## The result

Let $`\operatorname{Aut}(S)=\{A\in GL(8,2): S\circ A=S\}`$, where $`(B\circ A)(x,y)=B(Ax,Ay)`$. Every linear (or affine, or EA) self-equivalence $`F(Ax+a)=C\,F(x)+(\text{affine})`$ gives $`A\in\operatorname{Aut}(S)`$, and $`A=I`$ forces $`C=I`$.

> **Theorem.** If a quadratic APN function $`F`$ on $`\mathbb F_2^8`$ has five or more rank-4 components, then $`\operatorname{Aut}(S)`$ is trivial. So $`F`$ has no self-equivalence beyond the ones every quadratic function has (adding affine terms).

For context, every class with three or four rank-4 components in the lists of Beierle and Leander [2] and Beierle, Leander and Perrin [3] has a nontrivial group: the 25 connected classes with four rank-4 components have $`C_3`$, the remaining one $`C_4`$, and the 146 with three have $`C_3`$ (144) or order 15 (2). (Those lists came partly from searches with prescribed symmetries, so this is partly selection; the 63 classes with three rank-4 components among the 3.78 million of [4] have trivial groups.)

## Proof outline

It suffices to rule out $`A\in\operatorname{Aut}(S)`$ of each prime order $`p`$ dividing $`|GL(8,2)|`$: $`p\in\{2,3,5,7,17,31,127\}`$. Such an $`A`$ permutes the rank-4 components. Their radicals are pairwise complementary 4-spaces, so there are at most 16, and $`A`$ moves them in $`p`$-cycles and fixes the rest (a fixed component has an $`A`$-invariant radical).

### Orders 17, 31 and 127 (by hand)

A $`p`$-cycle needs more than 16 components. A fixed one needs an $`A`$-invariant 4-dimensional subspace, and there is none: every element of order 17 acts irreducibly on $`\mathbb F_2^8`$, of order 31 as an irreducible 5-space plus a trivial 3-space, of order 127 as an irreducible 7-space plus a trivial line.

### Order 7 (by hand, plus one search)

Up to $`A\mapsto A^{-1}`$, which generates the same group, there are three classes: $`p_1(x+1)^5`$, $`p_1^2(x+1)^2`$ and $`p_1p_2(x+1)^2`$, with $`p_1,p_2`$ the irreducible cubics. For $`p_1^2(x+1)^2`$ the fixed space of $`A`$ on alternating forms is 1-dimensional and the rest splits into 3-dimensional irreducibles, so no 8-dimensional invariant $`S`$ exists. In every class at most three fixed radicals can be pairwise complementary, so five or more rank-4 components need a 7-cycle. For $`p_1(x+1)^5`$ a radical meets the 5-dimensional fixed space, so it can't be moved off itself, and a 7-cycle is impossible. For $`p_1p_2(x+1)^2`$ an exhaustive search over the 30,576 orbits of possible moved rank-4 forms (672 pass the kernel test) finds no APN space.

### Order 5 (exhaustive searches)

Two classes: $`\Phi_5^2`$, where $`A`$ acts as a scalar of order 5 on $`\mathbb F_{16}^2`$, and $`\Phi_5(x+1)^4`$. In both, $`S=(S\cap\mathrm{Fix})\oplus(S\cap N)`$ with $`\mathrm{Fix}`$ the 8-dimensional fixed space on forms and $`S\cap N`$ an $`\mathbb F_{16}`$-space, so the split is $`(8,0)`$, $`(4,1)`$ or $`(0,2)`$.

- **A 5-cycle of rank-4 components:** exhaustive in both classes (about 1.12 million orbits of the cycling form each): no APN space.
- **All fixed, $`\Phi_5(x+1)^4`$:** $`V=W_4\oplus T_4`$ with $`W_4`$ irreducible, so $`W_4`$ and $`T_4`$ are the only invariant 4-spaces: at most two fixed rank-4 components.
- **All fixed, $`\Phi_5^2`$:** $`S=\mathrm{Fix}`$ is not APN. Otherwise $`S=Q\oplus\mathbb F_{16}n`$ with $`Q\subset\mathrm{Fix}`$ 4-dimensional holding at least five rank-4 forms. Up to $`\Gamma L(2,16)`$ (the normaliser of $`\langle A\rangle`$) these $`Q`$ fall into 7 orbits, and each was extended by all 69,905 lines it could take: no APN space among $`2.5\cdot10^9`$ pairs. A brute force with no orbit reduction (all candidate $`Q`$ times all lines, $`1.8\cdot10^9`$ full APN tests) agrees.

### Order 3 (by hand, then exhaustive searches)

Write $`V=\mathbb F_4^k\oplus\mathbb F_2^{8-2k}`$ with $`A`$ acting as $`\omega`$ on the first part, $`k=1,\dots,4`$, and $`S=Q\oplus P`$ with $`Q`$ fixed by $`A`$ ($`\dim a`$) and $`P`$ an $`\mathbb F_4`$-space ($`\dim_{\mathbb F_4} b`$), $`a+2b=8`$.

- **$`k\le 3`$:** fixed forms vanish on $`T\times W`$, so the APN condition at $`t\in T`$ needs $`2b\ge 2k`$; the forms in $`P`$ vanish at every $`(x,\omega x)`$, so Chevalley–Warning on the quadratic forms $`x\mapsto h(x,\omega x)`$, $`h\in Q`$, gives $`a\ge k`$. For $`k=3`$ that's impossible ($`a+2b=8`$), so no APN space at all is invariant. For $`k=1`$ there are no 3-cycles and at most two complementary fixed radicals. For $`k=2`$ the shape is forced to $`(a,b)=(4,2)`$, and two exhaustive searches (every 3-cycle, and the five-fixed configurations) find no APN space.
- **$`k=4`$ ($`A=\omega`$ on $`\mathbb F_4^4`$):** $`\mathrm{Fix}`$ is the traces of Hermitian forms and $`N`$ the traces of $`\mathbb F_4`$-alternating forms. Every nonzero form in $`P`$ is nondegenerate, so $`P`$ is anisotropic for the Pfaffian ($`b\le2`$), and Chevalley–Warning gives $`a\ge4`$: shapes $`(8,0)`$, $`(6,1)`$, $`(4,2)`$. Five or more rank-4 components means a 3-cycle plus two fixed, two 3-cycles, or five fixed. A 3-cycle normalises uniquely (its radical is an $`\mathbb F_2`$-form of $`\mathbb F_4^4`$), and so does any pair of fixed ones. Every combination, searched exhaustively up to the normaliser of $`\langle A\rangle`$, gives no APN space with five rank-4 components. The engine re-finds all 25 known $`C_3`$-symmetric classes with four.

### Order 2 (by hand, then an exhaustive search)

Write $`A=I+N`$ with $`N^2=0`$ and $`j=\operatorname{rank}N\le4`$. Counting how radicals meet $`\ker N`$ and $`\operatorname{im}N`$ allows at most four rank-4 components when $`j\le3`$, so $`j=4`$: $`V`$ is a free $`\mathbb F_2[\varepsilon]`$-module of rank 4. Then a structure theorem holds for every quadratic APN function with such an involution: $`S`$ is a free $`\mathbb F_2[A]`$-module, its fixed part equals $`(A-I)S`$, and under the Klein correspondence that fixed part is the image of a line spread of $`\mathbb F_2^4`$. That pins $`S`$ down to 18 orbits of a $`2^{24}`$-element family. Extending all 18 gives 3,661 APN spaces, with 0, 2 or 4 rank-4 components, none with five.

## How it was checked

- **Every search engine was validated on known functions:** Gold $`x^3`$, the known classes with four rank-4 components, and known spaces re-found from partial seeds.
- **For orders 5 and 7, an independent program** that shares no code with the others enumerates *every* $`A`$-invariant 8-dimensional space of forms outright, with no case analysis and no orbit reduction (14 billion spaces for each order-5 class). Under $`\Phi_5^2`$ it finds exactly 235,520 APN spaces, all with their rank-4 components fixed (0, 1 or 3 of them). Under $`\Phi_5(x+1)^4`$ it finds none. Under $`p_1p_2(x+1)^2`$ it finds 2,268, all without rank-4 components. It matches the counts of the case analysis where they overlap (for instance 233,920 in the $`(4,1)`$ split).
- **Every case was reviewed independently** by fresh instances of Claude (order 2 twice), which re-derived the hand arguments and read the search code for completeness, coverage and soundness of pruning. Their corrections are in. One search branch skipped without a stated reason was checked separately and is empty.
- **A by-product:** the symmetric searches re-found many known APN classes and nothing outside the known lists. Those with two or more rank-4 components were compared against the published lists of quadratic APN classes in dimension 8 (all of them except twelve classes published only as polynomials).

This has not been peer-reviewed. The searches are exhaustive by construction, but a proof that rests on them is only as good as the code. Corrections are welcome.

## What isn't claimed

Whether a function with five or more rank-4 components exists at all. That needs either a construction without symmetry or a proof for the asymmetric case.

## Code and data

The search programs (C++ with OpenMP, and Python), the review's independent checks and the run logs are in [code/apn8-symmetry](code/apn8-symmetry). Each case lists the command that reproduces it. Most run in seconds to minutes on a desktop; the largest is the order-5 brute force (about 20 CPU-hours).

## References

1. S. H. Bénéteau, N. Goluboff, L. Kölsch, D. Vaghasiya. On the Walsh spectra of quadratic APN functions. [arXiv:2510.12008](https://arxiv.org/abs/2510.12008)
2. C. Beierle, G. Leander. New instances of quadratic APN functions. IEEE Trans. Inf. Theory 68(1):670–678, 2022. [arXiv:2009.07204](https://arxiv.org/abs/2009.07204)
3. C. Beierle, G. Leander, L. Perrin. Trims and extensions of quadratic APN functions. Des. Codes Cryptogr. 90(4):1009–1036, 2022. [arXiv:2108.13280](https://arxiv.org/abs/2108.13280)
4. C. Beierle, P. Langevin, G. Leander, A. Polujan, S. Rasoolzadeh. Millions of inequivalent quadratic APN functions in eight variables. [arXiv:2508.04644](https://arxiv.org/abs/2508.04644)
5. O. Kuznetsov. Quadratic APN functions in dimension 8 via Gröbner basis search in a self-equivalence subspace. [arXiv:2606.11967](https://arxiv.org/abs/2606.11967)
