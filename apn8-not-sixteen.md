# No quadratic APN function in dimension 8 has fifteen or sixteen rank-4 components

*APN functions · dimension 8 · 2026-10-24 · settles two cases of an open problem*

## In plain words

> A 2026 paper listed the shapes a quadratic APN function on 8 bits could still have, and asked which of them actually occur. It singled out one as particularly interesting: sixteen components of a special low-rank kind ("rank-4"), the most the list allows. It would make a certain set as small as such sets can be. This note shows that shape doesn't occur, and neither does the next one down, fifteen. Every quadratic APN function in dimension 8 has at most fourteen rank-4 components.

A short argument shows that sixteen or fifteen rank-4 components would force a very rigid structure, a *spread*: a way of cutting $`\mathbb F_2^8`$ into seventeen 4-dimensional pieces. That makes a complete computer search small. Two searches, written separately and organised differently, find nothing in either case.

## The question

Notation as in the [companion note](apn8-no-symmetry.md): for a quadratic APN $`F`$ on $`\mathbb F_2^8`$, $`S=\{B_b\}`$ is the 8-dimensional space of the alternating forms $`B_b(x,y)=b\cdot(F(x+y)+F(x)+F(y)+F(0))`$, and the radicals of its nonzero forms partition $`\mathbb F_2^8\setminus 0`$.

Bénéteau, Goluboff, Kölsch and Vaghasiya [1, Thm IV.5] showed the amplitude distribution is $`[0^{190},2^{64},6^1]`$ or $`[0^{238-4i},2^{5i},4^{17-i}]`$ with $`1\le i\le 17`$, and asked in Open Problem VI.1 whether $`1\le i\le 12`$ occur. Of these they single out $`i=1`$, $`[0^{234},2^5,4^{16}]`$, where the set of non-bent components would be a smallest possible nontrivial blocking set.

> **Theorem.** No quadratic APN function $`F:\mathbb F_2^8\to\mathbb F_2^8`$ has amplitude distribution $`[0^{234},2^5,4^{16}]`$ or $`[0^{230},2^{10},4^{15}]`$. Together with [1, Thm IV.1] (which rules out 17), every quadratic APN function on $`\mathbb F_2^8`$ has at most 14 rank-4 components.

## Sixteen: the reduction (by hand)

Suppose $`F`$ has 16 rank-4 components, with 4-dimensional radicals $`R_1,\dots,R_{16}`$, and 5 of rank 6, with 2-dimensional radicals. The radicals partition the 255 nonzero points ($`16\cdot15+5\cdot3=255`$), so the $`R_i`$ are pairwise complementary.

**The leftover points form a 4-space.** Let $`T`$ be the 15 points outside $`\bigcup R_i`$, plus $`0`$. A hyperplane $`H`$ contains some $`k`$ of the $`R_i`$ and meets the others in 3-spaces, so

```math
15k+7(16-k)+|H\cap T\setminus 0|=127,\qquad |H\cap T\setminus0|=15-8k\in\{15,7\}.
```

So the character sum $`\sum_{t\in T}(-1)^{u\cdot t}`$ is 16 or 0 for every $`u\neq 0`$. By Parseval, exactly 15 nonzero $`u`$ give 16, i.e. $`|T^\perp|=16`$, so $`T`$ spans only 4 dimensions and, having 16 points, is a 4-space $`R_{17}`$. Then $`\{R_1,\dots,R_{17}\}`$ is a spread of $`\mathbb F_2^8`$, and the five rank-6 radicals are a line spread of $`R_{17}`$.

**Normalising.** $`GL(8,2)`$ is transitive on ordered triples of pairwise complementary 4-spaces, so we may take $`R_{17}=Y=\{(0,y)\}`$, $`R_1=X=\{(x,0)\}`$ and $`R_2=\{(x,x)\}`$. Every 4-space complementary to $`Y`$ is a graph $`\{(x,Nx)\}`$, so the spread is $`Y`$ together with the graphs of a *spread set* $`\Sigma`$: sixteen $`4\times4`$ matrices, containing $`0`$ and $`I`$, whose pairwise differences are invertible.

**The forms.** A form with radical exactly $`\{(x,Nx)\}`$ is $`d(Nx_1+x_2,\,Ny_1+y_2)`$ for one of the 28 nondegenerate alternating forms $`d`$ on $`\mathbb F_2^4`$. The maps $`\mathrm{diag}(g,g)`$ fix $`X`$, $`Y`$ and $`R_2`$ and permute the spread sets, so if every spread set is searched, the form on $`X`$ can be fixed.

What remains: for each $`\Sigma`$, choose one of 28 forms for each $`N\in\Sigma`$ so that the forms span an 8-dimensional space in which every nonzero form is nondegenerate, or has rank 4 with radical a graph of some $`N\in\Sigma`$, or has rank 6 with radical inside $`Y`$, with all radicals disjoint.

## Sixteen: the searches

**Search A** (`sixteen.cpp`) enumerates every spread set: all 14-cliques among the 5,824 matrices $`N`$ with $`N`$ and $`N+I`$ invertible, joined by invertible differences. That gives **1,128,864 spread sets** containing $`0`$ and $`I`$. For each one, a depth-first search assigns forms, taking the lowest uncovered radical next and pruning on every element of the growing span. Result: **no solutions**, and no case where the forms cover all points in fewer than 8 dimensions. About an hour on one desktop.

**Search B** (`sixteen_check.cpp`, written separately without reading A) never enumerates spread sets. It normalises by the same facts in a different order, which leaves three starting pairs of forms. It then grows the space of forms directly: at each step it takes the first point not yet in a radical and branches over every legal form whose radical contains it (rank-4 forms with a graph radical, or rank-6 forms whose radical is a line of $`Y`$). If all points are covered early, it completes the space by brute force over all-bent extensions. The spread emerges from the radicals instead of being chosen first. Result: **no solutions**, in three runs with different trees (orbit representatives; reversed point order; all 28 starting forms with no symmetry reduction). No branch ever got past dimension 6.

## Fifteen

Now 15 rank-4 radicals $`R_1,\dots,R_{15}`$ (again pairwise complementary) and 10 rank-6 radicals covering the remaining set $`D`$ of 30 points.

**The leftover points are two 4-spaces.** A hyperplane contains at most one $`R_i`$ (two complementary 4-spaces span everything), so it holds $`127-15k-7(15-k)=22-8k\in\{22,14\}`$ points of $`D`$, and $`W(u)=\sum_{d\in D}(-1)^{u\cdot d}`$ is $`14`$ or $`-2`$ for $`u\ne0`$. Since $`\sum_{u\ne0}W(u)=-30`$, exactly 30 of these values are $`14`$. They are the eigenvalues of the Cayley graph $`\Gamma=\mathrm{Cay}(\mathbb F_2^8,D)`$, so $`\Gamma`$ is strongly regular with parameters $`(256,30,14,2)`$. Only two consequences at the vertex $`0`$ are needed:

- On $`D`$ put $`x\sim y`$ when $`x+y\in D`$. This graph is 14-regular ($`\lambda=14`$).
- If $`y\neq z`$ are neighbours of $`x`$ with $`y\not\sim z`$, then $`0`$, $`x`$ and $`x+y+z`$ are common neighbours of $`y`$ and $`z`$ in $`\Gamma`$. As $`\mu=2`$, $`x+y+z\in\{0,x\}`$; the first gives $`y+z=x\in D`$, the second $`y=z`$, both contradictions. So each point and its 14 neighbours form a 15-clique, and $`D`$ splits into two of them.
- Within a clique $`C`$, $`x\ne y`$ gives $`x+y\in D`$ adjacent to $`x`$, so $`x+y\in C`$. Hence $`C\cup0`$ is a 4-space, and $`D=(R\cup R')\setminus0`$ with $`R,R'`$ complementary.

(Our first version used Shrikhande's theorem that the lattice graph $`L_2(16)`$ is determined by these parameters; the author of search B found the direct argument above.)

A rank-6 radical $`\{a,b,a+b\}`$ with $`a\in R`$, $`b\in R'`$ would have $`a+b\notin R\cup R'`$, so each lies inside $`R`$ or inside $`R'`$. $`R`$, $`R'`$ and the $`R_i`$ form a spread, so in the coordinates above ($`R=Y`$, $`R'=X`$, $`R_1=\{(x,x)\}`$) the $`R_i`$ are the graphs of the 15 nonzero members of a spread set $`\Sigma`$, and the rank-6 radicals are lines of $`X`$ or of $`Y`$.

**Search A** runs over the same 1,128,864 spread sets, assigning a form to each of the 15 graphs, with rank-6 radicals allowed in $`X`$ or $`Y`$. **No solutions.** In 36 cases the fifteen rank-4 forms cover all their radicals within a 6-dimensional span; each of those was extended exhaustively, to every APN 8-space containing it, and **none extends**.

**Search B** uses the extra symmetry of this case: the maps $`\mathrm{diag}(h,h)`$ with $`h`$ preserving the (fixed) form on $`R_1`$, and the swap $`(x,y)\mapsto(y,x)`$, a group of order 1440. Under it the 163,072 choices for a second rank-4 form fall into 200 orbits, and each orbit representative is searched as before. **No solutions**, in three runs with different trees (orbit representatives; reversed point order; the minimal-orbit pruning switched off). It meets the same near misses as search A, six of them up to symmetry, all 6-dimensional, and not one of them accepts even a single rank-6 form. Positive controls: many known APN functions have two complementary 4-spaces that miss their rank-4 radicals and are each cut into five rank-6 radicals; moved into these coordinates, with rank-6 radicals allowed elsewhere too, they are re-found from 5- or 6-dimensional seeds, and nothing else is.

## How it was checked

- **The reduction** was re-derived by hand by the author of search B.
- **Counts predicted by theory:** of the spread sets in search A, exactly **336** are fields, matching the 336 subfields $`\mathbb F_{16}\subset M_4(\mathbb F_2)`$ ($`|GL(4,2)|/60`$). Search B's self-tests reproduce the number of nondegenerate forms (28 and 13,888), the subspace counts of $`\mathbb F_2^8`$, and a brute force over all $`2^{28}`$ alternating forms.
- **Positive controls through the real code path:** search B was run with one rule relaxed (rank-6 radicals allowed outside $`Y`$) on known APN functions with three or four rank-4 components, moved into the same coordinates. Seeded with 6 or 5 of the function's dimensions, it re-finds the known function and nothing else. Seeded with only the two starting forms, and with the extra rank-6 candidates limited to the function's own, it still climbs to the known function, so the search carries a real solution all the way from its starting point.
- **An independent recount** of the first two levels of search B in separate code (numpy, radicals by brute force) matches it exactly.

The two programs share no code. Both are exhaustive by construction, and a result resting on them is only as good as the code, which is published below. This has not been peer-reviewed; corrections are welcome.

A side observation from the controls: 7 of the 9 known functions tried have a 4-space that misses all their rank-4 radicals and is cut exactly into five rank-6 radicals. That is the structure $`R_{17}`$ would have, so known functions come close to the forbidden configuration without reaching it.

## What isn't claimed

The other open cases of Problem VI.1, with 5 to 14 rank-4 components ($`3\le i\le12`$). For 14 the partial spread still extends to a spread [2], so the leftover points are three 4-spaces, but a rank-6 radical may cross all three, and that case needs its own search (running). Below 14 the partial spread needn't extend, and the reductions above don't apply. Their status is unchanged, apart from the companion result that such a function would have no symmetry.

## Code and data

Both searches (each with its sixteen and fifteen modes), the extension program for the 6-dimensional near misses, the recount scripts and the run logs are in [code/apn8-sixteen](code/apn8-sixteen). Each search runs in a few hours or less on a desktop.

## References

1. S. H. Bénéteau, N. Goluboff, L. Kölsch, D. Vaghasiya. On the Walsh spectra of quadratic APN functions. [arXiv:2510.12008](https://arxiv.org/abs/2510.12008)
2. S. Kurz. Divisible codes. [arXiv:2112.11763](https://arxiv.org/abs/2112.11763), Table 9.3 (maximal partial spreads of solids in PG(7,2)).
