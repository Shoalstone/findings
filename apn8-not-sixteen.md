# No quadratic APN function in dimension 8 has more than twelve rank-4 components

*APN functions · dimension 8 · 2026-10-25 · settles four cases of an open problem*

## In plain words

> A 2026 paper listed the shapes a quadratic APN function on 8 bits could still have, and asked which of them actually occur. It singled out one as particularly interesting: sixteen components of a special low-rank kind ("rank-4"), the most the list allows. It would make a certain set as small as such sets can be. This note shows that shape doesn't occur, and neither do the next three down: fifteen, fourteen or thirteen. Every quadratic APN function in dimension 8 has at most twelve rank-4 components.

The radicals of the rank-4 components are 4-dimensional pieces of $`\mathbb F_2^8`$ that pairwise meet only in $`0`$, a *partial spread*. With thirteen or more of them the structure is rigid enough to classify completely, which makes a complete computer search small. For each case, two searches, written separately and organised differently, find nothing.

*First version 2026-10-24 (sixteen and fifteen); fourteen and thirteen added 2026-10-25.*

## The question

Notation as in the [companion note](apn8-no-symmetry.md): for a quadratic APN $`F`$ on $`\mathbb F_2^8`$, $`S=\{B_b\}`$ is the 8-dimensional space of the alternating forms $`B_b(x,y)=b\cdot(F(x+y)+F(x)+F(y)+F(0))`$, and the radicals of its nonzero forms partition $`\mathbb F_2^8\setminus 0`$.

Bénéteau, Goluboff, Kölsch and Vaghasiya [1, Thm IV.5] showed the amplitude distribution is $`[0^{190},2^{64},6^1]`$ or $`[0^{238-4i},2^{5i},4^{17-i}]`$ with $`1\le i\le 17`$, and asked in Open Problem VI.1 whether $`1\le i\le 12`$ occur. Of these they single out $`i=1`$, $`[0^{234},2^5,4^{16}]`$, where the set of non-bent components would be a smallest possible nontrivial blocking set.

> **Theorem.** No quadratic APN function $`F:\mathbb F_2^8\to\mathbb F_2^8`$ has amplitude distribution $`[0^{238-4i},2^{5i},4^{17-i}]`$ for $`i=1,2,3,4`$, that is $`[0^{234},2^5,4^{16}]`$, $`[0^{230},2^{10},4^{15}]`$, $`[0^{226},2^{15},4^{14}]`$ or $`[0^{222},2^{20},4^{13}]`$. Together with [1, Thm IV.1] (which rules out 17), every quadratic APN function on $`\mathbb F_2^8`$ has at most 12 rank-4 components.

In every case the $`17-i`$ rank-4 radicals are pairwise complementary 4-spaces, since the radicals of the nonzero forms partition the nonzero points and two radicals meeting in a point would share it. The $`5i`$ rank-6 radicals are lines covering the $`15i`$ points left over, the *holes*. A hyperplane contains at most one of the 4-spaces (two complementary ones span everything), and counting points in hyperplanes gives the holes' character values: with $`D`$ the set of holes and $`W(u)=\sum_{d\in D}(-1)^{u\cdot d}`$, $`W(u)=16-i`$ for the $`15i`$ hyperplanes $`u^\perp`$ that contain no radical and $`W(u)=-i`$ for the rest. Everything turns on what the hole set can be.

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

## Fourteen

Now 14 rank-4 radicals and 15 rank-6 radicals covering $`D`$, the 45 holes. The character values ($`13`$ forty-five times, $`-3`$ otherwise) make $`D`$ a projective two-weight set, but unlike fifteen they don't pin it down: there are two such sets [2, §11].

**The partial spread extends.** Kurz's computer classification of maximal partial spreads of solids in $`PG(7,2)`$ [2, Table 9.3] has none of sizes 14, 15 or 16, so the 14 radicals extend to a spread and $`D`$ is three pairwise complementary 4-spaces $`H_1,H_2,H_3`$. Search B's author re-proved this by computation rather than citing it: put three of the 14 radicals at $`X`$, $`Y`$, $`\{(x,x)\}`$; the other 11 are graphs of an 11-clique in the graph on the 5,824 matrices above. Conjugation by $`GL(4,2)`$ moves one member of the smallest conjugacy class present to a fixed representative, which leaves 863,186 *conjugation-reduced* cliques (each standing for an exactly known number of cliques). Every one has three more matrices completing a spread set. As a count check, the weighted number of 11-cliques is $`410{,}906{,}496=\binom{14}{11}\cdot 1{,}128{,}864`$, so each lies in exactly one spread.

In search B's coordinates, $`H_1,H_2,H_3=Y,X,\{(x,x)\}`$, the lines inside $`D`$ are the 105 lines inside one of the three 4-spaces and 15 transversals $`\{(x,0),(0,x),(x,x)\}`$, one through each hole (in search A's coordinates below they are $`\{(x,0),(0,Nx),(x,Nx)\}`$). So unlike fifteen, a rank-6 radical can cross all three holes, and both searches allow every one of the 120 lines.

**Search A** (`sixteen.cpp`, mode `i3`) puts two holes at $`Y`$ and $`X`$ and one rank-4 radical at $`\{(x,x)\}`$, with its form fixed as in fifteen. For every one of the 1,128,864 spread sets it then tries each of the 14 remaining graphs as the third hole, 15.8 million cases in all, with no further symmetry used. **No solutions**, and no case where the 14 rank-4 forms cover their radicals in fewer than 8 dimensions. 15 hours, sharing the machine.

**Search B** (`fourteen_check.cpp`) classifies the configurations first: a spread with an unordered triple of hole solids, up to $`GL(8,2)`$. The residual group (conjugation by $`GL(4,2)`$ times the $`S_3`$ permuting the three holes, order 120,960) leaves **73 configurations**, whose orbit sizes sum to 1,128,864. On each, a depth-first search branches on the first uncovered point over every form that could cover it: the 28 forms on the radical through it, or, for a hole, the 13,888 rank-6 forms on each line of $`D`$ through it. **No solutions**, in four runs with different trees (ascending, descending and random point order, and all 2,873 conjugation-reduced spreads without the $`S_3`$ reduction). No span ever reached dimension 7, and no node ever had all 210 points outside $`D`$ covered, so the rank-6 part of the problem never came into play.

## Thirteen

Thirteen rank-4 radicals, 20 rank-6 radicals, 60 holes, with character values $`12`$ sixty times and $`-4`$ otherwise: a projective $`[60,8,\{24,32\}]`$ two-weight set. Here the partial spread need not extend. Kurz reports, by computer [2, §9 after Table 9.3], that exactly six hole sets occur. Three are unions of four disjoint solids (the partial spread extends to a spread), and three, printed there, are holes of maximal partial spreads of 13 solids. Search B's author proved this classification by computation as well, from all 2,358,714 conjugation-reduced 10-cliques (three radicals placed at $`X`$, $`Y`$, $`\{(x,x)\}`$ as before). The extendable ones number exactly $`\binom{14}{10}\cdot 1{,}128{,}864`$ when weighted, each lying in exactly one spread. Kurz's three sets have automorphism groups of orders 576, 14,400 and 96, and their complements split into 13 solids in 552, 10,440 and 12 ways, so they account for $`\sum 1716\cdot 20160\cdot\mathrm{covers}/|\mathrm{Aut}\,H| = 33{,}153{,}120+25{,}081{,}056+4{,}324{,}320=62{,}558{,}496`$ non-extendable cliques, exactly the weighted count found; and every one of the 141,165 conjugation-reduced non-extendable cliques canonicalises onto one of Kurz's covers.

**Search A** (`cover.cpp`, inputs from `cover_prep.py`) takes the two kinds in turn. *Non-extendable:* Kurz's three hole sets (checked by `i4_holes.py`), and all 11,004 ways of cutting their complements into 13 solids. *Extendable:* holes $`Y`$, $`X`$, $`\{(x,x)\}`$ and $`\{(x,Tx)\}`$. Any four disjoint solids can be moved there with $`T`$ up to conjugation, which leaves one $`T`$ for each of the 5 conjugacy classes with $`T`$ and $`T+I`$ invertible. The radicals are the other 13 graphs of each spread set through $`0`$, $`I`$ and $`T`$: 25,850 cases, which, weighted by class size, total $`14\cdot 1{,}128{,}864`$ as they should. The search branches on points as search B does, with no symmetry beyond this (all 28 forms on the first radical). **No solutions** on either side, with nothing short of dimension 8 covering every point. 1 hour and 2 hours.

**Search B** (`thirteen_check.cpp`) works from its own classification: **190 configurations** up to $`GL(8,2)`$ (159 extendable, 31 not), found by canonical forms over frames of the partial spread itself, whose orbit sizes sum exactly to the weighted counts above. **No solutions** in five runs: ascending, descending and random point order, then two that bypass it, one per side: the extendable side from the 73 configurations for fourteen with each of their 14 radicals in turn as the fourth hole, and the non-extendable side from all 11,004 of Kurz's covers in Kurz's own coordinates. Nothing reached dimension 7. Four 6-dimensional spans, in two configurations whose holes are the most symmetric of the three four-solid types, do hold all 13 rank-4 forms; every rank-6 form on every hole line through the first uncovered hole is illegal there (486,080 forms on 35 lines), so they die at once.

The two searches meet on the non-extendable side. Search B's unreduced run over Kurz's 11,004 covers visits 374,909,300 nodes, exactly search A's 374,898,296 plus one root node for each cover: two independently written programs walking the same tree.

## How it was checked

- **Who wrote what.** Search A is Rowen's throughout. Each search B was written by a separate Claude Opus 5.5 instance that didn't read search A until its own verdict was in (one instance for sixteen and fifteen, a fresh one for fourteen, another for thirteen, which reused the fourteen checker's utilities). The two sides share no code.
- **The reductions** were re-derived by hand by the author of search B for sixteen and fifteen. For fourteen and thirteen the classification of the hole sets, which search A takes from Kurz, was re-proved by search B's authors by computation, with the count identities above. Kurz's three printed hole sets for thirteen were checked from scratch (60 points, two weights, no solid inside, automorphism group orders 576, 14,400 and 96, as he gives them).
- **Counts predicted by theory:** of the spread sets in search A, exactly **336** are fields, matching the 336 subfields $`\mathbb F_{16}\subset M_4(\mathbb F_2)`$ ($`|GL(4,2)|/60`$), in each of the three runs over all spread sets. Search B for fourteen reaches the same count, 1,128,864, independently. The self-tests reproduce the number of nondegenerate forms (28 and 13,888), the subspace counts of $`\mathbb F_2^8`$, the 120 and 200 lines inside the hole sets, and a brute force over all $`2^{28}`$ alternating forms.
- **Positive controls through the real code path:** every search B, and `cover.cpp`, was run with one rule relaxed (a known function's own rank-6 radicals allowed outside the holes) on known APN functions moved into its coordinates. For sixteen and fifteen, seeded with 6 or 5 of the function's dimensions, or with only the two starting forms, search B re-finds the known function and nothing else. For fourteen and thirteen, search B was run from the empty span on every function with at least two rank-4 components in five lists [3–7] (7,732 functions, several placements each): **13,086 of 13,086** runs re-found the function for fourteen and **18,674 of 18,674** for thirteen, with no other space ever found. Between them the thirteen controls passed through every one of the 190 configurations. `cover.cpp` re-finds a known function from one rank-4 form plus 4, 3 or 2 of its dimensions, branching through one and two levels of rank-6 forms. `sixteen.cpp` (search A for sixteen, fifteen and fourteen) places only rank-4 forms, and known functions have at most four, so it has no such control; its checks are the counts above, and in fifteen search B independently reaches exactly the 6-dimensional near misses that `sixteen.cpp` reports.
- **Independent recounts** of the first levels of each search B in separate Python code (radicals by brute force) match exactly.
- **Cross-checks between the cases:** `cover.cpp` run on the 73 configurations for fourteen also finds nothing (a check of `cover.cpp`, not an independent one), and on the non-extendable side of thirteen the two searches' node counts agree exactly.

All searches are exhaustive by construction, and a result resting on them is only as good as the code, which is published below. This has not been peer-reviewed; corrections are welcome.

A side observation from the controls: 7 of the 9 known functions tried for sixteen have a 4-space that misses all their rank-4 radicals and is cut exactly into five rank-6 radicals. That is the structure $`R_{17}`$ would have, so known functions come close to the forbidden configuration without reaching it. The near misses tell the same story from the other side: in fifteen and thirteen, a few 6-dimensional spans carry every rank-4 form, but none can cover its first uncovered hole with any rank-6 form; in fourteen, nothing gets even that far.

## What isn't claimed

The other open cases of Problem VI.1, with 5 to 12 rank-4 components ($`5\le i\le12`$). For twelve the hole classification is still small: about 5.3 million conjugation-reduced 9-cliques, and roughly a thousand configurations expected. But some partial spreads of 12 solids extend to a spread in more than one way, and with fewer rank-4 forms the rank-6 branching will dominate the search. Their status is unchanged, apart from the companion result that such a function would have no symmetry.

## Code and data

All the searches, the extension program for the 6-dimensional near misses, the recount scripts, Kurz's hole sets and the run logs are in [code/apn8-sixteen](code/apn8-sixteen), with a table of commands and expected outputs. Each search B runs in minutes and each search A in hours on a desktop.

## References

1. S. H. Bénéteau, N. Goluboff, L. Kölsch, D. Vaghasiya. On the Walsh spectra of quadratic APN functions. [arXiv:2510.12008](https://arxiv.org/abs/2510.12008)
2. S. Kurz. Divisible codes. [arXiv:2112.11763](https://arxiv.org/abs/2112.11763), §9 and Table 9.3 (maximal partial spreads of solids in PG(7,2)), §11 (the two-weight codes).
3. C. Beierle, G. Leander. New instances of quadratic APN functions. IEEE Trans. Inf. Theory 68(1):670–678, 2022. [arXiv:2009.07204](https://arxiv.org/abs/2009.07204) (list: Zenodo 4738942)
4. C. Beierle, G. Leander, L. Perrin. Trims and extensions of quadratic APN functions. Des. Codes Cryptogr. 90(4):1009–1036, 2022. [arXiv:2108.13280](https://arxiv.org/abs/2108.13280) (list: Zenodo 5821018)
5. C. Beierle, P. Langevin, G. Leander, A. Polujan, S. Rasoolzadeh. Millions of inequivalent quadratic APN functions in eight variables. [arXiv:2508.04644](https://arxiv.org/abs/2508.04644) (list: Zenodo 16752428)
6. Y. Yu, L. Perrin. Constructing more quadratic APN functions with the QAM method. Cryptogr. Commun. 14:1359–1369, 2022 (list as collected in sboxU's apnDB)
7. Y. Yu, M. Wang, Y. Li. A matrix approach for constructing quadratic APN functions. Des. Codes Cryptogr. 73:587–600, 2014 (list as collected in sboxU's apnDB, with Edel and Pott's functions)
