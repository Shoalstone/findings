"""Test the 'SHORT dim 6' spans in Rowen's runs/a/fifteen.txt against sixteen_check's legality, coordinate-free.

Decode each 28-bit code as an alternating form (bit p = pair (i, j), i < j, lexicographic, as in sixteen_check;
other orderings are tried if this one doesn't give 15 rank-4 forms). For the 63 nonzero elements: ranks, radicals,
whether the radicals are pairwise disjoint, whether every rank-6 radical lies in the hole set H (the points outside
the 15 rank-4 radicals), and whether H ∪ 0 is two complementary 4-spaces.
"""
import itertools
import re
import sys

PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]


def decode(code, pairs):
    rows = [0] * 8
    for p, (i, j) in enumerate(pairs):
        if code >> p & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return rows


def apply(rows, v):
    r = 0
    for i in range(8):
        if v >> i & 1:
            r ^= rows[i]
    return r


def radical(rows):
    return frozenset(v for v in range(1, 256) if apply(rows, v) == 0)


def analyse(codes, pairs):
    forms = [decode(c, pairs) for c in codes]
    out = []
    for m in range(1, 64):
        rows = [0] * 8
        for k in range(6):
            if m >> k & 1:
                rows = [a ^ b for a, b in zip(rows, forms[k])]
        out.append(radical(rows))
    return out


def is_subspace(s):
    s = set(s) | {0}
    return all((a ^ b) in s for a in s for b in s)


text = open(sys.argv[1], encoding="utf-8").read()
orders = {"lex i<j": PAIRS, "colex": sorted(PAIRS, key=lambda p: (p[1], p[0]))}
for line in text.splitlines():
    if not line.startswith("SHORT"):
        continue
    codes = [int(x, 16) for x in line.split("basis")[1].split()]
    for name, pairs in orders.items():
        rads = analyse(codes, pairs)
        sizes = [len(r) for r in rads]
        n4 = sizes.count(15)
        if n4 != 15:
            continue
        r4 = [r for r in rads if len(r) == 15]
        union4 = set().union(*r4)
        holes = set(range(1, 256)) - union4
        disjoint = sum(sizes) == len(set().union(*rads))
        r6 = [r for r in rads if len(r) == 3]
        bad = [s for s in sizes if s not in (0, 3, 15)]
        r6_in_holes = all(r <= holes for r in r6)
        # holes as two complementary 4-spaces?
        two = False
        for a, b in itertools.combinations(sorted(holes), 2):
            pass
        comps = []
        rest = set(holes)
        while rest:
            x = min(rest)
            comp = {x} | {y for y in holes if y != x and (x ^ y) in holes}
            comps.append(comp)
            rest -= comp
        two = len(comps) == 2 and all(len(c) == 15 and is_subspace(c) for c in comps)
        print(f"{name}: rank4 {n4}, rank6 {len(r6)}, bent {sizes.count(0)}, other {len(bad)}, radicals disjoint {disjoint},"
              f" rank-6 radicals in holes {r6_in_holes}, holes = two 4-spaces {two}, holes covered by rank-6 {sum(len(r) for r in r6)}/30")
        break
    else:
        print("no ordering gives 15 rank-4 forms:", line[:60])
