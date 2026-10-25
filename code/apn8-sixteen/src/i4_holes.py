"""Hole sets for 13 rank-4 components (i = 4): check Kurz's three maximal-partial-spread hole configurations.

Kurz, Divisible codes (arXiv 2112.11763) §9 prints three 8x60 generator matrices; columns are the holes of a maximal
partial spread of 13 solids in PG(7,2). For each: 60 distinct nonzero points, every hyperplane meets it in 36 or 28
points (weights 24/32), the complement (195 points) has how many solids inside it and how many exact covers by 13
of them, and whether the holes contain a solid / a line partition (needed: the 20 rank-6 radicals).
usage: python i4_holes.py runs/i4/kurz13.txt
Rowen (Claude Opus 5.5), 2026-10-24.
"""
import sys
from itertools import combinations

def load(path):
    rows = [l.strip() for l in open(path) if l.strip()]
    assert len(rows) % 8 == 0
    mats = [rows[k:k + 8] for k in range(0, len(rows), 8)]
    out = []
    for m in mats:
        pts = []
        for c in range(60):
            v = 0
            for r in range(8):
                if m[r][c] == '1':
                    v |= 1 << r
            pts.append(v)
        out.append(pts)
    return out

def parity(x):
    return bin(x).count('1') & 1

def solids_in(pointset):
    """All 4-dim subspaces whose 15 nonzero points lie in pointset (a set of ints)."""
    pts = sorted(pointset)
    found = set()
    for a in pts:
        for b in pts:
            if b <= a or (a ^ b) not in pointset:
                continue
            plane = {a, b, a ^ b}
            for c in pts:
                if c in plane or c < b:
                    continue
                if not all((c ^ p) in pointset for p in plane):
                    continue
                s3 = plane | {c ^ p for p in plane} | {c}
                for d in pts:
                    if d in s3 or d < c:
                        continue
                    if all((d ^ p) in pointset for p in s3):
                        s4 = frozenset(s3 | {d ^ p for p in s3} | {d})
                        found.add(s4)
    return found

def exact_covers(universe, sets):
    """Count exact covers of universe by members of sets (Algorithm X, small)."""
    by_pt = {p: [s for s in sets if p in s] for p in universe}
    count = 0
    def go(left):
        nonlocal count
        if not left:
            count += 1
            return
        p = min(left, key=lambda q: sum(1 for s in by_pt[q] if s <= left))
        for s in by_pt[p]:
            if s <= left:
                go(left - s)
    go(frozenset(universe))
    return count

def lines_in(pointset):
    return {frozenset((a, b, a ^ b)) for a, b in combinations(sorted(pointset), 2) if (a ^ b) in pointset}

def line_partitions(pointset, cap=10**6):
    ls = list(lines_in(pointset))
    by_pt = {p: [l for l in ls if p in l] for p in pointset}
    count = 0
    def go(left):
        nonlocal count
        if count >= cap:
            return
        if not left:
            count += 1
            return
        p = min(left, key=lambda q: sum(1 for l in by_pt[q] if l <= left))
        for l in by_pt[p]:
            if l <= left:
                go(left - l)
    go(frozenset(pointset))
    return count, len(ls)

def main():
    for k, pts in enumerate(load(sys.argv[1])):
        H = set(pts)
        print(f"config {k}: {len(pts)} columns, {len(H)} distinct, zero in it: {0 in H}")
        ws = {}
        for u in range(1, 256):
            n = sum(1 for p in H if parity(u & p) == 0)
            ws[n] = ws.get(n, 0) + 1
        print("  points per hyperplane:", dict(sorted(ws.items())))
        comp = set(range(1, 256)) - H
        sol = solids_in(comp)
        print(f"  solids inside the complement: {len(sol)}")
        print(f"  exact covers of the complement by 13 solids: {exact_covers(comp, sol)}")
        print(f"  solids inside the holes: {len(solids_in(H))}")
        n, nl = line_partitions(H)
        print(f"  lines inside the holes: {nl}; partitions of the holes into 20 lines: {n}")

if __name__ == '__main__':
    main()
