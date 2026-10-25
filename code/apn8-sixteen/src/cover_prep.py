"""Inputs for cover.cpp: covers of V minus a hole set by solids (the 4-dim radicals), one per line.

Line format: <tag> <k> <4k hex bytes: a basis of each solid> [seed <28-bit coef hex> ...]
The holes are always the nonzero points outside the solids; rank-6 radicals must lie in them.

  python cover_prep.py ext OUT    i = 4, extendable case: holes Y, X, graph(I), graph(T), one T per conjugacy class
                                  (T, T + I invertible; 5 classes by minimal polynomial; diag(g, g) fixes X, Y,
                                  graph(I) and conjugates T); covers = the other 13 graphs of each spread set
                                  containing 0, I, T (13-cliques).
  python cover_prep.py max OUT [HOLES]   i = 4, the three non-extendable hole sets of Kurz (HOLES, default
                                  runs/i4/kurz13.txt): every exact cover of the 195 other points by 13 solids.
  python cover_prep.py ctl OUT FILE LINE [NR4 [KEEP]]   positive control: known function LINE (0-based) of FILE (a LUT
                                  list); cover = its rank-4 radicals, holes = everything else; seeds = NR4 (default
                                  1) of its rank-4 forms plus forms completing the span of all its rank-4 forms to a
                                  basis (only the first KEEP of those), so cover.cpp has to find the remaining
                                  rank-4 forms (and, with KEEP short, rank-6 forms) itself.
Coordinates: point x = x1 | x2 << 4 (x1 = X-coordinates), graph(N) = {(x1, N x1)}, Y = {(0, y)}.
Forms as 28-bit coefficient vectors, bit k = B(e_i, e_j) for the k-th pair i < j in lexicographic order (as in
sixteen.cpp's to_coef).
Rowen (Claude Opus 5.5), 2026-10-24.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ID = 0x8421

def mv(m, x):
    y = 0
    for i in range(4):
        if bin((m >> (4 * i)) & 15 & x).count('1') & 1:
            y |= 1 << i
    return y

def rank_vecs(vecs):
    basis = []  # (pivot bit, vector), reduced
    for v in vecs:
        for p, b in basis:
            if v >> p & 1:
                v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [(q, c ^ v if c >> p & 1 else c) for q, c in basis]
            basis.append((p, v))
    return len(basis)

def inv4(m):
    return rank_vecs([(m >> (4 * i)) & 15 for i in range(4)]) == 4

def mmul(a, b):
    c = 0
    for j in range(4):
        col = mv(a, mv(b, 1 << j))
        for i in range(4):
            if col >> i & 1:
                c |= 1 << (4 * i + j)
    return c

def minpoly(T):
    pw = [ID]
    for _ in range(4):
        pw.append(mmul(pw[-1], T))
    for d in range(1, 5):
        for mask in range(1 << d):
            acc = pw[d]
            for k in range(d):
                if mask >> k & 1:
                    acc ^= pw[k]
            if acc == 0:
                return mask | (1 << d)
    raise ValueError

def graph_basis(N):
    return [(1 << j) | (mv(N, 1 << j) << 4) for j in range(4)]

def fmt(tag, solids, seeds=()):
    s = f"{tag} {len(solids)} " + " ".join(f"{b:02x}" for sol in solids for b in sol)
    if seeds:
        s += " seed " + " ".join(f"{c:07x}" for c in seeds)
    return s

def ext(out):
    cand = [m for m in range(65536) if inv4(m) and inv4(m ^ ID) and m != ID]
    reps = {}
    for T in cand:
        reps.setdefault(minpoly(T), T)
    total = 0
    with open(out, 'w') as f:
        for mp, T in sorted(reps.items()):
            C = [m for m in cand if m != T and inv4(m ^ T)]
            adj = []
            for a, ma in enumerate(C):
                bits = 0
                for b, mb in enumerate(C):
                    if b > a and inv4(ma ^ mb):
                        bits |= 1 << b
                adj.append(bits)
            found = []
            pick = []
            def rec(P):
                if len(pick) == 13:
                    found.append(list(pick)); return
                if bin(P).count('1') < 13 - len(pick):
                    return
                while P:
                    v = (P & -P).bit_length() - 1
                    P &= P - 1
                    pick.append(v)
                    rec(P & adj[v])
                    pick.pop()
            rec((1 << len(C)) - 1)
            print(f"T {T:04x} (min poly {mp:b}): {len(C)} candidates, {len(found)} spread sets containing 0, I, T",
                  flush=True)
            for n, pk in enumerate(found):
                f.write(fmt(f"ext-{T:04x}-{n}", [graph_basis(C[v]) for v in pk]) + "\n")
            total += len(found)
    print("lines", total)

def span_points(basis):
    pts = {0}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return pts

def basis_of(points):
    basis, sp = [], {0}
    for p in sorted(points):
        if p not in sp:
            basis.append(p)
            sp |= {q ^ p for q in sp}
    return basis

def maxi(out, holes=None):
    from i4_holes import load, solids_in
    total = 0
    with open(out, 'w') as f:
        holes = holes or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs', 'i4', 'kurz13.txt')
        for k, pts in enumerate(load(holes)):
            H = set(pts)
            comp = frozenset(range(1, 256)) - H
            sols = list(solids_in(comp))
            by_pt = {p: [s for s in sols if p in s] for p in comp}
            covers = []
            def go(left, chosen):
                if not left:
                    covers.append(list(chosen)); return
                p = min(left)
                for s in by_pt[p]:
                    if s <= left:
                        chosen.append(s); go(left - s, chosen); chosen.pop()
            go(comp, [])
            print(f"config {k}: {len(covers)} covers", flush=True)
            for n, cv in enumerate(covers):
                assert len(cv) == 13
                f.write(fmt(f"max{k}-{n}", [basis_of(s) for s in cv]) + "\n")
            total += len(covers)
    print("lines", total)

PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]

def form_coef(F, b):
    def B(x, y):
        return bin(b & (F[x ^ y] ^ F[x] ^ F[y] ^ F[0])).count('1') & 1
    return sum(B(1 << i, 1 << j) << k for k, (i, j) in enumerate(PAIRS))

def radical_of(F, b):
    D = [F[x] ^ F[0] for x in range(256)]
    rad = []
    for x in range(256):
        if all((bin(b & (F[x ^ y] ^ F[x] ^ F[y] ^ F[0])).count('1') & 1) == 0 for y in (1 << j for j in range(8))):
            rad.append(x)
    return rad

def ctl(out, path, line, nr4=1, keep=99):
    with open(path) as fh:
        for i, l in enumerate(fh):
            if i == line:
                F = [int(t) for t in l.replace('{', '').replace('}', '').replace(',', ' ').split()]
                break
    assert len(F) == 256
    r4 = []  # (coef, radical points)
    for b in range(1, 256):
        rad = radical_of(F, b)
        if len(rad) == 16:
            r4.append((form_coef(F, b), set(rad) - {0}))
    print(f"line {line}: {len(r4)} rank-4 forms")
    solids = [basis_of(r) for _, r in r4]
    seeds = [c for c, _ in r4[:nr4]]
    r4span = [c for c, _ in r4]
    for b in range(1, 256):
        c = form_coef(F, b)
        if rank_vecs(r4span + [c]) > rank_vecs(r4span):
            r4span.append(c)
            if len(seeds) < nr4 + keep:
                seeds.append(c)
    print(f"  span of rank-4 forms {rank_vecs([c for c, _ in r4])}, seeds {len(seeds)} (dim {rank_vecs(seeds)})")
    with open(out, 'w') as f:
        f.write(fmt(f"ctl-{line}", solids, seeds) + "\n")

if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'ext':
        ext(sys.argv[2])
    elif m == 'max':
        maxi(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    elif m == 'ctl':
        ctl(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else 1,
            int(sys.argv[6]) if len(sys.argv) > 6 else 99)
