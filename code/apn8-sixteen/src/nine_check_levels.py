"""Independent recount of the first levels of nine_check's search (i = 8) for one configuration.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-28; adapted from ten_check_levels.py (9 members, 120 holes).
Shares no code with nine_check.cpp: a form with radical R is
B(x, y) = pi(x)^T d pi(y), pi = the last four coordinates in a basis (basis of R, completion), d over the 28
nondegenerate 4x4 alternating matrices (by their Pfaffian); radicals by brute force over all 256 vectors.
Usage: python nine_check_levels.py <configs.txt> <index> <maxdim>
Prints nodes by span dimension, ascending point order with off-D points first (as nine_check's `asc`).
"""
import sys
from itertools import product


def apply4(h, x):  # u16 matrix, row r in bits 4r..4r+3
    return sum(((bin(((h >> (4 * r)) & 15) & x).count("1") & 1) << r) for r in range(4))


def span(vs):
    s = {0}
    for v in vs:
        s |= {a ^ v for a in s}
    return s


def alt4_nondeg():
    out = []
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    for bits in product([0, 1], repeat=6):
        d = [[0] * 4 for _ in range(4)]
        for (i, j), b in zip(pairs, bits):
            d[i][j] = d[j][i] = b
        if (d[0][1] * d[2][3] + d[0][2] * d[1][3] + d[0][3] * d[1][2]) % 2:
            out.append(d)
    return out


def forms_with_radical(R, ds):
    pts = sorted(R)
    basis = []
    for p in pts:
        if p not in span(basis):
            basis.append(p)
    assert len(basis) == 4
    for v in range(1, 256):
        if len(basis) == 8:
            break
        if v not in span(basis):
            basis.append(v)
    coord = {}  # vector -> coordinate bits in this basis
    for z in range(256):
        v = 0
        for k in range(8):
            if z >> k & 1:
                v ^= basis[k]
        coord[v] = z
    pi = [coord[1 << i] >> 4 for i in range(8)]  # last four coordinates of e_i
    out = []
    for d in ds:
        rows = []
        for i in range(8):
            r = 0
            for j in range(8):
                s = 0
                for a in range(4):
                    if pi[i] >> a & 1:
                        for b in range(4):
                            if pi[j] >> b & 1:
                                s ^= d[a][b]
                if s:
                    r |= 1 << j
            rows.append(r)
        out.append(tuple(rows))
    return out


def radical(f):
    ker = []
    for v in range(1, 256):
        acc = 0
        for i in range(8):
            if v >> i & 1:
                acc ^= f[i]
        if acc == 0:
            ker.append(v)
    return frozenset(ker)


def main():
    path, idx, maxdim = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    Ns = None
    for line in open(path):
        if line.startswith("CONF") and int(line.split()[1]) == idx:
            Ns = [int(t, 16) for t in line.split(":")[1].split()]
    assert Ns and len(Ns) == 6
    X = frozenset(range(1, 16))
    Y = frozenset(x << 4 for x in range(1, 16))
    mem = [X, Y] + [frozenset(x | (apply4(N, x) << 4) for x in range(1, 16)) for N in [0x8421] + Ns]
    radof = {}
    for j, r in enumerate(mem):
        for v in r:
            assert v not in radof
            radof[v] = j
    D = frozenset(v for v in range(1, 256) if v not in radof)
    assert len(radof) == 135 and len(D) == 120
    ds = alt4_nondeg()
    assert len(ds) == 28
    F4 = [forms_with_radical(r, ds) for r in mem]
    for j in range(9):
        assert len(set(F4[j])) == 28 and all(radical(f) == mem[j] for f in F4[j])
    Rset = set(mem)
    order = sorted(radof)
    cache = {}

    def legal(f):
        if f in cache:
            return cache[f]
        k = radical(f)
        if not k:
            res = (True, k)
        elif len(k) == 15:
            res = (k in Rset, k)
        elif len(k) == 3:
            res = (k <= D, k)
        else:
            res = (False, k)
        if len(cache) < 2_000_000:
            cache[f] = res
        return res

    counts = [0] * (maxdim + 1)
    zero = (0,) * 8

    def rec(sp, cov, d):
        counts[d] += 1
        if d == maxdim:
            return
        a = next((v for v in order if v not in cov), None)
        assert a is not None, "all off-D points covered below maxdim: the recount doesn't branch in D"
        for f in F4[radof[a]]:
            ok, nc = True, set(cov)
            for s in sp:
                g = tuple(x ^ y for x, y in zip(f, s))
                if g == zero:
                    ok = False
                    break
                lg, k = legal(g)
                if not lg or (k & nc):
                    ok = False
                    break
                nc |= k
            if ok:
                rec(sp + [tuple(x ^ y for x, y in zip(f, s)) for s in sp], nc, d + 1)

    rec([zero], set(), 0)
    print(f"config {idx}: nodes by dim (ascending order, off-D points first):", counts)


main()
