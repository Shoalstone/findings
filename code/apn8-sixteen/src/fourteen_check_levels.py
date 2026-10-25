"""Independent recount of the first levels of fourteen_check's search (i = 3) for one spread configuration.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-24. Shares no code with fourteen_check.cpp: forms are built by the block
formula B = P^T d P with P = [N | I] (so rad B = graph(N)), radicals by brute force over all 256 vectors.
Usage: python fourteen_check_levels.py <spreads.txt> <index> <maxdim>
"""
import sys
from itertools import product


def mat_of_hex(h):  # u16, row r in bits 4r..4r+3 -> list of rows (bit c = column c)
    return [(h >> (4 * r)) & 15 for r in range(4)]


def apply4(M, x):
    return sum(((bin(M[r] & x).count("1") & 1) << r) for r in range(4))


def alt4_nondeg():
    out = []
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    for bits in product([0, 1], repeat=6):
        d = [[0] * 4 for _ in range(4)]
        for (i, j), b in zip(pairs, bits):
            d[i][j] = d[j][i] = b
        # Pfaffian of a 4x4 alternating matrix: d01 d23 + d02 d13 + d03 d12
        if (d[0][1] * d[2][3] + d[0][2] * d[1][3] + d[0][3] * d[1][2]) % 2:
            out.append(d)
    return out


def form_rows(N, d):
    # P: 4x8, column k of P = image of e_k: for k < 4, N e_k (x-part); for k >= 4, e_{k-4} (y-part)
    cols = [apply4(N, 1 << k) for k in range(4)] + [1 << k for k in range(4)]
    rows = []
    for i in range(8):
        r = 0
        for j in range(8):
            if i == j:
                continue
            s = 0
            for a in range(4):
                if cols[i] >> a & 1:
                    for b in range(4):
                        if cols[j] >> b & 1:
                            s ^= d[a][b]
            if s:
                r |= 1 << j
        rows.append(r)
    return tuple(rows)


def add(f, g):
    return tuple(a ^ b for a, b in zip(f, g))


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
        if line.startswith("SPREAD") and int(line.split()[1]) == idx:
            Ns = [int(t, 16) for t in line.split(":")[1].split()]
    assert Ns and len(Ns) == 14
    D = frozenset(v for v in range(1, 256) if (v & 15) == 0 or (v >> 4) == 0 or (v & 15) == (v >> 4))
    R = [frozenset((x | (apply4(mat_of_hex(N), x) << 4)) for x in range(1, 16)) for N in Ns]
    radof = {}
    for j, r in enumerate(R):
        for v in r:
            assert v not in radof and v not in D
            radof[v] = j
    assert len(radof) == 210 and len(D) == 45
    ds = alt4_nondeg()
    assert len(ds) == 28
    F4 = [[form_rows(mat_of_hex(N), d) for d in ds] for N in Ns]
    for j in range(14):
        assert len(set(F4[j])) == 28 and all(radical(f) == R[j] for f in F4[j])
    Rset = set(R)
    order = [v for v in range(1, 256) if v not in D]
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

    def rec(span, cov, d):
        counts[d] += 1
        if d == maxdim:
            return
        a = next((v for v in order if v not in cov), None)
        if a is None:
            return
        for f in F4[radof[a]]:
            ok, nc = True, set(cov)
            for s in span:
                g = add(f, s)
                if g == zero:
                    ok = False
                    break
                lg, k = legal(g)
                if not lg or (k & nc):
                    ok = False
                    break
                nc |= k
            if ok:
                rec(span + [add(f, s) for s in span], nc, d + 1)

    rec([zero], set(), 0)
    print(f"spread {idx}: nodes by dim (ascending order, off-D points first):", counts)


main()
