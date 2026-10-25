"""Random EA-equivalent copies of apn8 basis lines, to test that classify8's invariant really is invariant.
Input side: x -> Ax (alternating polar form B -> A^T B A); output side: a random invertible mix of the 8 basis
forms. Linear/constant terms don't appear in basis lines, so that covers EA up to the affine parts.
Usage: eacopies.py <in> <out> <copies per line> [seed]. Rowen (Claude Opus 5.5), 2026-10-24."""
import sys
import random

PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]


def to_mat(c):
    B = [[0] * 8 for _ in range(8)]
    for k, (i, j) in enumerate(PAIRS):
        if c >> k & 1:
            B[i][j] = B[j][i] = 1
    return B


def from_mat(B):
    return sum(1 << k for k, (i, j) in enumerate(PAIRS) if B[i][j])


def mul(X, Y):
    return [[sum(X[i][k] & Y[k][j] for k in range(8)) & 1 for j in range(8)] for i in range(8)]


def rank(rows):
    rows = list(rows); r = 0
    for bit in range(8):
        p = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def rand_gl(rng):
    while True:
        rows = [rng.randrange(256) for _ in range(8)]
        if rank(rows) == 8:
            return [[rows[i] >> j & 1 for j in range(8)] for i in range(8)]


def main():
    src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    out = []
    for ln, line in enumerate(open(src), 1):
        if "basis" not in line:
            continue
        head, tail = line.split("basis")
        basis = [int(h, 16) for h in tail.split()[:8]]
        for c in range(n):
            A = rand_gl(rng); At = [list(r) for r in zip(*A)]
            forms = [from_mat(mul(mul(At, to_mat(b)), A)) for b in basis]
            G = rand_gl(rng)
            mixed = []
            for s in range(8):
                v = 0
                for t in range(8):
                    if G[s][t]:
                        v ^= forms[t]
                mixed.append(v)
            out.append(f"{head.strip()} copy{ln}.{c} basis " + " ".join(f"{v:07x}" for v in mixed) + "\n")
    open(dst, "w").writelines(out)


main()
