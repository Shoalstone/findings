"""Independent recount of sixteen_check's first two levels (numpy; block formula for the forms, brute-force radicals).

B_(N,M)(v, w) = (N x_v + y_v)^T M (N x_w + y_w), i.e. the 8x8 matrix [[N^T M N, N^T M], [M N, M]].
B0 = (N=0, M=J). B1 = (N=I, M). Then at the first uncovered point a2 (non-Y points first, ascending), count the
(N, M') with N x2 = y2 whose form keeps every element of the 3-dim span legal with pairwise disjoint radicals.
sixteen_check main printed 5316 / 5280 / 5480 for M codes 12 / 13 / 33.
"""
import itertools
import numpy as np

PAIRS = [(k, l) for k in range(4) for l in range(k + 1, 4)]
V = np.array([[(v >> i) & 1 for i in range(8)] for v in range(256)], dtype=np.int64)  # row v = point v


def alt4(code):
    M = np.zeros((4, 4), dtype=np.int64)
    for p, (k, l) in enumerate(PAIRS):
        if code >> p & 1:
            M[k, l] = M[l, k] = 1
    return M


def form(N, M):
    top = np.concatenate([N.T @ M @ N, N.T @ M], axis=1)
    bot = np.concatenate([M @ N, M], axis=1)
    return np.concatenate([top, bot], axis=0) % 2


def radical(B):
    """set of nonzero points v with B v = 0"""
    z = ~((V @ B) % 2).any(axis=1)
    return frozenset(int(v) for v in np.nonzero(z)[0] if v)


def in_Y(v):
    return v & 15 == 0


def legal_radical(r):
    if len(r) == 0:
        return True
    if len(r) == 15:
        return not any(in_Y(v) for v in r)
    if len(r) == 3:
        return all(in_Y(v) for v in r)
    return False


def gf2_rank(A):
    A = A.copy() % 2
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        piv = [i for i in range(r, rows) if A[i, c]]
        if not piv:
            continue
        A[[r, piv[0]]] = A[[piv[0], r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r


nondeg = [c for c in range(64) if gf2_rank(alt4(c)) == 4]
assert len(nondeg) == 28
mats = [np.array(b, dtype=np.int64).reshape(4, 4) for b in itertools.product([0, 1], repeat=16)]
GL = [N for N in mats if gf2_rank(N) == 4]
assert len(GL) == 20160
J = alt4(1 | 1 << 5)
I4 = np.eye(4, dtype=np.int64)
Z4 = np.zeros((4, 4), dtype=np.int64)
order = [v for v in range(1, 256) if not in_Y(v)] + [v for v in range(1, 256) if in_Y(v)]

B0 = form(Z4, J)
for code in (12, 13, 33):
    B1 = form(I4, alt4(code))
    span = [B0, B1, (B0 + B1) % 2]
    rads = [radical(B) for B in span]
    assert all(legal_radical(r) for r in rads)
    cov = set()
    for r in rads:
        assert not (cov & r)
        cov |= r
    a2 = next(v for v in order if v not in cov)
    x2 = np.array([(a2 >> i) & 1 for i in range(4)])
    y2 = np.array([(a2 >> (4 + i)) & 1 for i in range(4)])
    count = 0
    for N in GL:
        if not ((N @ x2) % 2 == y2).all():
            continue
        for c2 in nondeg:
            B = form(N, alt4(c2))
            ok = True
            used = set(cov)
            for S in [np.zeros((8, 8), dtype=np.int64)] + span:
                r = radical((B + S) % 2)
                if not legal_radical(r) or (used & r):
                    ok = False
                    break
                used |= r
            count += ok
    print(f"M code {code}: second point {a2:02x}, {count} dim-3 children")
