"""Independent recount of sixteen_check's first branching for the case i = 2 (holes X and Y).

Forms from the block formula [[N^T M N, N^T M], [M N, M]] (radical graph(N)), radicals by brute force.
B1 = (I, J); B2 = (N, M) from `sixteen_check i2count j`. Rules: rank-4 radicals miss X and Y, rank-6 radicals lie
inside X or inside Y, radicals pairwise disjoint. Count the (N', M') covering the first uncovered point outside
X ∪ Y (ascending) that keep the 3-dim span legal. Usage: python sixteen_check_level2_i2.py <N hex> <M code> ...
"""
import itertools
import sys
import numpy as np

PAIRS = [(k, l) for k in range(4) for l in range(k + 1, 4)]
V = np.array([[(v >> i) & 1 for i in range(8)] for v in range(256)], dtype=np.int64)


def alt4(code):
    M = np.zeros((4, 4), dtype=np.int64)
    for p, (k, l) in enumerate(PAIRS):
        if code >> p & 1:
            M[k, l] = M[l, k] = 1
    return M


def mat_from_hex(h):
    return np.array([[(h >> (4 * r + c)) & 1 for c in range(4)] for r in range(4)], dtype=np.int64)


def form(N, M):
    return np.block([[N.T @ M @ N, N.T @ M], [M @ N, M]]) % 2


def radical(B):
    z = ~((V @ B) % 2).any(axis=1)
    return frozenset(int(v) for v in np.nonzero(z)[0] if v)


def in_X(v):
    return v >> 4 == 0


def in_Y(v):
    return v & 15 == 0


def legal(r):
    if len(r) == 0:
        return True
    if len(r) == 15:
        return not any(in_X(v) or in_Y(v) for v in r)
    if len(r) == 3:
        return all(in_X(v) for v in r) or all(in_Y(v) for v in r)
    return False


def gf2_rank(A):
    A = A.copy() % 2
    r = 0
    for c in range(A.shape[1]):
        piv = [i for i in range(r, A.shape[0]) if A[i, c]]
        if not piv:
            continue
        A[[r, piv[0]]] = A[[piv[0], r]]
        for i in range(A.shape[0]):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r


nondeg = [c for c in range(64) if gf2_rank(alt4(c)) == 4]
GL = [np.array(b, dtype=np.int64).reshape(4, 4) for b in itertools.product([0, 1], repeat=16)]
GL = [N for N in GL if gf2_rank(N) == 4]
I4 = np.eye(4, dtype=np.int64)
J = alt4(1 | 1 << 5)
order = [v for v in range(1, 256) if not in_X(v) and not in_Y(v)]

args = sys.argv[1:]
for k in range(0, len(args), 2):
    N2, M2 = mat_from_hex(int(args[k], 16)), alt4(int(args[k + 1]))
    span = [form(I4, J), form(N2, M2)]
    span.append((span[0] + span[1]) % 2)
    cov = set()
    for B in span:
        r = radical(B)
        assert legal(r) and not (cov & r)
        cov |= r
    a = next(v for v in order if v not in cov)
    x = np.array([(a >> i) & 1 for i in range(4)])
    y = np.array([(a >> (4 + i)) & 1 for i in range(4)])
    count = 0
    for N in GL:
        if not ((N @ x) % 2 == y).all():
            continue
        for c in nondeg:
            B = form(N, alt4(c))
            used, ok = set(cov), True
            for S in [np.zeros((8, 8), dtype=np.int64)] + span:
                r = radical((B + S) % 2)
                if not legal(r) or (used & r):
                    ok = False
                    break
                used |= r
            count += ok
    print(f"N {args[k]} M {args[k + 1]}: point {a:02x}, {count} dim-3 children")
