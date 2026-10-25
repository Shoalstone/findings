"""Reviewer check (Claude Opus 5.5, 2026-10-24), independent of c3k4.cpp's code.

c3k4 search2, (4,2) 3+3: the second 3-cycle h' + F4* g' is enumerated with h' != h1 (`if (c == h1) continue;`).
Question: can h1 + g be a rank-4 form whose A-orbit radicals avoid B1's, for g in N outside F4 n1, with
P = F4 n1 + F4 g anisotropic? If so, that case (seed P + h1, r = 3) is not covered by search2.

Coordinates as in c3k4: V = F_2^8, x = u + w v (u = bits 0..3, v = bits 4..7), A(u, v) = (v, u + v).
Forms: 8 row bitmasks of a symmetric zero-diagonal matrix. (B o g)(x, y) = B(gx, gy).
"""
import itertools

def A(x):
    u, v = x & 15, x >> 4
    return v | ((u ^ v) << 4)

def form_from_pairs(pairs):
    m = [0] * 8
    for i, j in pairs:
        m[i] ^= 1 << j; m[j] ^= 1 << i
    return tuple(m)

def bval(m, x, y):
    r = 0
    for i in range(8):
        if x >> i & 1:
            r ^= m[i]
    return bin(r & y).count("1") & 1

def compose(m, g):  # B o g
    cols = [g(1 << i) for i in range(8)]
    out = [0] * 8
    for i in range(8):
        for j in range(8):
            if bval(m, cols[i], cols[j]):
                out[i] |= 1 << j
    return tuple(out)

def add(a, b):
    return tuple(x ^ y for x, y in zip(a, b))

ZERO = (0,) * 8

def rank(m):
    rows = list(m); r = 0
    for bit in range(8):
        p = next((i for i in range(r, 8) if rows[i] >> bit & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(8):
            if i != r and rows[i] >> bit & 1:
                rows[i] ^= rows[r]
        r += 1
    return r

def radical(m):
    return frozenset(x for x in range(1, 256) if all(bval(m, x, 1 << j) == 0 for j in range(8)))

def act(m):
    return compose(m, A)

def fixpart(m):
    a = act(m)
    return add(add(m, a), act(a))

# all alternating forms as basis
basis = [form_from_pairs([(i, j)]) for i in range(8) for j in range(i + 1, 8)]
B1 = form_from_pairs([(4, 5), (6, 7)])
h1 = fixpart(B1)
n1 = add(B1, h1)
assert rank(B1) == 4 and rank(n1) == 8
RB = [radical(B1), radical(act(B1)), radical(act(act(B1)))]
assert all(len(r) == 15 for r in RB)

# N as the image of m -> m + fixpart(m); enumerate all its elements
def span(vecs):
    S = {ZERO}
    for v in vecs:
        S |= {add(s, v) for s in S}
    return S

Ngen = [add(b, fixpart(b)) for b in basis]
N = span(Ngen)
Fix = span([fixpart(b) for b in basis])
print("dim N", len(N).bit_length() - 1, "dim Fix", len(Fix).bit_length() - 1)
F4n1 = {ZERO, n1, act(n1), add(n1, act(n1))}

def f4line(g):
    return {ZERO, g, act(g), add(g, act(g))}

bad = []
for g in N:
    if g in F4n1:
        continue
    b = add(h1, g)
    if rank(b) != 4:
        continue
    ok = True
    bb = b
    for _ in range(3):
        R = radical(bb)
        if any(R & r for r in RB):
            ok = False; break
        bb = act(bb)
    if not ok:
        continue
    P = span([n1, act(n1), g, act(g)])
    aniso = all(rank(p) == 8 for p in P if p != ZERO)
    bad.append((g, aniso))
print("g in N \\ F4 n1 with h1 + g rank 4, orbit radicals avoiding B1's:", len(bad),
      "of which P anisotropic:", sum(1 for _, a in bad if a))
# also: how many rank-4 h1 + g at all, and for the record the rank of h1
print("rank h1", rank(h1))
cnt = sum(1 for g in N if g not in F4n1 and rank(add(h1, g)) == 4)
print("g outside F4 n1 with rank(h1 + g) = 4:", cnt)
