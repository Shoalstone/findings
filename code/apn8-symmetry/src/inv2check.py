# inv2check: checks the order-2 structure claims (README "p = 2") on a quadratic APN function with an involution A.
# python apn/inv2check.py <seed line file> <A cols hex x8>   (A e_i as in aut8's AUT8_PRINT output)
# Independent of the C++ code: plain Python on the look-up table.
import sys, itertools

def load_lut(path):
    s = open(path).read()
    s = s[s.index('{') + 1:s.index('}')]
    return [int(t) for t in s.split(',')]

def apply(cols, x):
    y = 0
    for i in range(8):
        if x >> i & 1: y ^= cols[i]
    return y

def rank_gf2(vecs):
    piv = {}
    r = 0
    for v in vecs:
        for b in sorted(piv, reverse=True):
            if v >> b & 1: v ^= piv[b]
        if v:
            piv[v.bit_length() - 1] = v; r += 1
            # keep reduced: not needed for rank
    return r

def span(vecs):
    out = {0}
    for v in vecs:
        out |= {x ^ v for x in out}
    return out

lut = load_lut(sys.argv[1])
A = [int(t, 16) for t in sys.argv[2:10]]
pairs = [(i, j) for i in range(8) for j in range(i + 1, 8)]

def comp(c):  # component c.F as a function table
    return [bin(c & lut[x]).count('1') & 1 for x in range(256)]

def form_of(f):  # alternating form as 28-bit vector over pairs (e_i, e_j)
    v = 0
    for k, (i, j) in enumerate(pairs):
        x, y = 1 << i, 1 << j
        if f[x ^ y] ^ f[x] ^ f[y] ^ f[0]: v |= 1 << k
    return v

def B(v, x, y):  # form value from its 28-bit vector
    s = 0
    for k, (i, j) in enumerate(pairs):
        if v >> k & 1:
            s ^= ((x >> i & 1) & (y >> j & 1)) ^ ((x >> j & 1) & (y >> i & 1))
    return s

def act(v):  # (A.B)(x, y) = B(Ax, Ay)
    w = 0
    for k, (i, j) in enumerate(pairs):
        if B(v, A[i], A[j]): w |= 1 << k
    return w

def radical(v):
    return [x for x in range(1, 256) if all(B(v, x, 1 << t) == 0 for t in range(8))]

forms = [form_of(comp(1 << i)) for i in range(8)]
S8 = span(forms)
assert len(S8) == 256, "components not independent"
assert all(act(v) in S8 for v in forms), "A doesn't preserve S8"
Kv = [x for x in range(1, 256) if apply(A, x) == x]
N = [A[i] ^ (1 << i) for i in range(8)]
j = rank_gf2(N)
print("j = rank N =", j, " dim K =", rank_gf2(Kv))
F = [v for v in S8 if v and act(v) == v]
D = {v ^ act(v) for v in S8}
print("dim S8^A =", (len(F) + 1).bit_length() - 1, " dim D(S8) =", (len(D)).bit_length() - 1, " D(S8) == S8^A:", D == set(F) | {0})
# b_C(u, v) = C(eps u, v) on W = complement of K (generators); check the radical sets in K
Kset = set(Kv)
rad_in_K = {}
for C in F:
    rk = frozenset(x for x in radical(C) if x in Kset)
    rad_in_K[C] = rk
sizes = sorted(len(r) for r in rad_in_K.values())
print("fixed forms by |rad & K \\ 0|:", {s: sizes.count(s) for s in set(sizes)})
cover = [x for r in rad_in_K.values() for x in r]
print("K\\0 covered exactly once by the fixed radicals:", sorted(cover) == sorted(Kv))
for C in F:
    r = radical(C)
    print("  fixed form rank", 8 - (len(r) + 1).bit_length() + 1, " |rad & K| =", len(rad_in_K[C]))
# rank-4 components and their orbits
r4 = [v for v in S8 if v and len(radical(v)) == 15]
print("rank-4 components:", len(r4), " fixed:", sum(act(v) == v for v in r4))
