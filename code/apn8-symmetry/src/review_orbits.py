"""Independent reviewer check of inv2's orbit reduction (subagent, Claude Opus 5.5, 2026-10-24).

Works directly with 8x8 alternating matrices, not inv2's tables. V = F_2^8, bit i = e_i, W = <e0..e3>,
eps e_i = e_{4+i}, A = 1 + eps. g in C(A) = (g0, g1): e_i -> g0 e_i + eps g1 e_i, e_{4+i} -> eps g0 e_i.
Checks:
  1. the 56 line spreads of PG(3,2) form one GL(4,2) orbit; |Stab(standard spread)| = 360;
     every 3-dim subspace of Alt(F_2^4) holds a nonzero degenerate form (no anisotropic 3-space).
  2. E (span of the spread's Klein points) has exactly 5 degenerate nonzero forms; inv2's printed E basis.
  3. F(L) decoded from inv2's L encoding matches the first 4 basis forms written in inv2-all.txt.
  4. b-part of C o g depends on g0 only; shears act on L by translations T (dim T); H = Stab x {g1 = 0}
     acts linearly and preserves T.
  5. the 18 reps' G-orbits, computed as unions of H-images mod T, are disjoint and cover all of Hom(E, Alt W),
     with the sizes inv2 reported.
  6. every line of inv2-all.txt: A-invariant, D(S8) = S8^A = F(L) as written, APN, and its rank distribution.
"""
import itertools, sys, random

ROOT = ""  # run from the folder that holds results/
PAIRS8 = [(i, j) for i in range(8) for j in range(i + 1, 8)]
PAIRS4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
par = lambda x: bin(x).count("1") & 1


# ---- 8-dim forms: list of 8 row masks ----
def coef_to_rows(c):
    r = [0] * 8
    for k, (i, j) in enumerate(PAIRS8):
        if c >> k & 1:
            r[i] ^= 1 << j; r[j] ^= 1 << i
    return r

def rows_to_coef(r):
    return sum(1 << k for k, (i, j) in enumerate(PAIRS8) if r[i] >> j & 1)

def B(r, x, y):
    s = 0
    for i in range(8):
        if x >> i & 1:
            s ^= par(r[i] & y)
    return s

def rvec(r, a):
    v = 0
    for i in range(8):
        if a >> i & 1:
            v ^= r[i]
    return v

def rank_vecs(vs):
    piv = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; break
    return len(piv)

def rank_form(r): return rank_vecs(r)

def compose(r, gimg):  # (B o g)(x, y) = B(g x, g y); gimg[i] = g e_i
    out = [0] * 8
    for i in range(8):
        for j in range(8):
            if i != j and B(r, gimg[i], gimg[j]): out[i] |= 1 << j
    return out

def apply_map(gimg, x):
    y = 0
    for i in range(8):
        if x >> i & 1: y ^= gimg[i]
    return y

A_IMG = [(1 << i) | (1 << (4 + i)) for i in range(4)] + [1 << (4 + i) for i in range(4)]

def cmap(g0, g1):  # g0, g1: 4 nibble columns
    return [g0[i] | (g1[i] << 4) for i in range(4)] + [g0[i] << 4 for i in range(4)]

# ---- 4-dim alternating forms (6 bits over PAIRS4) ----
def b4(b, x, y):
    s = 0
    for p, (i, j) in enumerate(PAIRS4):
        if b >> p & 1:
            s ^= ((x >> i & 1) & (y >> j & 1)) ^ ((x >> j & 1) & (y >> i & 1))
    return s

def rad4(b): return frozenset(x for x in range(1, 16) if all(b4(b, x, 1 << i) == 0 for i in range(4)))

def bpart(r): return sum(1 << p for p, (i, j) in enumerate(PAIRS4) if B(r, 1 << (4 + i), 1 << j))
def cpart(r): return sum(1 << p for p, (i, j) in enumerate(PAIRS4) if B(r, 1 << i, 1 << j))

def dm_form(b, c):  # the form in D(M) with b_C = b, c_C = c
    r = [0] * 8
    def set_(x, y):
        r[x] ^= 1 << y; r[y] ^= 1 << x
    for p, (i, j) in enumerate(PAIRS4):
        if b >> p & 1:
            set_(4 + i, j); set_(4 + j, i)
        if c >> p & 1:
            set_(i, j)
    return r

def gl4():
    for cols in itertools.product(range(1, 16), repeat=4):
        if rank_vecs(list(cols)) == 4: yield cols

def img4(g0, x):
    y = 0
    for i in range(4):
        if x >> i & 1: y ^= g0[i]
    return y


def main():
    random.seed(1)
    # ---- 1. spreads ----
    lines = set()
    for a in range(1, 16):
        for b in range(a + 1, 16):
            lines.add(frozenset({a, b, a ^ b}))
    lines = list(lines); assert len(lines) == 35
    spreads = [s for s in itertools.combinations(lines, 5) if len(frozenset().union(*s)) == 15]
    print("lines", len(lines), "spreads", len(spreads))
    std_masks = [0x000e, 0x1110, 0x8420, 0x4280, 0x2840]  # inv2's printed spread
    std = [frozenset(x for x in range(16) if m >> x & 1) for m in std_masks]
    assert frozenset().union(*std) == frozenset(range(1, 16)) and all(len(l) == 3 and l in lines for l in std)
    stdset = frozenset(std)
    G4 = list(gl4()); print("|GL(4,2)|", len(G4))
    orbit = set(); stab = []
    for g in G4:
        im = frozenset(frozenset(img4(g, x) for x in l) for l in std)
        orbit.add(im)
        if im == stdset: stab.append(g)
    print("spread orbit size", len(orbit), "(56 = one orbit)", "|Stab|", len(stab))
    assert len(orbit) == 56 == len(spreads) and len(stab) == 360

    # anisotropic 3-spaces of Alt(F_2^4)?
    deg = [b for b in range(64) if rad4(b)]  # b = 0 has full radical; rank 2 have radical size 3
    degset = set(b for b in range(1, 64) if len(rad4(b)) > 0)
    aniso = 0
    seen = set()
    for x, y, z in itertools.combinations(range(1, 64), 3):
        sp = frozenset({x, y, z, x ^ y, x ^ z, y ^ z, x ^ y ^ z})
        if 0 in sp or len(sp) != 7 or sp in seen: continue
        seen.add(sp)
        if not (sp & degset): aniso += 1
    print("3-dim subspaces of Alt(F_2^4):", len(seen), "anisotropic:", aniso)
    assert len(seen) == 1395 and aniso == 0

    # ---- 2. E ----
    bL = []
    for l in std:
        cands = [b for b in range(1, 64) if rad4(b) == l]
        assert len(cands) == 1; bL.append(cands[0])
    # inv2: first 4 independent in line order
    Eb = []
    for b in bL:
        if rank_vecs(Eb + [b]) > len(Eb): Eb.append(b)
        if len(Eb) == 4: break
    print("E basis %s (inv2 printed 20 01 2d 3b)" % " ".join("%02x" % b for b in Eb))
    assert Eb == [0x20, 0x01, 0x2d, 0x3b]
    Eel = {}
    for c in range(16):
        v = 0
        for i in range(4):
            if c >> i & 1: v ^= Eb[i]
        Eel[v] = c
    ndeg = sum(1 for v in Eel if v and rad4(v)); print("E: degenerate nonzero", ndeg)
    assert ndeg == 5

    def F_of(L):
        return [dm_form(Eb[k], (L >> (6 * k)) & 63) for k in range(4)]

    # D(M) check of dm_form: C = B o A + B for some B? (solve: D is linear; build image of the 28 basis forms)
    Dimg = []
    for k in range(28):
        r = coef_to_rows(1 << k); ra = compose(r, A_IMG)
        Dimg.append(rows_to_coef([ra[i] ^ r[i] for i in range(8)]))
    assert rank_vecs(Dimg) == 12
    for _ in range(50):
        c = rows_to_coef(dm_form(random.randrange(64), random.randrange(64)))
        assert rank_vecs(Dimg + [c]) == 12
    print("dm_form lands in D(M): ok; rank D = 12")

    # ---- 3. encoding vs inv2-all.txt ----
    allines = [l for l in open(ROOT + "results/inv2-all.txt") if "basis" in l]
    checkedL = set()
    for l in allines:
        t = l.split(); L = int(t[t.index("F") + 1], 16); bas = [int(h, 16) for h in t[t.index("basis") + 1:]]
        if L in checkedL: continue
        checkedL.add(L)
        mine = [rows_to_coef(r) for r in F_of(L)]
        print("rep %06x: F_of matches inv2-all's first 4 forms: %s" % (L, mine == bas[:4]))
        assert mine == bas[:4]

    # ---- 4. group action on L ----
    def act(gimg, L):
        X = [compose(r, gimg) for r in F_of(L)]
        coords = []
        for r in X:
            b = bpart(r); assert b in Eel, "g does not preserve E"
            # still in D(M) shape: K x K zero, K x W symmetric
            for i in range(4):
                for j in range(4):
                    assert B(r, 1 << (4 + i), 1 << (4 + j)) == 0
                    assert B(r, 1 << (4 + i), 1 << j) == B(r, 1 << (4 + j), 1 << i)
            coords.append(Eel[b])
        # combination of X with b-part Eb[m]
        Lp = 0
        for m in range(4):
            hit = None
            for comb in range(1, 16):
                v = 0
                for k in range(4):
                    if comb >> k & 1: v ^= Eb_vec(coords[k])
                if v == Eb[m]:
                    hit = comb; break
            Y = [0] * 8
            for k in range(4):
                if hit >> k & 1: Y = [Y[i] ^ X[k][i] for i in range(8)]
            assert bpart(Y) == Eb[m]
            Lp |= cpart(Y) << (6 * m)
        return Lp

    def Eb_vec(c):
        v = 0
        for i in range(4):
            if c >> i & 1: v ^= Eb[i]
        return v

    # every g in C(A) commutes with A
    for _ in range(20):
        g0 = random.choice(G4); g1 = [random.randrange(16) for _ in range(4)]
        gi = cmap(g0, g1)
        for x in range(256):
            assert apply_map(gi, apply_map(A_IMG, x)) == apply_map(A_IMG, apply_map(gi, x))
    # b-part depends on g0 only
    for _ in range(20):
        g0 = random.choice(G4); g1 = [random.randrange(16) for _ in range(4)]
        b, c = random.randrange(64), random.randrange(64)
        r1 = compose(dm_form(b, c), cmap(g0, g1)); r2 = compose(dm_form(b, c), cmap(g0, [0] * 4))
        assert bpart(r1) == bpart(r2)
    I4 = (1, 2, 4, 8)
    T = []
    for i in range(4):
        for j in range(4):
            g1 = [0] * 4; g1[i] = 1 << j
            gi = cmap(I4, g1)
            t0 = act(gi, 0)
            for _ in range(10):
                L = random.randrange(1 << 24)
                assert act(gi, L) == L ^ t0, "shear not a translation"
            T.append(t0)
    dimT = rank_vecs(T); print("shears act by translations; dim T =", dimT)
    piv = {}
    for v in T:
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; break
    def red(L):
        for h in sorted(piv, reverse=True):
            if L >> h & 1: L ^= piv[h]
        return L
    Hm = [cmap(g, [0] * 4) for g in stab]
    for gi in Hm[:60]:
        for t in T: assert red(act(gi, t)) == 0, "H doesn't preserve T"
        for _ in range(3):
            L1, L2 = random.randrange(1 << 24), random.randrange(1 << 24)
            assert act(gi, L1 ^ L2) == act(gi, L1) ^ act(gi, L2), "H not linear on L"
    print("H linear on L and preserves T (checked on 60 of 360)")

    # ---- 5. orbits of the reps ----
    reps = [(int(a, 16), int(b)) for a, b in (l.split() for l in open(ROOT + "results/inv2-reps.txt"))]
    union = set(); ok = True
    for L, sz in reps:
        orb = set(red(act(gi, L)) for gi in Hm)
        size = len(orb) << dimT
        dis = not (orb & union); union |= orb
        print("rep %06x: orbit %d (inv2 %d) disjoint %s" % (L, size, sz, dis))
        ok &= dis and size == sz
    print("cosets covered %d of %d; all reps consistent %s" % (len(union), 1 << (24 - dimT), ok))
    assert ok and len(union) == 1 << (24 - dimT)

    # ---- 6. outputs ----
    nbad = 0; dist = {}
    for l in allines:
        t = l.split(); L = int(t[t.index("F") + 1], 16); bas = [int(h, 16) for h in t[t.index("basis") + 1:]]
        R = [coef_to_rows(c) for c in bas]
        span = []  # all 255 nonzero forms as coef
        if rank_vecs(bas) != 8: nbad += 1; continue
        # A-invariant
        Aim = [rows_to_coef(compose(r, A_IMG)) for r in R]
        if rank_vecs(bas + Aim) != 8: nbad += 1; continue
        # D(S8) = span F(L) and S8^A = span F(L)
        Dv = [a ^ b for a, b in zip(Aim, bas)]
        Fc = bas[:4]
        if rank_vecs(Dv) != 4 or rank_vecs(Dv + Fc) != 4: nbad += 1; continue
        if any(rows_to_coef(compose(coef_to_rows(c), A_IMG)) != c for c in Fc): nbad += 1; continue
        # APN: for every a, the 8 vectors B_i(a, .) have rank 7
        if any(rank_vecs([rvec(r, a) for r in R]) != 7 for a in range(1, 256)): nbad += 1; continue
        cnt = {}
        for s in range(1, 256):
            r = [0] * 8
            for i in range(8):
                if s >> i & 1: r = [r[q] ^ R[i][q] for q in range(8)]
            k = rank_form(r); cnt[k] = cnt.get(k, 0) + 1
        claimed = {int(x.split(":")[0]): int(x.split(":")[1]) for x in t[1:5]}
        if any(cnt.get(k, 0) != claimed.get(k, 0) for k in (8, 6, 4, 2)): nbad += 1; continue
        key = tuple(cnt.get(k, 0) for k in (8, 6, 4, 2)); dist[key] = dist.get(key, 0) + 1
    print("inv2-all: %d lines, bad %d, distribution %s" % (len(allines), nbad, dist))


if __name__ == "__main__":
    main()
