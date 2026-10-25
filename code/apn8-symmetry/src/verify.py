"""Independent check of apn8 output lines: rebuilds each function's look-up table from the 8 coefficient vectors
(28 bits each, bit k = pair (i,j), i<j in lexicographic order) and computes, by brute force and the Walsh transform
rather than the engine's linear algebra, the differential uniformity and every component's rank.
Usage: verify.py <file> [--limit N]   (lines: 'ranks ... basis h1 .. h8')."""
import sys
import numpy as np

PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]
X = np.arange(256)
BITS = (X[:, None] >> np.arange(8)) & 1  # 256 x 8


def lut(basis):
    F = np.zeros(256, dtype=np.int64)
    for s, c in enumerate(basis):
        f = np.zeros(256, dtype=np.int64)
        for k, (i, j) in enumerate(PAIRS):
            if c >> k & 1:
                f ^= BITS[:, i] & BITS[:, j]
        F |= f << s
    return F


def delta(F):
    best = 0
    for a in range(1, 256):
        best = max(best, np.bincount(F ^ F[X ^ a], minlength=256).max())
    return int(best)


H = np.array([[1]])
for _ in range(8):
    H = np.block([[H, H], [H, -H]])  # Sylvester Hadamard, H[u,x] = (-1)^{u.x}


def rank_dist(F):
    cnt = {}
    for b in range(1, 256):
        comp = np.bitwise_count(F & b) & 1 if hasattr(np, "bitwise_count") else np.array([bin(v).count("1") & 1 for v in (F & b)])
        W = H @ (1 - 2 * comp.astype(np.int64))
        m = int(np.abs(W).max())          # quadratic of rank r: |W| in {0, 2^(8 - r/2)}
        r = 2 * (8 - (m.bit_length() - 1))
        cnt[r] = cnt.get(r, 0) + 1
    return cnt


def main():
    path = sys.argv[1]
    limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else 10**9
    n = bad = 0
    for line in open(path):
        if "basis" not in line:
            continue
        head, tail = line.split("basis")
        basis = [int(h, 16) for h in tail.split()]
        claimed = dict((int(t.split(":")[0]), int(t.split(":")[1])) for t in head.split()[1:] if ":" in t)
        F = lut(basis)
        d = delta(F)
        rd = rank_dist(F)
        ok = d == 2 and all(rd.get(r, 0) == claimed.get(r, 0) for r in (8, 6, 4, 2))
        n += 1
        bad += not ok
        if not ok or n <= 3:
            print(("OK  " if ok else "BAD ") + f"delta={d} ranks={dict(sorted(rd.items(), reverse=True))} claimed={claimed}")
        if n >= limit:
            break
    print(f"checked {n}, bad {bad}")


if __name__ == "__main__":
    main()
