# holes_tally.py <pspread3 holes out> <configs.txt> <Langevin-Hou fix-K.txt> [Kurz point sets] [Kurz maximal]
# Hole-set classes from `pspread3 holes`: configurations i, j share a hole set up to GL(8,2) iff some exact cover of
# C(P_i) lies on j. Checks every line's |Aut H| = count_j * |Aut P_j| agrees for all j it lists, that the classes are
# consistent (each member lists the same members, with the same |Aut H|), and compares the |Aut H| distribution and the
# number of point sets with Langevin-Hou class by class. Rowen (Claude Opus 5.5), 2026-10-30.
import sys, re, collections
from fractions import Fraction

out, confs, fix = sys.argv[1:4]
kurz = int(sys.argv[4]) if len(sys.argv) > 4 else None
kurzmax = int(sys.argv[5]) if len(sys.argv) > 5 else None

GL = 1
for i in range(8): GL *= 2**8 - 2**i

tags = [l.split()[0] for l in open(confs) if l.strip()]
w = [int(re.search(r'-w(\d+)$', t).group(1)) for t in tags]
n = len(tags)

rows = {}
bad = 0
for l in open(out):
    t = l.split()
    if not t: continue
    i, autp, nsol, ncov = int(t[1]), int(t[2]), int(t[3]), int(t[4])
    if 'MISSING' in t: print('MISSING', l.strip()); bad += 1; continue
    hits = {int(a): int(b) for a, b in (x.split(':') for x in t[5:])}
    rows[i] = (autp, ncov, hits)
assert len(rows) == n - bad, (len(rows), n)
autp = {i: r[0] for i, r in rows.items()}

autH = {}
for i, (ap, ncov, hits) in rows.items():
    if i not in hits: print('SELF NOT HIT', tags[i]); bad += 1
    if sum(hits.values()) != ncov: print('COUNT', tags[i]); bad += 1
    vals = {c * autp[j] for j, c in hits.items()}
    if len(vals) != 1: print('AUT H DISAGREES', tags[i], vals); bad += 1
    autH[i] = vals.pop()

classes = {}
for i, (ap, ncov, hits) in rows.items():
    key = frozenset(hits)
    for j in hits:
        if frozenset(rows[j][2]) != key or autH[j] != autH[i]: print('INCONSISTENT', tags[i], tags[j]); bad += 1
    classes[key] = autH[i]
print(f'configurations {n}, hole-set classes {len(classes)}' + (f' (Kurz {kurz})' if kurz else ''))

nmax = 0
for key in classes:
    ws = {w[j] for j in key}
    if len(ws & {0}) and len(ws) > 1: print('MIXED MAXIMAL', [tags[j] for j in key]); bad += 1
    nmax += ws == {0}
print(f'maximal (hole set contains no solid, w = 0): {nmax}' + (f' (Kurz {kurzmax})' if kurzmax else ''))

sizes = collections.Counter(len(k) for k in classes)
print('configurations per hole set:', dict(sorted(sizes.items())))

ours = collections.Counter(classes.values())
lh = {}
txt = open(fix).read()
for m in re.finditer(r'fix=(\d+)\s+number of class (\d+)\s+qte=(\d+)', txt):
    f, c, q = map(int, m.groups())
    lh[f] = c
    if q != c * GL // f or (c * GL) % f: print('LH qte inconsistent', f); bad += 1
lhtotal = int(re.search(r'Total=(\d+)', txt).group(1))
diff = [(f, ours.get(f, 0), lh.get(f, 0)) for f in sorted(set(ours) | set(lh)) if ours.get(f, 0) != lh.get(f, 0)]
print('|Aut H| distribution vs Langevin-Hou:', 'identical' if not diff else diff, f'({len(lh)} orders)')
total = sum(Fraction(GL, a) for a in classes.values())
print(f'point sets: {total} vs Langevin-Hou {lhtotal}:', 'MATCH' if total == lhtotal else 'DIFFER')
# labelled covers: each class contributes |GL| / |Aut P_j| partial spreads for each j on it
spreads = sum(Fraction(GL, autp[j]) for j in autp)
print(f'labelled partial spreads (sum |GL| / |Aut P|): {spreads}')
print('problems:', bad)
sys.exit(1 if bad or diff or total != lhtotal else 0)
