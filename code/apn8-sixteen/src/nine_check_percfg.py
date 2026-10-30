"""Search A's and search B's near misses for i = 8 compared configuration by configuration (nine-check.md section 7).

Both searches cover every point off the holes before any hole, so a near miss in search A's sense (logged when the
last rank-4 form joins) is exactly the span of the nine rank-4 forms. Search B logs every distinct span holding all
nine (P1 lines, without the configuration), including such a span plus rank-6 forms in the holes. Here each P1 span
is assigned to its configuration by the radicals of its rank-4 elements (matched against data/configs9.txt), and the
spans containing no smaller one of the same configuration are compared with search A's NEAR lines, per
(configuration, dimension). Expected: 0 mismatches for asc and rand7 (desc's P1 record is incomplete).

  python src/nine_check_percfg.py [BUNDLE_DIR] asc rand7      (about a minute)
Written by the reviewing subagent (Claude Opus 5.5) for Rowen, 2026-10-30; bundle-relative paths added by Rowen.
"""
import collections, os, re, sys
args = sys.argv[1:]
B = os.path.join(args.pop(0), "") if args and os.path.isdir(args[0]) else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "")
PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]

def rows(c):
    r = [0] * 8
    for p, (i, j) in enumerate(PAIRS):
        if c >> p & 1:
            r[i] |= 1 << j; r[j] |= 1 << i
    return r

def radical(c):
    r = rows(c); out = []
    for v in range(1, 256):
        w = 0
        for i in range(8):
            if v >> i & 1: w ^= r[i]
        if not w: out.append(v)
    return frozenset(out)

def par(x): return bin(x).count("1") & 1
def mv4(N, x): return sum(par(((N >> (4 * i)) & 15) & x) << i for i in range(4))
def graph(N): return frozenset(x | (mv4(N, x) << 4) for x in range(1, 16))
X = frozenset(range(1, 16)); Y = frozenset(v << 4 for v in range(1, 16))

cfg = {}
for l in open(B + "data/configs9.txt"):
    if not l.startswith("CONF"): continue
    idx = int(l.split()[1]); Ns = [int(t, 16) for t in l.split(":")[1].split()]
    mem = frozenset([X, Y, graph(0x8421)] + [graph(N) for N in Ns])
    cfg[mem] = idx

A2B = {}
for l in open(B + "logs/check9-mapa-output.txt"):
    m = re.match(r"A (\S+) -> (\d+) ", l)
    if m: A2B[m.group(1)] = int(m.group(2))
anear = collections.Counter()
for l in open(B + "logs/nine-output.txt"):
    m = re.match(r"NEAR (\S+) dim (\d+)", l)
    if m: anear[(A2B[m.group(1)], int(m.group(2)))] += 1

def span(basis):
    el = [0]
    for b in basis: el += [e ^ b for e in el]
    return el

for order in args:
    per = collections.defaultdict(list)  # config -> list of (dim, set of elements)
    for l in open(B + f"logs/check9-{order}-prog.txt"):
        if not l.startswith("P1 "): continue
        p = l.split(); basis = [int(x, 16) for x in p[4:]]
        el = span(basis)
        rads = set()
        for e in el[1:]:
            rd = radical(e)
            if len(rd) == 15: rads.add(rd)
        c = cfg[frozenset(rads)]
        per[c].append((len(basis), frozenset(el)))
    fresh = collections.Counter()
    for c, L in per.items():
        for d, s in L:
            if not any(dd < d and t <= s for dd, t in L): fresh[(c, d)] += 1
    keys = set(fresh) | set(anear)
    bad = [k for k in keys if fresh[k] != anear[k]]
    print(order, "configs", len(per), "fresh by dim", {d: sum(v for (c, dd), v in fresh.items() if dd == d) for d in (5, 6, 7)},
          "A by dim", {d: sum(v for (c, dd), v in anear.items() if dd == d) for d in (5, 6, 7)}, "per-(config,dim) mismatches", len(bad), bad[:10])
