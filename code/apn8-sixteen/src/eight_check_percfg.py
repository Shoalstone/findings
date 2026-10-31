"""Search A's and search B's near misses for i = 9 (eight rank-4 components) compared configuration by configuration.

nine_check_percfg.py (written by the reviewing subagent for nine) carried over to eight: both searches cover every
point off the holes before any hole, so a near miss in search A's sense (logged when the last rank-4 form joins) is
exactly the span of the eight rank-4 forms. Search B logs every distinct span holding all eight (P1 lines, without the
configuration), including such a span plus rank-6 forms in the holes. Here each P1 span is assigned to its
configuration by the radicals of its rank-4 elements (matched against data/configs8.txt), and the spans containing no
smaller one of the same configuration are compared with search A's NEAR lines (logs/eight-output-c*.txt, mapped through
logs/check8-mapa-output.txt), per (configuration, dimension). Expected: 0 mismatches in all three orders (search B's P1
record is complete in each, unlike desc's at nine).

  python src/eight_check_percfg.py [BUNDLE_DIR] [orders...]      (default: asc desc rand7; a few minutes each)
Written for this folder by Claude Opus 5.5, 2026-10-30, from nine_check_percfg.py (file names changed, the radical of a
form cached, A's NEAR lines read from the eight chunk outputs).
"""
import collections, functools, glob, os, re, sys
args = sys.argv[1:]
B = os.path.join(args.pop(0), "") if args and os.path.isdir(args[0]) else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "")
orders = args or ["asc", "desc", "rand7"]
PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]

def rows(c):
    r = [0] * 8
    for p, (i, j) in enumerate(PAIRS):
        if c >> p & 1:
            r[i] |= 1 << j; r[j] |= 1 << i
    return r

@functools.lru_cache(maxsize=None)
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
for l in open(B + "data/configs8.txt"):
    if not l.startswith("CONF"): continue
    idx = int(l.split()[1]); Ns = [int(t, 16) for t in l.split(":")[1].split()]
    mem = frozenset([X, Y, graph(0x8421)] + [graph(N) for N in Ns])
    cfg[mem] = idx

A2B = {}
for l in open(B + "logs/check8-mapa-output.txt"):
    m = re.match(r"A (\S+) -> (\d+) ", l)
    if m: A2B[m.group(1)] = int(m.group(2))
anear = collections.Counter()
for f in sorted(glob.glob(B + "logs/eight-output-c*.txt")):
    for l in open(f):
        m = re.match(r"NEAR (\S+) dim (\d+)", l)
        if m: anear[(A2B[m.group(1)], int(m.group(2)))] += 1

def span(basis):
    el = [0]
    for b in basis: el += [e ^ b for e in el]
    return el

for order in orders:
    per = collections.defaultdict(list)  # config -> list of (dim, set of elements)
    for l in open(B + f"logs/check8-{order}-prog.txt"):
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
          "A by dim", {d: sum(v for (c, dd), v in anear.items() if dd == d) for d in (5, 6, 7)}, "per-(config,dim) mismatches", len(bad), bad[:10], flush=True)
