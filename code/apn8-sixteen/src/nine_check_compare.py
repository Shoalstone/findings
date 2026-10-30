"""Section 7 of nine-check.md: search B's near-miss records against search A's (comparison only, after the verdict).

  python src/nine_check_compare.py [logs dir]      (default: logs/ beside src/)
Reads logs/check9-mapa-output.txt (A's config -> ours), logs/nine-output.txt (A's NEAR lines), and
logs/check9-{asc,desc,rand7}-prog.txt (our R lines: p1[d] = spans of dimension d holding all 9 rank-4 forms; P1 lines:
those spans' bases, without config). A logs NEAR only at the step where the last rank-4 form joins; ours logs every
distinct span holding all 9, so ours also counts a lower near-miss extended by a rank-6 form in D. Those are separated
by containment.
Rowen (Claude Opus 5.5), 2026-10-29. This folder's copy reads the files here (the original read Rowen's working tree:
runs/check9/mapa.out, runs/i8/cover9.out, runs/check9/<order>.prog); otherwise unchanged.
"""
import collections, os, re, sys

LOGS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'logs')

A2B = {}
for l in open(os.path.join(LOGS, "check9-mapa-output.txt")):
    m = re.match(r"A (\S+) -> (\d+) ", l)
    if m:
        A2B[m.group(1)] = int(m.group(2))
anear = collections.Counter()
for l in open(os.path.join(LOGS, "nine-output.txt")):
    m = re.match(r"NEAR (\S+) dim (\d+)", l)
    if m:
        anear[(A2B[m.group(1)], int(m.group(2)))] += 1


def rank(vs):
    rows, r = list(vs), 0
    for bit in range(28):
        piv = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


for order in ("asc", "desc", "rand7"):
    bnear = collections.Counter(); spans = collections.defaultdict(list)
    for l in open(os.path.join(LOGS, f"check9-{order}-prog.txt")):
        if l.startswith("R "):
            v = list(map(int, l.split()[1:]))
            for d in range(9):
                if v[19 + d]:
                    bnear[(v[0], d)] = v[19 + d]
        elif l.startswith("P1 "):
            p = l.split()
            spans[int(p[2])].append([int(x, 16) for x in p[4:]])
    print(f"== {order}")
    for d in (5, 6, 7):
        ca = {c for (c, dd) in anear if dd == d}; cb = {c for (c, dd) in bnear if dd == d}
        same = sum(anear[(c, d)] == bnear[(c, d)] for c in ca | cb)
        print(f"dim {d}: A {sum(anear[(c, d)] for c in ca)} on {len(ca)} configs, B {sum(bnear[(c, d)] for c in cb)} on "
              f"{len(cb)}; configs with equal counts {same} of {len(ca | cb)}")
    low5 = {c for (c, dd) in bnear if dd == 5}
    d6 = [c for (c, dd) in anear if dd == 6]
    off = [c for c in set(d6) | {c for (c, dd) in bnear if dd == 6} if c not in low5 and anear[(c, 6)] != bnear[(c, 6)]]
    print(f"dim 6 off the dim-5 configs: unequal on {len(off)}")
    if sum(len(x) for x in spans.values()) == sum(bnear.values()):
        for d in (6, 7):
            lower = [s for dd in range(5, d) for s in spans[dd]]
            fresh = 0
            for s in spans[d]:
                r = rank(s)
                if not any(rank(s + t) == r for t in lower):
                    fresh += 1
            print(f"dim {d} spans not containing a lower near-miss: {fresh} of {len(spans[d])}")
    else:
        print(f"P1 spans incomplete ({sum(len(x) for x in spans.values())} of {sum(bnear.values())}): no containment count")
