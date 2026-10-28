"""The near misses for i = 7 (ten), search A against search B, from the logs in this folder.

usage: python src/near_misses10.py [logs dir]      (default: logs/ beside src/)

Reads search A's NEAR lines (cover.cpp on data/p10.txt, logs/ten-output.txt), the map of A's configurations onto B's
(logs/check10-mapa-output.txt), and for each of B's three point orders its per-configuration counts (the p1 column of
logs/check10-<order>-output.txt) and its near-miss spans (logs/check10-<order>-p1.txt). Prints, per dimension, A's and
B's near misses configuration by configuration, and how many 7-dimensional near misses on each side contain a
6-dimensional one of the same side (A's cover.cpp logs NEAR only at the step where the last rank-4 form joins; B logs
every span holding all ten). Forms are 28-bit masks in hex, as both programs print them.
Written for this folder (Claude Opus 5.5, 2026-10-28).
"""
import os
import re
import sys
from collections import Counter

LOGS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'logs')


def rank(vectors):
    basis = []  # reduced rows, distinct leading bits
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def contains(big, small):
    return rank(big + small) == rank(big)


def spans(lines, prefix):
    out = []
    for l in lines:
        m = re.match(prefix + r'.*?dim (\d+) .*basis ((?:[0-9a-f]+ ?)+)$', l.strip())
        if m:
            out.append((int(m.group(1)), [int(x, 16) for x in m.group(2).split()]))
    return out


def nested(sp):
    six = [b for d, b in sp if d == 6]
    sev = [b for d, b in sp if d == 7]
    return len(sev), sum(any(contains(b7, b6) for b6 in six) for b7 in sev)


tag_to_b = {}
for l in open(os.path.join(LOGS, 'check10-mapa-output.txt')):
    m = re.match(r'A (\S+) -> (\d+) ', l)
    if m:
        tag_to_b[m.group(1)] = int(m.group(2))

a_lines = [l for l in open(os.path.join(LOGS, 'ten-output.txt')) if l.startswith('NEAR ')]
a = {6: Counter(), 7: Counter()}
for l in a_lines:
    tag, dim = l.split()[1], int(l.split()[3])
    a[dim][tag_to_b[tag]] += 1
a_sp = spans(a_lines, 'NEAR')
n7, n7in = nested(a_sp)
print(f'A: {sum(a[6].values())} near misses at dim 6 in {len(a[6])} configurations, {sum(a[7].values())} at dim 7 in '
      f'{len(a[7])}; configurations (B numbering) {len(set(a[6]) | set(a[7]))}; dim-7 ones containing a dim-6 one: {n7in}')

for order in ('asc', 'desc', 'rand7'):
    b = {6: Counter(), 7: Counter()}
    for l in open(os.path.join(LOGS, f'check10-{order}-output.txt')):
        m = re.match(r'config (\d+) .* p1((?: \d+){9}) ', l)
        if m:
            p1 = [int(x) for x in m.group(2).split()]
            for d in (6, 7):
                if p1[d]:
                    b[d][int(m.group(1))] = p1[d]
    sp = spans(open(os.path.join(LOGS, f'check10-{order}-p1.txt')), 'P1')
    n7, n7in = nested(sp)
    same6 = b[6] == a[6]
    agree7 = sum(b[7].get(c, 0) == a[7].get(c, 0) for c in set(a[6]) | set(a[7]))
    extra7 = sorted(c for c in b[7] if c not in a[7])
    print(f'B {order}: dim 6 {sum(b[6].values())} in {len(b[6])} configurations, the same counts as A on every '
          f'configuration: {"yes" if same6 else "NO"}; dim 7 {sum(b[7].values())} ({n7} spans in the p1 file), '
          f'{n7in} of them containing a dim-6 near miss, leaving {n7 - n7in}; dim-7 counts equal to A\'s on {agree7} of '
          f'A\'s {len(set(a[6]) | set(a[7]))} configurations; dim-7 near misses where A has none: {extra7}')
