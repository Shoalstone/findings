# cover_tally.py <cover3 .out>...: one summary over several cover3 outputs (chunks of one input).
# Sums the done: lines and per-dimension node counts; near misses by dimension, by w (from the tag), and the busiest
# configurations. Rowen (Claude Opus 5.5), 2026-10-30.
import sys, re, collections

tot = collections.Counter(); dims = None; secs = 0; ndone = 0
near = []
for fn in sys.argv[1:]:
    done = False
    for l in open(fn):
        if l.startswith('done: covers'):
            done = True; ndone += 1
            for part in l[6:].strip().split(', ')[:-1]:  # last part is the time
                k, v = part.rsplit(' ', 1); tot[k] += int(v)
            secs += int(re.search(r'(\d+) s$', l.strip()).group(1))
        elif l.startswith('nodes by span dimension'):
            v = list(map(int, l.split(':')[1].split()))
            dims = v if dims is None else [a + b for a, b in zip(dims, v)]
        elif l.startswith('NEAR '):
            t = l.split(); near.append((t[1], int(t[3])))
        elif l.startswith('SOLUTION') or l.startswith('FULL'):
            print('!!', fn, l.strip())
    if not done: print('NOT DONE', fn)
print(f'{ndone}/{len(sys.argv) - 1} outputs done, {secs} s')
for k, v in tot.items(): print(f'  {k}: {v:,}')
print('nodes by span dimension 1..8:', ' / '.join(f'{x:,}' for x in dims or []))
print('near misses', len(near), 'by dim', dict(sorted(collections.Counter(d for _, d in near).items())))
cfg = collections.Counter(t for t, _ in near)
print('configurations with near misses', len(cfg), 'by w',
      dict(sorted(collections.Counter(re.search(r'-w(\d+)$', t).group(1) for t in cfg).items())))
print('near misses by w', dict(sorted(collections.Counter(re.search(r'-w(\d+)$', t).group(1) for t, _ in near).items())))
first = min((d for _, d in near), default=None)
print('lowest near-miss dim', first, [t for t, d in near if d == first][:10])
print('most:', cfg.most_common(5))
