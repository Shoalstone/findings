"""Compare the bundle re-runs (this folder) with the originals. Ran from Rowen's working tree: B = the bundle,
C8 = search B's run folder (its outputs are in logs/ here under other names), runs/i9 = search A's. Output: compare.log.
Rowen (Claude Opus 5.5), 2026-10-30."""
import glob, os, re

RR = os.path.dirname(os.path.abspath(__file__)) + '/'  # the re-run outputs (this folder)
B = 'findings/code/apn8-sixteen/'
C8 = 'runs/check8/'


def lines(p):
    return open(p, errors='replace').read().splitlines()


def same_file(a, b, label, sort=False):
    la, lb = lines(a), lines(b)
    if sort:
        la, lb = sorted(la), sorted(lb)
    print(f"{label}: {'IDENTICAL' if la == lb else 'DIFFER'} ({len(la)} vs {len(lb)} lines)")


# Part A slices and their maxc
for c, s in ((4, 0), (4, 1), (3, 3)):
    same_file(RR + f'parta-c{c}-s{s}.txt', B + f'data/parta8-c{c}-s{s}.txt', f'parta c{c} s{s}')
    same_file(RR + f'maxc-c{c}-s{s}.txt', C8 + f'maxc-c{c}-s{s}.txt', f'maxc c{c} s{s}')
same_file(RR + 'orbits.txt', B + 'data/orbits8.txt', 'merge orbits')
same_file(RR + 'configs8.txt', B + 'data/configs8.txt', 'configs list')


def strip_time(ls):
    return [re.sub(r'\[?\d+(\.\d+)? ?s\]?|\d+ ms', 'T', l) for l in ls if not l.startswith('#')]


for name, orig in (('e8-merge.log', C8 + 'merge.out'), ('e8-configs.log', C8 + 'configs.out'), ('e8-mapa.log', C8 + 'mapa.out')):
    a, b = strip_time(lines(RR + name)), strip_time(lines(orig))
    diff = [x for x in a if x not in set(b)]
    print(f"{name} vs {orig}: {len(a)} vs {len(b)} lines, lines of the re-run not in the original: {len(diff)}", diff[:3])

# cover chunks: node line and NEAR set
for k in (0, 7):
    a, b = lines(RR + f'cover-p8-c{k}.log'), lines(f'runs/i9/cover/c{k}.out')
    na = [l for l in a if l.startswith('nodes by span')]; nb = [l for l in b if l.startswith('nodes by span')]
    ea = sorted(l for l in a if l.startswith('NEAR')); eb = sorted(l for l in b if l.startswith('NEAR'))
    da = [l for l in a if l.startswith('done:')]; db = [l for l in b if l.startswith('done:')]
    print(f"cover c{k}: nodes line equal {na == nb}; NEAR {len(ea)} vs {len(eb)} equal as sets {ea == eb}; done line (sans time) equal "
          f"{[re.sub(r', \d+ s$', '', x) for x in da] == [re.sub(r', \d+ s$', '', x) for x in db]}")


def rl(files):
    out = {}
    for f in files:
        for l in lines(f):
            if l.startswith('R '):
                v = l.split()
                out[int(v[1])] = v[2:]
    return out


for name, pat in (('asc', 'asc-r*.prog'), ('desc', 'desc-r*.prog'), ('rand7', 'rand7-r*.prog'), ('maxs', 'maxs-r*.prog'), ('s9', 's9-r*.prog')):
    a = rl([RR + name + '.prog']); b = rl(glob.glob(C8 + pat))
    common = [i for i in a if i in b]
    bad = [i for i in common if a[i] != b[i]]
    print(f"{name}: re-run R lines {len(a)}, in original {len(common)}, unequal {len(bad)}", bad[:5])

# ctlrot reps in original
orig = [l for f in glob.glob(C8 + 'ctlrot-*.out') for l in lines(f) if re.match(r'rep (91|92|9364|9365) ', l)]
new = [l for l in lines(RR + 'e8-ctlrot.log') if l.startswith('rep ')]
print('ctlrot reps:', 'IDENTICAL' if sorted(orig) == sorted(new) else 'DIFFER', orig, new)
# numpy 9191
o = [l for f in glob.glob(C8 + 'numpy-*.out') for l in lines(f) if l.startswith('config 9191:')]
n = [l for l in lines(RR + 'e8-numpy-9191.log') if l.startswith('config 9191:')]
strip = lambda x: re.sub(r' \[\d+ s\]', '', x)
print('numpy 9191:', [strip(x) for x in o] == [strip(x) for x in n], o[:1])
o = [strip(l) for l in lines(C8 + 'numpy-ctl-bl22.out') if l.startswith('control line 1') and l.split()[2] in ('15:', '16:')]
n = [strip(l) for l in lines(RR + 'e8-numpy-ctl.log') if l.startswith('control line')]
print('numpy ctl lines 15,16:', o == n, len(o), len(n))
# nond control summary
a = [l for l in lines(RR + 'e8-ctlnond.log') if not l.startswith('#')][-6:]
b = lines(C8 + 'ctl-ind-bl22.out')[-6:]
print('ctl nond tail equal:', strip_time(a) == strip_time(b))
same_file(RR + 'ctl-ind-bl22.txt', C8 + 'ctl-ind-bl22.txt', 'ctl nond txt')
