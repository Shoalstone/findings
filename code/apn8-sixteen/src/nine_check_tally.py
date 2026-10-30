"""Tallies of nine_check search runs (search B for i = 8), from the .prog files in runs/check9.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-28.
  python apn/nine_check_tally.py prog <configs.txt or -> <expected count> <prog>...
      sums the R lines (each configuration once; duplicates flagged), by type when a configs.txt is given; SOLUTION lines;
      distinct P1 spans over all the files.
R line: id, nodes[0..8], tries[0..8], p1[0..8], bentext, full, triesD[0..8].
"""
import sys


def conf_types(path):
    t = {}
    for line in open(path):
        if line.startswith("CONF"):
            p = line.split()
            t[int(p[1])] = p[2]
    return t


def prog(cpath, expect, files):
    types = conf_types(cpath) if cpath != "-" else {}
    seen, dup, sols, p1 = {}, 0, [], set()
    tot = {}
    for f in files:
        for line in open(f):
            if line.startswith("R "):
                v = list(map(int, line.split()[1:]))
                assert len(v) == 39, line
                i = v[0]
                if i in seen:
                    dup += 1
                    if seen[i] != v:
                        print(f"config {i}: two different R lines")
                    continue
                seen[i] = v
            elif line.startswith("SOLUTION"):
                sols.append(line.strip())
            elif line.startswith("P1 "):
                p1.add(line.split("basis", 1)[1].strip())
    for i, v in seen.items():
        for key in ("all", types.get(i, "all")):
            a = tot.setdefault(key, [0] * 38)
            for k in range(38):
                a[k] += v[1 + k]
            if key == "all" and types.get(i, "all") == "all":
                break
    n = len(seen)
    print(f"configurations done {n} (expect {expect}){'' if n == expect else '  INCOMPLETE'}; duplicate R lines {dup}; "
          f"ids outside 0..{expect - 1}: {sum(1 for i in seen if not 0 <= i < expect)}")
    count = {}
    for i in seen:
        count[types.get(i, "all")] = count.get(types.get(i, "all"), 0) + 1
    for key, a in tot.items():
        nodes, tries, p1n, bent, full, triesD = a[0:9], a[9:18], a[18:27], a[27], a[28], a[29:38]
        c = n if key == "all" else count[key]
        print(f"{key} ({c}): nodes {' / '.join(f'{x:,}' for x in nodes)}")
        print(f"   near-misses (p1) by dim {' / '.join(f'{x:,}' for x in p1n)}; tries in D {sum(triesD):,} "
              f"(dim 6 {triesD[6]:,}, dim 7 {triesD[7]:,}); bent extension {bent}; full spans {full}")
    print(f"distinct P1 spans {len(p1)}; SOLUTION lines {len(sols)}")
    for s in sols:
        print("  " + s)


if __name__ == "__main__":
    if sys.argv[1] == "prog":
        prog(sys.argv[2], int(sys.argv[3]), sys.argv[4:])
