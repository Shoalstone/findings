"""Tallies for ten-check (i = 7), from outputs in runs/check10.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-27.
  python apn/ten_check_tally.py reps <dir>    configurations no control reached (REPS lines of every ctl-*.txt there)
  python apn/ten_check_tally.py maxall <dir>  per max orbit, orbitcov tallies summed over maxall-orbits.out against n in orbits.txt
"""
import glob
import os
import re
import sys


def conf_types(d):
    t = {}
    for line in open(os.path.join(d, "configs.txt")):
        if line.startswith("CONF"):
            p = line.split()
            t[int(p[1])] = p[2]
    return t


def reps(d):
    types = conf_types(d)
    hit = set()
    for f in glob.glob(os.path.join(d, "ctl-*.txt")):
        for line in open(f):
            if line.startswith("REPS"):
                hit |= {int(t.split(":")[0]) for t in line.split()[2:]}
    miss = sorted(set(types) - hit)
    by = {}
    for i in miss:
        by[types[i]] = by.get(types[i], 0) + 1
    print(f"reps hit {len(hit & set(types))} of {len(types)}; not reached {len(miss)} {by}")
    print(" ".join(map(str, miss)))


def maxall(d):
    n = {}
    for line in open(os.path.join(d, "orbits.txt")):
        m = re.match(r"ORB (\d+) max .* n (\d+) ", line)
        if m:
            n[int(m.group(1))] = int(m.group(2))
    got, runs, bad = {}, 0, 0
    for line in open(os.path.join(d, "maxall-orbits.out")):
        m = re.match(r"orbitcov: (\d+) cliques, (\d+) on no listed orbit, (\d+) on a non-max orbit", line)
        if m:
            runs += 1
            bad += int(m.group(2)) + int(m.group(3))
        m = re.match(r"\s+orbit (\d+) \(config -?\d+\): (\d+) sampled", line)
        if m:
            got[int(m.group(1))] = got.get(int(m.group(1)), 0) + int(m.group(2))
    diff = [o for o in n if got.get(o, 0) != n[o]]
    print(f"orbitcov runs {runs} (expect 17); cliques off the max orbits {bad}; max orbits {len(n)}, "
          f"cliques {sum(got.values())} (Part A {sum(n.values())}); orbits whose tally differs from Part A's n: {len(diff)}")
    if diff:
        print(" ".join(f"{o}:{got.get(o, 0)}/{n[o]}" for o in diff))


{"reps": reps, "maxall": maxall}[sys.argv[1]](sys.argv[2])
