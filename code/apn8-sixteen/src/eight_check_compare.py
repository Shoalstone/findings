"""Cross-checks of eight_check's search outputs (search B, 8 rank-4 components), from runs/check8.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-30.
  python src/eight_check_compare.py orders [logs dir]   per configuration, near-misses across asc / desc / rand 7:
                                             dim-5 counts equal everywhere; dim-6 counts equal on configurations without
                                             a dim-5 one
  python src/eight_check_compare.py numpy [logs dir]    eight_check_numpy.py's results against asc's R lines: nodes,
                                             near-misses and tries in D by dimension, configuration by configuration
This folder's copy (made for it by Claude Opus 5.5) reads the files in logs/ (default: beside src/):
check8-{asc,desc,rand7}-prog.txt (each order's 8 ranges, concatenated) and check8-numpy-output.txt (the 32 numpy-<i>.out
concatenated). The original read runs/check8/<order>-r*.prog and runs/check8/numpy-<i>.out; otherwise unchanged.
"""
import glob
import os
import re
import sys

D = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs")


def rlines(order):
    out = {}
    for f in [os.path.join(D, f"check8-{order}-prog.txt")]:
        for line in open(f):
            if line.startswith("R "):
                v = list(map(int, line.split()[1:]))
                out[v[0]] = {"nodes": v[1:10], "p1": v[19:28], "triesD": v[30:39]}
    return out


def orders():
    o = {k: rlines(k) for k in ("asc", "desc", "rand7")}
    n5 = sum(1 for i in o["asc"] if len({o[k][i]["p1"][5] for k in o}) != 1)
    with5 = [i for i in o["asc"] if o["asc"][i]["p1"][5]]
    n6 = sum(1 for i in o["asc"] if i not in with5 and len({o[k][i]["p1"][6] for k in o}) != 1)
    s6 = {k: sum(o[k][i]["p1"][6] for i in o[k] if i not in with5) for k in o}
    conf6 = {k: {i for i in o[k] if o[k][i]["p1"][6]} for k in o}
    print(f"configurations {len(o['asc'])}; dim-5 near-misses differ across orders on {n5}; {len(with5)} configurations "
          f"with a dim-5 near-miss ({sum(o['asc'][i]['p1'][5] for i in with5)} in all); dim-6 counts differ on {n6} of the "
          f"others; dim-6 near-misses on configurations without a dim-5 one: {s6}; configurations with a dim-6 "
          f"near-miss: {', '.join(f'{k} {len(v)}' for k, v in conf6.items())}, same set: {len({frozenset(v) for v in conf6.values()}) == 1}")


def numpy():
    asc = rlines("asc")
    agree = bad = 0
    for m in re.finditer(r"^config (\d+): nodes \[(.*?)\] p1 \[(.*?)\] triesD \[(.*?)\] full (\d+) allcovered-below-8 (\d+)",
                         open(os.path.join(D, "check8-numpy-output.txt")).read(), re.M):
        i = int(m[1])
        got = {"nodes": list(map(int, m[2].split(","))), "p1": list(map(int, m[3].split(","))), "triesD": list(map(int, m[4].split(",")))}
        ok = got == asc[i] and m[5] == "0" and m[6] == "0"
        agree += ok
        bad += not ok
        print(f"config {i}: {'agrees' if ok else 'DIFFERS'}  nodes {got['nodes']} p1 {got['p1']} triesD {sum(got['triesD']):,}"
              + ("" if ok else f"  C++ {asc[i]}"))
    print(f"numpy re-search: {agree} configurations agree with asc in every count, {bad} differ")


{"orders": orders, "numpy": numpy}[sys.argv[1]]()
