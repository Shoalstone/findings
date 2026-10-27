"""Full-depth recount of twelve_check's search tree (asc order) for all 376 configurations, with
twelve_check_levels.py (which shares no code with twelve_check). twelve_check never branches at points of D for
twelve, so the recount to dimension 8 covers the whole tree. Resumable: one result file per configuration in OUTDIR.
Ends by comparing every count with the per-configuration lines of twelve_check's asc output.
Usage: python src/levels_all12.py data/configs12.txt logs/check12-asc-output.txt OUTDIR [workers]
"""
import ast
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

CONFIGS, ASC, OUT = sys.argv[1:4]
SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "twelve_check_levels.py")
N = 376


def one(i):
    path = os.path.join(OUT, f"{i}.txt")
    if os.path.exists(path):
        return
    r = subprocess.run([sys.executable, SCRIPT, CONFIGS, str(i), "8"], capture_output=True, text=True)
    if r.returncode != 0 or "nodes by dim" not in r.stdout:
        print(f"config {i}: FAILED rc {r.returncode}: {r.stderr.strip()[-300:]}", flush=True)
        return
    with open(path + ".tmp", "w") as f:
        f.write(r.stdout)
    os.replace(path + ".tmp", path)
    print(r.stdout.strip(), flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, range(N)))
    want = {}
    for line in open(ASC):
        m = re.match(r"config (\d+) \(\w+, \d+ lines in D\): nodes ([\d ]+?)\s+p1", line)
        if m:
            want[int(m.group(1))] = [int(x) for x in m.group(2).split()]
    bad = missing = 0
    for i in range(N):
        path = os.path.join(OUT, f"{i}.txt")
        if not os.path.exists(path):
            missing += 1
            continue
        got = ast.literal_eval(open(path).read().split(":", 2)[2].strip())
        w = want.get(i)
        # the recount stops at dim 8; pad both with zeros to compare
        L = max(len(got), len(w or []))
        if w is None or got + [0] * (L - len(got)) != w + [0] * (L - len(w)):
            bad += 1
            print(f"MISMATCH config {i}: python {got} vs C++ {w}")
    print(f"levels recount: {N - missing} of {N} configurations, {bad} mismatches, {missing} missing")
    print("total nodes (python):", sum(sum(ast.literal_eval(open(os.path.join(OUT, f'{i}.txt')).read().split(':', 2)[2].strip()))
                                       for i in range(N) if os.path.exists(os.path.join(OUT, f"{i}.txt"))))
    sys.exit(1 if bad or missing else 0)


main()
