"""Resumable driver for eight_check (search B for 8 rank-4 components), for the share runner.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-30; adapted from nine_check_slices.py. Every stage writes per-slice or
per-range outputs and skips finished ones, so an eviction loses at most the slice or range in progress.
  python apn/eight_check_slices.py count5     5-cliques (and 6-cliques) of V5824 by plain enumeration -> count5.out
  python apn/eight_check_slices.py parta      Part A: class 0 in 16 slices, 1 in 8, 2 and 3 in 4, 4 in 2
  python apn/eight_check_slices.py merge      merge of all 34 Part A slices -> orbits.txt, merge.out
  python apn/eight_check_slices.py configs    hole sets, covers, outside families -> configs.txt, configs.out
  python apn/eight_check_slices.py selftest
  python apn/eight_check_slices.py search     asc, desc, rand 7 on configs.txt, in 8 ranges each
  python apn/eight_check_slices.py cross <name>   cross-runs on inputs that don't use this classification (CROSS)
  python apn/eight_check_slices.py maxsample [step] / orbitcov / controls / ctlreps [n]
Outputs in runs/check8.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.join(ROOT, "runs", "bin", "eight_check.exe")
OUT = os.path.join(ROOT, "runs", "check8")
R = os.path.join(ROOT, "runs")
SPLIT = {0: 16, 1: 8, 2: 4, 3: 4, 4: 2}


def done(path):
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            n = f.tell()
            f.seek(max(0, n - 16))
            return f.read().strip().endswith(b"END")
    except OSError:
        return False


def run_logged(log, args, mode="w"):
    with open(os.path.join(OUT, log), mode) as lf:
        r = subprocess.run([EXE] + args, stdout=lf, stderr=subprocess.STDOUT)
    print(f"{args[0]} exit {r.returncode}", flush=True)
    return r.returncode


def count5():
    sys.exit(run_logged("count5.out", ["count5"]))


def parta():
    for c, ns in SPLIT.items():
        for s in range(ns):
            txt = os.path.join(OUT, f"parta-c{c}-s{s}.txt")
            if done(txt):
                continue
            mx = os.path.join(OUT, f"maxc-c{c}-s{s}.txt")
            print(f"class {c} slice {s} of {ns}", flush=True)
            rc = run_logged(f"parta-c{c}-s{s}.out", ["parta", str(c), txt, mx, str(s), str(ns)])
            if rc != 0 or not done(txt):
                print(f"class {c} slice {s}: exit {rc}", flush=True)
                sys.exit(1)
    print("parta: all slices done", flush=True)


def parts():
    return [os.path.join(OUT, f"parta-c{c}-s{s}.txt") for c, ns in SPLIT.items() for s in range(ns)]


def merge():
    assert all(done(p) for p in parts()), "Part A incomplete"
    sys.exit(run_logged("merge.out", ["merge", os.path.join(OUT, "orbits.txt"), os.path.join(R, "check9", "configs.txt")] + parts()))


def configs():
    sys.exit(run_logged("configs.out", ["configs", os.path.join(OUT, "orbits.txt"), os.path.join(R, "check9", "configs.txt"),
                                        os.path.join(R, "check10", "configs.txt"), os.path.join(R, "check11", "configs.txt"),
                                        os.path.join(R, "check12", "configs.txt"), os.path.join(R, "check14", "spreads.txt"),
                                        os.path.join(R, "check13", "configs.txt"), os.path.join(OUT, "configs.txt")]))


def selftest():
    sys.exit(run_logged("selftest.out", ["selftest", os.path.join(OUT, "configs.txt")]))


def prog_ids(path):
    ids = set()
    if os.path.exists(path):
        for line in open(path):
            if line.startswith("R "):
                ids.add(int(line.split()[1]))
    return ids


def ranged(tag, mode, inp, n, nr, extra):
    size = -(-n // nr)
    for k in range(nr):
        lo, hi = k * size, min(n, (k + 1) * size) - 1
        if lo > hi:
            continue
        t = tag if nr == 1 else f"{tag}-r{k}"
        pf = os.path.join(OUT, f"{t}.prog")
        if len({i for i in prog_ids(pf) if lo <= i <= hi}) == hi - lo + 1:
            continue
        print(f"{t}: configurations {lo}..{hi}", flush=True)
        rc = run_logged(f"{t}.out", [mode, inp, pf] + extra + [str(lo), str(hi)], "a")
        if rc or len({i for i in prog_ids(pf) if lo <= i <= hi}) != hi - lo + 1:
            print(f"{t}: exit {rc}, incomplete", flush=True)
            sys.exit(1)
    print(f"{tag}: all {n} done", flush=True)


def search():
    C = os.path.join(OUT, "configs.txt")
    n = sum(1 for line in open(C) if line.startswith("CONF"))
    for name, extra in (("asc", ["asc"]), ("desc", ["desc"]), ("rand7", ["rand", "7"])):
        ranged(name, "search", C, n, 8, extra)


# cross-runs on inputs that don't use this classification: mode, input, count (None: count MAXC lines), ranges
CROSS = {
    "s9": ("search9", os.path.join("check9", "configs.txt"), 5663 * 9, 8),
    "s13": ("search13", os.path.join("check13", "configs.txt"), 31 * 1287, 8),
    "ssp": ("searchsp", os.path.join("check14", "spreads.txt"), 73 * 3003, 32),
    "maxs": ("searchmax", os.path.join("check8", "maxc-sample.txt"), None, 4),
}


def cross(name):
    mode, inp, n, nr = CROSS[name]
    inp = os.path.join(R, inp)
    if n is None:
        n = sum(1 for line in open(inp) if line.startswith("MAXC"))
    ranged(name, mode, inp, n, nr, ["asc"])


def maxsample(step=100):
    # every step-th maximal clique of Part A, in slice order
    step = int(step)
    k = 0
    with open(os.path.join(OUT, "maxc-sample.txt"), "w") as of:
        for c, ns in SPLIT.items():
            for s in range(ns):
                for line in open(os.path.join(OUT, f"maxc-c{c}-s{s}.txt")):
                    if line.startswith("MAXC"):
                        if k % step == 0:
                            of.write(line)
                        k += 1
    print(f"{k} maximal cliques, sample of {(k + step - 1) // step} (every {step}th)", flush=True)


def orbitcov():
    # every maxc file onto the orbits
    log = os.path.join(OUT, "maxall-orbits.out")
    if os.path.exists(log):
        os.remove(log)
    for c, ns in SPLIT.items():
        for s in range(ns):
            rc = run_logged("maxall-orbits.out", ["orbitcov", os.path.join(OUT, "orbits.txt"), os.path.join(OUT, "configs.txt"),
                                                  os.path.join(OUT, f"maxc-c{c}-s{s}.txt")], "a")
            if rc:
                sys.exit(rc)


# controls as nine-check §6: name, LUT file, choices per path, extra flag
CTL = [("ind-bl22", "apn_8bit.txt", 3, ["nond"]), ("bl22", "apn_8bit.txt", 3, []), ("blp22", "blp22_8bit.txt", 6, []),
       ("yuperrin", "yuperrin_qam.txt", 6, []), ("yuwangli", "yuwangli_qam.txt", 6, []), ("new", "new_apns.txt", 1, [])]


def controls():
    # a batch isn't resumable: one with fewer than two SUMMARY lines is deleted and rerun
    worst = 0
    for name, lut, nch, extra in CTL:
        txt = os.path.join(OUT, f"ctl-{name}.txt")
        if os.path.exists(txt) and sum(1 for line in open(txt) if line.startswith("SUMMARY")) == 2:
            continue
        if os.path.exists(txt):
            os.remove(txt)
        rc = run_logged(f"ctl-{name}.out", ["ctlbatch", os.path.join(OUT, "configs.txt"), os.path.join(R, "known", lut),
                                            "2", "1000000", str(nch), txt] + extra)
        worst = max(worst, rc)
    sys.exit(worst)


def ctlreps(nfun="40"):
    # BL22 functions with k >= 3 embedded into every configuration no control batch reached (second line of
    # reps-unreached.txt, from `ten_check_tally.py reps runs/check8`), in chunks of 200 reps (ctlreps-<k>.out)
    reps = open(os.path.join(OUT, "reps-unreached.txt")).read().splitlines()[1].split()
    worst = 0
    for k in range(0, len(reps), 200):
        log = f"ctlreps-{k // 200}.out"
        p = os.path.join(OUT, log)
        if os.path.exists(p) and "exit" not in open(p).read() and sum(1 for line in open(p) if line.startswith("rep ")) == len(reps[k:k + 200]):
            continue
        worst = max(worst, run_logged(log, ["ctlreps", os.path.join(OUT, "configs.txt"), os.path.join(R, "known", "apn_8bit.txt"),
                                            "1", nfun] + reps[k:k + 200]))
    sys.exit(worst)


def ctlrot(pool="40"):
    # one control per configuration no batch reached (reps-unreached.txt), BL22's first `pool` functions with k >= 3 in
    # rotation, parallel over the configurations, chunks of 400 (ctlrot-<k>.out). Uses eight_check_rot.exe (the same
    # source, built while eight_check.exe was busy with s9).
    exe = os.path.join(ROOT, "runs", "bin", "eight_check_rot.exe")
    reps = open(os.path.join(OUT, "reps-unreached.txt")).read().splitlines()[1].split()
    worst = 0
    for k in range(0, len(reps), 400):
        p = os.path.join(OUT, f"ctlrot-{k // 400}.out")
        if os.path.exists(p) and sum(1 for line in open(p) if line.startswith("rep ")) == len(reps[k:k + 400]):
            continue
        with open(p, "w") as lf:
            r = subprocess.run([exe, "ctlrot", os.path.join(OUT, "configs.txt"), os.path.join(R, "known", "apn_8bit.txt"), "1", pool]
                               + reps[k:k + 400], stdout=lf, stderr=subprocess.STDOUT)
        print(f"ctlrot chunk {k // 400}: exit {r.returncode}", flush=True)
        worst = max(worst, r.returncode)
    sys.exit(worst)


# the numpy re-search sample: the first configuration of each type, every configuration with a dimension-5 near-miss
# in asc (15), the three heaviest others by tries in D, and 7 drawn at random (seed 8)
NUMPY_SAMPLE = [0, 584, 1107, 1457, 3103, 6984, 9191,
                1214, 552, 941, 497, 291, 241, 905, 240, 328, 309, 1377, 522, 555, 2453, 114,
                583, 465, 557]


def numpy_sample(nproc="8"):
    import random
    ids = list(NUMPY_SAMPLE)
    rng = random.Random(8)
    while len(ids) < len(NUMPY_SAMPLE) + 7:
        i = rng.randrange(9372)
        if i not in ids:
            ids.append(i)
    py = os.path.join(R, "venv", "Scripts", "python.exe")
    todo = [i for i in ids if not (os.path.exists(os.path.join(OUT, f"numpy-{i}.out")) and
                                    open(os.path.join(OUT, f"numpy-{i}.out")).read().startswith(f"config {i}:"))]
    print(f"numpy sample: {len(ids)} configurations, {len(todo)} to do", flush=True)
    import time
    todo.sort(key=lambda i: i in (114, 2453))  # the two heaviest last in the list, so popped first
    running = []
    while todo or running:
        while todo and len(running) < int(nproc):
            i = todo.pop()
            f = open(os.path.join(OUT, f"numpy-{i}.out"), "w")
            running.append((subprocess.Popen([py, os.path.join(ROOT, "apn", "eight_check_numpy.py"), "search",
                                              os.path.join(OUT, "configs.txt"), str(i)], stdout=f, stderr=subprocess.STDOUT), f))
        for p, f in list(running):
            if p.poll() is not None:
                f.close()
                running.remove((p, f))
        time.sleep(5)
    print("numpy sample done", flush=True)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    {"count5": count5, "parta": parta, "merge": merge, "configs": configs, "selftest": selftest, "search": search,
     "cross": lambda: cross(sys.argv[2]), "maxsample": lambda: maxsample(*sys.argv[2:3]), "orbitcov": orbitcov,
     "controls": controls, "ctlreps": lambda: ctlreps(*sys.argv[2:3]), "numpy": lambda: numpy_sample(*sys.argv[2:3]), "ctlrot": lambda: ctlrot(*sys.argv[2:3])}[sys.argv[1]]()
