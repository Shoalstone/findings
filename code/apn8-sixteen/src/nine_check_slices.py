"""Resumable driver for nine_check's sliced modes (search B for i = 8), for the share runner.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-28. Runs every slice whose output doesn't yet end in END, one at a time,
so an eviction loses at most the slice in progress; rerunning the job skips finished slices.
  python apn/nine_check_slices.py parta     Part A: class 0 in 16 slices, 1 in 8, 2 and 3 in 4, 4 in 2 (runs/check9/parta-c<c>-s<s>.txt)
  python apn/nine_check_slices.py merge     merge of all 34 Part A slices -> runs/check9/orbits.txt, merge.out
"""
import glob
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.join(ROOT, "runs", "bin", "nine_check.exe")
OUT = os.path.join(ROOT, "runs", "check9")
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


def parta():
    for c, ns in SPLIT.items():
        for s in range(ns):
            txt = os.path.join(OUT, f"parta-c{c}-s{s}.txt")
            if done(txt):
                continue
            mx = os.path.join(OUT, f"maxc-c{c}-s{s}.txt")
            log = os.path.join(OUT, f"parta-c{c}-s{s}.out")
            print(f"class {c} slice {s} of {ns}", flush=True)
            with open(log, "w") as lf:
                r = subprocess.run([EXE, "parta", str(c), txt, mx, str(s), str(ns)], stdout=lf, stderr=subprocess.STDOUT)
            if r.returncode != 0 or not done(txt):
                print(f"class {c} slice {s}: exit {r.returncode}", flush=True)
                sys.exit(1)
    print("parta: all slices done", flush=True)


def merge():
    parts = [os.path.join(OUT, f"parta-c{c}-s{s}.txt") for c, ns in SPLIT.items() for s in range(ns)]
    assert all(done(p) for p in parts), "Part A incomplete"
    with open(os.path.join(OUT, "merge.out"), "w") as lf:
        r = subprocess.run([EXE, "merge", os.path.join(OUT, "orbits.txt"), os.path.join(ROOT, "runs", "check10", "configs.txt")] + parts,
                           stdout=lf, stderr=subprocess.STDOUT)
    print(f"merge exit {r.returncode}", flush=True)
    sys.exit(r.returncode)


def run_logged(log, args, mode="w"):
    with open(os.path.join(OUT, log), mode) as lf:
        r = subprocess.run([EXE] + args, stdout=lf, stderr=subprocess.STDOUT)
    print(f"{args[0]} exit {r.returncode}", flush=True)
    return r.returncode


def configs():
    R = os.path.join(ROOT, "runs")
    sys.exit(run_logged("configs.out", ["configs", os.path.join(OUT, "orbits.txt"), os.path.join(R, "check10", "configs.txt"),
                                        os.path.join(R, "check11", "configs.txt"), os.path.join(R, "check12", "configs.txt"),
                                        os.path.join(R, "check14", "spreads.txt"), os.path.join(R, "check13", "configs.txt"),
                                        os.path.join(OUT, "configs.txt")]))


def search():
    # the three orders on configs.txt, one after another; each resumes from its .prog
    C = os.path.join(OUT, "configs.txt")
    for name, extra in (("asc", ["asc"]), ("desc", ["desc"]), ("rand7", ["rand", "7"])):
        rc = run_logged(f"{name}.out", ["search", C, os.path.join(OUT, f"{name}.prog")] + extra, "a")
        if rc:
            sys.exit(rc)


def prog_ids(path):
    ids = set()
    if os.path.exists(path):
        for line in open(path):
            if line.startswith("R "):
                ids.add(int(line.split()[1]))
    return ids


# cross-runs on inputs that don't use this classification: mode, input, count, ranges
CROSS = {
    "ssp": ("searchsp", os.path.join("check14", "spreads.txt"), 73 * 2002, 8),
    "s13": ("search13", os.path.join("check13", "configs.txt"), 31 * 715, 1),
    "s10": ("search10", os.path.join("check10", "configs.txt"), 1777 * 10, 1),
    "maxs1000": ("searchmax", os.path.join("check9", "maxc-s1000.txt"), None, 1),
}


def cross(name):
    mode, inp, n, nr = CROSS[name]
    inp = os.path.join(ROOT, "runs", inp)
    if n is None:
        n = sum(1 for line in open(inp) if line.startswith("MAXC"))
    size = -(-n // nr)
    for k in range(nr):
        lo, hi = k * size, min(n, (k + 1) * size) - 1
        tag = name if nr == 1 else f"{name}-r{k}"
        pf = os.path.join(OUT, f"{tag}.prog")
        if len({i for i in prog_ids(pf) if lo <= i <= hi}) == hi - lo + 1:
            continue
        print(f"{tag}: configurations {lo}..{hi}", flush=True)
        rc = run_logged(f"{tag}.out", [mode, inp, pf, "asc", str(lo), str(hi)], "a")
        if rc or len({i for i in prog_ids(pf) if lo <= i <= hi}) != hi - lo + 1:
            print(f"{tag}: exit {rc}, incomplete", flush=True)
            sys.exit(1)
    print(f"{name}: all {n} done", flush=True)


def maxsample(step=1000):
    # every step-th maximal clique of Part A, in slice order
    k = 0
    with open(os.path.join(OUT, f"maxc-s{step}.txt"), "w") as of:
        for c, ns in SPLIT.items():
            for s in range(ns):
                for line in open(os.path.join(OUT, f"maxc-c{c}-s{s}.txt")):
                    if line.startswith("MAXC"):
                        if k % step == 0:
                            of.write(line)
                        k += 1
    print(f"{k} maximal cliques, sample of {(k + step - 1) // step}", flush=True)


def orbitcov():
    # every maxc file onto the orbits; appends to maxall-orbits.out
    log = os.path.join(OUT, "maxall-orbits.out")
    if os.path.exists(log):
        os.remove(log)
    for c, ns in SPLIT.items():
        for s in range(ns):
            rc = run_logged("maxall-orbits.out", ["orbitcov", os.path.join(OUT, "orbits.txt"), os.path.join(OUT, "configs.txt"),
                                                  os.path.join(OUT, f"maxc-c{c}-s{s}.txt")], "a")
            if rc:
                sys.exit(rc)


# controls as ten-check §6: name, LUT file, choices per path, extra flag
CTL = [("ind-bl22", "apn_8bit.txt", 3, ["nond"]), ("bl22", "apn_8bit.txt", 3, []), ("blp22", "blp22_8bit.txt", 6, []),
       ("yuperrin", "yuperrin_qam.txt", 6, []), ("yuwangli", "yuwangli_qam.txt", 6, []), ("new", "new_apns.txt", 1, [])]


def controls():
    # each control batch isn't resumable: one with fewer than two SUMMARY lines is deleted and rerun
    worst = 0
    for name, lut, nch, extra in CTL:
        txt = os.path.join(OUT, f"ctl-{name}.txt")
        if os.path.exists(txt) and sum(1 for line in open(txt) if line.startswith("SUMMARY")) == 2:
            continue
        if os.path.exists(txt):
            os.remove(txt)
        rc = run_logged(f"ctl-{name}.out", ["ctlbatch", os.path.join(OUT, "configs.txt"), os.path.join(ROOT, "runs", "known", lut),
                                            "2", "1000000", str(nch), txt] + extra)
        worst = max(worst, rc)
    sys.exit(worst)


def ctlreps(nfun="40"):
    # BL22 functions with k >= 3 embedded into every configuration no control batch reached (second line of
    # reps-unreached.txt, from `ten_check_tally.py reps runs/check9`); not resumable, rerun whole if evicted
    reps = open(os.path.join(OUT, "reps-unreached.txt")).read().splitlines()[1].split()
    sys.exit(run_logged("ctlreps.out", ["ctlreps", os.path.join(OUT, "configs.txt"), os.path.join(ROOT, "runs", "known", "apn_8bit.txt"),
                                        "1", nfun] + reps))


{"ctlreps": lambda: ctlreps(*sys.argv[2:3]), "parta": parta, "merge": merge, "configs": configs, "search": search, "maxsample": maxsample, "orbitcov": orbitcov,
 "cross": lambda: cross(sys.argv[2]), "controls": controls}[sys.argv[1]]()
