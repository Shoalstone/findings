"""Re-runs for eight from the bundle's own files (cwd = bundle root), with fresh builds of its sources in rr/.
Each step: rr/out/<name>.log, times in rr/out/times.txt. Skips steps whose log ends with an 'exit' line (resumable)."""
import os
import subprocess
import sys
import time

R = os.path.dirname(os.path.abspath(__file__))
O = os.path.join(R, "out")
B = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # the bundle
K = os.environ.get("KNOWN", "known")  # folder holding the known-function lists (apn_8bit.txt), downloaded
VPY = os.environ.get("VPY", "python")  # a Python with numpy
E = os.path.join(R, "eight_check.exe")
P = os.path.join(R, "pspread3.exe")
C = os.path.join(R, "cover.exe")
os.makedirs(O, exist_ok=True)
os.chdir(B)

parts = [f"data/parta8-c{c}-s{s}.txt" for c, ns in {0: 16, 1: 8, 2: 4, 3: 4, 4: 2}.items() for s in range(ns)]
o = lambda f: os.path.join(O, f)
STEPS = [
    # search A
    ("pspread3-selftest8", [P, "selftest", "8", "data/p8.txt"]),
    ("pspread3-holes8", [P, "holes", "8", "data/p8.txt", o("holes8.txt")]),
    ("pspread3-p17", [P, "enum", "17", o("p17.txt")]),
    ("cover-p8-c0", [C, o("c0.txt")]),
    ("cover-p8-c7", [C, o("c7.txt")]),
    # search B
    ("e8-parta4s0", [E, "parta", "4", o("parta-c4-s0.txt"), o("maxc-c4-s0.txt"), "0", "2"]),
    ("e8-parta4s1", [E, "parta", "4", o("parta-c4-s1.txt"), o("maxc-c4-s1.txt"), "1", "2"]),
    ("e8-parta3s3", [E, "parta", "3", o("parta-c3-s3.txt"), o("maxc-c3-s3.txt"), "3", "4"]),
    ("e8-merge", [E, "merge", o("orbits.txt"), "data/configs9.txt"] + parts),
    ("e8-configs", [E, "configs", o("orbits.txt"), "data/configs9.txt", "data/configs10.txt", "data/configs11.txt",
                    "data/configs12.txt", "data/spreads.txt", "data/configs.txt", o("configs8.txt")]),
    ("e8-mapa", [E, "mapa", "data/configs8.txt", "data/p8.txt"]),
    ("e8-orbitcov-s100", [E, "orbitcov", "data/orbits8.txt", "data/configs8.txt", "data/maxc8-s100.txt"]),
    ("e8-ctlnond", [E, "ctlbatch", "data/configs8.txt", os.path.join(K, "apn_8bit.txt"), "2", "1000000", "3",
                    o("ctl-ind-bl22.txt"), "nond"]),
    ("e8-ctlrot", [E, "ctlrot", "data/configs8.txt", os.path.join(K, "apn_8bit.txt"), "1", "40", "91", "92", "9364", "9365"]),
    ("e8-asc-8204-9371", [E, "search", "data/configs8.txt", o("asc.prog"), "asc", "8204", "9371"]),
    ("e8-desc-2344-2843", [E, "search", "data/configs8.txt", o("desc.prog"), "desc", "2344", "2843"]),
    ("e8-rand7-2344-2843", [E, "search", "data/configs8.txt", o("rand7.prog"), "rand", "7", "2344", "2843"]),
    ("e8-maxs-0-499", [E, "searchmax", "data/maxc8-s100.txt", o("maxs.prog"), "asc", "0", "499"]),
    ("e8-s9-0-99", [E, "search9", "data/configs9.txt", o("s9.prog"), "asc", "0", "99"]),
    ("e8-numpy-9191", [VPY, "src/eight_check_numpy.py", "search", "data/configs8.txt", "9191"]),
    ("e8-numpy-ctl", [VPY, "src/eight_check_numpy.py", "control", os.path.join(K, "apn_8bit.txt"), "1", "2", "1"]),
]

# the cover chunks: data/p8.txt cut as the run cut it (8 chunks of 1,172 lines, the last 1,168)
lines = open("data/p8.txt").read().splitlines(True)
for k in (0, 7):
    with open(o(f"c{k}.txt"), "w", newline="") as f:
        f.writelines(lines[k * 1172:(k + 1) * 1172])

only = set(sys.argv[1:])
for name, cmd in STEPS:
    if only and name not in only:
        continue
    log = o(name + ".log")
    if os.path.exists(log) and open(log, errors="replace").read().rstrip().splitlines()[-1:] and \
            open(log, errors="replace").read().rstrip().splitlines()[-1].startswith("# exit"):
        continue
    t0 = time.time()
    with open(log, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    with open(log, "a") as lf:
        lf.write(f"# exit {rc}\n")
    with open(o("times.txt"), "a") as tf:
        tf.write(f"{name}: exit {rc}, {time.time() - t0:.0f} s\n")
    print(f"{name}: exit {rc}, {time.time() - t0:.0f} s", flush=True)
