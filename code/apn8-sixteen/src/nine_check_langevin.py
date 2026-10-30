"""Search B's 5,442 hole sets for i = 8 against Langevin and Hou's count of partial spread functions in eight variables.

The hole set of a configuration is the complement of the union of its nine solids, so the classes of hole sets up to
GL(8,2) are the classes of 9-solid partial spreads by the point set they cover, which Langevin and Hou counted
(Counting partial spread functions in eight variables, IEEE Trans. Inf. Theory 57(4), 2011; data at
https://langevin.univ-tln.fr/project/spread/psp.html). Their fix-9.txt gives the number of classes for each stabiliser
order. This compares it with the |Aut H| on the HOLE lines of data/configs9.txt, and the total number of point sets.

  python src/nine_check_langevin.py PATH/TO/fix-9.txt [data/configs9.txt]
Expected: 5442 classes on both sides, identical distribution, total 24520650576127040978944.
Rowen (Claude Opus 5.5), 2026-10-30.
"""
import collections, os, re, sys

GL82 = 5348063769211699200
fix = sys.argv[1]
conf = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "configs9.txt")
ours = collections.Counter(int(re.search(r"\baut (\d+)", l).group(1)) for l in open(conf) if l.startswith("HOLE"))
theirs, total = collections.Counter(), None
for l in open(fix):
    m = re.match(r"fix=(\d+)\s+number of class (\d+)", l)
    if m:
        theirs[int(m.group(1))] = int(m.group(2))
    m = re.match(r"Total=(\d+)", l)
    if m:
        total = int(m.group(1))
print(f"classes: ours {sum(ours.values())}, Langevin-Hou {sum(theirs.values())}")
diff = {k: (ours[k], theirs[k]) for k in set(ours) | set(theirs) if ours[k] != theirs[k]}
print("stabiliser orders: identical distribution" if not diff else f"stabiliser orders differ: {diff}")
mine = sum(GL82 // k * v for k, v in ours.items())
print(f"point sets: ours {mine}, Langevin-Hou {total}, {'equal' if mine == total else 'DIFFERENT'}")
