"""fourteen_check's spread representatives (spreads.txt) as cover.cpp input: for i = 3 the 14 rank-4 radicals are the
graphs of the listed N, the holes the rest (X, Y, graph(I)). Both programs store N row-major (row r in bits 4r..4r+3)
and put x in the low nibble, so graph(N) has basis (e_j, N e_j).
usage: python spreads_to_cover.py spreads.txt > i3-reps.txt
Written for this bundle (Claude Opus 5.5, 2026-10-25); it reproduces the i3-reps.txt that cover.cpp was run on.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cover_prep import graph_basis, fmt

for line in open(sys.argv[1]):
    if not line.startswith('SPREAD'):
        continue
    head, _, mats = line.partition(':')
    k = int(head.split()[1])
    print(fmt(f"i3rep-{k}", [graph_basis(int(h, 16)) for h in mats.split()]))
