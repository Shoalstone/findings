#!/bin/sh
# build.sh [program ...]   Linux or macOS, g++ or clang++ with OpenMP. Builds every program (or the ones named)
# into bin/ (or $BIN). CXX picks the compiler (default g++); CXXFLAGS adds flags (e.g. -march=native).
# macOS: Apple's clang has no OpenMP. Either use Homebrew's (brew install llvm libomp; CXX=$(brew --prefix llvm)/bin/clang++)
# or Homebrew gcc (CXX=g++-14). compat/portable.h maps the MSVC intrinsics the sources use onto gcc/clang builtins.
set -e
cd "$(dirname "$0")"
CXX=${CXX:-g++}
BIN=${BIN:-bin}
mkdir -p "$BIN"
sed -n 14,205p src/c3k4.cpp > src/c3k4_engine.inc   # review_odd_engine includes c3k4.cpp's Engine from this file
PROGS=${*:-orbit5 fix5check symp c3k4 inv2 aut8 classify8 review_engine review_inv2 review_odd_brute review_odd_engine}
for p in $PROGS; do
    echo "$p"
    $CXX -O2 -fopenmp -std=c++17 $CXXFLAGS -include compat/portable.h -Icompat "src/$p.cpp" -o "$BIN/$p"
done
