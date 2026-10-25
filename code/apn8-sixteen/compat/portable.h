// portable.h: lets the sources (written for MSVC) build with gcc or clang. Force-included by the gcc/clang build
// line (-include compat/portable.h); MSVC builds don't use it. It changes no logic: each shim reproduces the MSVC
// intrinsic's result exactly, including at 0.
//   __popcnt(x), __lzcnt(x), _tzcnt_u32(x)   MSVC bit intrinsics on unsigned int (lzcnt/tzcnt of 0 = 32)
//   __argc, __argv                           MSVC CRT globals, read by review_engine.cpp and review_inv2.cpp
//                                            from a static initialiser (so before main)
#pragma once
#if !defined(_MSC_VER)

// On x86, pull in the compiler's own intrinsic headers first, so their declarations (some compilers declare
// __lzcnt or _tzcnt_u32 there) come before the macros below and aren't rewritten by them.
#if defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
#endif

static inline unsigned quire_popcnt(unsigned x) { return (unsigned)__builtin_popcount(x); }
static inline unsigned quire_lzcnt(unsigned x) { return x ? (unsigned)__builtin_clz(x) : 32u; }
static inline unsigned quire_tzcnt(unsigned x) { return x ? (unsigned)__builtin_ctz(x) : 32u; }
#undef __popcnt
#undef __lzcnt
#undef _tzcnt_u32
#define __popcnt(x) quire_popcnt((unsigned)(x))
#define __lzcnt(x) quire_lzcnt((unsigned)(x))
#define _tzcnt_u32(x) quire_tzcnt((unsigned)(x))

#if defined(__APPLE__)
#include <crt_externs.h>
#define __argc (*_NSGetArgc())
#define __argv (*_NSGetArgv())
#elif defined(__linux__)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
// The command line from /proc/self/cmdline, read on first use (works during static initialisation).
struct QuireArgs {
    int argc = 0; char** argv = nullptr;
    QuireArgs() {
        static std::vector<char> buf; static std::vector<char*> ptrs;
        FILE* f = fopen("/proc/self/cmdline", "rb"); if (!f) return;
        char tmp[4096]; size_t n;
        while ((n = fread(tmp, 1, sizeof tmp, f)) > 0) buf.insert(buf.end(), tmp, tmp + n);
        fclose(f);
        for (size_t i = 0; i < buf.size(); i += strlen(&buf[i]) + 1) ptrs.push_back(&buf[i]);
        ptrs.push_back(nullptr);
        argc = (int)ptrs.size() - 1; argv = ptrs.data();
    }
};
static inline QuireArgs& quire_args() { static QuireArgs a; return a; }
#define __argc (quire_args().argc)
#define __argv (quire_args().argv)
#endif

#endif // !_MSC_VER
