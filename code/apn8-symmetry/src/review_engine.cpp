// review_engine: reviewer's checks of inv2's lift engine against direct kernel computations (subagent, Claude Opus
// 5.5, 2026-10-24). Includes inv2.cpp unchanged; a static initialiser runs review_main() and exits before inv2's main.
//   review_engine lvl <rep> [ntriples]   level_ok(j) == [max_kernel(F + B_0..B_{j-1}) <= 5 - j] on every g at each
//                                        level along a few sampled paths; solve() vs brute force over all 4096 g3
//   review_engine sols <every>           every <every>-th line of inv2-all.txt: recover (g0..g3), check the DFS path
//                                        accepts it and solve() == brute force over g3 on that triple
// Single-threaded.
#include <cstdlib>
#include <random>
#include <set>
#include "inv2.cpp"

static const char* REPS = "results/inv2-reps.txt";
static const char* ALL = "results/inv2-all.txt";

static std::vector<u32> load_reps() { std::vector<u32> r; FILE* f = fopen(REPS, "r"); unsigned L; long sz; while (fscanf(f, "%x %ld", &L, &sz) == 2) r.push_back(L); fclose(f); return r; }

struct Rec { std::set<u32>* s; void operator()(u32* b8, int*) { s->insert(b8[7]); } };

static int crit_kernel(const Tables& T, const u32* g, int j) {
    u64 m[8]; for (int k = 0; k < 4; k++) m[k] = coef_to_mat(T.S[k]);
    for (int t = 0; t < j; t++) m[4 + t] = coef_to_mat(T.lift(t, g[t]));
    return max_kernel(m, 4 + j);
}
static std::set<u32> brute_g3(const Tables& T, const u32* g) {
    std::set<u32> s; u64 m[8]; for (int k = 0; k < 4; k++) m[k] = coef_to_mat(T.S[k]);
    for (int t = 0; t < 3; t++) m[4 + t] = coef_to_mat(T.lift(t, g[t]));
    for (u32 x = 0; x < 4096; x++) { u32 l = T.lift(3, x); m[7] = coef_to_mat(l); if (max_kernel(m, 8) == 1) s.insert(l); }
    return s;
}

static int run_lvl(int ri, int ntrip) {
    auto reps = load_reps(); u32 C[4]; F_of(reps[ri], C);
    Tables* T = new Tables; if (!T->init(C)) { printf("rep %d init failed\n", ri); return 1; }
    std::set<u32> got; Rec rc{&got}; Dfs<Rec>* d = new Dfs<Rec>(*T, rc);
    std::mt19937 rng(ri * 7919 + 1);
    long long mism[4] = {0, 0, 0, 0}, tested[4] = {0, 0, 0, 0}, pass[4] = {0, 0, 0, 0};
    u32 g[4];
    std::vector<u32> p1;
    for (u32 a = 0; a < 4096; a++) { g[0] = a; d->place(0, a); bool lk = d->level_ok(1), cr = crit_kernel(*T, g, 1) <= 4; tested[1]++; pass[1] += lk; if (lk != cr) mism[1]++; if (lk) p1.push_back(a); }
    std::shuffle(p1.begin(), p1.end(), rng); if (p1.size() > 3) p1.resize(3);
    std::vector<std::array<u32, 2>> p2;
    for (u32 a : p1) { g[0] = a; d->place(0, a); for (u32 b = 0; b < 4096; b++) { g[1] = b; d->place(1, b); bool lk = d->level_ok(2), cr = crit_kernel(*T, g, 2) <= 3; tested[2]++; pass[2] += lk; if (lk != cr) mism[2]++; if (lk) p2.push_back({a, b}); } }
    std::shuffle(p2.begin(), p2.end(), rng); if (p2.size() > 3) p2.resize(3);
    std::vector<std::array<u32, 3>> p3, f3;
    for (auto& ab : p2) { g[0] = ab[0]; g[1] = ab[1]; d->place(0, g[0]); d->place(1, g[1]);
        for (u32 c = 0; c < 4096; c++) { g[2] = c; d->place(2, c); bool lk = d->level_ok(3), cr = crit_kernel(*T, g, 3) <= 2; tested[3]++; pass[3] += lk; if (lk != cr) mism[3]++; if (lk) p3.push_back({ab[0], ab[1], c}); else f3.push_back({ab[0], ab[1], c}); } }
    std::shuffle(p3.begin(), p3.end(), rng); if ((int)p3.size() > ntrip) p3.resize(ntrip);
    std::shuffle(f3.begin(), f3.end(), rng); if (f3.size() > 3) f3.resize(3);
    int smis = 0, nsol = 0, ntr = 0;
    auto cmp = [&](const std::array<u32, 3>& t, bool passed) {
        for (int i = 0; i < 3; i++) { g[i] = t[i]; d->place(i, t[i]); }
        got.clear(); if (passed) d->solve();
        auto bf = brute_g3(*T, g); ntr++; nsol += (int)bf.size(); if (bf != got) smis++;
    };
    for (auto& t : p3) cmp(t, true);
    for (auto& t : f3) cmp(t, false);
    printf("rep %2d %06x: level 1 %lld/%lld pass, %lld mismatch; level 2 %lld/%lld, %lld; level 3 %lld/%lld, %lld; solve vs brute on %d triples (%d solutions): %d mismatch\n",
           ri, reps[ri], pass[1], tested[1], mism[1], pass[2], tested[2], mism[2], pass[3], tested[3], mism[3], ntr, nsol, smis);
    return (mism[1] || mism[2] || mism[3] || smis) ? 1 : 0;
}

static bool recover(const Tables& T, int t, u32 B, u32& g) {
    // B ^ lam_t in span(Y): echelon with tracking
    u32 piv[28] = {0}, src[28] = {0};
    for (int i = 0; i < 12; i++) { u32 v = T.Y[i], s = 1u << i; for (int b = 27; b >= 0 && v; b--) if ((v >> b & 1) && piv[b]) { v ^= piv[b]; s ^= src[b]; } if (v) { int b = 31 - (int)__lzcnt(v); piv[b] = v; src[b] = s; } }
    u32 v = B ^ T.lam[t], s = 0; for (int b = 27; b >= 0 && v; b--) if ((v >> b & 1) && piv[b]) { v ^= piv[b]; s ^= src[b]; }
    if (v) return false; g = s; return true;
}
static int run_sols(int every) {
    auto reps = load_reps(); FILE* f = fopen(ALL, "r"); char line[512]; int n = 0, tested = 0, bad = 0, total_sol = 0;
    std::map<u32, Tables*> tabs;
    while (fgets(line, sizeof line, f)) {
        if (!strstr(line, "basis")) continue; if (n++ % every) continue;
        unsigned L; sscanf(strstr(line, " F ") + 3, "%x", &L);
        u32 b[8]; const char* p = strstr(line, "basis") + 5; for (int i = 0; i < 8; i++) { unsigned v; int c; sscanf(p, "%x%n", &v, &c); b[i] = v; p += c; }
        if (!tabs.count(L)) { u32 C[4]; F_of(L, C); Tables* T = new Tables; if (!T->init(C)) { printf("init failed\n"); return 1; } tabs[L] = T; }
        Tables& T = *tabs[L];
        u32 g[4]; bool ok = true; for (int t = 0; t < 4; t++) ok &= recover(T, t, b[4 + t], g[t]);
        for (int k = 0; k < 4; k++) ok &= (b[k] == T.S[k]);
        std::set<u32> got; Rec rc{&got}; Dfs<Rec> d(T, rc);
        d.place(0, g[0]); bool l1 = d.level_ok(1); d.place(1, g[1]); bool l2 = d.level_ok(2); d.place(2, g[2]); bool l3 = d.level_ok(3);
        if (l1 && l2 && l3) d.solve();
        auto bf = brute_g3(T, g); total_sol += (int)bf.size();
        bool in = got.count(b[7]) > 0;
        tested++; if (!ok || !l1 || !l2 || !l3 || !in || bf != got) { bad++; printf("  line %d: recover %d level %d%d%d in-solve %d solve==brute %d\n", n, ok, l1, l2, l3, in, bf == got); }
    }
    fclose(f);
    printf("sols: tested %d lines, bad %d, brute solutions on those triples %d\n", tested, bad, total_sol);
    return bad ? 1 : 0;
}

static int review_main() {
    if (setup(false)) { fprintf(stderr, "setup failed\n"); exit(1); }
    make_spread();
    int argc = __argc; char** argv = __argv; int r = 2;
    if (argc >= 3 && !strcmp(argv[1], "lvl")) {
        int nt = argc > 3 ? atoi(argv[3]) : 20;
        if (!strcmp(argv[2], "all")) { r = 0; for (int i = 0; i < 18; i++) r |= run_lvl(i, nt); } else r = run_lvl(atoi(argv[2]), nt);
    }
    else if (argc >= 3 && !strcmp(argv[1], "sols")) r = run_sols(atoi(argv[2]));
    else fprintf(stderr, "usage: review_engine lvl <rep|all> [ntriples] | sols <every>\n");
    fflush(stdout); exit(r);
}
static int review_run = review_main();
