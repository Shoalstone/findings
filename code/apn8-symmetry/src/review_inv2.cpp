// review_inv2: reviewer's checks on inv2 (subagent, Claude Opus 5.5, 2026-10-24). Includes inv2.cpp unchanged;
// inv2's own main is bypassed by a static initialiser that runs review_main() and exits.
//   review_inv2 gold <ntrials>        Gold x^3 with Frobenius x^16: random C(A)-conjugate -> normalise -> BFS to a rep
//                                     -> check the engine's DFS path accepts it
//   review_inv2 known <ntrials>       same, starting from the first type-5 line of results/inv2-search-12.txt
//   review_inv2 brute <rep> <npairs>  per (g0,g1) pair: engine completions vs unpruned brute force over (g2,g3)
//   review_inv2 actl <n>              act_L vs the direct subspace F(L)∘g for random L, g in C(A)
#include <cstdlib>
#include <random>
#include "inv2.cpp"

static const char* REPS = "results/inv2-reps.txt";
static const char* KNOWN = "results/inv2-search-12.txt";
static std::mt19937_64 rng(20261024);

static std::vector<u32> load_reps() { std::vector<u32> r; FILE* f = fopen(REPS, "r"); unsigned L; long sz; while (fscanf(f, "%x %ld", &L, &sz) == 2) r.push_back(L); fclose(f); return r; }
static void rand_g0(u8* g) { for (;;) { for (int i = 0; i < 4; i++) g[i] = (u8)(rng() & 15); if (rank8(g, 4) == 4) return; } }
static bool same_space(const u32* a, int na, const u32* b, int nb) {
    Ech28 ea, eb; for (int i = 0; i < na; i++) ea.add(a[i]); for (int i = 0; i < nb; i++) eb.add(b[i]);
    if (ea.n != eb.n) return false; for (int i = 0; i < nb; i++) if (ea.reduce(b[i])) return false; return true;
}
static bool is_apn(const u32* b) { u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(b[i]); return max_kernel(m, 8) == 1; }
static bool a_invariant(const u32* b) { Ech28 e; for (int i = 0; i < 8; i++) e.add(b[i]); for (int i = 0; i < 8; i++) if (e.reduce(act(b[i]))) return false; return true; }

// generators in group_gens' order, as Maps
static std::vector<Map> gen_maps() {
    std::vector<Map> g; u8 I4[4] = {1, 2, 4, 8}, Z[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { u8 g1[4] = {0, 0, 0, 0}; g1[i] = (u8)(1 << j); g.push_back(cmap(I4, g1)); }
    for (auto& s : stab_elements()) g.push_back(cmap(s.data(), Z));
    return g;
}

// standard-coordinate A-invariant APN space b (8 coef vectors): random C(A) conjugate, normalise, BFS to a rep,
// then check the engine's DFS path through its (g0..g3) reaches it. Returns true if found.
static bool pipeline(const u32* bin, const std::vector<u32>& reps, const std::vector<Map>& gm, const std::vector<Aff>& ga, bool randomise) {
    u32 b[8]; memcpy(b, bin, sizeof b);
    if (randomise) { u8 g0[4], g1[4]; rand_g0(g0); for (int i = 0; i < 4; i++) g1[i] = (u8)(rng() & 15); Map h = cmap(g0, g1); for (int i = 0; i < 8; i++) b[i] = apply(h, b[i]); }
    if (!a_invariant(b) || !is_apn(b)) { printf("  conjugate not A-invariant APN\n"); return false; }
    u32 C[8]; int nc = 0; { Ech28 e; u32 c = 0; for (unsigned g = 1; g < 256; g++) { c ^= b[_tzcnt_u32(g)]; u32 d = Dop(c); if (d && e.add(d)) C[nc++] = d; } }
    if (nc != 4) { printf("  dim D(S8) = %d\n", nc); return false; }
    // g0 in GL(4,2) taking b(F) to E (search all of GL(4,2))
    u8 Z[4] = {0, 0, 0, 0}; Map G; bool got = false;
    for (unsigned q = 0; q < 65536 && !got; q++) {
        u8 g0[4] = {(u8)(q & 15), (u8)(q >> 4 & 15), (u8)(q >> 8 & 15), (u8)(q >> 12)}; if (rank8(g0, 4) != 4) continue;
        Map t = cmap(g0, Z); bool ok = true; for (int k = 0; k < 4 && ok; k++) if (Ebas_inv[bpart(apply(t, C[k]))] == 0xff) ok = false;
        if (ok) { G = t; got = true; }
    }
    if (!got) { printf("  b(F) is not GL(4,2)-conjugate to E\n"); return false; }
    for (int i = 0; i < 8; i++) b[i] = apply(G, b[i]);
    for (int k = 0; k < 4; k++) C[k] = apply(G, C[k]);
    u32 L = 0;
    for (int m = 0; m < 4; m++) { bool f = false; for (unsigned c = 1; c < 16 && !f; c++) { u32 Y = 0; for (int k = 0; k < 4; k++) if (c >> k & 1) Y ^= C[k]; if (bpart(Y) == Eb[m]) { L |= (u32)cpart(Y) << (6 * m); f = true; } } if (!f) return false; }
    { u32 FL[4]; F_of(L, FL); if (!same_space(FL, 4, C, 4)) { printf("  F_of(L) != F\n"); return false; } }
    // BFS with parents until a rep
    std::set<u32> repset(reps.begin(), reps.end());
    std::vector<u32> par(1u << 24, 0xffffffffu); std::vector<uint16_t> via(1u << 24, 0);
    std::vector<u32> q{L}; par[L] = L; u32 hit = 0xffffffffu;
    if (repset.count(L)) hit = L;
    for (size_t h = 0; h < q.size() && hit == 0xffffffffu; h++) for (size_t gi = 0; gi < ga.size(); gi++) {
        u32 y = aff_apply(ga[gi], q[h]); if (par[y] != 0xffffffffu) continue; par[y] = q[h]; via[y] = (uint16_t)gi; q.push_back(y);
        if (repset.count(y)) { hit = y; break; }
    }
    if (hit == 0xffffffffu) { printf("  L %06x reaches no rep\n", L); return false; }
    std::vector<int> path; for (u32 y = hit; y != L; y = par[y]) path.push_back(via[y]);
    std::reverse(path.begin(), path.end());
    for (int gi : path) for (int i = 0; i < 8; i++) b[i] = apply(gm[gi], b[i]);
    int ri = (int)(std::find(reps.begin(), reps.end(), hit) - reps.begin());
    u32 CR[4]; F_of(hit, CR);
    { u32 D8[8]; int nd = 0; Ech28 e; u32 c = 0; for (unsigned g = 1; g < 256; g++) { c ^= b[_tzcnt_u32(g)]; u32 d = Dop(c); if (d && e.add(d)) D8[nd++] = d; } if (!same_space(D8, nd, CR, 4)) { printf("  D(S8) != F(rep) after the path\n"); return false; } }
    Tables* T = new Tables; if (!T->init(CR)) { printf("  init failed on rep %d\n", ri); return false; }
    Ech28 Fe; for (int k = 0; k < 4; k++) Fe.add(CR[k]);
    u32 gs[4];
    for (int t = 0; t < 4; t++) {
        u32 X = 0; bool f = false; u32 c = 0; for (unsigned g = 1; g < 256 && !f; g++) { c ^= b[_tzcnt_u32(g)]; if (Dop(c) == T->S[t]) { X = c; f = true; } }
        if (!f) { printf("  no preimage of C_%d in S8\n", t); return false; }
        f = false; for (u32 g = 0; g < 4096 && !f; g++) if (!Fe.reduce(T->lift(t, g) ^ X)) { gs[t] = g; f = true; }
        if (!f) { printf("  lift %d not in lambda + Y\n", t); return false; }
    }
    u32 target[8]; memcpy(target, b, sizeof target); canon8(target);
    bool found = false; long long nsol = 0;
    auto fnd = [&](u32* b8, int*) { u32 c[8]; memcpy(c, b8, sizeof c); canon8(c); nsol++; if (!memcmp(c, target, sizeof c)) found = true; };
    Dfs<decltype(fnd)>* d = new Dfs<decltype(fnd)>(*T, fnd);
    d->place(0, gs[0]); bool l1 = d->level_ok(1); d->place(1, gs[1]); bool l2 = d->level_ok(2); d->place(2, gs[2]); bool l3 = d->level_ok(3);
    if (l1 && l2 && l3) d->solve();
    printf("  L %06x -> rep %d (%06x) in %zu steps; g = %03x %03x %03x %03x; level_ok %d %d %d; solve gave %lld, target found %d\n",
           L, ri, hit, path.size(), gs[0], gs[1], gs[2], gs[3], l1, l2, l3, nsol, found);
    delete d; delete T;
    return found;
}

static u8 gmul(u8 a, u8 b) { u8 r = 0; while (b) { if (b & 1) r ^= a; a = (u8)((a << 1) ^ ((a & 0x80) ? 0x1d : 0)); b >>= 1; } return r; }

static int run_gold(int ntr) {
    auto reps = load_reps(); auto gm = gen_maps(); auto ga = group_gens(false);
    u8 F[256], fro[256]; for (unsigned x = 0; x < 256; x++) { u8 y = (u8)x; F[x] = gmul(gmul(y, y), y); u8 z = y; for (int k = 0; k < 4; k++) z = gmul(z, z); fro[x] = z; }
    u32 b[8]; for (int s = 0; s < 8; s++) { b[s] = 0; for (int k = 0; k < 28; k++) { unsigned x = 1u << KI[k], y = 1u << KJ[k]; if ((F[x ^ y] ^ F[x] ^ F[y] ^ F[0]) >> s & 1) b[s] |= 1u << k; } }
    { u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(b[i]); int cnt[5]; rank_dist(m, 8, cnt); printf("gold: APN %d, ranks 8:%d 6:%d 4:%d 2:%d\n", is_apn(b), cnt[4], cnt[3], cnt[2], cnt[1]); }
    // standard coordinates for the Frobenius involution (as inv2 check)
    u8 w[4]; int n = 0; { u8 piv[8] = {0}; for (unsigned x = 1; x < 256; x++) if (fro[x] == x) { u8 y = (u8)x; for (int bb = 7; bb >= 0 && y; bb--) if (y >> bb & 1) { if (piv[bb]) y ^= piv[bb]; else { piv[bb] = y; y = 0; } } }
        for (unsigned x = 1; x < 256 && n < 4; x++) { u8 y = (u8)x; for (int bb = 7; bb >= 0 && y; bb--) if (y >> bb & 1) { if (piv[bb]) y ^= piv[bb]; else { piv[bb] = y; w[n++] = (u8)x; y = 0; } } } }
    Map M; for (unsigned x = 0; x < 256; x++) { u8 y = 0; for (int i = 0; i < 4; i++) { if (x >> i & 1) y ^= w[i]; if (x >> (4 + i) & 1) y ^= (u8)(w[i] ^ fro[w[i]]); } M.f[x] = y; }
    make_cols(M); for (int i = 0; i < 8; i++) b[i] = apply(M, b[i]);
    printf("gold in standard coords: A-invariant %d\n", a_invariant(b));
    int ok = 0; for (int t = 0; t < ntr; t++) ok += pipeline(b, reps, gm, ga, t > 0);
    printf("gold: found %d / %d\n", ok, ntr); return 0;
}
static int run_known(int ntr) {
    auto reps = load_reps(); auto gm = gen_maps(); auto ga = group_gens(false);
    FILE* f = fopen(KNOWN, "r"); char line[512]; if (!fgets(line, sizeof line, f)) return 1; fclose(f);
    u32 b[8]; const char* p = strstr(line, "basis") + 5; for (int i = 0; i < 8; i++) { unsigned v; int c; sscanf(p, "%x%n", &v, &c); b[i] = v; p += c; }
    int ok = 0; for (int t = 0; t < ntr; t++) ok += pipeline(b, reps, gm, ga, true);
    printf("known: found %d / %d\n", ok, ntr); return 0;
}

// unpruned check: completions over (g2, g3) of a fixed (g0, g1), direct APN test on 8 forms
static long long brute_pair(const Tables& T, u32 g0, u32 g1, std::vector<u32>* out) {
    u64 base[6]; for (int k = 0; k < 4; k++) base[k] = coef_to_mat(T.S[k]); base[4] = coef_to_mat(T.lift(0, g0)); base[5] = coef_to_mat(T.lift(1, g1));
    long long cnt = 0; std::mutex mu;
    #pragma omp parallel for schedule(dynamic, 16) reduction(+:cnt) num_threads(4)
    for (int g2 = 0; g2 < 4096; g2++) {
        u64 m2 = coef_to_mat(T.lift(2, (u32)g2));
        for (u32 g3 = 0; g3 < 4096; g3++) {
            u64 m3 = coef_to_mat(T.lift(3, g3)); bool apn = true;
            for (unsigned a = 1; a < 256 && apn; a++) { u8 r[8]; for (int k = 0; k < 6; k++) r[k] = rvec(base[k], a); r[6] = rvec(m2, a); r[7] = rvec(m3, a); if (rank8(r, 8) != 7) apn = false; }
            if (apn) { cnt++; if (out) { std::lock_guard<std::mutex> lk(mu); out->push_back((u32)g2 << 12 | g3); } }
        }
    }
    return cnt;
}
static int run_brute(int ri, int npairs) {
    auto reps = load_reps(); u32 C[4]; F_of(reps[ri], C);
    Tables* T = new Tables; if (!T->init(C)) { printf("init failed\n"); return 1; }
    // pairs: some (g0,g1) with engine solutions, some pruned at level 2, some passing level 2 with no solutions
    std::vector<std::pair<u32, u32>> with, pruned, empty;
    long long cur = 0; auto fnd = [&](u32*, int*) { cur++; };
    Dfs<decltype(fnd)>* d = new Dfs<decltype(fnd)>(*T, fnd);
    for (u32 g0 = 0; g0 < 4096 && (with.size() < 3 || pruned.size() < 2 || empty.size() < 2); g0 = (g0 + 977) % 4096) {
        d->place(0, g0); if (!d->level_ok(1)) { if (pruned.size() < 1) pruned.push_back({g0, (u32)(rng() & 4095)}); continue; }
        for (u32 g1 = 0; g1 < 4096; g1++) {
            d->place(1, g1);
            if (!d->level_ok(2)) { if (pruned.size() < 2 && (rng() & 1023) == 0) pruned.push_back({g0, g1}); continue; }
            cur = 0; for (u32 g2 = 0; g2 < 4096; g2++) { d->place(2, g2); d->rec(3); }
            if (cur && with.size() < 3) with.push_back({g0, g1}); else if (!cur && empty.size() < 2) empty.push_back({g0, g1});
            if (with.size() >= 3 && empty.size() >= 2) break;
        }
        if (g0 == 0 && false) break;
        if ((int)(with.size() + pruned.size() + empty.size()) >= npairs) break;
    }
    auto test = [&](const char* kind, u32 g0, u32 g1) {
        cur = 0; d->place(0, g0); bool l1 = d->level_ok(1); d->place(1, g1); bool l2 = d->level_ok(2);
        if (l1 && l2) for (u32 g2 = 0; g2 < 4096; g2++) { d->place(2, g2); d->rec(3); }
        auto t0 = std::chrono::steady_clock::now(); long long bf = brute_pair(*T, g0, g1, nullptr);
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        printf("  %-7s g0 %03x g1 %03x: level_ok %d %d, engine %lld, brute force %lld (%.0f s) %s\n", kind, g0, g1, l1, l2, cur, bf, el, cur == bf ? "agree" : "DISAGREE");
        fflush(stdout);
    };
    printf("rep %d (%06x)\n", ri, reps[ri]);
    for (auto& p : with) test("with", p.first, p.second);
    for (auto& p : empty) test("empty", p.first, p.second);
    for (auto& p : pruned) test("pruned", p.first, p.second);
    return 0;
}

static int run_actl(int n) {
    int bad = 0;
    for (int t = 0; t < n; t++) {
        u8 g0[4], g1[4], Z[4] = {0, 0, 0, 0}; bool st = (t & 1);
        if (st) { auto s = stab_elements(); auto& e = s[rng() % s.size()]; memcpy(g0, e.data(), 4); } else { g0[0] = 1; g0[1] = 2; g0[2] = 4; g0[3] = 8; }
        for (int i = 0; i < 4; i++) g1[i] = (u8)(rng() & 15);
        Map g = cmap(g0, g1); u32 L = (u32)(rng() & 0xffffff);
        u32 C[4], X[4], C2[4]; F_of(L, C); for (int k = 0; k < 4; k++) X[k] = apply(g, C[k]);
        u32 L2 = act_L(g, L); F_of(L2, C2);
        if (L2 == 0xffffffffu || !same_space(X, 4, C2, 4)) bad++;
        // and a non-stabilising g0 must be rejected
        u8 h0[4]; rand_g0(h0); if (!stabilises_spread(h0)) { Map h = cmap(h0, Z); if (act_L(h, L) != 0xffffffffu) bad++; }
    }
    printf("actl: %d bad of %d\n", bad, n); return 0;
}

static int run_init() {
    auto reps = load_reps(); int ok = 0;
    for (size_t r = 0; r < reps.size(); r++) { u32 C[4]; F_of(reps[r], C); Tables* T = new Tables; bool b = T->init(C); int mx = 0; if (b) for (unsigned a = 1; a < 256; a++) mx = std::max(mx, T->D[a]); printf("  rep %zu %06x init %d max D %d\n", r, reps[r], b, mx); ok += b; delete T; }
    printf("init ok %d / %zu; |Stab(spread)| %zu\n", ok, reps.size(), stab_elements().size()); return 0;
}
static int review_main() {
    if (setup(false)) { fprintf(stderr, "setup failed\n"); exit(1); }
    make_spread();
    int argc = __argc; char** argv = __argv; int r = 2;
    if (argc >= 3 && !strcmp(argv[1], "gold")) r = run_gold(atoi(argv[2]));
    else if (argc >= 3 && !strcmp(argv[1], "known")) r = run_known(atoi(argv[2]));
    else if (argc >= 4 && !strcmp(argv[1], "brute")) r = run_brute(atoi(argv[2]), atoi(argv[3]));
    else if (argc >= 2 && !strcmp(argv[1], "init")) r = run_init();
    else if (argc >= 3 && !strcmp(argv[1], "actl")) r = run_actl(atoi(argv[2]));
    else fprintf(stderr, "usage: review_inv2 gold|known <n> | brute <rep> <npairs> | actl <n>\n");
    fflush(stdout); exit(r);
}
static int review_run = review_main();
