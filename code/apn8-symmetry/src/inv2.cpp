// inv2: quadratic APN spaces S8 on F_2^8 with an involution A in Aut(S8). README "p = 2" shows A has four Jordan
// blocks, so in these coordinates V = W + eps W (W = <e0..e3>, eps e_i = e_{4+i}, K = eps W) and A = 1 + eps.
// Also shown there: S8 is a free F_2[A]-module and F = S8^A = D(S8) (D = A - 1 on forms), 4-dim. For C in D(M):
// C vanishes on K x K, b_C(u, v) = C(eps u, v) is alternating on W, and c_C = C|W x W; C <-> (b_C, c_C) is a
// bijection D(M) <-> Alt(W)^2 (12 dims). b(F) must be the Klein image E of a line spread of W (all 56 spreads are
// one GL(4,2)-orbit). So F = {(b, L b) : b in E}, L in Hom(E, Alt W) (24 bits), up to C(A) restricted to Stab(E):
// shears (1, S) and (g0, 0) with g0 in Stab(spread), 360. Then S8 = F + <B1..B4>, D B_t = C_t (basis of F), each
// B_t in lambda_t + Fix, taken mod F: 12 bits each. B1..B3 enumerated with kernel-signature pruning, B4 solved.
// Modes: inv2 check <seedfile> <line>   structure checks + re-find a known function's S8 from its own F
//        inv2 forbits <repfile>         orbits of L; reps passing the partition test -> file
//        inv2 search <repfile> <out> [first last [min4, default 4]]
// Rowen (Claude Opus 5.5), 2026-10-24.
#define main apn8_main
#include "apn8.cpp"
#undef main
#include <set>
#include <mutex>
#include <map>
#include <array>
#include <algorithm>
#include <string>
#include <atomic>
#include <chrono>

static inline int bform(u64 m, unsigned x, unsigned y) { return __popcnt(rvec(m, x) & y) & 1; }
struct Set256 { u64 w[4] = {0, 0, 0, 0}; void add(unsigned x) { w[x >> 6] |= 1ull << (x & 63); } bool meets(const Set256& o) const { return (w[0] & o.w[0] & ~1ull) | (w[1] & o.w[1]) | (w[2] & o.w[2]) | (w[3] & o.w[3]); } };
static void radical_set(u64 m, Set256& S) { S = Set256(); S.add(0); for (unsigned x = 1; x < 256; x++) if (rvec(m, x) == 0) S.add(x); }
static bool partition_ok(const u32* fs, int n) {
    std::vector<std::pair<u64, Set256>> el; std::vector<int> rk; Ech28 e;
    for (int i = 0; i < n; i++) {
        if (!e.add(fs[i])) return false;
        u64 fm = coef_to_mat(fs[i]); size_t old = el.size();
        std::vector<u64> news; news.push_back(fm); for (size_t j = 0; j < old; j++) news.push_back(fm ^ el[j].first);
        for (u64 m : news) {
            int r = rank_mat(m); if (r == 0) return false;
            Set256 R; if (r < 8) { radical_set(m, R); for (size_t j = 0; j < el.size(); j++) if (rk[j] < 8 && R.meets(el[j].second)) return false; }
            el.push_back({m, R}); rk.push_back(r);
        }
    }
    return true;
}
struct Map { u8 f[256]; u32 col[28]; };
static void make_cols(Map& g) { for (int k = 0; k < 28; k++) { u64 m = coef_to_mat(1u << k); u32 r = 0; for (int q = 0; q < 28; q++) if (bform(m, g.f[1u << KI[q]], g.f[1u << KJ[q]])) r |= 1u << q; g.col[k] = r; } }
static inline u32 apply(const Map& g, u32 c) { u32 r = 0; while (c) { r ^= g.col[_tzcnt_u32(c)]; c &= c - 1; } return r; }

static Map Am;
static inline u32 act(u32 c) { return apply(Am, c); }
static inline u32 Dop(u32 c) { return c ^ act(c); }

// g = (g0, g1) in C(A) = GL(4, F_2[eps]): w + eps u -> g0 w + eps (g1 w + g0 u); g0, g1 as 4 nibble columns
static Map cmap(const u8* g0, const u8* g1) {
    Map g;
    for (unsigned x = 0; x < 256; x++) {
        unsigned w = x & 15, u = x >> 4, a = 0, b = 0;
        for (int i = 0; i < 4; i++) { if (w >> i & 1) { a ^= g0[i]; b ^= g1[i]; } if (u >> i & 1) b ^= g0[i]; }
        g.f[x] = (u8)(a | b << 4);
    }
    make_cols(g); return g;
}

// alternating forms on W: 6 bits over pairs
static const int P4I[6] = {0, 0, 0, 1, 1, 2}, P4J[6] = {1, 2, 3, 2, 3, 3};
static inline int b4(u8 b, unsigned x, unsigned y) { int s = 0; for (int p = 0; p < 6; p++) if (b >> p & 1) s ^= ((x >> P4I[p] & 1) & (y >> P4J[p] & 1)) ^ ((x >> P4J[p] & 1) & (y >> P4I[p] & 1)); return s; }
static unsigned rad4(u8 b) { unsigned r = 0; for (unsigned x = 1; x < 16; x++) { bool z = true; for (int i = 0; i < 4; i++) if (b4(b, x, 1u << i)) z = false; if (z) r |= 1u << x; } return r; } // bitmask of nonzero radical vectors
static u8 bpart(u32 C) { u64 m = coef_to_mat(C); u8 b = 0; for (int p = 0; p < 6; p++) if (bform(m, 1u << (4 + P4I[p]), 1u << P4J[p])) b |= 1 << p; return b; }
static u8 cpart(u32 C) { u64 m = coef_to_mat(C); u8 c = 0; for (int p = 0; p < 6; p++) if (bform(m, 1u << P4I[p], 1u << P4J[p])) c |= 1 << p; return c; }

// linear algebra over 28-bit forms with tracking
static std::vector<u32> fixb;            // Fix basis (16)
static u32 dmform[4096];                 // (b | c << 6) -> form in D(M)
static u32 dpiv[28], dsrc[28];           // echelon of D(M) with preimages: dpiv[b] = D(dsrc[b])
static u32 preimage(u32 C) { u32 x = 0; for (int b = 27; b >= 0; b--) if (C >> b & 1) { if (!dpiv[b]) return 0xffffffffu; C ^= dpiv[b]; x ^= dsrc[b]; } return x; }

static int setup(bool verbose) {
    init_index();
    for (unsigned x = 0; x < 256; x++) Am.f[x] = (u8)(x ^ ((x & 15) << 4));
    make_cols(Am);
    memset(dpiv, 0, sizeof dpiv); memset(dsrc, 0, sizeof dsrc);
    Ech28 fx; fixb.clear();
    // D on the 28 basis forms: echelon of images with sources; kernel from tracking
    std::vector<std::pair<u32, u32>> rows; for (int k = 0; k < 28; k++) rows.push_back({Dop(1u << k), 1u << k});
    for (auto& rw : rows) {
        u32 v = rw.first, s = rw.second;
        for (int b = 27; b >= 0 && v; b--) if ((v >> b & 1) && dpiv[b]) { v ^= dpiv[b]; s ^= dsrc[b]; }
        if (v) { int b = 31 - (int)__lzcnt(v); dpiv[b] = v; dsrc[b] = s; }
        else if (fx.add(s)) fixb.push_back(s);
    }
    int dimDM = 0; for (int b = 0; b < 28; b++) if (dpiv[b]) dimDM++;
    // D(M) <-> (b, c): every form in D(M) vanishes on K x K, b alternating; key map bijective
    std::vector<u32> dmb; for (int b = 0; b < 28; b++) if (dpiv[b]) dmb.push_back(dpiv[b]);
    bool kk = true, alt = true;
    for (u32 C : dmb) {
        u64 m = coef_to_mat(C);
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { if (bform(m, 1u << (4 + i), 1u << (4 + j))) kk = false; if (bform(m, 1u << (4 + i), 1u << j) != bform(m, 1u << (4 + j), 1u << i)) alt = false; }
        for (int i = 0; i < 4; i++) if (bform(m, 1u << (4 + i), 1u << i)) alt = false;
    }
    // key echelon
    u32 kp[12] = {0}, kf[12] = {0};
    for (u32 C : dmb) { u32 key = bpart(C) | (u32)cpart(C) << 6, f = C; for (int b = 11; b >= 0 && key; b--) if ((key >> b & 1) && kp[b]) { key ^= kp[b]; f ^= kf[b]; } if (key) { int b = 31 - (int)__lzcnt(key); kp[b] = key; kf[b] = f; } }
    int krank = 0; for (int b = 0; b < 12; b++) if (kp[b]) krank++;
    for (u32 key = 0; key < 4096; key++) { u32 k = key, f = 0; for (int b = 11; b >= 0 && k; b--) if (k >> b & 1) { k ^= kp[b]; f ^= kf[b]; } dmform[key] = f; }
    if (verbose) printf("setup: dim Fix %zu, dim D(M) %d, D(M) vanishes on KxK %d, b alternating %d, (b,c) rank %d\n", fixb.size(), dimDM, kk, alt, krank);
    return (fixb.size() == 16 && dimDM == 12 && kk && alt && krank == 12) ? 0 : 1;
}

// ---------- the standard spread and E ----------
static unsigned lines[5]; static u8 bL[5]; static u8 Eb[4]; static u8 Ebas_inv[64]; // Ebas_inv: element of E -> coords (4 bits), 0xff outside
static u8 mulw(u8 z) { return (u8)(((z >> 1) & 1) | ((((z & 1) ^ (z >> 1)) & 1) << 1)); } // (x0 + x1 w) w
static void make_spread() {
    // W = F_4^2: z1 = bits 0,1; z2 = bits 2,3
    u8 vs[5][2] = {{1, 0}, {0, 1}, {1, 1}, {1, 2}, {1, 3}};
    for (int l = 0; l < 5; l++) {
        unsigned m = 0; u8 a = vs[l][0], b = vs[l][1];
        for (int k = 0; k < 3; k++) { m |= 1u << (a | b << 2); a = mulw(a); b = mulw(b); }
        lines[l] = m;
        for (u8 bb = 1; bb < 64; bb++) if (rad4(bb) == m) bL[l] = bb;
    }
    // basis of E from the bL, echelon on 6 bits
    u8 pv[6] = {0}; int n = 0;
    for (int l = 0; l < 5; l++) { u8 x = bL[l]; for (int b = 5; b >= 0 && x; b--) if (x >> b & 1) { if (pv[b]) x ^= pv[b]; else { pv[b] = x; Eb[n++] = bL[l]; x = 0; } } }
    memset(Ebas_inv, 0xff, sizeof Ebas_inv);
    for (unsigned c = 0; c < 16; c++) { u8 v = 0; for (int i = 0; i < 4; i++) if (c >> i & 1) v ^= Eb[i]; Ebas_inv[v] = (u8)c; }
}
static bool stabilises_spread(const u8* g0) {
    for (int l = 0; l < 5; l++) { unsigned im = 0; for (unsigned x = 1; x < 16; x++) if (lines[l] >> x & 1) { unsigned y = 0; for (int i = 0; i < 4; i++) if (x >> i & 1) y ^= g0[i]; im |= 1u << y; } bool ok = false; for (int q = 0; q < 5; q++) if (lines[q] == im) ok = true; if (!ok) return false; }
    return true;
}

// F(L) basis: C_k = dmform(Eb[k] | L_k << 6)
static void F_of(u32 L, u32* C) { for (int k = 0; k < 4; k++) C[k] = dmform[Eb[k] | ((L >> (6 * k)) & 63) << 6]; }
// g . F(L) as L' (or 0xffffffff if g doesn't preserve E)
static u32 act_L(const Map& g, u32 L) {
    u32 C[4]; F_of(L, C);
    u32 X[4]; u8 coord[4];
    for (int k = 0; k < 4; k++) { X[k] = apply(g, C[k]); u8 bb = bpart(X[k]); if (Ebas_inv[bb] == 0xff) return 0xffffffffu; coord[k] = Ebas_inv[bb]; }
    // invert the 4x4 coordinate matrix: find, for each m, combination of X with b-part Eb[m]
    u8 rowc[4]; u8 comb[4]; for (int k = 0; k < 4; k++) { rowc[k] = coord[k]; comb[k] = (u8)(1 << k); }
    for (int b = 0; b < 4; b++) { int p = -1; for (int i = b; i < 4; i++) if (rowc[i] >> b & 1) { p = i; break; } if (p < 0) return 0xffffffffu; std::swap(rowc[p], rowc[b]); std::swap(comb[p], comb[b]); for (int i = 0; i < 4; i++) if (i != b && (rowc[i] >> b & 1)) { rowc[i] ^= rowc[b]; comb[i] ^= comb[b]; } }
    u32 Lp = 0;
    for (int m = 0; m < 4; m++) { u32 Y = 0; for (int k = 0; k < 4; k++) if (comb[m] >> k & 1) Y ^= X[k]; if (bpart(Y) != Eb[m]) return 0xffffffffu; Lp |= (u32)cpart(Y) << (6 * m); }
    return Lp;
}
struct Aff { u32 col[24]; u32 t; };
static bool make_aff(const Map& g, Aff& A) {
    A.t = act_L(g, 0); if (A.t == 0xffffffffu) return false;
    for (int i = 0; i < 24; i++) { u32 v = act_L(g, 1u << i); if (v == 0xffffffffu) return false; A.col[i] = v ^ A.t; }
    // affinity check on random points
    u64 s = 0x1234567ull;
    for (int r = 0; r < 50; r++) { s = s * 6364136223846793005ull + 1442695040888963407ull; u32 L = (u32)(s >> 40) & 0xffffff; u32 v = A.t; for (int i = 0; i < 24; i++) if (L >> i & 1) v ^= A.col[i]; if (v != act_L(g, L)) return false; }
    return true;
}
static inline u32 aff_apply(const Aff& A, u32 L) { u32 v = A.t; while (L) { v ^= A.col[_tzcnt_u32(L)]; L &= L - 1; } return v; }

static std::vector<std::array<u8, 4>> stab_elements() {
    std::vector<std::array<u8, 4>> res;
    for (unsigned a = 1; a < 16; a++) for (unsigned b = 1; b < 16; b++) for (unsigned c = 1; c < 16; c++) for (unsigned d = 1; d < 16; d++) {
        u8 g[4] = {(u8)a, (u8)b, (u8)c, (u8)d}; if (rank8(g, 4) != 4) continue;
        if (stabilises_spread(g)) res.push_back({g[0], g[1], g[2], g[3]});
    }
    return res;
}
static std::vector<Aff> group_gens(bool verbose) {
    std::vector<Aff> gens;
    u8 I4[4] = {1, 2, 4, 8}, Z[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { u8 g1[4] = {0, 0, 0, 0}; g1[i] = (u8)(1 << j); Map g = cmap(I4, g1); Aff A; if (!make_aff(g, A)) { fprintf(stderr, "shear not affine on L\n"); exit(1); } gens.push_back(A); }
    auto st = stab_elements();
    if (verbose) printf("Stab(spread) in GL(4,2): %zu\n", st.size());
    for (auto& s : st) { Map g = cmap(s.data(), Z); Aff A; if (!make_aff(g, A)) { fprintf(stderr, "stab element not affine on L\n"); exit(1); } gens.push_back(A); }
    return gens;
}

// ---------- the lift engine ----------
struct Tables {
    u32 S[4], lam[4], Y[12];
    int D[256]; u8 ybas[256][4];
    u8 sigY[12][256], sigL[4][256];
    std::vector<std::array<u64, 16>> plY; std::array<u64, 16> plL[4];
    u64 maskD[5][4];
    bool init(const u32* C) {
        for (int k = 0; k < 4; k++) { S[k] = C[k]; lam[k] = preimage(C[k]); if (lam[k] == 0xffffffffu) return false; if (Dop(lam[k]) != C[k]) return false; }
        { Ech28 e; for (int k = 0; k < 4; k++) if (!e.add(C[k])) return false; int n = 0; for (u32 f : fixb) if (e.add(f)) Y[n++] = f; if (n != 12) return false; }
        memset(maskD, 0, sizeof maskD);
        u64 sm[4]; for (int k = 0; k < 4; k++) sm[k] = coef_to_mat(S[k]);
        for (unsigned a = 1; a < 256; a++) {
            u8 rr[4]; for (int k = 0; k < 4; k++) rr[k] = rvec(sm[k], a);
            u8 piv[8] = {0}; std::vector<u8> bas;
            auto addv = [&](u8 x) { u8 y = x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; bas.push_back(x); return true; } } return false; };
            addv((u8)a);
            for (unsigned y = 1; y < 256; y++) { bool z = true; for (int i = 0; i < 4; i++) if (__popcnt(rr[i] & y) & 1) { z = false; break; } if (z) addv((u8)y); }
            D[a] = (int)bas.size() - 1;
            if (D[a] > 4) return false;
            for (int k = 0; k < D[a]; k++) ybas[a][k] = bas[k + 1];
            maskD[D[a]][a >> 6] |= 1ull << (a & 63);
            for (int i = 0; i < 12; i++) { u64 cm = coef_to_mat(Y[i]); u8 sg = 0; for (int k = 0; k < D[a]; k++) if (bform(cm, a, ybas[a][k])) sg |= 1 << k; sigY[i][a] = sg; }
            for (int i = 0; i < 4; i++) { u64 cm = coef_to_mat(lam[i]); u8 sg = 0; for (int k = 0; k < D[a]; k++) if (bform(cm, a, ybas[a][k])) sg |= 1 << k; sigL[i][a] = sg; }
        }
        auto plane = [&](const u8* sg, std::array<u64, 16>& P) { P.fill(0); for (unsigned a = 1; a < 256; a++) for (int k = 0; k < 4; k++) if (sg[a] >> k & 1) P[k * 4 + (a >> 6)] |= 1ull << (a & 63); };
        plY.assign(4096, {});
        for (u32 g = 1; g < 4096; g++) { int i = _tzcnt_u32(g); std::array<u64, 16> Pi; plane(sigY[i], Pi); auto& P = plY[g]; P = plY[g & (g - 1)]; for (int k = 0; k < 16; k++) P[k] ^= Pi[k]; }
        for (int i = 0; i < 4; i++) plane(sigL[i], plL[i]);
        return true;
    }
    u32 lift(int t, u32 g) const { u32 v = lam[t]; for (int i = 0; i < 12; i++) if (g >> i & 1) v ^= Y[i]; return v; }
};

struct IStats { long long n1 = 0, n2 = 0, leaves = 0, passed = 0, sols = 0, bugs = 0; std::map<std::string, long long> types; };

template <class Found>
struct Dfs {
    const Tables& T; Found& found; IStats st;
    std::array<u64, 16> combo[16]; u32 gs[4];
    Dfs(const Tables& t, Found& f) : T(t), found(f) {}
    bool level_ok(int j) const {
        const int r = 4; int nc = (1 << j) - 1;
        for (int w = 0; w < 4; w++) {
            u64 nz[16]; for (int c = 1; c <= nc; c++) nz[c] = combo[c][w] | combo[c][4 + w] | combo[c][8 + w] | combo[c][12 + w];
            u64 all = ~0ull, any = 0, pairzero = 0;
            for (int c = 1; c <= nc; c++) all &= nz[c];
            for (int i = 0; i < j; i++) any |= nz[1 << i];
            for (int c = 1; c <= nc; c++) for (int c2 = c + 1; c2 <= nc; c2++) pairzero |= ~nz[c] & ~nz[c2];
            for (int d = 1; d <= 4; d++) {
                u64 mk = T.maskD[d][w]; if (!mk) continue;
                int R = d - (r - j); if (R <= 0) continue;
                u64 ok; if (R == j) ok = all; else if (R == j - 1) ok = ~pairzero; else ok = any;
                if (mk & ~ok) return false;
            }
        }
        return true;
    }
    void place(int j, u32 g) {
        gs[j] = g; int base = 1 << j;
        for (int k = 0; k < 16; k++) combo[base][k] = T.plL[j][k] ^ T.plY[g][k];
        for (int c = 1; c < base; c++) for (int k = 0; k < 16; k++) combo[base + c][k] = combo[c][k] ^ combo[base][k];
    }
    u8 sigat(int i, unsigned a) const { u8 s = 0; for (int k = 0; k < 4; k++) if (combo[1 << i][k * 4 + (a >> 6)] >> (a & 63) & 1) s |= 1 << k; return s; }
    void rec(int j) {
        if (j == 3) { st.leaves++; if (level_ok(3)) { st.passed++; solve(); } return; }
        if (j > 0 && !level_ok(j)) return;
        if (j == 1) st.n1++; if (j == 2) st.n2++;
        for (u32 g = 0; g < 4096; g++) { place(j, g); rec(j + 1); }
    }
    void run_from(u32 g0) { place(0, g0); rec(1); }
    void solve() {
        u32 rows[300]; int nr = 0;
        for (unsigned a = 1; a < 256; a++) {
            int d = T.D[a]; if (!d) continue;
            u8 v[3]; for (int i = 0; i < 3; i++) v[i] = sigat(i, a);
            u8 pv[4] = {0, 0, 0, 0}; int rk = 0;
            for (int i = 0; i < 3; i++) { u8 x = v[i]; for (int b = 3; b >= 0 && x; b--) if (x >> b & 1) { if (pv[b]) x ^= pv[b]; else { pv[b] = x; rk++; x = 0; } } }
            if (rk == d) continue;
            if (rk < d - 1) return;
            u8 z = 0; for (u8 c = 1; c < (1 << d); c++) { bool ok = true; for (int i = 0; i < 3; i++) if (__popcnt((unsigned)(v[i] & c)) & 1) { ok = false; break; } if (ok) { z = c; break; } }
            u32 row = 0; for (int i = 0; i < 12; i++) if (__popcnt((unsigned)(T.sigY[i][a] & z)) & 1) row |= 1u << i;
            int rhs = 1 ^ (__popcnt((unsigned)(T.sigL[3][a] & z)) & 1);
            rows[nr++] = row | (u32)rhs << 12;
        }
        int piv[12]; for (int i = 0; i < 12; i++) piv[i] = -1; int rk = 0;
        for (int b = 11; b >= 0; b--) {
            int p = -1; for (int i = rk; i < nr; i++) if (rows[i] >> b & 1) { p = i; break; }
            if (p < 0) continue;
            std::swap(rows[p], rows[rk]);
            for (int i = 0; i < nr; i++) if (i != rk && (rows[i] >> b & 1)) rows[i] ^= rows[rk];
            piv[b] = rk; rk++;
        }
        for (int i = rk; i < nr; i++) if (rows[i]) return;
        std::vector<int> fr; for (int b = 0; b < 12; b++) if (piv[b] < 0) fr.push_back(b);
        for (u32 f = 0; f < (1u << fr.size()); f++) {
            u32 x = 0; for (size_t q = 0; q < fr.size(); q++) if (f >> q & 1) x |= 1u << fr[q];
            for (int b = 0; b < 12; b++) if (piv[b] >= 0) { u32 rw = rows[piv[b]]; int val = (rw >> 12) & 1; val ^= __popcnt(rw & x & 0xfff & ~(1u << b)) & 1; if (val) x |= 1u << b; }
            u32 b8[8]; for (int k = 0; k < 4; k++) b8[k] = T.S[k]; for (int t = 0; t < 3; t++) b8[4 + t] = T.lift(t, gs[t]); b8[7] = T.lift(3, x);
            st.sols++;
            u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(b8[i]);
            if (max_kernel(m, 8) != 1) { st.bugs++; continue; }
            int cnt[5]; rank_dist(m, 8, cnt);
            char t[64]; snprintf(t, sizeof t, "ranks 8:%d 6:%d 4:%d 2:%d", cnt[4], cnt[3], cnt[2], cnt[1]); st.types[t]++;
            found(b8, cnt);
        }
    }
};

static void merge(IStats& a, const IStats& b) { a.n1 += b.n1; a.n2 += b.n2; a.leaves += b.leaves; a.passed += b.passed; a.sols += b.sols; a.bugs += b.bugs; for (auto& kv : b.types) a.types[kv.first] += kv.second; }

// run all lifts for one F, in parallel over B1; found(b8, cnt) must be thread-safe
template <class Found>
static IStats run_F(const u32* C, Found& found, bool& ok) {
    Tables T; ok = T.init(C); IStats tot; if (!ok) return tot;
    std::mutex mu;
    #pragma omp parallel
    {
        Dfs<Found> d(T, found);
        #pragma omp for schedule(dynamic, 8)
        for (int g0 = 0; g0 < 4096; g0++) d.run_from((u32)g0);
        std::lock_guard<std::mutex> lk(mu); merge(tot, d.st);
    }
    return tot;
}
static void print_stats(const IStats& s) {
    printf("  B1 kept %lld, (B1,B2) kept %lld, leaves %lld, passed %lld, solutions %lld, bugs %lld\n", s.n1, s.n2, s.leaves, s.passed, s.sols, s.bugs);
    for (auto& kv : s.types) printf("    %s: %lld\n", kv.first.c_str(), kv.second);
}

// ---------- check mode: a known function with an involution ----------
static int check(const char* path, int line) {
    std::vector<Found> fs; if (!load_seeds(path, fs) || line >= (int)fs.size()) { fprintf(stderr, "can't read %s:%d\n", path, line); return 1; }
    Found& Fn = fs[line];
    // its involutions (as aut8, rank-4 radical triples), then into standard coordinates
    u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(Fn.basis[i]);
    Ech28 S; for (int i = 0; i < 8; i++) S.add(Fn.basis[i]);
    // brute force over involutions is too big; use the aut8 route
    struct Rad { std::vector<u8> el; u8 mem[256]; };
    auto radical_of = [](u64 mm) { Rad r; memset(r.mem, 0, 256); r.mem[0] = 1; for (unsigned x = 1; x < 256; x++) if (rvec(mm, x) == 0) { r.el.push_back((u8)x); r.mem[x] = 1; } return r; };
    auto graph_map = [](const Rad& R1, const Rad& R2, const Rad& R3, u8* Tm) { u8 p1[256], p2[256]; memset(p1, 0, 256); memset(p2, 0, 256); std::vector<u8> a1(R1.el), a2(R2.el); a1.push_back(0); a2.push_back(0); for (u8 a : a1) for (u8 b : a2) { p1[a ^ b] = a; p2[a ^ b] = b; } memset(Tm, 0, 256); for (u8 r : R3.el) { if (!p1[r] || !p2[r]) return false; Tm[p1[r]] = p2[r]; } return true; };
    std::vector<Rad> R; { u64 cur = 0; for (unsigned g = 1; g < 256; g++) { cur ^= m[_tzcnt_u32(g)]; if (rank_mat(cur) == 4) R.push_back(radical_of(cur)); } }
    int k = (int)R.size(); if (k < 3) { printf("fewer than 3 rank-4 components\n"); return 1; }
    u8 e[4]; { u8 piv[8] = {0}; int n = 0; for (u8 x : R[0].el) { u8 y = x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; y = 0; e[n++] = x; } } if (n == 4) break; } }
    u8 Tm[256]; graph_map(R[0], R[1], R[2], Tm);
    std::vector<std::array<u8, 256>> invs;
    for (int s1 = 0; s1 < k; s1++) for (int s2 = 0; s2 < k; s2++) for (int s3 = 0; s3 < k; s3++) {
        if (s1 == s2 || s1 == s3 || s2 == s3) continue;
        u8 Tp[256]; if (!graph_map(R[s1], R[s2], R[s3], Tp)) continue;
        for (u8 u0 : R[s1].el) for (u8 u1 : R[s1].el) for (u8 u2 : R[s1].el) for (u8 u3 : R[s1].el) {
            u8 u[4] = {u0, u1, u2, u3}; if (rank8(u, 4) != 4) continue;
            u8 b[8], im[8]; for (int i = 0; i < 4; i++) { b[i] = e[i]; im[i] = u[i]; b[4 + i] = Tm[e[i]]; im[4 + i] = Tp[u[i]]; }
            std::array<u8, 256> Af; { u8 x = 0, y = 0; Af[0] = 0; for (unsigned g = 1; g < 256; g++) { int t = _tzcnt_u32(g); x ^= b[t]; y ^= im[t]; Af[x] = y; } }
            bool inv = true; for (unsigned x = 1; x < 256; x++) if (Af[Af[x]] != x) inv = false;
            bool id = true; for (unsigned x = 1; x < 256; x++) if (Af[x] != x) id = false;
            if (!inv || id) continue;
            bool ok = true; for (int t = 0; t < 8 && ok; t++) { u32 c = 0; for (int q = 0; q < 28; q++) if (bform(m[t], Af[1u << KI[q]], Af[1u << KJ[q]])) c |= 1u << q; if (S.reduce(c)) ok = false; }
            if (ok) invs.push_back(Af);
        }
    }
    printf("involutions in Aut(S8): %zu\n", invs.size());
    int refound = 0;
    for (auto& Ap : invs) {
        // N = Ap + 1; generators w_i: a basis of a complement of ker N; M: e_i -> w_i, e_{4+i} -> N w_i
        int fixdim = 0; for (unsigned x = 1; x < 256; x++) if (Ap[x] == x) fixdim++;
        u8 w[4]; int n = 0; { u8 piv[8] = {0}; std::vector<u8> kb; for (unsigned x = 1; x < 256; x++) if (Ap[x] == x) { u8 y = (u8)x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; y = 0; } } }
            for (unsigned x = 1; x < 256 && n < 4; x++) { u8 y = (u8)x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; w[n++] = (u8)x; y = 0; } } } }
        printf("  involution: fixed vectors %d (2^4 - 1 = 15 expected), complement basis %d\n", fixdim, n);
        if (fixdim != 15 || n != 4) continue;
        Map M; for (unsigned x = 0; x < 256; x++) { u8 y = 0; for (int i = 0; i < 4; i++) { if (x >> i & 1) y ^= w[i]; if (x >> (4 + i) & 1) y ^= (u8)(w[i] ^ Ap[w[i]]); } M.f[x] = y; }
        make_cols(M);
        u32 b2[8]; for (int i = 0; i < 8; i++) b2[i] = apply(M, Fn.basis[i]);
        Ech28 S2; for (int i = 0; i < 8; i++) S2.add(b2[i]);
        bool inv = true; for (int i = 0; i < 8; i++) if (S2.reduce(act(b2[i]))) inv = false;
        // F = D(S8), S8^A
        Ech28 Fe, Xe; u32 C[4]; int nc = 0, nfix = 0;
        { u32 c = 0; for (unsigned g = 1; g < 256; g++) { c ^= b2[_tzcnt_u32(g)]; if (act(c) == c) { nfix++; } u32 d = Dop(c); if (d && Fe.add(d)) { if (nc < 4) C[nc] = d; nc++; } } }
        bool fisfix = true; for (int i = 0; i < nc && i < 4; i++) if (act(C[i]) != C[i]) fisfix = false;
        // b-parts: radicals form a spread of W?
        std::vector<unsigned> rads; int r2 = 0, r0 = 0, other = 0;
        { u32 c = 0; for (unsigned g = 1; g < 16; g++) { c ^= C[_tzcnt_u32(g)]; u8 bb = bpart(c); unsigned rd = rad4(bb); int sz = __popcnt(rd); if (sz == 3) { r2++; rads.push_back(rd); } else if (sz == 0) r0++; else other++; } }
        unsigned cover = 0; bool disj = true; for (unsigned rd : rads) { if (cover & rd) disj = false; cover |= rd; }
        printf("  invariant %d, dim S8^A %d, dim D(S8) %d, D(S8) fixed %d, b-forms: rank-2 %d (radicals disjoint %d, cover W %d), nondegenerate %d, other %d\n",
               inv, 31 - (int)__lzcnt((unsigned)nfix + 1), nc, fisfix, r2, disj, cover == 0xfffe, r0, other);
        if (!inv || nc != 4) continue;
        u32 target[8]; memcpy(target, b2, sizeof target); canon8(target);
        bool hit = false; std::mutex mu; long long cnt4 = 0;
        auto found = [&](u32* b8, int* cnt) { u32 c[8]; memcpy(c, b8, sizeof c); canon8(c); std::lock_guard<std::mutex> g(mu); if (!memcmp(c, target, sizeof c)) hit = true; if (cnt[2] >= 4) cnt4++; };
        bool ok; auto t0 = std::chrono::steady_clock::now();
        IStats s = run_F(C, found, ok);
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        printf("  engine on its own F (%.1f s): init %d, re-found %d, APN spaces with >= 4 rank-4: %lld\n", el, ok, hit, cnt4);
        print_stats(s);
        refound += hit;
        break; // one involution is enough
    }
    printf("check: %s\n", refound ? "re-found" : "NOT re-found");
    return refound ? 0 : 1;
}

// ---------- orbits of F ----------
static int forbits(const char* outp) {
    auto gens = group_gens(true);
    std::vector<u64> seen((1u << 24) / 64, 0);
    std::vector<u32> q; long orbits = 0, kept = 0; std::map<long, long> hist;
    FILE* f = fopen(outp, "w");
    auto t0 = std::chrono::steady_clock::now();
    for (u32 L = 0; L < (1u << 24); L++) {
        if (seen[L >> 6] >> (L & 63) & 1) continue;
        orbits++; q.clear(); q.push_back(L); seen[L >> 6] |= 1ull << (L & 63);
        for (size_t h = 0; h < q.size(); h++) for (auto& A : gens) { u32 y = aff_apply(A, q[h]); if (!(seen[y >> 6] >> (y & 63) & 1)) { seen[y >> 6] |= 1ull << (y & 63); q.push_back(y); } }
        hist[(long)q.size()]++;
        u32 C[4]; F_of(L, C);
        u64 mm[4]; for (int i = 0; i < 4; i++) mm[i] = coef_to_mat(C[i]);
        bool part = partition_ok(C, 4), ker = max_kernel(mm, 4) <= 5;
        if (part && ker) { kept++; fprintf(f, "%06x %zu\n", L, q.size()); }
    }
    fclose(f);
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("orbits of L: %ld (%.1f s), passing partition + kernel tests: %ld\n  sizes:", orbits, el, kept);
    for (auto& kv : hist) printf(" %ld^%ld", kv.first, kv.second);
    printf("\n");
    return 0;
}

static int search(const char* repp, const char* outp, long first, long last, int min4) {
    std::vector<std::pair<u32, long>> reps; { FILE* f = fopen(repp, "r"); if (!f) return 1; unsigned L; long sz; while (fscanf(f, "%x %ld", &L, &sz) == 2) reps.push_back({L, sz}); fclose(f); }
    if (last > (long)reps.size()) last = (long)reps.size();
    FILE* out = fopen(outp, "a"); std::mutex mu; IStats all; long long found5 = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (long r = first; r < last; r++) {
        u32 C[4]; F_of(reps[r].first, C);
        auto found = [&](u32* b8, int* cnt) {
            if (cnt[2] < min4) return;
            std::lock_guard<std::mutex> g(mu);
            fprintf(out, "ranks 8:%d 6:%d 4:%d 2:%d inv2 F %06x basis", cnt[4], cnt[3], cnt[2], cnt[1], reps[r].first); for (int i = 0; i < 8; i++) fprintf(out, " %07x", b8[i]); fprintf(out, "\n"); fflush(out);
            if (cnt[2] >= 5) { found5++; printf("FOUND ranks 8:%d 6:%d 4:%d\n", cnt[4], cnt[3], cnt[2]); fflush(stdout); }
        };
        bool ok; IStats s = run_F(C, found, ok);
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        printf("F %ld (%06x, orbit %ld) [%.0f s]: init %d", r, reps[r].first, reps[r].second, el, ok); if (ok) { printf("\n"); print_stats(s); } else printf("\n");
        fflush(stdout);
        merge(all, s);
    }
    printf("inv2 search reps %ld..%ld: total\n", first, last); print_stats(all);
    printf("with >= 5 rank-4: %lld\n", found5);
    fclose(out);
    return 0;
}

int main(int argc, char** argv) {
    if (setup(true)) { fprintf(stderr, "setup checks failed\n"); return 1; }
    make_spread();
    printf("spread lines:"); for (int l = 0; l < 5; l++) printf(" %04x", lines[l]); printf("; E basis %02x %02x %02x %02x\n", Eb[0], Eb[1], Eb[2], Eb[3]);
    { int deg = 0; for (unsigned c = 1; c < 16; c++) { u8 v = 0; for (int i = 0; i < 4; i++) if (c >> i & 1) v ^= Eb[i]; if (rad4(v)) deg++; } printf("E: degenerate nonzero forms %d (want 5)\n", deg); }
    fflush(stdout);
    if (argc >= 4 && !strcmp(argv[1], "check")) return check(argv[2], atoi(argv[3]));
    if (argc >= 3 && !strcmp(argv[1], "forbits")) return forbits(argv[2]);
    if (argc >= 4 && !strcmp(argv[1], "search")) return search(argv[2], argv[3], argc > 5 ? atol(argv[4]) : 0, argc > 5 ? atol(argv[5]) : 0x7fffffffL, argc > 6 ? atoi(argv[6]) : 4);
    fprintf(stderr, "usage: inv2 check <seedfile> <line> | forbits <repfile> | search <repfile> <out> [first last]\n");
    return 2;
}
