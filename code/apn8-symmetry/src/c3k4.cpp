// c3k4: quadratic APN spaces S8 on F_2^8 invariant under a fixed-point-free order-3 map, with >= 5 rank-4
// components. V = F_4^4, x = u + w v (u = low nibble, v = high nibble), A = multiplication by w: (u, v) -> (v, u + v).
// Alt(V) = Fix (16, Hermitian forms) + N (12, F_4-alternating forms); S8 = Q + P with Q = S8 & Fix (F_2-dim a) and
// P = S8 & N (F_4-dim b). Shown in README: (a, b) in {(8,0), (6,1), (4,2)}.
// This file: case (6,1), pattern 3+1+1: a rank-4 3-cycle B1 (normalised to radical u-space and form c on v, which
// fixes P) plus two fixed rank-4 components h4, h5, pairs up to G1 = stabiliser (Sp(4,2) x Frobenius); then
// Q = <h1, h4, h5> + 3 more Hermitian forms (two enumerated, one solved linearly). Anything APN found has >= 5 rank-4.
// Modes:  c3k4 search <out> [first last]     the search (seed range for splitting)
//         c3k4 validate <seedfile>           transform known C3 classes into these coordinates and re-find them
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
#include <unordered_map>
#include <functional>

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

// a linear (or semilinear) map of V as a LUT; form transform B -> B o g as a 28x28 column table
struct Map { u8 f[256]; u32 col[28]; };
static void make_cols(Map& g) { for (int k = 0; k < 28; k++) { u64 m = coef_to_mat(1u << k); u32 r = 0; for (int q = 0; q < 28; q++) if (bform(m, g.f[1u << KI[q]], g.f[1u << KJ[q]])) r |= 1u << q; g.col[k] = r; } }
static inline u32 apply(const Map& g, u32 c) { u32 r = 0; while (c) { r ^= g.col[_tzcnt_u32(c)]; c &= c - 1; } return r; }
static Map Am;
static inline u32 act(u32 c) { return apply(Am, c); }
static inline u32 fixpart(u32 c) { u32 a = act(c); return c ^ a ^ act(a); }

static std::vector<u32> fixb; // basis of Fix (16)
static u32 B1, h1, n1, n1a; // standard 3-cycle component and its parts

static void setup() {
    for (unsigned x = 0; x < 256; x++) { unsigned u = x & 15, v = x >> 4; Am.f[x] = (u8)(v | ((u ^ v) << 4)); }
    make_cols(Am);
    for (unsigned x = 1; x < 256; x++) { unsigned y = Am.f[Am.f[Am.f[x]]]; if (y != x || Am.f[x] == x) { fprintf(stderr, "A\n"); exit(1); } }
    u32 rows[28], comb[28]; for (int k = 0; k < 28; k++) { rows[k] = act(1u << k) ^ (1u << k); comb[k] = 1u << k; }
    int rk = 0;
    for (int bit = 27; bit >= 0; bit--) {
        int p = -1; for (int i = rk; i < 28; i++) if (rows[i] >> bit & 1) { p = i; break; }
        if (p < 0) continue;
        std::swap(rows[p], rows[rk]); std::swap(comb[p], comb[rk]);
        for (int i = 0; i < 28; i++) if (i != rk && (rows[i] >> bit & 1)) { rows[i] ^= rows[rk]; comb[i] ^= comb[rk]; }
        rk++;
    }
    for (int i = rk; i < 28; i++) fixb.push_back(comb[i]);
    if (fixb.size() != 16) { fprintf(stderr, "dim Fix %zu\n", fixb.size()); exit(1); }
    B1 = (1u << K[4][5]) | (1u << K[6][7]);
    h1 = fixpart(B1); n1 = B1 ^ h1; n1a = act(n1);
    if (rank_mat(coef_to_mat(B1)) != 4) { fprintf(stderr, "B1\n"); exit(1); }
    if (!n1 || n1a == n1 || rank_mat(coef_to_mat(n1)) != 8 || rank_mat(coef_to_mat(n1a)) != 8 || rank_mat(coef_to_mat(n1 ^ n1a)) != 8) { fprintf(stderr, "P degenerate\n"); exit(1); }
}

// G1: elements of GL(4,2) acting on u and v alike, possibly times Frobenius (u, v) -> (u + v, v), that preserve
// Fix, P = <n1, n1a> and B1.
static std::vector<Map> G1;
static void make_G1() {
    std::vector<unsigned> cols4;
    for (unsigned c0 = 1; c0 < 16; c0++) for (unsigned c1 = 1; c1 < 16; c1++) for (unsigned c2 = 1; c2 < 16; c2++) for (unsigned c3 = 1; c3 < 16; c3++) {
        u8 vv[4] = {(u8)c0, (u8)c1, (u8)c2, (u8)c3}; if (rank8(vv, 4) != 4) continue;
        unsigned M[16]; for (unsigned x = 0; x < 16; x++) { unsigned y = 0; if (x & 1) y ^= c0; if (x & 2) y ^= c1; if (x & 4) y ^= c2; if (x & 8) y ^= c3; M[x] = y; }
        for (int fr = 0; fr < 2; fr++) {
            Map g; for (unsigned x = 0; x < 256; x++) { unsigned u = M[x & 15], v = M[x >> 4]; if (fr) u ^= v; g.f[x] = (u8)(u | (v << 4)); }
            make_cols(g);
            if (apply(g, B1) != B1) continue;
            u32 p1 = apply(g, n1), p2 = apply(g, n1a);
            auto inP = [&](u32 c) { return c == n1 || c == n1a || c == (n1 ^ n1a); };
            if (!inP(p1) || !inP(p2)) continue;
            bool ok = true; Ech28 F; for (u32 b : fixb) F.add(b);
            for (u32 b : fixb) if (F.reduce(apply(g, b))) { ok = false; break; }
            if (ok) G1.push_back(g);
        }
    }
}

// ---------- the extension engine: seed S + r forms from Fix (r <= 4): r-1 enumerated, the last solved ----------
// K_a(S) = {y : B(a, y) = 0 for B in S}, D_a = dim K_a(S)/<a>. The final S8 is APN iff for every a the r new forms,
// restricted to K_a/<a>, span its dual: rank D_a. So after j enumerated forms the rank must be >= D_a - (r - j), and
// at the leaf a deficient a (rank D_a - 1) gives one linear equation t_r(a, y) = 1 (y spanning the common kernel).
struct Engine {
    std::vector<u32> S, C; int m = 0, r = 0;
    int D[256]; u8 sig[16][256]; u8 ybas[256][4];
    std::vector<std::array<u64, 16>> pl; // per combination of C: plane k (< 4) x 4 words, bit a = k-th sig bit
    u64 maskD[5][4];
    long long leaves = 0, passed = 0, sols = 0;
    bool init(const std::vector<u32>& seed, int r_) {
        r = r_; S = seed; Ech28 e; for (u32 x : S) if (!e.add(x)) return false;
        if ((int)S.size() + r != 8) return false;
        C.clear(); for (u32 b : fixb) { if (e.add(b)) C.push_back(b); }
        m = (int)C.size();
        memset(maskD, 0, sizeof maskD);
        std::vector<u64> sm; for (u32 x : S) sm.push_back(coef_to_mat(x));
        for (unsigned a = 1; a < 256; a++) {
            u8 rr[8]; int n = 0; for (u64 x : sm) rr[n++] = rvec(x, a);
            u8 piv[8] = {0}; std::vector<u8> bas;
            auto addv = [&](u8 x) { u8 y = x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; bas.push_back(x); return true; } } return false; };
            addv((u8)a);
            for (unsigned y = 1; y < 256; y++) { bool z = true; for (int i = 0; i < n; i++) if (__popcnt(rr[i] & y) & 1) { z = false; break; } if (z) addv((u8)y); }
            D[a] = (int)bas.size() - 1;
            if (D[a] > r) return false;
            for (int k = 0; k < D[a]; k++) ybas[a][k] = bas[k + 1];
            maskD[D[a]][a >> 6] |= 1ull << (a & 63);
            for (int i = 0; i < m; i++) { u64 cm = coef_to_mat(C[i]); u8 sg = 0; for (int k = 0; k < D[a]; k++) if (bform(cm, a, ybas[a][k])) sg |= 1 << k; sig[i][a] = sg; }
        }
        pl.assign((size_t)1 << m, {});
        for (u32 g = 1; g < (1u << m); g++) {
            int i = _tzcnt_u32(g); auto& P = pl[g]; P = pl[g & (g - 1)];
            for (unsigned a = 1; a < 256; a++) for (int k = 0; k < 4; k++) if (sig[i][a] >> k & 1) P[k * 4 + (a >> 6)] ^= 1ull << (a & 63);
        }
        return true;
    }
    u8 sigof(u32 g, unsigned a) const { u8 s = 0; for (int k = 0; k < 4; k++) if (pl[g][k * 4 + (a >> 6)] >> (a & 63) & 1) s |= 1 << k; return s; }
    // combos[j]: planes of all nonzero combinations of the first j enumerated forms (index = bitmask of members)
    std::array<u64, 16> combo[16];
    u32 gs[4];
    bool level_ok(int j) const {
        // j forms chosen; combos 1..2^j-1 in combo[]
        int nc = (1 << j) - 1;
        for (int w = 0; w < 4; w++) {
            u64 nz[16]; for (int c = 1; c <= nc; c++) nz[c] = combo[c][w] | combo[c][4 + w] | combo[c][8 + w] | combo[c][12 + w];
            u64 all = ~0ull, any = 0, pairzero = 0;
            for (int c = 1; c <= nc; c++) { all &= nz[c]; }
            for (int i = 0; i < j; i++) any |= nz[1 << i];
            for (int c = 1; c <= nc; c++) for (int c2 = c + 1; c2 <= nc; c2++) pairzero |= ~nz[c] & ~nz[c2];
            for (int d = 1; d <= 4; d++) {
                u64 mk = maskD[d][w]; if (!mk) continue;
                int R = d - (r - j);
                if (R <= 0) continue;
                u64 ok;
                if (R == j) ok = all; else if (R == j - 1) ok = ~pairzero; else ok = any; // R == 1 < j - 1
                if (mk & ~ok) return false;
            }
        }
        return true;
    }
    template <class F> void rec(int j, int pmin, unsigned piv, F& found) {
        if (j == r - 1) { leaves++; if (level_ok(j)) { passed++; solve(j, found); } return; }
        if (j > 0 && !level_ok(j)) return;
        int left = (r - 1) - j;
        for (int p = pmin; p <= m - left; p++) {
            unsigned M = ((1u << p) - 1) & ~piv;
            for (unsigned s = M;; s = (s - 1) & M) {
                u32 g = (1u << p) | s; gs[j] = g;
                int base = 1 << j; combo[base] = pl[g];
                for (int c = 1; c < base; c++) for (int k = 0; k < 16; k++) combo[base + c][k] = combo[c][k] ^ pl[g][k];
                rec(j + 1, p + 1, piv | 1u << p, found);
                if (!s) break;
            }
        }
    }
    template <class F> void run(F&& found) { rec(0, 0, 0, found); }
    template <class F> void solve(int j, F& found) {
        u32 rowsv[300]; int nr = 0;
        for (unsigned a = 1; a < 256; a++) {
            int d = D[a]; if (!d) continue;
            u8 v[4]; for (int i = 0; i < j; i++) v[i] = sigof(gs[i], a);
            u8 pv[4] = {0, 0, 0, 0}; int rk = 0;
            for (int i = 0; i < j; i++) { u8 x = v[i]; for (int b = 3; b >= 0 && x; b--) if (x >> b & 1) { if (pv[b]) x ^= pv[b]; else { pv[b] = x; rk++; x = 0; } } }
            if (rk == d) continue;
            if (rk < d - 1) return; // can't happen after level_ok
            u8 z = 0; for (u8 c = 1; c < (1 << d); c++) { bool ok = true; for (int i = 0; i < j; i++) if (__popcnt((unsigned)(v[i] & c)) & 1) { ok = false; break; } if (ok) { z = c; break; } }
            u32 row = 1u << m; for (int i = 0; i < m; i++) if (__popcnt((unsigned)(sig[i][a] & z)) & 1) row |= 1u << i;
            rowsv[nr++] = row;
        }
        int piv[20]; for (int i = 0; i < 20; i++) piv[i] = -1; int rk = 0;
        for (int b = m - 1; b >= 0; b--) {
            int p = -1; for (int i = rk; i < nr; i++) if (rowsv[i] >> b & 1) { p = i; break; }
            if (p < 0) continue;
            std::swap(rowsv[p], rowsv[rk]);
            for (int i = 0; i < nr; i++) if (i != rk && (rowsv[i] >> b & 1)) rowsv[i] ^= rowsv[rk];
            piv[b] = rk; rk++;
        }
        for (int i = rk; i < nr; i++) if (rowsv[i]) return;
        std::vector<int> fr; for (int b = 0; b < m; b++) if (piv[b] < 0) fr.push_back(b);
        auto comb = [&](u32 g) { u32 v = 0; for (int i = 0; i < m; i++) if (g >> i & 1) v ^= C[i]; return v; };
        for (u32 f = 0; f < (1u << fr.size()); f++) {
            u32 x = 0; for (size_t q = 0; q < fr.size(); q++) if (f >> q & 1) x |= 1u << fr[q];
            for (int b = 0; b < m; b++) if (piv[b] >= 0) { u32 rw = rowsv[piv[b]]; int val = (rw >> m) & 1; val ^= __popcnt(rw & x & ((1u << m) - 1) & ~(1u << b)) & 1; if (val) x |= 1u << b; }
            u32 b8[8]; int n = 0; for (u32 y : S) b8[n++] = y;
            for (int i = 0; i < j; i++) b8[n++] = comb(gs[i]);
            b8[n++] = comb(x);
            sols++;
            found(b8);
        }
    }
};

struct Out {
    std::mutex mu; std::set<std::vector<u32>> seen; FILE* f = nullptr; long apn = 0; std::map<std::string, long> types;
    void report(u32* b8, const char* how) {
        Ech28 e; for (int i = 0; i < 8; i++) if (!e.add(b8[i])) return;
        u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(b8[i]);
        if (max_kernel(m, 8) != 1) { std::lock_guard<std::mutex> g(mu); fprintf(stderr, "engine produced a non-APN space\n"); return; }
        int cnt[5]; rank_dist(m, 8, cnt);
        canon8(b8); std::vector<u32> key(b8, b8 + 8);
        std::lock_guard<std::mutex> g(mu);
        if (!seen.insert(key).second) return;
        apn++; char t[64]; snprintf(t, sizeof t, "ranks 8:%d 6:%d 4:%d 2:%d", cnt[4], cnt[3], cnt[2], cnt[1]); types[t]++;
        if (f) { fprintf(f, "%s %s basis", t, how); for (int i = 0; i < 8; i++) fprintf(f, " %07x", b8[i]); fprintf(f, "\n"); fflush(f); }
        if (cnt[2] >= 5) { printf("FOUND %s\n", t); fflush(stdout); }
    }
};

// ---------- validation: bring a known class with an order-3 automorphism into these coordinates ----------
struct Rad { std::vector<u8> el; u8 mem[256]; };
static Rad radical_of(u64 m) { Rad r; memset(r.mem, 0, 256); r.mem[0] = 1; for (unsigned x = 1; x < 256; x++) if (rvec(m, x) == 0) { r.el.push_back((u8)x); r.mem[x] = 1; } return r; }
static bool graph_map(const Rad& R1, const Rad& R2, const Rad& R3, u8* T) {
    u8 p1[256], p2[256]; memset(p1, 0, 256); memset(p2, 0, 256);
    std::vector<u8> a1(R1.el), a2(R2.el); a1.push_back(0); a2.push_back(0);
    for (u8 a : a1) for (u8 b : a2) { p1[a ^ b] = a; p2[a ^ b] = b; }
    memset(T, 0, 256); for (u8 r : R3.el) { if (!p1[r] || !p2[r]) return false; T[p1[r]] = p2[r]; }
    return true;
}
// order-3 automorphisms of S8 (as in aut8)
static std::vector<std::array<u8, 256>> order3_auts(const u32* basis) {
    std::vector<std::array<u8, 256>> res;
    u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(basis[i]);
    Ech28 S; for (int i = 0; i < 8; i++) S.add(basis[i]);
    std::vector<Rad> R; u64 cur = 0;
    for (unsigned g = 1; g < 256; g++) { cur ^= m[_tzcnt_u32(g)]; if (rank_mat(cur) == 4) R.push_back(radical_of(cur)); }
    int k = (int)R.size(); if (k < 3) return res;
    u8 e[4]; { u8 piv[8] = {0}; int n = 0; for (u8 x : R[0].el) { u8 y = x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; y = 0; e[n++] = x; } } if (n == 4) break; } }
    u8 T[256]; graph_map(R[0], R[1], R[2], T);
    for (int s1 = 0; s1 < k; s1++) for (int s2 = 0; s2 < k; s2++) for (int s3 = 0; s3 < k; s3++) {
        if (s1 == s2 || s1 == s3 || s2 == s3) continue;
        u8 Tp[256]; if (!graph_map(R[s1], R[s2], R[s3], Tp)) continue;
        const auto& E = R[s1].el;
        for (u8 u0 : E) for (u8 u1 : E) for (u8 u2 : E) for (u8 u3 : E) {
            u8 u[4] = {u0, u1, u2, u3}; if (rank8(u, 4) != 4) continue;
            u8 b[8], im[8]; for (int i = 0; i < 4; i++) { b[i] = e[i]; im[i] = u[i]; b[4 + i] = T[e[i]]; im[4 + i] = Tp[u[i]]; }
            std::array<u8, 256> Af; { u8 x = 0, y = 0; Af[0] = 0; for (unsigned g = 1; g < 256; g++) { int t = _tzcnt_u32(g); x ^= b[t]; y ^= im[t]; Af[x] = y; } }
            bool ok = true;
            for (int t = 0; t < 8 && ok; t++) { u32 c = 0; for (int q = 0; q < 28; q++) if (bform(m[t], Af[1u << KI[q]], Af[1u << KJ[q]])) c |= 1u << q; if (S.reduce(c)) ok = false; }
            if (!ok) continue;
            bool o3 = true, fpf = true; for (unsigned x = 1; x < 256; x++) { if (Af[Af[Af[x]]] != x) o3 = false; if (Af[x] == x) fpf = false; }
            if (o3 && fpf) res.push_back(Af);
        }
    }
    return res;
}

static int validate(const char* path) {
    std::vector<Found> fs; if (!load_seeds(path, fs)) { fprintf(stderr, "can't read %s\n", path); return 1; }
    int tried = 0, refound = 0;
    for (size_t idx = 0; idx < fs.size(); idx++) {
        Found& F = fs[idx]; if (F.cnt[2] < 4 || F.cnt[1]) continue;
        auto auts = order3_auts(F.basis); if (auts.empty()) continue;
        auto& Ap = auts[0];
        u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(F.basis[i]);
        // a rank-4 component moved by Ap
        u64 cur = 0, mv = 0; bool got = false;
        for (unsigned g = 1; g < 256 && !got; g++) {
            cur ^= m[_tzcnt_u32(g)]; if (rank_mat(cur) != 4) continue;
            Rad R = radical_of(cur); if (!R.mem[Ap[R.el[0]]] || !R.mem[Ap[R.el[1]]] || !R.mem[Ap[R.el[2]]]) { mv = cur; got = true; }
        }
        if (!got) continue;
        Rad R = radical_of(mv);
        u8 r[4]; { u8 piv[8] = {0}; int n = 0; for (u8 x : R.el) { u8 y = x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; y = 0; r[n++] = x; } } if (n == 4) break; } }
        // M: (u, v) -> sum u_i r_i + sum v_i Ap(r_i); then my A = M^-1 Ap M
        Map Mm; for (unsigned x = 0; x < 256; x++) { u8 y = 0; for (int i = 0; i < 4; i++) { if (x >> i & 1) y ^= r[i]; if (x >> (4 + i) & 1) y ^= Ap[r[i]]; } Mm.f[x] = y; }
        make_cols(Mm);
        u32 b2[8]; for (int i = 0; i < 8; i++) b2[i] = apply(Mm, F.basis[i]);
        u32 mvc = apply(Mm, mat_to_coef(mv));
        // normalise the moved component to B1: N in GL(4,2) on u and v alike with mvc o N = B1
        bool norm = false;
        for (unsigned c0 = 1; c0 < 16 && !norm; c0++) for (unsigned c1 = 1; c1 < 16 && !norm; c1++) for (unsigned c2 = 1; c2 < 16 && !norm; c2++) for (unsigned c3 = 1; c3 < 16 && !norm; c3++) {
            u8 vv[4] = {(u8)c0, (u8)c1, (u8)c2, (u8)c3}; if (rank8(vv, 4) != 4) continue;
            unsigned Mx[16]; for (unsigned x = 0; x < 16; x++) { unsigned y = 0; if (x & 1) y ^= c0; if (x & 2) y ^= c1; if (x & 4) y ^= c2; if (x & 8) y ^= c3; Mx[x] = y; }
            Map g; for (unsigned x = 0; x < 256; x++) g.f[x] = (u8)(Mx[x & 15] | (Mx[x >> 4] << 4)); make_cols(g);
            if (apply(g, mvc) != B1) continue;
            for (int i = 0; i < 8; i++) b2[i] = apply(g, b2[i]);
            norm = true;
        }
        if (!norm) { printf("%zu: could not normalise B1\n", idx); continue; }
        // check: invariant under A, contains B1, b = 1
        Ech28 S; for (int i = 0; i < 8; i++) S.add(b2[i]);
        bool inv = true; for (int i = 0; i < 8; i++) if (S.reduce(act(b2[i]))) inv = false;
        int nfix = 0; { u32 c = 0; for (unsigned g = 1; g < 256; g++) { c ^= b2[_tzcnt_u32(g)]; if (act(c) == c) nfix++; } }
        bool hasB1 = !S.reduce(B1), hasP = !S.reduce(n1) && !S.reduce(n1a);
        // Q and a fixed rank-4 component h4 in it
        std::vector<u32> Q; { Ech28 q; u32 c = 0; for (unsigned g = 1; g < 256; g++) { c ^= b2[_tzcnt_u32(g)]; if (act(c) == c && q.add(c)) Q.push_back(c); } }
        u32 h4 = 0; { u32 c = 0; for (unsigned g = 1; g < 256; g++) { c ^= b2[_tzcnt_u32(g)]; if (act(c) == c && rank_mat(coef_to_mat(c)) == 4) { h4 = c; break; } } }
        tried++;
        printf("%zu: invariant %d, contains B1 %d, contains P %d, dim Q %d, fixed rank-4 %s", idx, inv, hasB1, hasP, (int)Q.size(), h4 ? "yes" : "no");
        if (!inv || !hasB1 || !hasP || Q.size() != 6 || !h4) { printf("\n"); continue; }
        // seeds: P + h1 + h4 + q for q in Q outside <h1, h4>; the engine must give back b2's space
        u32 target[8]; memcpy(target, b2, sizeof target); canon8(target);
        bool hit = false; int seeds = 0;
        for (u32 q : Q) {
            Ech28 t; t.add(h1); t.add(h4); if (!t.add(q)) continue;
            Engine E; std::vector<u32> seed = {n1, n1a, h1, h4, q};
            if (!E.init(seed, 3)) { continue; }
            seeds++;
            E.run([&](u32* b8) { u32 c[8]; memcpy(c, b8, sizeof c); canon8(c); if (!memcmp(c, target, sizeof c)) hit = true; });
            if (hit) break;
        }
        printf(", seeds %d, refound %d\n", seeds, hit); fflush(stdout);
        refound += hit;
    }
    printf("validate: %d classes tried, %d refound\n", tried, refound);
    return 0;
}

// run the engine on a list of seeds, all with the same r
static void run_seeds(const std::vector<std::vector<u32>>& seeds, int r, const char* outp, const char* tag, long first, long last) {
    if (last > (long)seeds.size()) last = (long)seeds.size();
    Out out; out.f = fopen(outp, "a");
    std::atomic<long long> leaves(0), passed(0), sols(0), bad(0), done(0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (long s = first; s < last; s++) {
        Engine E; if (!E.init(seeds[s], r)) { bad++; continue; }
        E.run([&](u32* b8) { out.report(b8, tag); });
        leaves += E.leaves; passed += E.passed; sols += E.sols;
        long long d = ++done; if (d % 20 == 0) { printf("  %lld seeds done\n", d); fflush(stdout); }
    }
    printf("%s seeds %ld..%ld (r = %d): leaves %lld, passed %lld, solutions %lld, dead at init %lld, APN spaces %ld\n", tag, first, last, r, leaves.load(), passed.load(), sols.load(), bad.load(), out.apn);
    for (auto& kv : out.types) printf("  %s: %ld\n", kv.first.c_str(), kv.second);
    fclose(out.f);
}

// F_4 = {0, 1, w, w^2} coded 0..3; g0(x, y) in F_4 from B_n1 = Tr g0: (Tr z, Tr wz) = (B_n1(x,y), B_n1(x, wy))
static const u8 F4MUL[4][4] = {{0, 0, 0, 0}, {0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}};
static u8 smul(u8 z, u8 w) { return z == 0 ? 0 : z == 1 ? w : z == 2 ? Am.f[w] : Am.f[Am.f[w]]; }
static u8 g0val(u64 N1m, u8 x, u8 y) { int t = bform(N1m, x, y), tw = bform(N1m, x, Am.f[y]); return t ? (tw ? 2 : 3) : (tw ? 1 : 0); }

// generators of a subgroup of Stab(P) in GammaL(4,4): random symplectic transvections for g0, A, Frobenius
static std::vector<Map> stabP_generators(int ntrans, u64 seed) {
    std::vector<Map> gens; u64 N1m = coef_to_mat(n1); Rng R(seed);
    auto inP = [&](u32 c) { return c == n1 || c == n1a || c == (n1 ^ n1a); };
    Ech28 F; for (u32 b : fixb) F.add(b);
    auto valid = [&](Map& g) {
        for (unsigned x = 0; x < 256; x++) for (unsigned y = 0; y < 256; y += 37) if (g.f[x ^ y] != (g.f[x] ^ g.f[y])) return false;
        u8 seen[256] = {0}; for (unsigned x = 0; x < 256; x++) { if (seen[g.f[x]]) return false; seen[g.f[x]] = 1; }
        make_cols(g); if (!inP(apply(g, n1)) || !inP(apply(g, n1a))) return false;
        for (u32 b : fixb) if (F.reduce(apply(g, b))) return false;
        return true;
    };
    for (int t = 0; t < ntrans; ) {
        u8 w = (u8)(1 + R.below(255)); u8 lam = (u8)(1 + R.below(3));
        Map g; for (unsigned x = 0; x < 256; x++) g.f[x] = (u8)(x ^ smul(F4MUL[lam][g0val(N1m, (u8)x, w)], w));
        if (!valid(g)) { fprintf(stderr, "transvection not in Stab(P): g0 reconstruction wrong\n"); exit(1); }
        gens.push_back(g); t++;
    }
    gens.push_back(Am);
    Map fr; for (unsigned x = 0; x < 256; x++) { unsigned u = x & 15, v = x >> 4; fr.f[x] = (u8)((u ^ v) | (v << 4)); }
    if (valid(fr)) gens.push_back(fr); else fprintf(stderr, "(Frobenius doesn't fix P)\n");
    return gens;
}
struct UF { std::vector<int> p; UF(size_t n) : p(n) { for (size_t i = 0; i < n; i++) p[i] = (int)i; } int f(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; } void u(int a, int b) { a = f(a); b = f(b); if (a != b) { if (a < b) p[b] = a; else p[a] = b; } } };

int main(int argc, char** argv) {
    init_index(); setup(); make_G1();
    printf("c3k4: dim Fix 16, |G1| = %zu\n", G1.size()); fflush(stdout);
    if (argc >= 3 && !strcmp(argv[1], "validate")) return validate(argv[2]);
    long first = argc > 4 ? atol(argv[3]) : 0, lastarg = argc > 4 ? atol(argv[4]) : 0x7fffffffL;
    Set256 RB[3]; { u32 b = B1; for (int i = 0; i < 3; i++) { radical_set(coef_to_mat(b), RB[i]); b = act(b); } }
    if (argc >= 3 && !strcmp(argv[1], "search33")) {
        // (6,1), a second rank-4 3-cycle {h' + g : g in P \ 0}; seed h' up to G1
        std::vector<std::vector<u32>> seeds; long cand = 0;
        u32 c = 0;
        for (u32 g = 1; g < (1u << 16); g++) {
            c ^= fixb[_tzcnt_u32(g)];
            if (c == h1 || rank_mat(coef_to_mat(c ^ n1)) != 4) continue;
            bool okr = true; { u32 b = c ^ n1; for (int i = 0; i < 3 && okr; i++) { Set256 R; radical_set(coef_to_mat(b), R); if (R.meets(RB[0]) || R.meets(RB[1]) || R.meets(RB[2])) okr = false; b = act(b); } }
            if (!okr) continue;
            cand++;
            bool canon = true; for (auto& gg : G1) if (apply(gg, c) < c) { canon = false; break; }
            if (!canon) continue;
            u32 s4[4] = {n1, n1a, h1, c}; if (!partition_ok(s4, 4)) continue;
            seeds.push_back({n1, n1a, h1, c});
        }
        printf("3+3: %ld candidates h', %zu seeds up to G1 with partition ok\n", cand, seeds.size()); fflush(stdout);
        run_seeds(seeds, 4, argv[2], "c3k4-33", first, lastarg);
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "search2")) {
        // (4,2) with a rank-4 3-cycle B1 (standard): P = F4 n1 + F4 n2 anisotropic; 3+1+1 (r = 1) and 3+3 (r = 2)
        std::vector<u32> nb; { Ech28 e; e.add(n1); e.add(n1a); for (int k = 0; k < 28; k++) { u32 v = (1u << k) ^ fixpart(1u << k); if (e.add(v)) { nb.push_back(v); e.add(act(v)); } } }
        if (nb.size() != 5) { fprintf(stderr, "N basis %zu\n", nb.size()); return 1; }
        auto f4mul = [&](unsigned c, u32 v) { u32 r = 0; if (c & 1) r ^= v; if (c & 2) r ^= act(v); return r; }; // c0 + c1 A
        auto spankey = [&](u32 a, u32 b) { u32 bs[8] = {a, act(a), b, act(b), 0, 0, 0, 0}; Ech28 e; for (int i = 0; i < 4; i++) e.add(bs[i]); for (int x = 27; x >= 0; x--) if (e.piv[x]) for (int y = x + 1; y < 28; y++) if (e.piv[y] && (e.piv[y] >> x & 1)) e.piv[y] ^= e.piv[x]; std::array<u32, 4> k{}; int n = 0; for (int x = 27; x >= 0; x--) if (e.piv[x]) k[n++] = e.piv[x]; return k; };
        std::vector<u32> n2s;
        for (int lead = 0; lead < 5; lead++) for (u32 t = 0; t < (1u << (2 * (4 - lead))); t++) {
            u32 v = nb[lead]; for (int j = 0; j < 4 - lead; j++) v ^= f4mul(t >> (2 * j) & 3, nb[lead + 1 + j]);
            u32 bs[4] = {n1, n1a, v, act(v)}; bool aniso = true; u32 c = 0;
            for (unsigned g = 1; g < 16 && aniso; g++) { c ^= bs[_tzcnt_u32(g)]; if (rank_mat(coef_to_mat(c)) != 8) aniso = false; }
            if (aniso) n2s.push_back(v);
        }
        // orbits of P under G1 by its span key
        std::map<std::array<u32, 4>, int> keyidx; for (size_t i = 0; i < n2s.size(); i++) keyidx[spankey(n1, n2s[i])] = (int)i;
        std::vector<int> rep(n2s.size(), -1); std::vector<int> Preps;
        for (size_t i = 0; i < n2s.size(); i++) { if (rep[i] >= 0) continue; Preps.push_back((int)i); for (auto& g : G1) { auto it = keyidx.find(spankey(apply(g, n1), apply(g, n2s[i]))); if (it == keyidx.end()) { fprintf(stderr, "P orbit\n"); return 1; } rep[it->second] = (int)i; } }
        printf("(4,2): %zu anisotropic P through n1, %zu up to G1\n", n2s.size(), Preps.size()); fflush(stdout);
        std::vector<u32> H; { u32 c = 0; for (u32 g = 1; g < (1u << 16); g++) { c ^= fixb[_tzcnt_u32(g)]; if (rank_mat(coef_to_mat(c)) != 4) continue; Set256 R; radical_set(coef_to_mat(c), R); if (R.meets(RB[0]) || R.meets(RB[1]) || R.meets(RB[2])) continue; H.push_back(c); } }
        std::vector<std::vector<u32>> seeds311, seeds33;
        for (int pi : Preps) {
            u32 n2 = n2s[pi]; auto key = spankey(n1, n2);
            std::vector<const Map*> St; for (auto& g : G1) if (spankey(apply(g, n1), apply(g, n2)) == key) St.push_back(&g);
            std::vector<Set256> HR(H.size()); for (size_t i = 0; i < H.size(); i++) radical_set(coef_to_mat(H[i]), HR[i]);
            for (size_t i = 0; i < H.size(); i++) for (size_t j = i + 1; j < H.size(); j++) {
                if (HR[i].meets(HR[j])) continue;
                u32 lo = std::min(H[i], H[j]), hi = std::max(H[i], H[j]); bool canon = true;
                for (auto* g : St) { u32 a = apply(*g, H[i]), b = apply(*g, H[j]); if (a > b) std::swap(a, b); if (a < lo || (a == lo && b < hi)) { canon = false; break; } }
                if (!canon) continue;
                u32 s7[7] = {n1, n1a, n2, act(n2), h1, H[i], H[j]}; if (!partition_ok(s7, 7)) continue;
                seeds311.push_back({n1, n1a, n2, act(n2), h1, H[i], H[j]});
            }
            // 3+3: h' + g' rank 4, g' over the 5 F4-lines of P, h' up to St
            u32 lines[5] = {n1, n2, n1 ^ n2, n1 ^ act(n2), n1 ^ act(act(n2))};
            u32 c = 0;
            for (u32 g = 1; g < (1u << 16); g++) {
                c ^= fixb[_tzcnt_u32(g)]; if (c == h1) continue;
                for (int l = 0; l < 5; l++) {
                    u32 b = c ^ lines[l]; if (rank_mat(coef_to_mat(b)) != 4) continue;
                    bool okr = true; u32 bb = b; for (int i = 0; i < 3 && okr; i++) { Set256 R; radical_set(coef_to_mat(bb), R); if (R.meets(RB[0]) || R.meets(RB[1]) || R.meets(RB[2])) okr = false; bb = act(bb); }
                    if (!okr) continue;
                    // canonical: c minimal in its St-orbit, and the line l minimal among lines that map to it is not needed:
                    // the seed is (P, h1, c) whatever l is, so dedup by c only
                    bool canon = true; for (auto* gg : St) if (apply(*gg, c) < c) { canon = false; break; }
                    if (!canon) break;
                    u32 s6[6] = {n1, n1a, n2, act(n2), h1, c}; if (partition_ok(s6, 6)) seeds33.push_back({n1, n1a, n2, act(n2), h1, c});
                    break;
                }
            }
        }
        printf("(4,2) seeds: 3+1+1 %zu, 3+3 %zu\n", seeds311.size(), seeds33.size()); fflush(stdout);
        run_seeds(seeds311, 1, argv[2], "c3k4-311-b2", 0, 0x7fffffffL);
        run_seeds(seeds33, 2, argv[2], "c3k4-33-b2", 0, 0x7fffffffL);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "searchf")) {
        // >= 5 fixed rank-4 components, b = 0 or 2. Any two fixed ones are, up to GammaL(4,4), the standard pair
        // h4 = x3 x3' + x4 x4', h5 = x1 x1' + x2 x2' (F4 coordinates x_i = u_i + w v_i, ' = conjugate), with
        // stabiliser (U(2,2) x U(2,2)) : (swap x Frobenius), listed explicitly.
        int bsel = atoi(argv[2]); const char* outp = argv[3];
        auto co = [](unsigned x, int i) { return (u8)(((x >> i) & 1) | (((x >> (4 + i)) & 1) << 1)); };
        auto enc = [](const u8* c) { unsigned x = 0; for (int i = 0; i < 4; i++) { x |= (c[i] & 1u) << i; x |= ((c[i] >> 1) & 1u) << (4 + i); } return (u8)x; };
        auto conj = [](u8 z) { return (u8)(z == 2 ? 3 : z == 3 ? 2 : z); };
        auto herm = [&](const u8 H[4][4]) { // B = Tr h(x, y), h = sum H_ij x_i conj(y_j)
            u32 c = 0;
            for (int q = 0; q < 28; q++) {
                unsigned x = 1u << KI[q], y = 1u << KJ[q]; u8 s = 0;
                for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) s ^= F4MUL[F4MUL[H[i][j]][co(x, i)]][conj(co(y, j))];
                if (s >= 2) c |= 1u << q;
            }
            return c;
        };
        u8 H4[4][4] = {{0}}, H5[4][4] = {{0}}; H4[2][2] = H4[3][3] = 1; H5[0][0] = H5[1][1] = 1;
        u32 h4 = herm(H4), h5 = herm(H5);
        Set256 R4, R5; radical_set(coef_to_mat(h4), R4); radical_set(coef_to_mat(h5), R5);
        if (act(h4) != h4 || act(h5) != h5 || rank_mat(coef_to_mat(h4)) != 4 || rank_mat(coef_to_mat(h5)) != 4 || R4.meets(R5)) { fprintf(stderr, "standard pair\n"); return 1; }
        // U(2,2): 2x2 F4 matrices M with M^* M = I
        std::vector<std::array<u8, 4>> U2;
        for (unsigned t = 0; t < 256; t++) {
            u8 m[4] = {(u8)(t & 3), (u8)(t >> 2 & 3), (u8)(t >> 4 & 3), (u8)(t >> 6 & 3)}; // m[0] m[1] / m[2] m[3]
            // columns c0 = (m0, m2), c1 = (m1, m3): <c_i, c_j> = sum_k c_i[k] conj(c_j[k]) = delta
            auto ip = [&](int i, int j) { u8 a = m[i], b = m[2 + i], c = m[j], d = m[2 + j]; return (u8)(F4MUL[a][conj(c)] ^ F4MUL[b][conj(d)]); };
            if (ip(0, 0) == 1 && ip(1, 1) == 1 && ip(0, 1) == 0) U2.push_back({m[0], m[1], m[2], m[3]});
        }
        std::vector<Map> Gp;
        for (auto& M1 : U2) for (auto& M2 : U2) for (int sw = 0; sw < 2; sw++) for (int fr = 0; fr < 2; fr++) {
            Map g;
            for (unsigned x = 0; x < 256; x++) {
                u8 c[4], d[4]; for (int i = 0; i < 4; i++) { c[i] = co(x, i); if (fr) c[i] = conj(c[i]); }
                d[0] = F4MUL[M1[0]][c[0]] ^ F4MUL[M1[1]][c[1]]; d[1] = F4MUL[M1[2]][c[0]] ^ F4MUL[M1[3]][c[1]];
                d[2] = F4MUL[M2[0]][c[2]] ^ F4MUL[M2[1]][c[3]]; d[3] = F4MUL[M2[2]][c[2]] ^ F4MUL[M2[3]][c[3]];
                if (sw) { std::swap(d[0], d[2]); std::swap(d[1], d[3]); }
                g.f[x] = enc(d);
            }
            make_cols(g);
            u32 a = apply(g, h4), b = apply(g, h5);
            if (!((a == h4 && b == h5) || (a == h5 && b == h4))) { fprintf(stderr, "Gp element doesn't fix the pair\n"); return 1; }
            Ech28 F; for (u32 x : fixb) F.add(x); bool ok = true; for (u32 x : fixb) if (F.reduce(apply(g, x))) ok = false;
            if (!ok) { fprintf(stderr, "Gp element doesn't fix Fix\n"); return 1; }
            Gp.push_back(g);
        }
        printf("standard fixed pair ok; |U(2,2)| = %zu, |Gp| = %zu\n", U2.size(), Gp.size()); fflush(stdout);
        std::vector<std::vector<u32>> seeds; int r = 0;
        auto key4 = [&](std::vector<u32> v) { Ech28 e; for (u32 x : v) e.add(x); for (int x = 27; x >= 0; x--) if (e.piv[x]) for (int y = x + 1; y < 28; y++) if (e.piv[y] && (e.piv[y] >> x & 1)) e.piv[y] ^= e.piv[x]; std::vector<u32> k; for (int x = 27; x >= 0; x--) if (e.piv[x]) k.push_back(e.piv[x]); return k; };
        if (bsel == 2) {
            // anisotropic F4-lines of N, up to Gp
            std::vector<u32> nbv; { Ech28 e; for (int k = 0; k < 28; k++) { u32 v = (1u << k) ^ fixpart(1u << k); if (e.add(v)) { nbv.push_back(v); e.add(act(v)); } } }
            auto f4mul = [&](unsigned c, u32 v) { u32 rr = 0; if (c & 1) rr ^= v; if (c & 2) rr ^= act(v); return rr; };
            std::vector<u32> pts;
            for (int lead = 0; lead < 6; lead++) for (u32 t = 0; t < (1u << (2 * (5 - lead))); t++) { u32 v = nbv[lead]; for (int j = 0; j < 5 - lead; j++) v ^= f4mul(t >> (2 * j) & 3, nbv[lead + 1 + j]); if (rank_mat(coef_to_mat(v)) == 8) pts.push_back(v); }
            std::set<std::vector<u32>> lines;
            for (size_t i = 0; i < pts.size(); i++) for (size_t j = i + 1; j < pts.size(); j++) {
                u32 bs[4] = {pts[i], act(pts[i]), pts[j], act(pts[j])}; bool an = true; u32 c = 0;
                for (unsigned g = 1; g < 16 && an; g++) { c ^= bs[_tzcnt_u32(g)]; if (rank_mat(coef_to_mat(c)) != 8) an = false; }
                if (an) lines.insert(key4({bs[0], bs[1], bs[2], bs[3]}));
            }
            std::set<std::vector<u32>> done; long orb = 0;
            for (auto& L : lines) {
                if (done.count(L)) continue; orb++;
                for (auto& g : Gp) { std::vector<u32> im; for (u32 x : L) im.push_back(apply(g, x)); done.insert(key4(im)); }
                std::vector<u32> s6(L.begin(), L.end()); s6.push_back(h4); s6.push_back(h5);
                if (partition_ok(s6.data(), 6)) seeds.push_back(s6);
            }
            printf("nondegenerate points %zu, anisotropic lines %zu, up to Gp %ld, seeds (partition ok) %zu\n", pts.size(), lines.size(), orb, seeds.size()); fflush(stdout);
            r = 2;
        } else {
            // b = 0: three more fixed rank-4 forms, each radical complementary to the others, up to stabiliser chains
            std::vector<u32> H; { u32 c = 0; for (u32 g = 1; g < (1u << 16); g++) { c ^= fixb[_tzcnt_u32(g)]; if (rank_mat(coef_to_mat(c)) == 4) H.push_back(c); } }
            std::vector<Set256> HR(H.size()); for (size_t i = 0; i < H.size(); i++) radical_set(coef_to_mat(H[i]), HR[i]);
            std::function<void(std::vector<u32>&, std::vector<int>&, std::vector<const Map*>&)> grow = [&](std::vector<u32>& cur, std::vector<int>& used, std::vector<const Map*>& St) {
                if (cur.size() == 5) {
                    std::vector<u32> k = key4(cur);
                    if ((int)k.size() < 4) { printf("warning: 5 fixed forms span only %zu dims, skipped\n", k.size()); return; }
                    if (!partition_ok(k.data(), (int)k.size())) return;
                    seeds.push_back(k); return;
                }
                Set256 Rs[5]; for (size_t i = 0; i < cur.size(); i++) radical_set(coef_to_mat(cur[i]), Rs[i]);
                std::set<u32> seen;
                for (size_t i = 0; i < H.size(); i++) {
                    bool ok = true; for (size_t q = 0; q < cur.size(); q++) if (HR[i].meets(Rs[q])) { ok = false; break; }
                    if (!ok || seen.count(H[i])) continue;
                    std::vector<const Map*> St2;
                    for (auto* g : St) { u32 im = apply(*g, H[i]); seen.insert(im); if (im == H[i]) St2.push_back(g); }
                    std::vector<u32> nx = cur; nx.push_back(H[i]);
                    // keep only stabiliser elements that fix the new form too (they fix the earlier set already)
                    grow(nx, used, St2);
                }
            };
            std::vector<u32> cur = {h4, h5}; std::vector<int> used; std::vector<const Map*> St; for (auto& g : Gp) St.push_back(&g);
            grow(cur, used, St);
            std::set<std::vector<u32>> uniq(seeds.begin(), seeds.end()); seeds.assign(uniq.begin(), uniq.end());
            printf("(8,0) seeds (5 fixed rank-4 forms, partition ok): %zu\n", seeds.size()); fflush(stdout);
        }
        if (bsel == 2) run_seeds(seeds, r, outp, "c3k4-ff-b2", 0, 0x7fffffffL);
        else {
            std::vector<std::vector<u32>> s5, s4; for (auto& s : seeds) (s.size() == 5 ? s5 : s4).push_back(s);
            run_seeds(s5, 3, outp, "c3k4-f-b0", 0, 0x7fffffffL);
            if (!s4.empty()) run_seeds(s4, 4, outp, "c3k4-f-b0-dep", 0, 0x7fffffffL);
        }
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "searchff")) {
        // (6,1), two fixed rank-4 components h4, h5 (pairs up to a subgroup of Stab(P)); finds every (6,1) space
        // with >= 2 fixed rank-4 components
        auto gens = stabP_generators(12, 12345);
        std::vector<u32> H; std::unordered_map<u32, int> idx;
        { u32 c = 0; for (u32 g = 1; g < (1u << 16); g++) { c ^= fixb[_tzcnt_u32(g)]; if (rank_mat(coef_to_mat(c)) == 4) { idx[c] = (int)H.size(); H.push_back(c); } } }
        size_t n = H.size();
        std::vector<std::vector<int>> perm(gens.size(), std::vector<int>(n));
        for (size_t gi = 0; gi < gens.size(); gi++) for (size_t i = 0; i < n; i++) { auto it = idx.find(apply(gens[gi], H[i])); if (it == idx.end()) { fprintf(stderr, "rank-4 Hermitians not invariant\n"); return 1; } perm[gi][i] = it->second; }
        { UF u1(n); for (size_t gi = 0; gi < gens.size(); gi++) for (size_t i = 0; i < n; i++) u1.u((int)i, perm[gi][i]); long k = 0; for (size_t i = 0; i < n; i++) if (u1.f((int)i) == (int)i) k++; printf("rank-4 Hermitian forms: %zu, orbits under %zu generators: %ld\n", n, gens.size(), k); }
        std::vector<Set256> HR(n); for (size_t i = 0; i < n; i++) radical_set(coef_to_mat(H[i]), HR[i]);
        UF uf(n * n);
        for (size_t i = 0; i < n; i++) for (size_t j = i + 1; j < n; j++) {
            if (HR[i].meets(HR[j])) continue;
            for (size_t gi = 0; gi < gens.size(); gi++) { int a = perm[gi][i], b = perm[gi][j]; if (a > b) std::swap(a, b); uf.u((int)(i * n + j), (int)(a * n + b)); }
        }
        std::vector<std::vector<u32>> seeds; long orbits = 0;
        for (size_t i = 0; i < n; i++) for (size_t j = i + 1; j < n; j++) {
            if (HR[i].meets(HR[j])) continue;
            if (uf.f((int)(i * n + j)) != (int)(i * n + j)) continue;
            orbits++;
            u32 s4[4] = {n1, n1a, H[i], H[j]}; if (!partition_ok(s4, 4)) continue;
            seeds.push_back({n1, n1a, H[i], H[j]});
        }
        printf("fixed pairs: %ld orbits, %zu with partition ok\n", orbits, seeds.size()); fflush(stdout);
        run_seeds(seeds, 4, argv[2], "c3k4-ff", first, lastarg);
        return 0;
    }
    if (argc < 3 || strcmp(argv[1], "search")) { fprintf(stderr, "usage: c3k4 search <out> [first last] | c3k4 validate <seedfile>\n"); return 2; }
    // fixed rank-4 candidates: Hermitian rank-4 forms with radical missing the 3 radicals of B1's orbit
    std::vector<u32> H; std::unordered_map<u32, int> idx;
    { u32 c = 0; for (u32 g = 1; g < (1u << 16); g++) { c ^= fixb[_tzcnt_u32(g)]; if (rank_mat(coef_to_mat(c)) != 4) continue; Set256 R; radical_set(coef_to_mat(c), R); if (R.meets(RB[0]) || R.meets(RB[1]) || R.meets(RB[2])) continue; idx[c] = (int)H.size(); H.push_back(c); } }
    printf("fixed rank-4 candidates: %zu\n", H.size()); fflush(stdout);
    std::vector<std::vector<int>> perm(G1.size(), std::vector<int>(H.size()));
    for (size_t gi = 0; gi < G1.size(); gi++) for (size_t i = 0; i < H.size(); i++) { auto it = idx.find(apply(G1[gi], H[i])); if (it == idx.end()) { fprintf(stderr, "H not G1-invariant\n"); return 1; } perm[gi][i] = it->second; }
    std::vector<Set256> HR(H.size()); for (size_t i = 0; i < H.size(); i++) radical_set(coef_to_mat(H[i]), HR[i]);
    // canonical pairs {i < j}
    std::vector<std::pair<int, int>> seeds;
    #pragma omp parallel for schedule(dynamic, 16)
    for (long i = 0; i < (long)H.size(); i++) {
        for (size_t j = i + 1; j < H.size(); j++) {
            if (HR[i].meets(HR[j])) continue;
            bool canon = true;
            for (size_t gi = 0; gi < G1.size() && canon; gi++) { int a = perm[gi][i], b = perm[gi][j]; if (a > b) std::swap(a, b); if (a < i || (a == i && b < (int)j)) canon = false; }
            if (!canon) continue;
            u32 s5[5] = {n1, n1a, h1, H[i], H[(int)j]};
            if (!partition_ok(s5, 5)) continue;
            #pragma omp critical
            seeds.push_back({(int)i, (int)j});
        }
    }
    std::sort(seeds.begin(), seeds.end());
    printf("seed pairs (up to G1, partition ok): %zu\n", seeds.size()); fflush(stdout);
    long last = argc > 4 ? lastarg : (long)seeds.size();
    if (last > (long)seeds.size()) last = (long)seeds.size();
    Out out; out.f = fopen(argv[2], "a");
    std::atomic<long long> leaves(0), passed(0), bad(0), done(0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (long s = first; s < last; s++) {
        Engine E; std::vector<u32> seed = {n1, n1a, h1, H[seeds[s].first], H[seeds[s].second]};
        if (!E.init(seed, 3)) { bad++; continue; }
        E.run([&](u32* b8) { out.report(b8, "c3k4-311"); });
        leaves += E.leaves; passed += E.passed;
        long long d = ++done; if (d % 200 == 0) { printf("  %lld seeds done\n", d); fflush(stdout); }
    }
    printf("seeds %ld..%ld: leaves %lld, passed filter %lld, seeds dead at init %lld, APN spaces %ld\n", first, last, leaves.load(), passed.load(), bad.load(), out.apn);
    for (auto& kv : out.types) printf("  %s: %ld\n", kv.first.c_str(), kv.second);
    fclose(out.f);
    return 0;
}
