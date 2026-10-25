// fourteen_check.cpp -- independent check of Open Problem VI.1 (arXiv 2510.12008) for i = 3: [0^226, 2^15, 4^14].
// Subagent (Claude Opus 5.5) for Rowen, 2026-10-24. Written without reading sixteen.cpp or runs/a/fourteen*.txt.
// Utilities (forms, radicals, generators, verifier) adapted from sixteen_check.cpp (the earlier checker's).
//
// Coordinates: v = (x, y), x = low nibble, y = high nibble. X = {y = 0}, Y = {x = 0}, Z = graph(I) = {y = x}.
// A 4-space complementary to Y is graph(N) = {(x, Nx)}; complementary to X iff N invertible, to Z iff N + I is.
// V5824 = {N : N, N + I invertible}; two graphs are disjoint iff N + N' is invertible ("adjacent").
//
// Part A (mode holes): every partial spread of 14 solids of PG(7,2) extends to a spread. Three of its solids are
//   moved to Y, X, Z; the other 11 are an 11-clique of V5824. Each 11-clique (up to GL(4,2)-conjugation, by the
//   minimal-class trick) is checked to have a triangle in its common neighbourhood.
// Part B (mode spreads): 14-cliques of V5824 (the rank-4 radicals once the holes are X, Y, Z) up to
//   G = GL(4,2)-conjugation x S3 (S3 = <N -> N^-1, N -> N + I>, which permutes X, Y, Z), orbit reps written out.
// Part C (mode search): per orbit rep, a point-branching search over the forms.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <mutex>
#include <set>
#include <map>
#include <chrono>
#include <functional>
#include <array>
#include <unordered_map>
#include <climits>
#include <random>
#include <atomic>
#include <intrin.h>
#include <immintrin.h>
#include <omp.h>
using namespace std;
typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16; typedef uint8_t u8;

static inline int par(u32 v) { return __popcnt(v) & 1; }
static inline u8 rowof(u64 f, int i) { return (u8)(f >> (8 * i)); }
static inline u8 applyf(u64 f, u8 v) { u8 r = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) r ^= rowof(f, i); return r; }
static inline int evalf(u64 f, u8 u, u8 w) { return par(u & applyf(f, w)); }

static int rank8(u64 f) {
    u8 r[8]; for (int i = 0; i < 8; i++) r[i] = rowof(f, i);
    int rk = 0;
    for (int c = 0; c < 8; c++) {
        int p = -1; for (int i = rk; i < 8; i++) if (r[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        u8 t = r[p]; r[p] = r[rk]; r[rk] = t;
        for (int i = rk + 1; i < 8; i++) if (r[i] >> c & 1) r[i] ^= t;
        rk++;
    }
    return rk;
}

struct Mask { u64 w[4]; };
static inline void mclear(Mask& m) { m.w[0] = m.w[1] = m.w[2] = m.w[3] = 0; }
static inline void mset(Mask& m, int v) { m.w[v >> 6] |= 1ULL << (v & 63); }
static inline bool mget(const Mask& m, int v) { return m.w[v >> 6] >> (v & 63) & 1; }
static inline bool mdisj(const Mask& a, const Mask& b) { return ((a.w[0] & b.w[0]) | (a.w[1] & b.w[1]) | (a.w[2] & b.w[2]) | (a.w[3] & b.w[3])) == 0; }
static inline bool msub(const Mask& a, const Mask& b) { return ((a.w[0] & ~b.w[0]) | (a.w[1] & ~b.w[1]) | (a.w[2] & ~b.w[2]) | (a.w[3] & ~b.w[3])) == 0; }
static inline void mor(Mask& a, const Mask& b) { for (int i = 0; i < 4; i++) a.w[i] |= b.w[i]; }
static inline int mcount(const Mask& m) { return (int)(__popcnt64(m.w[0]) + __popcnt64(m.w[1]) + __popcnt64(m.w[2]) + __popcnt64(m.w[3])); }
static inline bool meq(const Mask& a, const Mask& b) { return !memcmp(&a, &b, sizeof(Mask)); }

static int radical(u64 f, Mask& m) {
    mclear(m); int cnt = 0; u8 cur = 0, v = 0;
    for (u32 i = 1; i < 256; i++) { int b = _tzcnt_u32(i); v ^= (u8)(1u << b); cur ^= rowof(f, b); if (!cur) { mset(m, v); cnt++; } }
    return cnt;
}

// ---- forms with a prescribed radical: pullbacks of nondegenerate forms along functionals u_1..u_m
static u64 pairform(u8 u, u8 w) {
    u64 f = 0;
    for (int i = 0; i < 8; i++) { u8 r = 0; if (u >> i & 1) r ^= w; if (w >> i & 1) r ^= u; f |= (u64)r << (8 * i); }
    return f;
}
static vector<u32> nondeg_codes(int m) {
    int np = m * (m - 1) / 2; vector<u32> out;
    for (u32 c = 0; c < (1u << np); c++) {
        u8 rows[8] = { 0 }; int p = 0;
        for (int k = 0; k < m; k++) for (int l = k + 1; l < m; l++) { if (c >> p & 1) { rows[k] |= 1 << l; rows[l] |= 1 << k; } p++; }
        int rk = 0;
        for (int col = 0; col < m; col++) {
            int piv = -1; for (int i = rk; i < m; i++) if (rows[i] >> col & 1) { piv = i; break; }
            if (piv < 0) continue;
            swap(rows[piv], rows[rk]);
            for (int i = 0; i < m; i++) if (i != rk && (rows[i] >> col & 1)) rows[i] ^= rows[rk];
            rk++;
        }
        if (rk == m) out.push_back(c);
    }
    return out;
}
static vector<u32> ND4, ND6;
static void forms_from(const u8* u, int m, const vector<u32>& codes, vector<u64>& out) {
    u64 Q[15]; int p = 0;
    for (int k = 0; k < m; k++) for (int l = k + 1; l < m; l++) Q[p++] = pairform(u[k], u[l]);
    for (u32 c : codes) { u64 f = 0; u32 cc = c; while (cc) { int b = _tzcnt_u32(cc); cc &= cc - 1; f ^= Q[b]; } out.push_back(f); }
}

// ---- 4x4 matrices: u16, row r in bits 4r..4r+3
static inline u8 mv4(u16 N, u8 x) { u8 r = 0; for (int i = 0; i < 4; i++) r |= par(((N >> (4 * i)) & 15) & x) << i; return r; }
static inline u16 mm4(u16 A, u16 B) {
    u16 C = 0;
    for (int r = 0; r < 4; r++) { u8 row = 0; for (int k = 0; k < 4; k++) if (A >> (4 * r + k) & 1) row ^= (B >> (4 * k)) & 15; C |= (u16)row << (4 * r); }
    return C;
}
static bool inv4(u16 N) { for (u8 x = 1; x < 16; x++) if (!mv4(N, x)) return false; return true; }
static u16 fromcols4(const u8 col[4]) { u16 N = 0; for (int i = 0; i < 4; i++) for (int r = 0; r < 4; r++) if (col[i] >> r & 1) N |= 1 << (4 * r + i); return N; }
static const u16 ID4 = 0x8421;

static vector<u16> GL4; static vector<u16> GL4inv;
static u16 inverse4(u16 N);

static int gl4idx[65536];
static u16 inverse4(u16 N) { return GL4inv[gl4idx[N]]; }
// ---- V5824 and its graph
static const int NW = 91; // 5824 / 64
struct BS { u64 w[NW]; };
static vector<u16> VN;              // vertex -> matrix
static vector<int> cls;             // vertex -> conjugacy class index
static vector<int> classSize, classRep;
static int vid[65536];
static vector<BS> adj;
static inline void bset(BS& b, int i) { b.w[i >> 6] |= 1ULL << (i & 63); }
static inline bool bget(const BS& b, int i) { return b.w[i >> 6] >> (i & 63) & 1; }
static inline int bcount(const BS& b) { int c = 0; for (int i = 0; i < NW; i++) c += (int)__popcnt64(b.w[i]); return c; }
static inline void band(BS& o, const BS& a, const BS& b) { for (int i = 0; i < NW; i++) o.w[i] = a.w[i] & b.w[i]; }
static inline bool bempty(const BS& b) { for (int i = 0; i < NW; i++) if (b.w[i]) return false; return true; }

static void init_V() {
    for (u32 N = 0; N < 65536; N++) if (inv4((u16)N)) GL4.push_back((u16)N);
    GL4inv.resize(GL4.size()); for (size_t i = 0; i < GL4.size(); i++) gl4idx[GL4[i]] = (int)i;
    for (size_t i = 0; i < GL4.size(); i++) { u16 g = GL4[i], p = g, prev = ID4; while (p != ID4) { prev = p; p = mm4(p, g); } GL4inv[i] = (g == ID4) ? ID4 : prev; if (mm4(GL4inv[i], g) != ID4) { fprintf(stderr, "inverse\n"); exit(1); } }
    vector<u16> raw; for (u16 N : GL4) if (inv4(N ^ ID4)) raw.push_back(N);
    // conjugacy classes, numbered by first appearance in increasing N
    unordered_map<u16, int> c0; int nc = 0;
    for (u16 N : raw) if (!c0.count(N)) { for (size_t i = 0; i < GL4.size(); i++) c0[mm4(mm4(GL4[i], N), GL4inv[i])] = nc; nc++; }
    // vertices sorted by (class, N)
    vector<pair<int, u16>> srt; for (u16 N : raw) srt.push_back({ c0[N], N }); sort(srt.begin(), srt.end());
    memset(vid, -1, sizeof vid);
    classSize.assign(nc, 0); classRep.assign(nc, -1);
    for (auto& p : srt) { int v = (int)VN.size(); VN.push_back(p.second); cls.push_back(p.first); vid[p.second] = v; classSize[p.first]++; if (classRep[p.first] < 0) classRep[p.first] = v; }
    adj.assign(VN.size(), BS{});
    for (size_t i = 0; i < VN.size(); i++) for (size_t j = 0; j < VN.size(); j++) if (i != j && inv4(VN[i] ^ VN[j])) bset(adj[i], (int)j);
}

// ---- G = conj x S3 acting on vertices
static vector<u16> s3perm[6]; // s3perm[k][v]
static vector<vector<u16>> conjTab; // conjTab[h][v] = h N h^-1
static void init_group(bool tables) {
    int n = (int)VN.size();
    vector<u16> s(n), t(n), id(n);
    for (int v = 0; v < n; v++) { id[v] = v; s[v] = (u16)vid[inverse4(VN[v])]; t[v] = (u16)vid[VN[v] ^ ID4]; }
    auto comp = [&](const vector<u16>& a, const vector<u16>& b) { vector<u16> c(n); for (int v = 0; v < n; v++) c[v] = a[b[v]]; return c; };
    s3perm[0] = id; s3perm[1] = s; s3perm[2] = t; s3perm[3] = comp(s, t); s3perm[4] = comp(t, s); s3perm[5] = comp(s, comp(t, s));
    if (tables) {
        conjTab.assign(GL4.size(), vector<u16>(n));
#pragma omp parallel for
        for (int h = 0; h < (int)GL4.size(); h++) for (int v = 0; v < n; v++) conjTab[h][v] = (u16)vid[mm4(mm4(GL4[h], VN[v]), GL4inv[h])];
    }
}

// ---- clique enumeration (each clique once: increasing vertex index)
// at depth `need` remaining, call leaf(clique)
struct CliqueStats { long long nodes[16] = { 0 }; };
template <class Leaf> static void cliques(vector<int>& cur, const BS& cand, int need, CliqueStats& st, Leaf& leaf) {
    st.nodes[cur.size()]++;
    if (need == 0) { leaf(cur); return; }
    // quick bound: not enough candidates
    if (bcount(cand) < need) return;
    for (int wi = 0; wi < NW; wi++) {
        u64 w = cand.w[wi];
        while (w) {
            int b = (int)_tzcnt_u64(w); w &= w - 1; int v = wi * 64 + b;
            BS nc; for (int i = 0; i < NW; i++) nc.w[i] = 0;
            // keep only vertices > v adjacent to v
            nc.w[wi] = cand.w[wi] & adj[v].w[wi] & (b == 63 ? 0 : (~0ULL << (b + 1)));
            for (int i = wi + 1; i < NW; i++) nc.w[i] = cand.w[i] & adj[v].w[i];
            cur.push_back(v);
            cliques(cur, nc, need - 1, st, leaf);
            cur.pop_back();
        }
    }
}
static bool has_triangle(const BS& Q) {
    for (int wi = 0; wi < NW; wi++) { u64 w = Q.w[wi]; while (w) { int u = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
        BS a; band(a, Q, adj[u]);
        for (int wj = 0; wj < NW; wj++) { u64 x = a.w[wj]; while (x) { int t = wj * 64 + (int)_tzcnt_u64(x); x &= x - 1; BS b; band(b, a, adj[t]); if (!bempty(b)) return true; } }
    } }
    return false;
}

// Part A. For class c: r = classRep[c], P = adj(r) ∩ {class >= c}; enumerate (k-1)-cliques K in P, so K + r is a
// k-clique whose minimal class is c, with r moved to the class rep (conjugation). Weighted count:
// #k-cliques = sum_c |class c| * sum_K 1 / m(K + r), m = number of class-c elements.
static int mode_holes(int k, const char* outfile) {
    FILE* of = fopen(outfile, "a");
    printf("Part A: %d-cliques of V5824 (= partial spreads of %d solids with Y, X, Z among them)\n", k, k + 3);
    double total = 0; long long allLeaves = 0, allBad = 0;
    for (int c = 0; c < (int)classSize.size(); c++) {
        int r = classRep[c];
        BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
        for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
        // parallel over the first vertex
        vector<int> firsts; for (int v = 0; v < (int)VN.size(); v++) if (bget(P, v)) firsts.push_back(v);
        long long leaves = 0, bad = 0; double wsum = 0; CliqueStats tot; mutex mu;
        auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1) reduction(+:leaves, bad, wsum)
        for (int fi = 0; fi < (int)firsts.size(); fi++) {
            int v = firsts[fi]; BS nc; for (int i = 0; i < NW; i++) nc.w[i] = 0;
            int wi = v >> 6, b = v & 63;
            nc.w[wi] = P.w[wi] & adj[v].w[wi] & (b == 63 ? 0 : (~0ULL << (b + 1)));
            for (int i = wi + 1; i < NW; i++) nc.w[i] = P.w[i] & adj[v].w[i];
            vector<int> cur = { v }; CliqueStats st;
            auto leaf = [&](const vector<int>& K) {
                leaves++;
                int m = 1; for (int u : K) m += cls[u] == c;
                wsum += 1.0 / m;
                BS Q = adj[r]; for (int u : K) { BS t; band(t, Q, adj[u]); Q = t; }
                if (!has_triangle(Q)) {
                    bad++;
                    lock_guard<mutex> g(mu);
                    fprintf(of, "NONEXT c%d r%04x", c, VN[r]); for (int u : K) fprintf(of, " %04x", VN[u]); fprintf(of, " common %d\n", bcount(Q)); fflush(of);
                }
            };
            cliques(cur, nc, k - 2, st, leaf);
            lock_guard<mutex> g(mu); for (int i = 0; i < 16; i++) tot.nodes[i] += st.nodes[i];
        }
        double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        printf("class %d (size %d, rep %04x, |P| = %d): %lld leaves, %lld without extension, weighted %.3f, %.0f s\n  nodes by clique size (incl. rep):",
            c, classSize[c], VN[r], bcount(P), leaves, bad, wsum, el);
        for (int i = 0; i < k; i++) printf(" %lld", tot.nodes[i]); printf("\n"); fflush(stdout);
        fprintf(of, "CLASS %d k %d leaves %lld bad %lld weighted %.6f\n", c, k, leaves, bad, wsum); fflush(of);
        total += classSize[c] * wsum; allLeaves += leaves; allBad += bad;
    }
    printf("total %d-cliques (weighted count): %.3f; reduced leaves %lld; without extension: %lld\n", k, total, allLeaves, allBad);
    fprintf(of, "TOTAL k %d cliques %.3f leaves %lld bad %lld\n", k, total, allLeaves, allBad);
    fclose(of);
    return allBad ? 2 : 0;
}

// Part B: 14-cliques up to G.
static int mode_spreads(const char* outfile) {
    init_group(true);
    printf("Part B: 14-cliques of V5824 up to G = conj x S3 (order %zu)\n", GL4.size() * 6);
    // check S3 closure
    { set<vector<u16>> S; for (int a = 0; a < 6; a++) S.insert(s3perm[a]); bool closed = S.size() == 6;
      for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) { vector<u16> c(VN.size()); for (size_t v = 0; v < VN.size(); v++) c[v] = s3perm[a][s3perm[b][v]]; closed &= S.count(c) > 0; }
      printf("S3: 6 distinct vertex permutations, closed under composition: %s\n", closed ? "yes" : "NO"); if (!closed) return 1; }
    vector<vector<int>> red; double total = 0;
    for (int c = 0; c < (int)classSize.size(); c++) {
        int r = classRep[c];
        BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
        for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
        vector<int> cur = { r }; CliqueStats st; double wsum = 0; size_t before = red.size();
        auto leaf = [&](const vector<int>& K) { int m = 0; for (int u : K) m += cls[u] == c; wsum += 1.0 / m; red.push_back(K); };
        cliques(cur, P, 13, st, leaf);
        printf("class %d: %zu reduced 14-cliques, weighted %.4f\n", c, red.size() - before, wsum); fflush(stdout);
        total += classSize[c] * wsum;
    }
    printf("total 14-cliques by the weighted count: %.3f\n", total);
    { string fn = string(outfile) + ".all"; FILE* af = fopen(fn.c_str(), "w");
      for (size_t i = 0; i < red.size(); i++) { fprintf(af, "SPREAD %zu conj-reduced :", i); for (int v : red[i]) fprintf(af, " %04x", VN[v]); fprintf(af, "\n"); }
      fclose(af); printf("all %zu conjugation-reduced 14-cliques written to %s\n", red.size(), fn.c_str()); }
    // canonical forms
    int nr = (int)red.size();
    vector<vector<u16>> canon(nr); vector<int> stab(nr);
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < nr; i++) {
        vector<u16> best(14, 0xffff), img(14); int st = 0;
        vector<u16> self(red[i].begin(), red[i].end()); sort(self.begin(), self.end());
        for (int a = 0; a < 6; a++) for (size_t h = 0; h < GL4.size(); h++) {
            for (int j = 0; j < 14; j++) img[j] = conjTab[h][s3perm[a][red[i][j]]];
            sort(img.begin(), img.end());
            if (img < best) best = img;
            if (img == self) st++;
        }
        canon[i] = best; stab[i] = st;
    }
    map<vector<u16>, int> orbit; // canon -> first index
    for (int i = 0; i < nr; i++) if (!orbit.count(canon[i])) orbit[canon[i]] = i;
    long long sum = 0; int fields = 0;
    FILE* of = fopen(outfile, "w");
    int idx = 0;
    for (auto& kv : orbit) {
        int i = kv.second; long long osz = (long long)GL4.size() * 6 / stab[i]; sum += osz;
        // is {0, I} + clique closed under multiplication (a field)?
        set<u16> F = { 0, ID4 }; for (u16 v : kv.first) F.insert(VN[v]);
        bool fld = true; for (u16 a : F) for (u16 b : F) if (!F.count(mm4(a, b))) fld = false;
        bool add = true; for (u16 a : F) for (u16 b : F) if (!F.count(a ^ b)) add = false;
        fields += fld;
        fprintf(of, "SPREAD %d stab %d orbit %lld additive %d field %d :", idx, stab[i], osz, (int)add, (int)fld);
        for (u16 v : kv.first) fprintf(of, " %04x", VN[v]); fprintf(of, "\n");
        idx++;
    }
    fclose(of);
    printf("%zu G-orbits; sum of orbit sizes %lld (weighted count %.3f); %d orbit reps are fields\n", orbit.size(), sum, total, fields);
    return 0;
}

// ================================================================================================ Part C: forms
static u32 code28(u64 f) { u32 c = 0; int p = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { if (rowof(f, i) >> j & 1) c |= 1u << p; p++; } return c; }
static u64 E28[28];
static u64 TLO[1 << 14], THI[1 << 14];
static void init_code28() {
    int p = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { E28[p++] = (1ULL << (8 * i + j)) | (1ULL << (8 * j + i)); }
    for (u32 c = 0; c < (1u << 14); c++) { u64 a = 0, b = 0; for (int k = 0; k < 14; k++) if (c >> k & 1) { a ^= E28[k]; b ^= E28[14 + k]; } TLO[c] = a; THI[c] = b; }
}
static inline u64 form28(u32 c) { return TLO[c & 0x3fff] ^ THI[c >> 14]; }
static vector<u32> rref28(vector<u32> v) {
    vector<u32> out;
    for (int bit = 27; bit >= 0; bit--) {
        int p = -1; for (size_t i = 0; i < v.size(); i++) if (v[i] >> bit & 1) { p = (int)i; break; }
        if (p < 0) continue;
        u32 r = v[p]; v.erase(v.begin() + p);
        for (auto& x : v) if (x >> bit & 1) x ^= r;
        for (auto& x : out) if (x >> bit & 1) x ^= r;
        out.push_back(r);
    }
    return out;
}
static string space_key(const vector<u64>& basis) {
    vector<u32> c; for (u64 f : basis) c.push_back(code28(f));
    auto r = rref28(c); sort(r.begin(), r.end());
    string s; char buf[16]; for (u32 x : r) { sprintf(buf, "%07x ", x); s += buf; } return s;
}
static bool independent_of(const vector<u64>& basis, u64 f) {
    vector<u32> c; for (u64 b : basis) c.push_back(code28(b)); size_t r = rref28(c).size(); c.push_back(code28(f)); return rref28(c).size() > r;
}
// from scratch: LUT, derivative counts, Walsh spectrum
struct Verdict { bool apn; int dist[9]; };
static Verdict verify_space(const vector<u64>& basis) {
    Verdict V; memset(&V, 0, sizeof V);
    int n = (int)basis.size(); u8 F[256];
    for (int x = 0; x < 256; x++) {
        u8 y = 0;
        for (int i = 0; i < n; i++) { int s = 0; for (int j = 0; j < 8; j++) for (int k = j + 1; k < 8; k++) if ((x >> j & 1) && (x >> k & 1) && (rowof(basis[i], j) >> k & 1)) s ^= 1; y |= s << i; }
        F[x] = y;
    }
    V.apn = (n == 8);
    for (int a = 1; a < 256 && V.apn; a++) { int cnt[256] = { 0 }; for (int x = 0; x < 256; x++) if (++cnt[F[x ^ a] ^ F[x]] > 2) { V.apn = false; break; } }
    for (int b = 1; b < (1 << n); b++) {
        int W[256]; for (int x = 0; x < 256; x++) W[x] = par(b & F[x]) ? -1 : 1;
        for (int h = 1; h < 256; h <<= 1) for (int i = 0; i < 256; i += 2 * h) for (int j = i; j < i + h; j++) { int a = W[j], c = W[j + h]; W[j] = a + c; W[j + h] = a - c; }
        int mx = 0; for (int i = 0; i < 256; i++) mx = max(mx, abs(W[i]));
        int k = 0; while ((1 << (8 + k)) < mx * mx) k++;
        V.dist[k]++;
    }
    return V;
}
static u64 transform(u64 f, const u8 p[8]) { // B'(e_i, e_j) = B(p_i, p_j)
    u64 g = 0; for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (i != j && evalf(f, p[i], p[j])) g |= 1ULL << (8 * i + j); return g;
}
static void pick_independent(const vector<u8>& cand, int need, u8* out) {
    bool inspan[256] = { 0 }; inspan[0] = true; int k = 0;
    for (u8 u : cand) {
        if (k == need) break;
        if (inspan[u]) continue;
        out[k++] = u;
        bool nw[256]; memcpy(nw, inspan, 256);
        for (int v = 0; v < 256; v++) if (inspan[v]) nw[v ^ u] = true;
        memcpy(inspan, nw, 256);
    }
    if (k != need) { fprintf(stderr, "pick_independent failed\n"); exit(1); }
}
static void annihilator(const Mask& m, int need, u8* out) {
    vector<u8> ann; for (int u = 1; u < 256; u++) { bool ok = true; for (int v = 1; v < 256 && ok; v++) if (mget(m, v) && par(u & v)) ok = false; if (ok) ann.push_back((u8)u); }
    pick_independent(ann, need, out);
}

// The problem: the 14 rank-4 radicals (graphs of a 14-clique), holes D = X ∪ Y ∪ Z.
// Relaxations for controls: r6any (rank-6 radicals anywhere), oracle (outside-D rank-6 candidates from a list).
struct Problem {
    Mask R[14]; int radOf[256]; Mask D;
    vector<u64> f4[14];
    vector<Mask> lineM; vector<vector<u64>> lineF; vector<int> linesAt[256];
    vector<int> order;
    bool r6any = false;                 // rule: rank-6 radical anywhere (control)
    vector<u64> oracle; vector<Mask> oracleM;
    bool genAllLines = false;           // generator: rank-6 forms on every line through the point (control)
};
static Mask graph_mask(u16 N) { Mask m; mclear(m); for (u8 x = 1; x < 16; x++) mset(m, x | (mv4(N, x) << 4)); return m; }
static void build_problem(Problem& P, const vector<u16>& Ns) {
    mclear(P.D);
    for (int v = 1; v < 256; v++) { int x = v & 15, y = v >> 4; if (x == 0 || y == 0 || x == y) mset(P.D, v); }
    for (int v = 0; v < 256; v++) P.radOf[v] = -1;
    for (int j = 0; j < 14; j++) {
        P.R[j] = graph_mask(Ns[j]);
        for (int v = 1; v < 256; v++) if (mget(P.R[j], v)) { if (P.radOf[v] >= 0 || mget(P.D, v)) { fprintf(stderr, "radicals overlap\n"); exit(1); } P.radOf[v] = j; }
        u8 u[4]; for (int k = 0; k < 4; k++) u[k] = (u8)(((Ns[j] >> (4 * k)) & 15) | (1 << (4 + k)));
        P.f4[j].clear(); forms_from(u, 4, ND4, P.f4[j]);
    }
    // lines inside D
    set<vector<u64>> seen;
    P.lineM.clear(); P.lineF.clear(); for (int v = 0; v < 256; v++) P.linesAt[v].clear();
    for (int a = 1; a < 256; a++) for (int b = a + 1; b < 256; b++) if (mget(P.D, a) && mget(P.D, b) && mget(P.D, a ^ b)) {
        Mask m; mclear(m); mset(m, a); mset(m, b); mset(m, a ^ b);
        vector<u64> key(m.w, m.w + 4); if (seen.count(key)) continue; seen.insert(key);
        int li = (int)P.lineM.size(); P.lineM.push_back(m);
        u8 ann[6]; annihilator(m, 6, ann); vector<u64> fs; forms_from(ann, 6, ND6, fs); P.lineF.push_back(fs);
        for (int v : { a, b, a ^ b }) P.linesAt[v].push_back(li);
    }
}
static vector<Mask> allLines; static vector<array<u8, 6>> allLinesAnn; static vector<int> allLinesAt[256];
static void init_all_lines() {
    set<vector<u64>> seen;
    for (int a = 1; a < 256; a++) for (int b = a + 1; b < 256; b++) {
        Mask m; mclear(m); mset(m, a); mset(m, b); mset(m, a ^ b);
        vector<u64> key(m.w, m.w + 4); if (seen.count(key)) continue; seen.insert(key);
        array<u8, 6> ann; annihilator(m, 6, ann.data());
        int li = (int)allLines.size(); allLines.push_back(m); allLinesAnn.push_back(ann);
        for (int v : { a, b, a ^ b }) allLinesAt[v].push_back(li);
    }
}

static inline bool legal(const Problem& P, u64 f, Mask& m, int& cnt) {
    if (rank8(f) == 8) { cnt = 0; return true; }
    cnt = radical(f, m);
    if (cnt == 15) { int v = 0; for (v = 1; v < 256; v++) if (mget(m, v)) break; int j = P.radOf[v]; return j >= 0 && meq(m, P.R[j]); }
    if (cnt == 3) {
        if (P.r6any || msub(m, P.D)) return true;
        for (auto& o : P.oracleM) if (meq(o, m)) return true;
        return false;
    }
    return false;
}
struct State { vector<u64> span; Mask cov; int d; vector<u64> basis; };
static State empty_state() { State S; S.span = { 0 }; mclear(S.cov); S.d = 0; return S; }
static bool try_add(const Problem& P, const State& S, u64 f, State& out) {
    Mask nc = S.cov; size_t n = S.span.size();
    for (size_t i = 0; i < n; i++) {
        Mask m; int c;
        if (!legal(P, f ^ S.span[i], m, c)) return false;
        if (c) { if (!mdisj(m, nc)) return false; mor(nc, m); }
    }
    out.span.resize(2 * n);
    for (size_t i = 0; i < n; i++) { out.span[i] = S.span[i]; out.span[n + i] = S.span[i] ^ f; }
    out.cov = nc; out.d = S.d + 1; out.basis = S.basis; out.basis.push_back(f);
    return true;
}
template <class CB> static void for_each_candidate(const Problem& P, int a, const Mask& cov, CB cb) {
    int j = P.radOf[a];
    if (j >= 0) for (u64 f : P.f4[j]) cb(f);
    if (P.genAllLines) {
        vector<u64> buf;
        for (int li : allLinesAt[a]) if (mdisj(allLines[li], cov)) { buf.clear(); forms_from(allLinesAnn[li].data(), 6, ND6, buf); for (u64 f : buf) cb(f); }
        return;
    }
    if (j < 0) { for (int li : P.linesAt[a]) if (mdisj(P.lineM[li], cov)) for (u64 f : P.lineF[li]) cb(f); }
    // control only: rank-6 forms whose radical is not inside D, from a known function's list
    for (size_t i = 0; i < P.oracle.size(); i++) if (mget(P.oracleM[i], a) && !msub(P.oracleM[i], P.D) && mdisj(P.oracleM[i], cov)) cb(P.oracle[i]);
}
struct Stats { long long nodes[10] = { 0 }, tries[10] = { 0 }, bentext = 0, full = 0, p1[10] = { 0 }; };
static void add_stats(Stats& a, const Stats& b) { for (int i = 0; i < 10; i++) { a.nodes[i] += b.nodes[i]; a.tries[i] += b.tries[i]; a.p1[i] += b.p1[i]; } a.bentext += b.bentext; a.full += b.full; }
static mutex solmu;
static map<string, vector<u64>> solutions;
static map<string, vector<u64>> p1spaces; // spans in which every point off D is covered (the "short" analogue)
static bool g_logp1 = true;
static thread_local map<string, vector<u64>>* tl_sol = nullptr; // batch controls: per-thread solution sets
static void record(const State& S, Stats& st) {
    st.full++;
    string k = space_key(S.basis);
    if (tl_sol) { (*tl_sol)[k] = S.basis; return; }
    lock_guard<mutex> g(solmu);
    if (!solutions.count(k)) { solutions[k] = S.basis; printf("SOLUTION %s\n", k.c_str()); fflush(stdout); }
}
static void bent_extend(const State& S, Stats& st) {
    st.bentext++;
    vector<u32> c; for (u64 f : S.basis) c.push_back(code28(f));
    auto B = rref28(c);
    u32 piv = 0; for (u32 b : B) { unsigned long hb; _BitScanReverse(&hb, b); piv |= 1u << hb; }
    u32 freemask = ((1u << 28) - 1) & ~piv; int nf = __popcnt(freemask);
    for (u64 t = 1; t < (1ULL << nf); t++) {
        u64 f = form28(_pdep_u32((u32)t, freemask));
        bool ok = true; for (u64 s : S.span) if (rank8(f ^ s) != 8) { ok = false; break; }
        if (!ok) continue;
        State T; T.span.resize(2 * S.span.size());
        for (size_t i = 0; i < S.span.size(); i++) { T.span[i] = S.span[i]; T.span[S.span.size() + i] = S.span[i] ^ f; }
        T.cov = S.cov; T.d = S.d + 1; T.basis = S.basis; T.basis.push_back(f);
        if (T.d == 8) record(T, st); else bent_extend(T, st);
    }
}
static int first_uncovered(const Problem& P, const State& S) { for (int v : P.order) if (!mget(S.cov, v)) return v; return -1; }
static void note_node(const Problem& P, const State& S, int a, Stats& st) {
    st.nodes[S.d]++;
    if (a < 0 || P.radOf[a] < 0) { // every point off D is covered
        st.p1[S.d]++;
        if (g_logp1) { lock_guard<mutex> g(solmu); if (p1spaces.size() < 200000) p1spaces[space_key(S.basis)] = S.basis; }
    }
}
static void dfs(const Problem& P, const State& S, Stats& st) {
    int a = first_uncovered(P, S);
    note_node(P, S, a, st);
    if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, st); return; }
    if (S.d == 8) return;
    State T;
    for_each_candidate(P, a, S.cov, [&](u64 f) { st.tries[S.d]++; if (try_add(P, S, f, T)) dfs(P, T, st); });
}
static void expand(const Problem& P, const State& S, Stats& st, vector<State>& out) {
    int a = first_uncovered(P, S);
    note_node(P, S, a, st);
    if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, st); return; }
    if (S.d == 8) return;
    State T;
    for_each_candidate(P, a, S.cov, [&](u64 f) { st.tries[S.d]++; if (try_add(P, S, f, T)) out.push_back(T); });
}
static Stats run_search(const Problem& P, const State& S0, int mintasks) {
    Stats tot; vector<State> tasks = { S0 };
    while (tasks.size() && (int)tasks.size() < mintasks) {
        vector<State> nxt; for (auto& s : tasks) expand(P, s, tot, nxt);
        tasks.swap(nxt);
        if (tasks.size() && tasks[0].d >= 7) break;
    }
    int nt = (int)tasks.size();
    vector<Stats> per(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nt; i++) dfs(P, tasks[i], per[omp_get_thread_num()]);
    for (auto& p : per) add_stats(tot, p);
    return tot;
}
static void make_order(Problem& P, int how, unsigned seed) {
    P.order.clear();
    vector<int> off, on; for (int v = 1; v < 256; v++) (P.radOf[v] >= 0 ? off : on).push_back(v);
    if (how == 1) { reverse(off.begin(), off.end()); reverse(on.begin(), on.end()); }
    if (how == 2) { mt19937 g(seed); shuffle(off.begin(), off.end(), g); shuffle(on.begin(), on.end(), g); }
    for (int v : off) P.order.push_back(v); for (int v : on) P.order.push_back(v);
}
static vector<pair<int, vector<u16>>> read_spreads(const char* file) {
    vector<pair<int, vector<u16>>> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "SPREAD", 6)) continue;
        int idx = atoi(line + 7); const char* p = strchr(line, ':'); if (!p) continue; p++;
        vector<u16> Ns; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
        if (Ns.size() != 14) { printf("bad spread line\n"); exit(1); }
        out.push_back({ idx, Ns });
    }
    fclose(fp); return out;
}
static void print_stats(const char* tag, const Stats& s) {
    printf("%s nodes by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.nodes[i]);
    printf("\n%s tries by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.tries[i]);
    printf("\n%s nodes with every off-D point covered, by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.p1[i]);
    printf("\n%s bent-extension nodes %lld, full-span records %lld\n", tag, s.bentext, s.full); fflush(stdout);
}
static void write_stats(FILE* pf, int j, const Stats& s) {
    fprintf(pf, "R %d", j); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.nodes[k]); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.tries[k]);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.p1[k]); fprintf(pf, " %lld %lld\n", s.bentext, s.full); fflush(pf);
}
static bool read_stats(const char* line, int& id, Stats& s) {
    const char* p = line; if (*p != 'R') return false; p++;
    char* e; id = (int)strtol(p, &e, 10); if (e == p) return false; p = e;
    long long v[29]; for (int i = 0; i < 29; i++) { v[i] = strtoll(p, &e, 10); if (e == p) return false; p = e; }
    for (int k = 0; k < 9; k++) { s.nodes[k] = v[k]; s.tries[k] = v[9 + k]; s.p1[k] = v[18 + k]; } s.bentext = v[27]; s.full = v[28]; return true;
}
// search over every spread rep (or a range); progress file for resuming
static int mode_search(const char* spreadfile, const char* progfile, int how, unsigned seed, int first, int last) {
    auto sp = read_spreads(spreadfile);
    printf("search: %zu spread configurations; point order %s\n", sp.size(), how == 0 ? "ascending" : how == 1 ? "descending" : "random");
    set<int> done; Stats tot;
    if (FILE* fp = fopen(progfile, "r")) {
        char line[1024];
        while (fgets(line, sizeof line, fp)) { int id; Stats s; if (read_stats(line, id, s) && !done.count(id)) { done.insert(id); add_stats(tot, s); } else if (!strncmp(line, "SOLUTION", 8)) printf("(previous run) %s", line); }
        fclose(fp);
    }
    FILE* pf = fopen(progfile, "a");
    auto t0 = chrono::steady_clock::now();
    for (auto& [idx, Ns] : sp) {
        if (idx < first || idx > last || done.count(idx)) continue;
        Problem P; build_problem(P, Ns); make_order(P, how, seed + idx);
        if (P.lineM.size() != 120) { printf("expected 120 lines in D, got %zu\n", P.lineM.size()); return 1; }
        Stats st = run_search(P, empty_state(), 400);
        { lock_guard<mutex> g(solmu); for (auto& kv : solutions) fprintf(pf, "SOLUTION %d %s\n", idx, kv.first.c_str()); }
        write_stats(pf, idx, st);
        add_stats(tot, st);
        printf("spread %d: nodes", idx); for (int k = 0; k <= 8; k++) printf(" %lld", st.nodes[k]);
        printf("  off-D-covered"); for (int k = 0; k <= 8; k++) printf(" %lld", st.p1[k]);
        printf("  full %lld bentext %lld  [%.0f s]\n", st.full, st.bentext, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    }
    // the off-D-covered spans, for comparison with Rowen's SHORT lines
    for (auto& kv : p1spaces) { fprintf(pf, "P1 dim %zu basis", kv.second.size()); for (u64 f : kv.second) fprintf(pf, " %07x", code28(f)); fprintf(pf, "\n"); }
    fclose(pf);
    print_stats("search", tot);
    printf("distinct off-D-covered spans this run: %zu\n", p1spaces.size());
    printf("solutions: %zu\n", solutions.size());
    for (auto& kv : solutions) { Verdict V = verify_space(kv.second); printf("  %s APN=%d dist:", kv.first.c_str(), V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("\n"); }
    return 0;
}

// ================================================================================ positive control: known functions
static bool read_lut(const string& line, u8 F[256]) {
    int n = 0; const char* p = line.c_str();
    while (*p && n < 256) { if (*p >= '0' && *p <= '9') { F[n++] = (u8)strtol(p, (char**)&p, 10); } else p++; }
    return n == 256;
}
static vector<u64> forms_of_lut(const u8 F[256]) {
    vector<u64> B(8, 0);
    for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) for (int k = 0; k < 8; k++) if (j != k) {
        u8 d = F[(1 << j) ^ (1 << k)] ^ F[1 << j] ^ F[1 << k] ^ F[0];
        if (d >> i & 1) B[i] |= 1ULL << (8 * j + k);
    }
    return B;
}
struct Profile { int n4 = 0, n6 = 0, n8 = 0, bad = 0; vector<u64> f4, f6; vector<Mask> r4, r6; };
static Profile profile(const vector<u64>& B) {
    Profile P;
    for (int c = 1; c < (1 << B.size()); c++) {
        u64 f = 0; for (size_t i = 0; i < B.size(); i++) if (c >> i & 1) f ^= B[i];
        int rk = rank8(f); Mask m;
        if (rk == 8) P.n8++;
        else if (rk == 6) { P.n6++; radical(f, m); P.f6.push_back(f); P.r6.push_back(m); }
        else if (rk == 4) { P.n4++; radical(f, m); P.f4.push_back(f); P.r4.push_back(m); }
        else P.bad++;
    }
    return P;
}
static bool load_line(const char* file, int line, vector<u64>& B) {
    FILE* fp = fopen(file, "r"); if (!fp) return false;
    static char buf[1 << 16]; int ln = 0;
    while (fgets(buf, sizeof buf, fp)) { if (++ln == line) { u8 F[256]; fclose(fp); if (!read_lut(buf, F)) return false; B = forms_of_lut(F); return true; } }
    fclose(fp); return false;
}
// linear maps of V by column images: (A v) = xor of col[i] over the bits of v
typedef array<u8, 8> Lin;
static inline u8 lin_apply(const Lin& A, u8 v) { u8 r = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) r ^= A[i]; return r; }
static Lin lin_comp(const Lin& A, const Lin& B) { Lin C; for (int i = 0; i < 8; i++) C[i] = lin_apply(A, B[i]); return C; } // A after B
static Lin lin_inv(const Lin& A) { Lin C; for (int i = 0; i < 8; i++) { int f = -1; for (int v = 0; v < 256; v++) if (lin_apply(A, (u8)v) == (1 << i)) { f = v; break; } if (f < 0) { fprintf(stderr, "singular\n"); exit(1); } C[i] = (u8)f; } return C; }
static Mask mask_image(const Lin& A, const Mask& m) { Mask r; mclear(r); for (int v = 1; v < 256; v++) if (mget(m, v)) mset(r, lin_apply(A, (u8)v)); return r; }
// new coordinates w = A v: forms become B o A^-1
static void change_coords(vector<u64>& T, const Lin& A) { Lin Ai = lin_inv(A); for (auto& f : T) f = transform(f, Ai.data()); }
static void basis_of_mask(const Mask& m, int dim, u8* out) { vector<u8> c; for (int v = 1; v < 256; v++) if (mget(m, v)) c.push_back((u8)v); pick_independent(c, dim, out); }
// N with graph(N) = m (m complementary to Y), or -1
static int matrix_of_graph(const Mask& m) {
    u8 col[4]; for (int i = 0; i < 4; i++) { col[i] = 0xff; for (int v = 1; v < 256; v++) if (mget(m, v) && (v & 15) == (1 << i)) col[i] = (u8)(v >> 4); if (col[i] == 0xff) return -1; }
    u16 N = fromcols4(col); return meq(graph_mask(N), m) ? N : -1;
}
static bool find_clique(vector<int>& cur, const BS& cand, int need) {
    if (need == 0) return true;
    for (int wi = 0; wi < NW; wi++) { u64 w = cand.w[wi]; while (w) { int v = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
        BS nc; band(nc, cand, adj[v]); for (int i = 0; i <= wi; i++) nc.w[i] &= (i < wi) ? 0 : ((v & 63) == 63 ? 0 : (~0ULL << ((v & 63) + 1)));
        cur.push_back(v); if (find_clique(cur, nc, need - 1)) return true; cur.pop_back(); } }
    return false;
}
// a random spread through given elements: randomised order at each level, first success
static bool find_clique_rand(vector<int>& cur, const BS& cand, int need, mt19937& rng) {
    if (need == 0) return true;
    vector<int> vs; for (int wi = 0; wi < NW; wi++) { u64 w = cand.w[wi]; while (w) { vs.push_back(wi * 64 + (int)_tzcnt_u64(w)); w &= w - 1; } }
    if ((int)vs.size() < need) return false;
    shuffle(vs.begin(), vs.end(), rng);
    int tries = 0;
    for (int v : vs) { BS nc; band(nc, cand, adj[v]); cur.push_back(v); if (find_clique_rand(cur, nc, need - 1, rng)) return true; cur.pop_back(); if (++tries > 64) break; }
    return false;
}
// Positive control: a known function with >= 2 rank-4 components. Its rank-4 radicals are put into a spread, three
// other spread elements are made the holes X, Y, Z, and the configuration is moved onto its orbit representative in
// spreads.txt (searching G for it). Then the main search on that representative, with one relaxation: rank-6
// radicals outside D are allowed if they are the known function's own (those forms are offered at their points;
// everything else is generated exactly as in the main search).
struct CtlCtx {
    vector<pair<int, vector<u16>>> sp; map<vector<u16>, int> repOf; vector<Lin> S3; vector<Problem> probs;
};
static CtlCtx ctl_ctx(const char* spreadfile) {
    CtlCtx X; X.sp = read_spreads(spreadfile);
    for (auto& [idx, Ns] : X.sp) { vector<u16> s = Ns; sort(s.begin(), s.end()); X.repOf[s] = idx; }
    Lin sw, ta, id; for (int i = 0; i < 4; i++) { sw[i] = (u8)(1 << (4 + i)); sw[4 + i] = (u8)(1 << i); ta[i] = (u8)((1 << i) | (1 << (4 + i))); ta[4 + i] = (u8)(1 << (4 + i)); id[i] = (u8)(1 << i); id[4 + i] = (u8)(1 << (4 + i)); }
    set<Lin> S3 = { id }; bool grew = true; while (grew) { grew = false; for (Lin a : vector<Lin>(S3.begin(), S3.end())) for (Lin g : { sw, ta }) if (S3.insert(lin_comp(g, a)).second) grew = true; }
    if (S3.size() != 6) { printf("S3 has %zu elements\n", S3.size()); exit(1); }
    X.S3.assign(S3.begin(), S3.end());
    X.probs.resize(X.sp.size()); for (size_t i = 0; i < X.sp.size(); i++) { if (X.sp[i].first != (int)i) { printf("spread file not indexed 0..n-1\n"); exit(1); } build_problem(X.probs[i], X.sp[i].second); }
    return X;
}
struct CtlResult { int status = 1; int rep = -1; int n4 = 0, inD = 0, nsol = 0, napn = 0; bool known = false; string why; Stats st; };
static CtlResult ctl_one(const CtlCtx& X, vector<u64> T, int choice, int how, int seedDim, bool verbose) {
    CtlResult R; Profile Q = profile(T); R.n4 = Q.n4;
    if (Q.n4 < 2 || Q.bad) { R.why = "fewer than 2 rank-4"; return R; }
    if (!verify_space(T).apn) { R.why = "forms not an APN space (line not quadratic APN)"; return R; }
    Mask Xm, Ym; mclear(Xm); mclear(Ym); for (int v = 1; v < 256; v++) { if (!(v >> 4)) mset(Xm, v); if (!(v & 15)) mset(Ym, v); }
    // Q1 -> X, Q2 -> Y
    Lin Ai; basis_of_mask(Q.r4[0], 4, Ai.data()); basis_of_mask(Q.r4[1], 4, Ai.data() + 4);
    change_coords(T, lin_inv(Ai)); Q = profile(T);
    // a third radical (if any) -> graph(I) by (x, y) -> (M x, y)
    vector<int> rest; for (size_t i = 0; i < Q.r4.size(); i++) if (!meq(Q.r4[i], Xm) && !meq(Q.r4[i], Ym)) rest.push_back((int)i);
    if (!rest.empty()) {
        int M = matrix_of_graph(Q.r4[rest[0]]); if (M < 0) { R.why = "not a graph"; return R; }
        Lin B; for (int i = 0; i < 4; i++) { B[i] = mv4((u16)M, (u8)(1 << i)); B[4 + i] = (u8)(1 << (4 + i)); }
        change_coords(T, B); Q = profile(T);
    }
    vector<int> qN; bool hasI = false;
    for (auto& m : Q.r4) { if (meq(m, Xm) || meq(m, Ym)) continue; int N = matrix_of_graph(m); if (N < 0 || (N != ID4 && vid[N] < 0)) { R.why = "radical not in V5824"; return R; } if (N == ID4) hasI = true; else qN.push_back(N); }
    if (!rest.empty() && !hasI) { R.why = "normalisation to graph(I) failed"; return R; }
    // a spread set {0, I} + 14-clique containing the qN
    BS cand; for (int i = 0; i < NW; i++) cand.w[i] = ~0ULL;
    vector<int> cur; for (int N : qN) { cur.push_back(vid[N]); BS t; band(t, cand, adj[vid[N]]); cand = t; }
    bool gotSpread = false;
    if (choice >= 1000) { // randomised spread through the radicals (to reach more representatives)
        mt19937 rng((unsigned)choice); vector<int> c0 = cur;
        for (int t = 0; t < 50 && !gotSpread; t++) { cur = c0; gotSpread = find_clique_rand(cur, cand, 14 - (int)cur.size(), rng); }
        if (!gotSpread) cur = c0;
    }
    if (!gotSpread && !find_clique(cur, cand, 14 - (int)cur.size())) { R.why = "no spread through the rank-4 radicals"; return R; }
    vector<Mask> others, spreadAll = { Xm, Ym, graph_mask(ID4) };
    for (int v : cur) spreadAll.push_back(graph_mask(VN[v]));
    for (auto& m : spreadAll) { bool isQ = false; for (auto& r : Q.r4) if (meq(r, m)) isQ = true; if (!isQ) others.push_back(m); }
    vector<array<int, 3>> triples; for (size_t a = 0; a < others.size(); a++) for (size_t b = a + 1; b < others.size(); b++) for (size_t c = b + 1; c < others.size(); c++) triples.push_back({ (int)a, (int)b, (int)c });
    auto tr = triples[choice % triples.size()];
    const Mask& H1 = others[tr[0]]; const Mask& H2 = others[tr[1]]; const Mask& H3 = others[tr[2]];
    // H1 -> X, H2 -> Y, H3 -> Z: basis a_i of H1, b_i = the H2-part of a_i in H2 + H3
    Lin Ci; basis_of_mask(H1, 4, Ci.data());
    for (int i = 0; i < 4; i++) { int bi = -1; for (int b = 1; b < 256; b++) if (mget(H2, b) && mget(H3, Ci[i] ^ b)) bi = b; if (bi < 0) { R.why = "holes not complementary"; return R; } Ci[4 + i] = (u8)bi; }
    Lin C = lin_inv(Ci); change_coords(T, C);
    vector<Mask> rad14; for (auto& m : spreadAll) if (!meq(m, H1) && !meq(m, H2) && !meq(m, H3)) rad14.push_back(mask_image(C, m));
    // onto an orbit representative: diag(h, h) o L over L in S3, h in GL(4,2)
    int found = -1; Lin Mfound{};
    for (const Lin& L : X.S3) {
        vector<u16> NL; bool ok = true;
        for (auto& m : rad14) { int N = matrix_of_graph(mask_image(L, m)); if (N < 0 || vid[N] < 0) { ok = false; break; } NL.push_back((u16)N); }
        if (!ok) { R.why = "S3 image not in V5824"; return R; }
        for (size_t hi = 0; hi < GL4.size() && found < 0; hi++) {
            u16 h = GL4[hi], hinv = GL4inv[hi]; vector<u16> img(14);
            for (int j = 0; j < 14; j++) img[j] = mm4(mm4(h, NL[j]), hinv);
            sort(img.begin(), img.end());
            auto it = X.repOf.find(img);
            if (it != X.repOf.end()) { found = it->second; Lin Dh; for (int i = 0; i < 4; i++) { u8 he = mv4(h, (u8)(1 << i)); Dh[i] = he; Dh[4 + i] = (u8)(he << 4); } Mfound = lin_comp(Dh, L); }
        }
        if (found >= 0) break;
    }
    if (found < 0) { R.status = 2; R.why = "configuration matches no orbit representative"; return R; }
    change_coords(T, Mfound); R.rep = found;
    string key = space_key(T);
    Problem P = X.probs[found]; make_order(P, how, 1);
    Q = profile(T);
    for (auto& m : Q.r4) { bool in = false; for (int j = 0; j < 14; j++) if (meq(m, P.R[j])) in = true; if (!in) { R.status = 2; R.why = "rank-4 radical not in the rep's spread"; return R; } }
    for (size_t i = 0; i < Q.f6.size(); i++) { if (msub(Q.r6[i], P.D)) R.inD++; else { P.oracle.push_back(Q.f6[i]); P.oracleM.push_back(Q.r6[i]); } }
    State S = empty_state(), S2;
    for (u64 f : Q.f4) { if (S.d >= seedDim) break; if (!independent_of(S.basis, f)) continue; if (!try_add(P, S, f, S2)) { R.status = 2; R.why = "seed fail"; return R; } S = S2; }
    for (u64 f : Q.f6) { if (S.d >= seedDim) break; if (!independent_of(S.basis, f)) continue; if (!try_add(P, S, f, S2)) { R.status = 2; R.why = "seed fail"; return R; } S = S2; }
    map<string, vector<u64>> sol; tl_sol = &sol;
    R.st = run_search(P, S, verbose ? 400 : 1);
    tl_sol = nullptr;
    R.known = sol.count(key) > 0; R.nsol = (int)sol.size();
    for (auto& kv : sol) { Verdict V = verify_space(kv.second); R.napn += V.apn;
        if (verbose) { printf("  APN=%d dist:", V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("%s\n", kv.first == key ? "  <- known" : ""); } }
    R.status = R.known ? 0 : 3;
    if (verbose) printf("hole triple %d of %zu -> spread rep %d; rank-6 radicals %d inside D, %zu oracle; found %d spaces (%d APN), known among them: %s\n",
        choice % (int)triples.size(), triples.size(), found, R.inD, P.oracle.size(), R.nsol, R.napn, R.known ? "YES" : "NO");
    return R;
}
static int mode_ctl(const char* spreadfile, const char* file, int line, int choice, int seedDim, int how) {
    CtlCtx X = ctl_ctx(spreadfile);
    vector<u64> T; if (!load_line(file, line, T)) { printf("cannot load line %d\n", line); return 1; }
    g_logp1 = false;
    CtlResult R = ctl_one(X, T, choice, how, seedDim, true);
    if (R.status == 1 || R.status == 2) printf("ctl: %s\n", R.why.c_str());
    print_stats("ctl", R.st);
    return R.status;
}
// every function in a file with >= minr4 rank-4 components (up to maxn of them), `nch` hole choices each
static bool g_randSpread = false;
static int mode_ctlbatch(const char* spreadfile, const char* file, int minr4, long long maxn, int nch, const char* outfile) {
    CtlCtx X = ctl_ctx(spreadfile);
    g_logp1 = false;
    FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); return 1; }
    FILE* of = fopen(outfile, "a");
    vector<pair<long long, vector<u64>>> batch; long long ln = 0, taken = 0;
    long long tot = 0, ok = 0, notfound = 0, nospread = 0, fail = 0, extra = 0, nonapn = 0; map<int, int> repHits; map<int, int> byN4;
    static char buf[1 << 16];
    auto flush = [&]() {
        int n = (int)batch.size(); vector<vector<CtlResult>> res(n, vector<CtlResult>(nch));
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n; i++) for (int c = 0; c < nch; c++) res[i][c] = ctl_one(X, batch[i].second, g_randSpread ? (int)(1000 + batch[i].first * 64 + c) : (int)(batch[i].first * 7 + c * 131), 0, 0, false);
        for (int i = 0; i < n; i++) for (int c = 0; c < nch; c++) {
            auto& R = res[i][c]; tot++;
            if (R.status == 0) { ok++; repHits[R.rep]++; byN4[R.n4]++; extra += R.nsol - 1; nonapn += R.nsol - R.napn; }
            else if (R.why == "no spread through the rank-4 radicals") nospread++;
            else if (R.status == 3) { notfound++; fprintf(of, "NOTFOUND line %lld choice %d rep %d\n", batch[i].first, c, R.rep); }
            else { fail++; fprintf(of, "FAIL line %lld choice %d: %s\n", batch[i].first, c, R.why.c_str()); }
            if (R.status == 0 && R.nsol > 1) fprintf(of, "EXTRA line %lld choice %d rep %d: %d spaces, %d APN\n", batch[i].first, c, R.rep, R.nsol, R.napn);
        }
        fflush(of); batch.clear();
        printf("  %lld lines read, %lld taken; controls %lld: re-found %lld, NOT found %lld, no spread %lld, setup failures %lld\n", ln, taken, tot, ok, notfound, nospread, fail); fflush(stdout);
    };
    while (fgets(buf, sizeof buf, fp)) {
        ln++; u8 F[256]; if (!read_lut(buf, F)) continue;
        vector<u64> B = forms_of_lut(F);
        // quick rank-4 count
        int n4 = 0; for (int c = 1; c < 256; c++) { u64 f = 0; for (int i = 0; i < 8; i++) if (c >> i & 1) f ^= B[i]; if (rank8(f) == 4) n4++; }
        if (n4 < minr4) continue;
        batch.push_back({ ln, B }); taken++;
        if ((int)batch.size() >= 256) flush();
        if (taken >= maxn) break;
    }
    if (!batch.empty()) flush();
    fclose(fp);
    printf("ctlbatch %s: %lld controls; re-found %lld; NOT found %lld; no spread through the radicals %lld; setup failures %lld; other spaces found %lld (%lld not APN)\n", file, tot, ok, notfound, nospread, fail, extra, nonapn);
    printf("  re-found by rank-4 count:"); for (auto& kv : byN4) printf(" %d:%d", kv.first, kv.second); printf("\n");
    printf("  spread reps hit: %zu of %zu\n", repHits.size(), X.sp.size());
    fprintf(of, "SUMMARY %s controls %lld refound %lld notfound %lld nospread %lld fail %lld extra %lld nonapn %lld reps %zu\n", file, tot, ok, notfound, nospread, fail, extra, nonapn, repHits.size());
    fprintf(of, "REPS"); for (auto& kv : repHits) fprintf(of, " %d:%d", kv.first, kv.second); fprintf(of, "\n");
    fclose(of);
    return (notfound || fail) ? 3 : 0;
}

// generators and problem tables, for every configuration in a spread file
static int mode_selftest(const char* spreadfile) {
    auto sp = read_spreads(spreadfile); int bad = 0;
    printf("selftest: ND4 %zu (28), ND6 %zu (13888), GL4 %zu (20160), V5824 %zu\n", ND4.size(), ND6.size(), GL4.size(), VN.size());
    bad += ND4.size() != 28 || ND6.size() != 13888 || GL4.size() != 20160 || VN.size() != 5824;
    long long checked = 0;
    for (auto& [idx, Ns] : sp) {
        Problem P; build_problem(P, Ns);
        bool ok = mcount(P.D) == 45 && P.lineM.size() == 120;
        Mask all = P.D; for (int j = 0; j < 14; j++) { ok &= mdisj(all, P.R[j]) && mcount(P.R[j]) == 15; mor(all, P.R[j]); }
        ok &= mcount(all) == 255;
        // pairwise adjacency of the clique (so the 14 graphs + X, Y, Z are pairwise complementary)
        for (int a = 0; a < 14; a++) for (int b = a + 1; b < 14; b++) ok &= inv4(Ns[a] ^ Ns[b]) && inv4(Ns[a]) && inv4(Ns[a] ^ ID4);
        int inSolid = 0, trans = 0;
        for (auto& m : P.lineM) {
            int cx = 0, cy = 0, cz = 0; for (int v = 1; v < 256; v++) if (mget(m, v)) { int x = v & 15, y = v >> 4; if (!y) cx++; else if (!x) cy++; else if (x == y) cz++; }
            if (cx == 3 || cy == 3 || cz == 3) inSolid++; else if (cx == 1 && cy == 1 && cz == 1) trans++;
        }
        ok &= inSolid == 105 && trans == 15;
        for (int v = 1; v < 256; v++) if (mget(P.D, v)) ok &= P.linesAt[v].size() == 8;
        for (int j = 0; j < 14; j++) { set<u64> s(P.f4[j].begin(), P.f4[j].end()); ok &= s.size() == 28; for (u64 f : P.f4[j]) { Mask m; ok &= radical(f, m) == 15 && meq(m, P.R[j]); checked++; } }
        if (idx % 8 == 0) for (size_t li = 0; li < P.lineM.size(); li++) {
            set<u64> s(P.lineF[li].begin(), P.lineF[li].end()); ok &= s.size() == 13888;
            for (u64 f : P.lineF[li]) { Mask m; ok &= rank8(f) == 6 && radical(f, m) == 3 && meq(m, P.lineM[li]); checked++; }
        }
        if (!ok) { printf("spread %d: FAIL\n", idx); bad++; }
    }
    // brute force: forms with radical exactly graph(N) for the first rep's N_0, among all 2^28
    { Problem P; build_problem(P, sp[0].second); long long cnt = 0;
#pragma omp parallel for reduction(+:cnt) schedule(static, 1 << 16)
      for (int c = 0; c < (1 << 28); c++) { u64 f = form28((u32)c); if (rank8(f) == 4) { Mask m; radical(f, m); if (meq(m, P.R[0])) cnt++; } }
      printf("brute force: %lld forms with radical exactly R_0 of spread 0 (expect 28)\n", cnt); bad += cnt != 28; }
    printf("selftest: %zu configurations, %lld generated forms checked, %d failures\n", sp.size(), checked, bad);
    return bad;
}

int main(int argc, char** argv) {
    if (argc < 2) { printf("modes: holes <k> <out> | spreads <out> | search <spreads> <prog> [asc|desc|rand <seed>] [first last] | ctl <spreads> <lutfile> <line> <choice> <seeddim> [desc]\n"); return 1; }
    string mode = argv[1];
    ND4 = nondeg_codes(4); ND6 = nondeg_codes(6);
    init_V();
    printf("V5824: %zu vertices, %zu classes:", VN.size(), classSize.size()); for (size_t c = 0; c < classSize.size(); c++) printf(" %d (rep %04x)", classSize[c], VN[classRep[c]]); printf("\n");
    { long long e = 0; for (auto& b : adj) e += bcount(b); printf("edges %lld\n", e / 2); }
    fflush(stdout);
    if (mode == "holes") return mode_holes(atoi(argv[2]), argv[3]);
    if (mode == "spreads") return mode_spreads(argv[2]);
    init_code28();
    if (mode == "search") {
        int how = 0; unsigned seed = 1; int first = 0, last = INT_MAX, ai = 4;
        if (argc > ai && !strcmp(argv[ai], "desc")) { how = 1; ai++; }
        else if (argc > ai && !strcmp(argv[ai], "rand")) { how = 2; seed = (unsigned)atoi(argv[ai + 1]); ai += 2; }
        else if (argc > ai && !strcmp(argv[ai], "asc")) ai++;
        if (argc > ai + 1) { first = atoi(argv[ai]); last = atoi(argv[ai + 1]); }
        return mode_search(argv[2], argv[3], how, seed, first, last);
    }
    if (mode == "ctl") return mode_ctl(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), (argc > 7 && !strcmp(argv[7], "desc")) ? 1 : 0);
    if (mode == "selftest") return mode_selftest(argv[2]);
    if (mode == "ctlbatch" && argc > 8 && !strcmp(argv[8], "randspread")) g_randSpread = true;
    if (mode == "ctlbatch") return mode_ctlbatch(argv[2], argv[3], atoi(argv[4]), atoll(argv[5]), atoi(argv[6]), argv[7]);
    return 1;
}
