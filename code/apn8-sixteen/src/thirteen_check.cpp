// thirteen_check.cpp -- independent check of Open Problem VI.1 (arXiv 2510.12008) for i = 4: [0^222, 2^20, 4^13].
// Subagent (Claude Opus 5.5) for Rowen, 2026-10-24. Written without reading cover.cpp, cover_prep.py, i4_holes.py,
// README "Thirteen" or runs/i4/. Utilities (forms, radicals, generators, V5824, clique enumeration, LUT/Walsh
// verifier, coordinate changes) copied from fourteen_check.cpp (the i = 3 checker's).
//
// Coordinates: v = (x, y), x = low nibble, y = high nibble. X = {y = 0}, Y = {x = 0}, Z = graph(I) = {y = x}.
// A 4-space complementary to Y is graph(N) = {(x, Nx)}; complementary to X iff N invertible, to Z iff N + I is.
// V5824 = {N : N, N + I invertible}; two graphs are disjoint iff N + N' is invertible ("adjacent").
//
// holes   : every 10-clique of V5824 (= 13-partial spread with X, Y, Z among its members), up to conjugation by the
//           minimal-class trick; extendable (4-clique in the common neighbourhood) or not; hole invariants; counts.
// configs : 13-partial spreads up to GL(8,2) by a canonical form over all frames (ordered member triple -> X, Y, Z,
//           then GL(4,2)-conjugation). Sources: own 14-cliques x 14 (extendable), the holes run's non-extendable
//           cliques, and exact covers of the complements of Kurz's three printed hole sets. Count identities.
// search  : per configuration, point-branching over the forms (radicals fixed).
// ctl / ctlbatch : known functions through the same search path.
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

struct Mask { u64 w[4]; bool operator<(const Mask& o) const { return memcmp(w, o.w, sizeof w) < 0; } bool operator==(const Mask& o) const { return !memcmp(w, o.w, sizeof w); } };
static inline void mclear(Mask& m) { m.w[0] = m.w[1] = m.w[2] = m.w[3] = 0; }
static inline void mset(Mask& m, int v) { m.w[v >> 6] |= 1ULL << (v & 63); }
static inline bool mget(const Mask& m, int v) { return m.w[v >> 6] >> (v & 63) & 1; }
static inline bool mdisj(const Mask& a, const Mask& b) { return ((a.w[0] & b.w[0]) | (a.w[1] & b.w[1]) | (a.w[2] & b.w[2]) | (a.w[3] & b.w[3])) == 0; }
static inline bool msub(const Mask& a, const Mask& b) { return ((a.w[0] & ~b.w[0]) | (a.w[1] & ~b.w[1]) | (a.w[2] & ~b.w[2]) | (a.w[3] & ~b.w[3])) == 0; }
static inline void mor(Mask& a, const Mask& b) { for (int i = 0; i < 4; i++) a.w[i] |= b.w[i]; }
static inline int mcount(const Mask& m) { return (int)(__popcnt64(m.w[0]) + __popcnt64(m.w[1]) + __popcnt64(m.w[2]) + __popcnt64(m.w[3])); }
static inline bool meq(const Mask& a, const Mask& b) { return !memcmp(&a, &b, sizeof(Mask)); }
static Mask mcompl(const Mask& m) { Mask r; for (int i = 0; i < 4; i++) r.w[i] = ~m.w[i]; r.w[0] &= ~1ULL; return r; } // nonzero points not in m
static vector<u8> mpoints(const Mask& m) { vector<u8> p; for (int v = 1; v < 256; v++) if (mget(m, v)) p.push_back((u8)v); return p; }

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
    memset(gl4idx, -1, sizeof gl4idx);
    for (u32 N = 0; N < 65536; N++) if (inv4((u16)N)) GL4.push_back((u16)N);
    GL4inv.resize(GL4.size()); for (size_t i = 0; i < GL4.size(); i++) gl4idx[GL4[i]] = (int)i;
    for (size_t i = 0; i < GL4.size(); i++) { u16 g = GL4[i], p = g, prev = ID4; while (p != ID4) { prev = p; p = mm4(p, g); } GL4inv[i] = (g == ID4) ? ID4 : prev; if (mm4(GL4inv[i], g) != ID4) { fprintf(stderr, "inverse\n"); exit(1); } }
    vector<u16> raw; for (u16 N : GL4) if (inv4(N ^ ID4)) raw.push_back(N);
    unordered_map<u16, int> c0; int nc = 0;
    for (u16 N : raw) if (!c0.count(N)) { for (size_t i = 0; i < GL4.size(); i++) c0[mm4(mm4(GL4[i], N), GL4inv[i])] = nc; nc++; }
    vector<pair<int, u16>> srt; for (u16 N : raw) srt.push_back({ c0[N], N }); sort(srt.begin(), srt.end());
    memset(vid, -1, sizeof vid);
    classSize.assign(nc, 0); classRep.assign(nc, -1);
    for (auto& p : srt) { int v = (int)VN.size(); VN.push_back(p.second); cls.push_back(p.first); vid[p.second] = v; classSize[p.first]++; if (classRep[p.first] < 0) classRep[p.first] = v; }
    adj.assign(VN.size(), BS{});
    for (size_t i = 0; i < VN.size(); i++) for (size_t j = 0; j < VN.size(); j++) if (i != j && inv4(VN[i] ^ VN[j])) bset(adj[i], (int)j);
}

// ---- conjugation tables: conjTab[h][v] = vertex of h N_v h^-1; toRep[v] = some h with h N_v h^-1 = class rep;
// cent[c] = the centraliser of the class-c rep (as GL4 indices)
static vector<vector<u16>> conjTab; static vector<int> toRep; static vector<vector<int>> cent;
static void init_conj() {
    int n = (int)VN.size();
    conjTab.assign(GL4.size(), vector<u16>(n));
#pragma omp parallel for
    for (int h = 0; h < (int)GL4.size(); h++) for (int v = 0; v < n; v++) conjTab[h][v] = (u16)vid[mm4(mm4(GL4[h], VN[v]), GL4inv[h])];
    toRep.assign(n, -1); cent.assign(classSize.size(), {});
    for (int h = 0; h < (int)GL4.size(); h++) for (int v = 0; v < n; v++) if (toRep[v] < 0 && conjTab[h][v] == classRep[cls[v]]) toRep[v] = h;
    for (int c = 0; c < (int)classSize.size(); c++) for (int h = 0; h < (int)GL4.size(); h++) if (conjTab[h][classRep[c]] == classRep[c]) cent[c].push_back(h);
    for (int v = 0; v < n; v++) if (toRep[v] < 0) { fprintf(stderr, "toRep\n"); exit(1); }
    for (int c = 0; c < (int)classSize.size(); c++) if ((long long)cent[c].size() * classSize[c] != (long long)GL4.size()) { fprintf(stderr, "centraliser size\n"); exit(1); }
}

// ---- clique enumeration (each clique once: increasing vertex index)
struct CliqueStats { long long nodes[16] = { 0 }; };
template <class Leaf> static void cliques(vector<int>& cur, const BS& cand, int need, CliqueStats& st, Leaf& leaf) {
    st.nodes[cur.size()]++;
    if (need == 0) { leaf(cur); return; }
    if (bcount(cand) < need) return;
    for (int wi = 0; wi < NW; wi++) {
        u64 w = cand.w[wi];
        while (w) {
            int b = (int)_tzcnt_u64(w); w &= w - 1; int v = wi * 64 + b;
            BS nc; for (int i = 0; i < NW; i++) nc.w[i] = 0;
            nc.w[wi] = cand.w[wi] & adj[v].w[wi] & (b == 63 ? 0 : (~0ULL << (b + 1)));
            for (int i = wi + 1; i < NW; i++) nc.w[i] = cand.w[i] & adj[v].w[i];
            cur.push_back(v);
            cliques(cur, nc, need - 1, st, leaf);
            cur.pop_back();
        }
    }
}
static inline void above(BS& o, const BS& a, const BS& b, int v) { // o = a & b & {> v}
    int wi = v >> 6, bb = v & 63;
    for (int i = 0; i < wi; i++) o.w[i] = 0;
    o.w[wi] = a.w[wi] & b.w[wi] & (bb == 63 ? 0 : (~0ULL << (bb + 1)));
    for (int i = wi + 1; i < NW; i++) o.w[i] = a.w[i] & b.w[i];
}
static long long count4(const BS& Q) {
    long long n = 0;
    for (int wi = 0; wi < NW; wi++) { u64 w = Q.w[wi]; while (w) { int u = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
        BS a; above(a, Q, adj[u], u);
        for (int wj = 0; wj < NW; wj++) { u64 x = a.w[wj]; while (x) { int t = wj * 64 + (int)_tzcnt_u64(x); x &= x - 1;
            BS b; above(b, a, adj[t], t);
            for (int wk = 0; wk < NW; wk++) { u64 y = b.w[wk]; while (y) { int s = wk * 64 + (int)_tzcnt_u64(y); y &= y - 1;
                BS c; above(c, b, adj[s], s); n += bcount(c); } }
        } }
    } }
    return n;
}
static Mask graph_mask(u16 N) { Mask m; mclear(m); for (u8 x = 1; x < 16; x++) mset(m, x | (mv4(N, x) << 4)); return m; }
static Mask XMASK() { Mask m; mclear(m); for (int v = 1; v < 16; v++) mset(m, v); return m; }
static Mask YMASK() { Mask m; mclear(m); for (int v = 1; v < 16; v++) mset(m, v << 4); return m; }
static vector<Mask> members_of(const vector<u16>& Ns) { vector<Mask> M = { XMASK(), YMASK(), graph_mask(ID4) }; for (u16 N : Ns) M.push_back(graph_mask(N)); return M; }

// ---- hole-set invariants
static void pair_degrees(const Mask& H, int deg[256]) { // deg[v] = #{p in H : p ^ v in H}, for every v
    vector<u8> p = mpoints(H); for (int v = 0; v < 256; v++) deg[v] = 0;
    for (size_t i = 0; i < p.size(); i++) for (size_t j = i + 1; j < p.size(); j++) deg[p[i] ^ p[j]] += 2;
}
// pc[v] = number of planes inside H through v (v in H), or number of planes inside H + v through v (v not in H)
static int plane_counts(const Mask& H, int pc[256]) {
    vector<u8> p = mpoints(H); set<Mask> planes; vector<array<u8, 2>> lines;
    for (int v = 0; v < 256; v++) pc[v] = 0;
    for (size_t i = 0; i < p.size(); i++) for (size_t j = i + 1; j < p.size(); j++) {
        u8 a = p[i], b = p[j]; if (!mget(H, a ^ b)) continue;
        if (a < (a ^ b) && b < (a ^ b)) lines.push_back({ a, b });
        for (u8 c : p) { if (c == a || c == b || c == (a ^ b)) continue;
            if (mget(H, c ^ a) && mget(H, c ^ b) && mget(H, c ^ a ^ b)) {
                Mask m; mclear(m); u8 q[7] = { a, b, (u8)(a ^ b), c, (u8)(c ^ a), (u8)(c ^ b), (u8)(c ^ a ^ b) }; for (u8 x : q) mset(m, x); planes.insert(m); } }
    }
    for (auto& m : planes) for (int v = 1; v < 256; v++) if (mget(m, v)) pc[v]++;
    // v not in H: lines L of H with the other 4 points of <L, v> in H (planes meeting H in all but v)
    for (int v = 1; v < 256; v++) if (!mget(H, v)) {
        int c = 0;
        for (auto& l : lines) { u8 a = l[0], b = l[1]; if (mget(H, a ^ v) && mget(H, b ^ v) && mget(H, a ^ b ^ v)) c++; }
        pc[v] = c;
    }
    return (int)planes.size();
}
static string hole_key(const Mask& H) {
    vector<u8> p = mpoints(H); int deg[256]; pair_degrees(H, deg);
    int pc[256]; int np = plane_counts(H, pc);
    int lines = 0; map<int, int> dh, ho, hy, ph, po;
    for (u8 v : p) { lines += deg[v]; dh[deg[v] / 2]++; ph[pc[v]]++; }
    for (int v = 1; v < 256; v++) if (!mget(H, v)) po[pc[v]]++;
    for (int v = 1; v < 256; v++) if (!mget(H, v)) ho[deg[v] / 2]++;
    for (int u = 1; u < 256; u++) { int c = 0; for (u8 v : p) c += !par(u & v); hy[c]++; }
    char buf[64]; string s; sprintf(buf, "n%zu L%d deg", p.size(), lines / 6); s += buf;
    for (auto& kv : dh) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    s += " out"; for (auto& kv : ho) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    s += " hyp"; for (auto& kv : hy) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    sprintf(buf, " planes %d pdeg", np); s += buf; for (auto& kv : ph) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    s += " pout"; for (auto& kv : po) { sprintf(buf, " %d:%d", kv.first, kv.second); s += buf; }
    return s;
}

// ---- linear maps of V by column images
typedef array<u8, 8> Lin;
static inline u8 lin_apply(const Lin& A, u8 v) { u8 r = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) r ^= A[i]; return r; }
static Lin lin_comp(const Lin& A, const Lin& B) { Lin C; for (int i = 0; i < 8; i++) C[i] = lin_apply(A, B[i]); return C; } // A after B
static bool lin_inv_ok(const Lin& A, Lin& out) {
    u8 vec[8], comb[8]; for (int i = 0; i < 8; i++) { vec[i] = A[i]; comb[i] = (u8)(1 << i); }
    for (int bit = 0; bit < 8; bit++) {
        int p = -1; for (int i = bit; i < 8; i++) if (vec[i] >> bit & 1) { p = i; break; }
        if (p < 0) return false;
        swap(vec[p], vec[bit]); swap(comb[p], comb[bit]);
        for (int i = 0; i < 8; i++) if (i != bit && (vec[i] >> bit & 1)) { vec[i] ^= vec[bit]; comb[i] ^= comb[bit]; }
    }
    for (int k = 0; k < 8; k++) out[k] = comb[k];
    return true;
}
static Lin lin_inv(const Lin& A) { Lin C; if (!lin_inv_ok(A, C)) { fprintf(stderr, "singular\n"); exit(1); } return C; }
static Lin lin_id() { Lin I; for (int i = 0; i < 8; i++) I[i] = (u8)(1 << i); return I; }
static Mask mask_image(const Lin& A, const Mask& m) { Mask r; mclear(r); for (int v = 1; v < 256; v++) if (mget(m, v)) mset(r, lin_apply(A, (u8)v)); return r; }
static Lin diag_hh(u16 h) { Lin D; for (int i = 0; i < 4; i++) { u8 he = mv4(h, (u8)(1 << i)); D[i] = he; D[4 + i] = (u8)(he << 4); } return D; }
static Lin diag_ab(u16 a, u16 b) { Lin D; for (int i = 0; i < 4; i++) { D[i] = mv4(a, (u8)(1 << i)); D[4 + i] = (u8)(mv4(b, (u8)(1 << i)) << 4); } return D; }
static Lin random_gl8(mt19937& g) { while (true) { Lin A; for (int i = 0; i < 8; i++) A[i] = (u8)(g() & 255); Lin B; if (lin_inv_ok(A, B)) return A; } }

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
static void basis_of_mask(const Mask& m, int dim, u8* out) { pick_independent(mpoints(m), dim, out); }
static void annihilator(const Mask& m, int need, u8* out) {
    vector<u8> ann; for (int u = 1; u < 256; u++) { bool ok = true; for (int v = 1; v < 256 && ok; v++) if (mget(m, v) && par(u & v)) ok = false; if (ok) ann.push_back((u8)u); }
    pick_independent(ann, need, out);
}
// N with graph(N) = m (m complementary to Y), or -1
static int matrix_of_graph(const Mask& m) {
    u8 col[4]; for (int i = 0; i < 4; i++) { col[i] = 0xff; for (int v = 1; v < 256; v++) if (mget(m, v) && (v & 15) == (1 << i)) col[i] = (u8)(v >> 4); if (col[i] == 0xff) return -1; }
    u16 N = fromcols4(col); return meq(graph_mask(N), m) ? N : -1;
}

// ================================================================================ Part A: 10-cliques
// For class c: r = classRep[c], P = adj(r) ∩ {class >= c}; 9-cliques K in P, so K + r is a 10-clique with minimal
// class c. Weighted count #10-cliques = sum_c |class c| * sum_K 1/m(K + r) (m = number of class-c members); kept as
// integers by the factor 2520 = lcm(1..10).
static int mode_holes(const char* outfile) {
    FILE* of = fopen(outfile, "w");
    printf("Part A: 10-cliques of V5824 (= 13-partial spreads of solids with Y, X, Z among them)\n");
    long long Wext = 0, W4 = 0, Wnon = 0, Wpart = 0, allLeaves = 0, allNon = 0, allPart = 0;
    map<string, long long> invW; map<string, long long> invN; map<int, long long> qW; map<long long, long long> n4W;
    for (int c = 0; c < (int)classSize.size(); c++) {
        int r = classRep[c];
        BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
        for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
        vector<int> firsts; for (int v = 0; v < (int)VN.size(); v++) if (bget(P, v)) firsts.push_back(v);
        long long leaves = 0, non = 0, part = 0, wext = 0, w4 = 0, wnon = 0, wpart = 0; CliqueStats tot; mutex mu;
        auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1) reduction(+:leaves, non, part, wext, w4, wnon, wpart)
        for (int fi = 0; fi < (int)firsts.size(); fi++) {
            int v = firsts[fi]; BS nc; above(nc, P, adj[v], v);
            vector<int> cur = { v }; CliqueStats st;
            map<string, long long> lW, lN; map<int, long long> lq; map<long long, long long> ln4;
            auto leaf = [&](const vector<int>& K) {
                leaves++;
                int m = 1; for (int u : K) m += cls[u] == c;
                long long w = 2520 / m;
                BS Q = adj[r]; for (int u : K) { BS t; band(t, Q, adj[u]); Q = t; }
                if (bempty(Q)) {
                    non++; wnon += w;
                    Mask cov = XMASK(); mor(cov, YMASK()); mor(cov, graph_mask(ID4)); mor(cov, graph_mask(VN[r])); for (int u : K) mor(cov, graph_mask(VN[u]));
                    Mask H = mcompl(cov);
                    string key = hole_key(H); lW[key] += w; lN[key]++;
                    lock_guard<mutex> g(mu);
                    fprintf(of, "NONEXT c%d m%d :", c, m); fprintf(of, " %04x", VN[r]); for (int u : K) fprintf(of, " %04x", VN[u]); fprintf(of, " | %s\n", key.c_str());
                    return;
                }
                long long n4 = count4(Q);
                lq[bcount(Q)] += w; ln4[n4] += w;
                if (n4 == 0) {
                    part++; wpart += w;
                    lock_guard<mutex> g(mu);
                    fprintf(of, "PARTIAL c%d :", c); fprintf(of, " %04x", VN[r]); for (int u : K) fprintf(of, " %04x", VN[u]); fprintf(of, " common %d\n", bcount(Q));
                    return;
                }
                wext += w; w4 += w * n4;
            };
            cliques(cur, nc, 8, st, leaf);
            lock_guard<mutex> g(mu);
            for (int i = 0; i < 16; i++) tot.nodes[i] += st.nodes[i];
            for (auto& kv : lW) invW[kv.first] += kv.second * classSize[c];
            for (auto& kv : lN) invN[kv.first] += kv.second;
            for (auto& kv : lq) qW[kv.first] += kv.second * classSize[c];
            for (auto& kv : ln4) n4W[kv.first] += kv.second * classSize[c];
        }
        double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        printf("class %d (size %d, rep %04x, |P| = %d): %lld leaves, %lld non-extendable, %lld partial; %.0f s\n  nodes by clique size (incl. rep):",
            c, classSize[c], VN[r], bcount(P), leaves, non, part, el);
        for (int i = 0; i < 10; i++) printf(" %lld", tot.nodes[i]); printf("\n"); fflush(stdout);
        fprintf(of, "CLASS %d leaves %lld non %lld part %lld wext %lld w4 %lld wnon %lld wpart %lld (x2520)\n", c, leaves, non, part, wext, w4, wnon, wpart); fflush(of);
        Wext += wext * classSize[c]; W4 += w4 * classSize[c]; Wnon += wnon * classSize[c]; Wpart += wpart * classSize[c];
        allLeaves += leaves; allNon += non; allPart += part;
    }
    auto chk = [](long long x) { if (x % 2520) printf("  (not divisible by 2520!)\n"); return x / 2520; };
    printf("reduced leaves %lld (non-extendable %lld, partial %lld)\n", allLeaves, allNon, allPart);
    printf("10-cliques, weighted: extendable %lld, non-extendable %lld, common nbhd nonempty but no 4-clique %lld\n", chk(Wext), chk(Wnon), chk(Wpart));
    printf("sum over extendable of #4-cliques in the common nbhd: %lld (C(14,10) * 1128864 = %lld)\n", chk(W4), 1001LL * 1128864);
    printf("extendable by |common nbhd|:"); for (auto& kv : qW) printf(" %d:%lld", kv.first, chk(kv.second)); printf("\n");
    printf("extendable by #4-cliques:"); for (auto& kv : n4W) printf(" %lld:%lld", kv.first, chk(kv.second)); printf("\n");
    printf("non-extendable hole invariants (weighted count, reduced leaves):\n");
    for (auto& kv : invW) printf("  %lld  %lld  %s\n", chk(kv.second), invN[kv.first], kv.first.c_str());
    fprintf(of, "TOTAL ext %lld non %lld partial %lld w4 %lld leaves %lld\n", chk(Wext), chk(Wnon), chk(Wpart), chk(W4), allLeaves);
    for (auto& kv : invW) fprintf(of, "INV %lld %lld | %s\n", chk(kv.second), invN[kv.first], kv.first.c_str());
    fclose(of);
    return 0;
}

// Sizing for other i (mode holesk k): k-cliques (k + 3 solids), split by the largest clique in the common
// neighbourhood (capped at 14 - k, a full spread), weighted counts. No hole classification.
static int maxclique(const BS& Q, int cap) {
    int best = 0;
    function<void(const BS&, int)> rec = [&](const BS& c, int d) {
        if (d > best) best = d;
        if (best >= cap) return;
        if (d + bcount(c) <= best) return;
        for (int wi = 0; wi < NW; wi++) { u64 w = c.w[wi]; while (w) { int v = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
            BS n; above(n, c, adj[v], v); rec(n, d + 1); if (best >= cap) return; } }
    };
    rec(Q, 0); return best;
}
static int mode_holesk(int k) {
    printf("sizing: %d-cliques of V5824 (= partial spreads of %d solids with Y, X, Z among them)\n", k, k + 3);
    int cap = 14 - k; vector<long long> W(cap + 1, 0), L(cap + 1, 0); long long leavesAll = 0;
    for (int c = 0; c < (int)classSize.size(); c++) {
        int r = classRep[c];
        BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
        for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
        vector<int> firsts; for (int v = 0; v < (int)VN.size(); v++) if (bget(P, v)) firsts.push_back(v);
        vector<long long> w(cap + 1, 0), l(cap + 1, 0); mutex mu;
#pragma omp parallel for schedule(dynamic, 1)
        for (int fi = 0; fi < (int)firsts.size(); fi++) {
            int v = firsts[fi]; BS nc; above(nc, P, adj[v], v);
            vector<int> cur = { v }; CliqueStats st; vector<long long> lw(cap + 1, 0), ll(cap + 1, 0);
            auto leaf = [&](const vector<int>& K) {
                int m = 1; for (int u : K) m += cls[u] == c;
                BS Q = adj[r]; for (int u : K) { BS t; band(t, Q, adj[u]); Q = t; }
                int om = maxclique(Q, cap); lw[om] += 2520 / m; ll[om]++;
            };
            cliques(cur, nc, k - 2, st, leaf);
            lock_guard<mutex> g(mu); for (int i = 0; i <= cap; i++) { w[i] += lw[i]; l[i] += ll[i]; }
        }
        for (int i = 0; i <= cap; i++) { W[i] += w[i] * classSize[c]; L[i] += l[i]; leavesAll += l[i]; }
        printf("class %d done\n", c); fflush(stdout);
    }
    printf("reduced leaves %lld\n", leavesAll);
    for (int i = 0; i <= cap; i++) printf("  extends by at most %d solids: weighted %lld%s, reduced leaves %lld\n", i, W[i] / 2520, W[i] % 2520 ? " (not integral!)" : "", L[i]);
    return 0;
}

// ================================================================================ isomorphism of point sets
// Linear maps g with g(H1) = H2, by backtracking over images of a basis inside H1. Pruning: on every element of
// the span built so far, v in H1 <=> g(v) in H2 and deg1(v) = deg2(g v) (deg = number of H-pairs summing to v).
struct IsoCtx { Mask H1, H2; int d1[256], d2[256]; u8 b[8]; u8 vs[256], im[256]; long long count; bool first; Lin found; vector<u8> cand2; };
static void iso_rec(IsoCtx& X, int k) {
    if (k == 8) { X.count++; if (X.count == 1) for (int i = 0; i < 8; i++) X.found[i] = X.im[1 << i]; return; }
    int n = 1 << k; u8 bk = X.b[k];
    for (u8 c : X.cand2) {
        if (X.d2[c] != X.d1[bk]) continue;
        bool ok = true;
        for (int s = 0; s < n && ok; s++) {
            u8 v = X.vs[s] ^ bk, w = X.im[s] ^ c;
            if (!w || mget(X.H1, v) != mget(X.H2, w) || X.d1[v] != X.d2[w]) ok = false;
        }
        if (!ok) continue;
        for (int s = 0; s < n; s++) { X.vs[n + s] = X.vs[s] ^ bk; X.im[n + s] = X.im[s] ^ c; }
        iso_rec(X, k + 1);
        if (X.first && X.count) return;
    }
}
static long long iso_count(const Mask& H1, const Mask& H2, bool first, Lin* out) {
    IsoCtx X; X.H1 = H1; X.H2 = H2; pair_degrees(H1, X.d1); pair_degrees(H2, X.d2); X.count = 0; X.first = first;
    if (mcount(H1) != mcount(H2)) return 0;
    { int p1[256], p2[256]; plane_counts(H1, p1); plane_counts(H2, p2); for (int v = 0; v < 256; v++) { X.d1[v] = X.d1[v] * 65536 + p1[v]; X.d2[v] = X.d2[v] * 65536 + p2[v]; } }
    // basis inside H1, rarest degree first
    vector<u8> p1 = mpoints(H1); map<int, int> freq; for (u8 v : p1) freq[X.d1[v]]++;
    stable_sort(p1.begin(), p1.end(), [&](u8 a, u8 b) { return freq[X.d1[a]] < freq[X.d1[b]]; });
    pick_independent(p1, 8, X.b);
    X.cand2 = mpoints(H2);
    X.vs[0] = 0; X.im[0] = 0;
    iso_rec(X, 0);
    if (out && X.count) { // found maps b_i -> im[1<<i]; as column images of e_i: g = M2 * M1^-1
        Lin B1; for (int i = 0; i < 8; i++) B1[i] = X.b[i];
        *out = lin_comp(X.found, lin_inv(B1));
        if (!meq(mask_image(*out, H1), H2)) { fprintf(stderr, "iso map check failed\n"); exit(1); }
    }
    return X.count;
}

// ================================================================================ canonical form of a 13-partial spread
// Frames: an ordered triple (A, B, C) of members, sent to (X, Y, Z) by Cmap (basis a_i of A -> e_i, the B-part b_i
// of a_i in B + C -> e_{4+i}); every g in GL(8,2) with g(A, B, C) = (X, Y, Z) is diag(h, h) o Cmap. The image of
// the other 10 members is a 10-clique of V5824. Key = the least sorted vertex list over the frames whose ordering
// has the least sorted conjugacy-class vector (an intrinsic condition on the ordered triple, so the key is an
// invariant of the GL(8,2)-orbit, and it's an image of P, so equal keys mean equivalent P). Within those, only h
// moving a least-class member onto its class rep can give the least list. Every frame giving the key lies in this
// set, and those frames correspond one-to-one to Aut(P) (setwise stabiliser of the 13 members).
struct Canon { vector<u16> key; long long aut = 0; Lin frame; };
static bool frame_of(const vector<array<u8, 4>>& bas, const vector<vector<u8>>& pts, const vector<Mask>& M, int a, int b, int c, Lin& Cm) {
    Lin Ci;
    for (int i = 0; i < 4; i++) {
        u8 ai = bas[a][i]; int bi = -1;
        for (u8 p : pts[b]) if (mget(M[c], ai ^ p)) { bi = p; break; }
        if (bi < 0) return false;
        Ci[i] = ai; Ci[4 + i] = (u8)bi;
    }
    return lin_inv_ok(Ci, Cm);
}
static int vertex_in_frame(const Lin& Cm, const array<u8, 4>& bd) {
    u8 x[4], y[4];
    for (int j = 0; j < 4; j++) { u8 w = lin_apply(Cm, bd[j]); x[j] = w & 15; y[j] = w >> 4; }
    u16 Xm = fromcols4(x), Ym = fromcols4(y);
    if (gl4idx[Xm] < 0) return -1;
    u16 N = mm4(Ym, inverse4(Xm));
    return vid[N];
}
static Canon canon13(const vector<Mask>& M) {
    int n = (int)M.size(); if (n != 13) { fprintf(stderr, "canon13: %d members\n", n); exit(1); }
    vector<array<u8, 4>> bas(n); vector<vector<u8>> pts(n);
    for (int i = 0; i < n; i++) { basis_of_mask(M[i], 4, bas[i].data()); pts[i] = mpoints(M[i]); }
    struct Ord { Lin Cm; array<u16, 10> v; };
    vector<Ord> keep; array<int, 10> bestCls; bestCls.fill(INT_MAX);
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) for (int c = 0; c < n; c++) {
        if (a == b || a == c || b == c) continue;
        Ord o; if (!frame_of(bas, pts, M, a, b, c, o.Cm)) { fprintf(stderr, "canon13: members not complementary\n"); exit(1); }
        int k = 0; array<int, 10> cv;
        for (int d = 0; d < n; d++) { if (d == a || d == b || d == c) continue; int v = vertex_in_frame(o.Cm, bas[d]); if (v < 0) { fprintf(stderr, "canon13: not in V5824\n"); exit(1); } o.v[k] = (u16)v; cv[k] = cls[v]; k++; }
        sort(cv.begin(), cv.end());
        if (cv < bestCls) { bestCls = cv; keep.clear(); }
        if (cv == bestCls) keep.push_back(o);
    }
    Canon R; R.key.assign(10, 0xffff); array<u16, 10> img;
    int c0 = bestCls[0], r0 = classRep[c0];
    for (auto& o : keep) for (int e = 0; e < 10; e++) {
        if (cls[o.v[e]] != c0) continue;
        u16 h0 = GL4[toRep[o.v[e]]];
        for (int z : cent[c0]) {
            int hi = gl4idx[mm4(GL4[z], h0)];
            for (int j = 0; j < 10; j++) img[j] = conjTab[hi][o.v[j]];
            sort(img.begin(), img.end());
            if (img[0] != r0) { fprintf(stderr, "canon13: rep\n"); exit(1); }
            int cmp = 0; for (int j = 0; j < 10 && !cmp; j++) cmp = img[j] < R.key[j] ? -1 : img[j] > R.key[j] ? 1 : 0;
            if (cmp < 0) { R.key.assign(img.begin(), img.end()); R.aut = 1; R.frame = lin_comp(diag_hh(GL4[hi]), o.Cm); }
            else if (cmp == 0) R.aut++;
        }
    }
    return R;
}
static vector<u16> key_matrices(const vector<u16>& key) { vector<u16> N; for (u16 v : key) N.push_back(VN[v]); return N; }

// ================================================================================ Kurz's hole sets, solids, covers
static vector<Mask> read_kurz(const char* file) {
    FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    vector<string> rows; char line[256];
    while (fgets(line, sizeof line, fp)) { string s(line); while (s.size() && (s.back() == '\n' || s.back() == '\r')) s.pop_back(); if (s.size() == 60) rows.push_back(s); }
    fclose(fp);
    if (rows.size() != 24) { printf("expected 24 rows of 60, got %zu\n", rows.size()); exit(1); }
    vector<Mask> out;
    for (int b = 0; b < 3; b++) {
        Mask m; mclear(m); int distinct = 0;
        for (int j = 0; j < 60; j++) { int p = 0; for (int i = 0; i < 8; i++) if (rows[8 * b + i][j] == '1') p |= 1 << i; if (p && !mget(m, p)) distinct++; mset(m, p); }
        if (distinct != 60 || mget(m, 0)) { printf("Kurz set %d: %d distinct nonzero columns\n", b + 1, distinct); exit(1); }
        out.push_back(m);
    }
    return out;
}
static vector<Mask> ALLSOLIDS;
static void init_solids() {
    for (int p0 = 0; p0 < 8; p0++) for (int p1 = p0 + 1; p1 < 8; p1++) for (int p2 = p1 + 1; p2 < 8; p2++) for (int p3 = p2 + 1; p3 < 8; p3++) {
        int piv[4] = { p0, p1, p2, p3 }; u8 pm = (u8)((1 << p0) | (1 << p1) | (1 << p2) | (1 << p3));
        vector<pair<int, int>> fr; for (int k = 0; k < 4; k++) for (int q = piv[k] + 1; q < 8; q++) if (!(pm >> q & 1)) fr.push_back({ k, q });
        for (u32 t = 0; t < (1u << fr.size()); t++) {
            u8 row[4]; for (int k = 0; k < 4; k++) row[k] = (u8)(1 << piv[k]);
            for (size_t f = 0; f < fr.size(); f++) if (t >> f & 1) row[fr[f].first] |= (u8)(1 << fr[f].second);
            Mask m; mclear(m); for (int s = 1; s < 16; s++) { u8 v = 0; for (int k = 0; k < 4; k++) if (s >> k & 1) v ^= row[k]; mset(m, v); }
            ALLSOLIDS.push_back(m);
        }
    }
}
static bool contains_solid(const Mask& H) { for (auto& s : ALLSOLIDS) if (msub(s, H)) return true; return false; }
static void covers_rec(const vector<Mask>& S, const vector<vector<int>>& at, const Mask& C, Mask& cov, vector<int>& ch, vector<vector<int>>& out) {
    int best = -1, bc = INT_MAX;
    for (int v = 1; v < 256; v++) if (mget(C, v) && !mget(cov, v)) {
        int c = 0; for (int s : at[v]) if (mdisj(S[s], cov)) c++;
        if (c < bc) { bc = c; best = v; if (!c) break; }
    }
    if (best < 0) { out.push_back(ch); return; }
    if (!bc) return;
    for (int s : at[best]) if (mdisj(S[s], cov)) { Mask save = cov; mor(cov, S[s]); ch.push_back(s); covers_rec(S, at, C, cov, ch, out); ch.pop_back(); cov = save; }
}
static vector<vector<Mask>> covers_of(const Mask& H, int& nsolids) {
    Mask C = mcompl(H); vector<Mask> S; for (auto& s : ALLSOLIDS) if (msub(s, C)) S.push_back(s);
    nsolids = (int)S.size();
    vector<vector<int>> at(256); for (int i = 0; i < (int)S.size(); i++) for (int v = 1; v < 256; v++) if (mget(S[i], v)) at[v].push_back(i);
    Mask cov; mclear(cov); vector<int> ch; vector<vector<int>> raw; covers_rec(S, at, C, cov, ch, raw);
    vector<vector<Mask>> out; for (auto& r : raw) { vector<Mask> M; for (int s : r) M.push_back(S[s]); out.push_back(M); }
    return out;
}

// ================================================================================ configs
static vector<vector<int>> spread_cliques(double& weighted) {
    vector<vector<int>> red; weighted = 0;
    for (int c = 0; c < (int)classSize.size(); c++) {
        int r = classRep[c];
        BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
        for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
        vector<int> cur = { r }; CliqueStats st; double wsum = 0;
        auto leaf = [&](const vector<int>& K) { int m = 0; for (int u : K) m += cls[u] == c; wsum += 1.0 / m; red.push_back(K); };
        cliques(cur, P, 13, st, leaf);
        weighted += classSize[c] * wsum;
    }
    return red;
}
struct ConfInfo { int type; int hole; long long aut; long long hits = 0; int srcExt = 0, srcPartA = 0, srcKurz = 0; };
static int mode_configs(const char* holesfile, const char* kurzfile, const char* outfile) {
    init_conj(); init_solids();
    printf("solids in PG(7,2): %zu (expect 200787)\n", ALLSOLIDS.size()); if (ALLSOLIDS.size() != 200787) return 1;
    const long long FR = 1716LL * 20160; // frames per partial spread / |GL(8,2)| times #10-cliques
    // weighted counts from Part A
    long long Aext = -1, Anon = -1;
    { FILE* fp = fopen(holesfile, "r"); if (!fp) { printf("no holes file\n"); return 1; } char line[4096];
      while (fgets(line, sizeof line, fp)) if (!strncmp(line, "TOTAL", 5)) sscanf(line, "TOTAL ext %lld non %lld", &Aext, &Anon);
      fclose(fp); }
    printf("Part A weighted counts: extendable %lld, non-extendable %lld\n", Aext, Anon);
    map<vector<u16>, ConfInfo> conf;
    // --- hole types: Kurz's three, and the extendable types found below
    vector<Mask> kurz = read_kurz(kurzfile); vector<Mask> holeRep; vector<string> holeKey; vector<long long> holeAut;
    auto hole_class = [&](const Mask& H, bool add) -> int {
        string k = hole_key(H);
        for (size_t i = 0; i < holeRep.size(); i++) if (holeKey[i] == k && iso_count(H, holeRep[i], true, nullptr)) return (int)i;
        if (!add) return -1;
        holeRep.push_back(H); holeKey.push_back(k); holeAut.push_back(iso_count(H, H, false, nullptr)); return (int)holeRep.size() - 1;
    };
    for (int i = 0; i < 3; i++) {
        Mask H = kurz[i]; Mask C = mcompl(H);
        int rk = 0; { vector<u8> p = mpoints(H); bool sp[256] = { 0 }; sp[0] = 1; for (u8 v : p) if (!sp[v]) { rk++; bool nw[256]; memcpy(nw, sp, 256); for (int w = 0; w < 256; w++) if (sp[w]) nw[w ^ v] = 1; memcpy(sp, nw, 256); } }
        int hc = hole_class(H, true);
        printf("Kurz set %d: rank %d, contains a solid: %s, |Aut| = %lld, class %d\n  %s\n", i + 1, rk, contains_solid(H) ? "YES" : "no", holeAut[hc], hc, holeKey[hc].c_str());
    }
    for (int i = 0; i < 3; i++) for (int j = 0; j < i; j++) if (iso_count(kurz[i], kurz[j], true, nullptr)) printf("  Kurz sets %d and %d are isomorphic!\n", i + 1, j + 1);
    if (holeRep.size() != 3) { printf("Kurz's sets are not pairwise non-isomorphic\n"); return 1; }
    // --- Kurz covers
    vector<long long> kurzCovers(3);
    for (int i = 0; i < 3; i++) {
        int ns; auto cv = covers_of(kurz[i], ns); kurzCovers[i] = (long long)cv.size();
        map<vector<u16>, long long> orb; long long sumOrb = 0;
        for (auto& M : cv) {
            Mask cov; mclear(cov); for (auto& m : M) { if (!mdisj(cov, m)) { printf("cover overlap\n"); return 1; } mor(cov, m); }
            if (!meq(mcompl(cov), kurz[i])) { printf("cover holes mismatch\n"); return 1; }
            Canon K = canon13(M); orb[K.key] = K.aut;
            auto& ci = conf[K.key]; ci.type = 1; ci.hole = i; ci.aut = K.aut; ci.srcKurz++;
        }
        for (auto& kv : orb) sumOrb += holeAut[i] / kv.second;
        long long w = 0; for (auto& kv : orb) w += FR / kv.second;
        printf("Kurz set %d: %d solids in the complement, %zu covers by 13 solids, %zu orbits; sum |Aut H|/|Aut P| = %lld (= covers?)%s; 10-cliques %lld (= 1716*20160*covers/|Aut H| = %lld)\n",
            i + 1, ns, cv.size(), orb.size(), sumOrb, sumOrb == (long long)cv.size() ? " yes" : " NO", w, FR * (long long)cv.size() / holeAut[i]);
        for (auto& kv : orb) printf("    orbit |Aut P| = %lld\n", kv.second);
    }
    // --- Part A non-extendable cliques
    { FILE* fp = fopen(holesfile, "r"); char line[4096]; long long nl = 0, bad = 0; map<int, long long> byClass;
      vector<vector<u16>> lines;
      while (fgets(line, sizeof line, fp)) {
          if (strncmp(line, "NONEXT", 6)) continue;
          const char* p = strchr(line, ':') + 1; vector<u16> Ns; char* e;
          while (true) { while (*p == ' ') p++; if (*p == '|') break; long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
          if (Ns.size() != 10) { printf("bad NONEXT line\n"); return 1; }
          lines.push_back(Ns);
      }
      fclose(fp);
      // canonical form of each; one whose orbit is not among the Kurz covers' gets an explicit isomorphism test
      int n = (int)lines.size(); vector<Canon> can(n);
#pragma omp parallel for schedule(dynamic, 16)
      for (int i = 0; i < n; i++) can[i] = canon13(members_of(lines[i]));
      long long sample = 0, sampleOk = 0;
      for (int i = 0; i < n; i++) {
          nl++;
          auto it = conf.find(can[i].key);
          vector<Mask> M = members_of(lines[i]); Mask cov; mclear(cov); for (auto& m : M) mor(cov, m); Mask H = mcompl(cov);
          if (it == conf.end() || it->second.srcKurz == 0) {
              int hc = -1; for (int j = 0; j < 3 && hc < 0; j++) if (iso_count(H, kurz[j], true, nullptr)) hc = j;
              printf("Part A non-extendable clique %lld: orbit not among the Kurz covers' (hole set iso to Kurz set %d)\n", nl, hc + 1); bad++;
              auto& ci = conf[can[i].key]; ci.type = 1; ci.hole = hc; ci.aut = can[i].aut; ci.srcPartA++; continue;
          }
          if (it->second.aut != can[i].aut) { printf("aut mismatch\n"); bad++; }
          byClass[it->second.hole]++; it->second.srcPartA++;
          if (i % 997 == 0) { sample++; sampleOk += iso_count(H, kurz[it->second.hole], true, nullptr) > 0; } // direct check on a sample
      }
      printf("Part A non-extendable cliques: %lld, orbit not among the Kurz covers': %lld; by Kurz set:", nl, bad); for (auto& kv : byClass) printf(" %d:%lld", kv.first + 1, kv.second); printf("\n");
      printf("  direct isomorphism test on a sample: %lld of %lld hole sets isomorphic to their orbit's Kurz set\n", sampleOk, sample);
      long long kurzOrbitsHit = 0, kurzOrbits = 0; for (auto& kv : conf) if (kv.second.srcKurz) { kurzOrbits++; kurzOrbitsHit += kv.second.srcPartA > 0; }
      printf("  Kurz-cover orbits reached by Part A: %lld of %lld\n", kurzOrbitsHit, kurzOrbits);
    }
    // --- extendable: own spread sets x 14 choices of the fourth hole
    double wsp; auto red = spread_cliques(wsp);
    printf("14-cliques: %zu conjugation-reduced, weighted %.3f (expect 1128864)\n", red.size(), wsp);
    { int n = (int)red.size(); vector<vector<Canon>> can(n, vector<Canon>(14)); vector<vector<int>> hcl(n, vector<int>(14));
      map<int, int> mClass; // hole type of X + Y + Z + graph(M), by the vertex M
      for (int i = 0; i < n; i++) for (int j = 0; j < 14; j++) if (!mClass.count(red[i][j])) {
          Mask H = XMASK(); mor(H, YMASK()); mor(H, graph_mask(ID4)); mor(H, graph_mask(VN[red[i][j]]));
          mClass[red[i][j]] = hole_class(H, true);
      }
      printf("extendable: %zu distinct fourth holes M, hole types so far %zu\n", mClass.size(), holeRep.size()); fflush(stdout);
#pragma omp parallel for schedule(dynamic, 2)
      for (int i = 0; i < n; i++) for (int j = 0; j < 14; j++) {
          vector<u16> Ns; for (int t = 0; t < 14; t++) if (t != j) Ns.push_back(VN[red[i][t]]);
          vector<Mask> M; for (u16 N : Ns) M.push_back(graph_mask(N));
          can[i][j] = canon13(M);
          hcl[i][j] = mClass.at(red[i][j]);
      }
      for (int i = 0; i < n; i++) for (int j = 0; j < 14; j++) {
          auto it = conf.find(can[i][j].key);
          if (it != conf.end() && it->second.type == 1) { printf("extendable config equals a non-extendable one!\n"); return 1; }
          auto& ci = conf[can[i][j].key]; if (!ci.srcExt) { ci.type = 0; ci.hole = hcl[i][j]; ci.aut = can[i][j].aut; } else if (ci.hole != hcl[i][j] || ci.aut != can[i][j].aut) { printf("ext mismatch\n"); return 1; }
          ci.srcExt++;
      }
    }
    printf("hole types: %zu\n", holeRep.size());
    for (size_t i = 0; i < holeRep.size(); i++) printf("  type %zu: |Aut| = %lld, contains a solid: %s  %s\n", i, holeAut[i], contains_solid(holeRep[i]) ? "yes" : "no", holeKey[i].c_str());
    // --- identities and output
    long long sExt = 0, sNon = 0, nExt = 0, nNon = 0; map<int, long long> perHole, perHoleN;
    for (auto& kv : conf) { long long w = FR / kv.second.aut; if (FR % kv.second.aut) printf("aut does not divide\n");
        if (kv.second.type == 0) { sExt += w; nExt++; } else { sNon += w; nNon++; } perHole[kv.second.hole] += w; perHoleN[kv.second.hole]++; }
    printf("configurations up to GL(8,2): %lld extendable, %lld non-extendable\n", nExt, nNon);
    printf("sum 1716*20160/|Aut P|: extendable %lld (Part A %lld) %s, non-extendable %lld (Part A %lld) %s\n", sExt, Aext, sExt == Aext ? "MATCH" : "DIFFER", sNon, Anon, sNon == Anon ? "MATCH" : "DIFFER");
    for (auto& kv : perHole) printf("  hole type %d: %lld configurations, 10-cliques %lld\n", kv.first, perHoleN[kv.first], kv.second);
    FILE* of = fopen(outfile, "w"); int idx = 0;
    for (int t = 0; t < 2; t++) for (auto& kv : conf) if (kv.second.type == t) {
        fprintf(of, "CONF %d %s hole %d aut %lld src %d %d %d :", idx++, t ? "nonext" : "ext", kv.second.hole, kv.second.aut, kv.second.srcExt, kv.second.srcPartA, kv.second.srcKurz);
        for (u16 N : key_matrices(kv.first)) fprintf(of, " %04x", N); fprintf(of, "\n");
    }
    fclose(of);
    printf("written %d configurations to %s\n", idx, outfile);
    return 0;
}

// ================================================================================ forms: the search
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
struct Verdict { bool apn; int dist[9]; };
static Verdict verify_space(const vector<u64>& basis) { // from scratch: LUT, derivative counts, Walsh spectrum
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
static void change_coords(vector<u64>& T, const Lin& A) { Lin Ai = lin_inv(A); for (auto& f : T) f = transform(f, Ai.data()); } // w = A v

// The problem: 13 rank-4 radicals (members), holes D = the rest. Relaxation for controls: oracle rank-6 forms.
struct Problem {
    int nr = 0; Mask R[13]; int radOf[256]; Mask D;
    vector<u64> f4[13];
    vector<Mask> lineM; vector<vector<u64>> lineF; vector<int> linesAt[256];
    vector<int> order;
    vector<u64> oracle; vector<Mask> oracleM;
};
static void build_problem(Problem& P, const vector<Mask>& M) {
    P.nr = (int)M.size();
    Mask cov; mclear(cov);
    for (int v = 0; v < 256; v++) P.radOf[v] = -1;
    for (int j = 0; j < P.nr; j++) {
        P.R[j] = M[j];
        if (mcount(M[j]) != 15 || !mdisj(cov, M[j])) { fprintf(stderr, "members not a partial spread\n"); exit(1); }
        mor(cov, M[j]);
        for (int v = 1; v < 256; v++) if (mget(M[j], v)) P.radOf[v] = j;
        u8 u[4]; annihilator(M[j], 4, u);
        P.f4[j].clear(); forms_from(u, 4, ND4, P.f4[j]);
    }
    P.D = mcompl(cov);
    set<Mask> seen;
    P.lineM.clear(); P.lineF.clear(); for (int v = 0; v < 256; v++) P.linesAt[v].clear();
    for (int a = 1; a < 256; a++) for (int b = a + 1; b < 256; b++) if (mget(P.D, a) && mget(P.D, b) && mget(P.D, a ^ b)) {
        Mask m; mclear(m); mset(m, a); mset(m, b); mset(m, a ^ b);
        if (!seen.insert(m).second) continue;
        int li = (int)P.lineM.size(); P.lineM.push_back(m);
        u8 ann[6]; annihilator(m, 6, ann); vector<u64> fs; forms_from(ann, 6, ND6, fs); P.lineF.push_back(fs);
        for (int v : { a, b, a ^ b }) P.linesAt[v].push_back(li);
    }
}
static inline bool legal(const Problem& P, u64 f, Mask& m, int& cnt) {
    if (rank8(f) == 8) { cnt = 0; return true; }
    cnt = radical(f, m);
    if (cnt == 15) { int v = 0; for (v = 1; v < 256; v++) if (mget(m, v)) break; int j = P.radOf[v]; return j >= 0 && meq(m, P.R[j]); }
    if (cnt == 3) {
        if (msub(m, P.D)) return true;
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
    else for (int li : P.linesAt[a]) if (mdisj(P.lineM[li], cov)) for (u64 f : P.lineF[li]) cb(f);
    for (size_t i = 0; i < P.oracle.size(); i++) if (mget(P.oracleM[i], a) && !msub(P.oracleM[i], P.D) && mdisj(P.oracleM[i], cov)) cb(P.oracle[i]);
}
struct Stats { long long nodes[10] = { 0 }, tries[10] = { 0 }, bentext = 0, full = 0, p1[10] = { 0 }; };
static void add_stats(Stats& a, const Stats& b) { for (int i = 0; i < 10; i++) { a.nodes[i] += b.nodes[i]; a.tries[i] += b.tries[i]; a.p1[i] += b.p1[i]; } a.bentext += b.bentext; a.full += b.full; }
static mutex solmu;
static map<string, vector<u64>> solutions;
static map<string, vector<u64>> p1spaces; // spans in which every point off D is covered (all 13 rank-4 forms present)
static bool g_logp1 = true;
static thread_local map<string, vector<u64>>* tl_sol = nullptr;
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
    if (a < 0 || P.radOf[a] < 0) {
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
    auto* sink = tl_sol; // the caller's solution sink (controls), for whichever thread runs a task
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nt; i++) { auto* save = tl_sol; tl_sol = sink; dfs(P, tasks[i], per[omp_get_thread_num()]); tl_sol = save; }
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
struct Conf { int idx; bool nonext; int hole; long long aut; vector<u16> Ns; vector<Mask> mem; };
// Cross-run input: lines "MEM i k : 13 x 4 basis vectors" (mode covers: every cover of Kurz set k, his coordinates)
static vector<Conf> read_mems(const char* file) {
    vector<Conf> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "MEM", 3)) continue;
        Conf c; c.nonext = true; c.aut = 0;
        if (sscanf(line, "MEM %d %d", &c.idx, &c.hole) != 2) continue;
        const char* p = strchr(line, ':') + 1; char* e; vector<u8> b;
        while (true) { long v = strtol(p, &e, 16); if (e == p) break; b.push_back((u8)v); p = e; }
        if (b.size() != 52) { printf("bad MEM line\n"); exit(1); }
        for (int j = 0; j < 13; j++) { Mask m; mclear(m); for (int s = 1; s < 16; s++) { u8 v = 0; for (int k = 0; k < 4; k++) if (s >> k & 1) v ^= b[4 * j + k]; mset(m, v); } c.mem.push_back(m); }
        out.push_back(c);
    }
    fclose(fp); return out;
}
static int mode_covers(const char* kurzfile, const char* outfile) {
    init_solids(); vector<Mask> kurz = read_kurz(kurzfile);
    FILE* of = fopen(outfile, "w"); int idx = 0;
    for (int k = 0; k < 3; k++) {
        int ns; auto cv = covers_of(kurz[k], ns);
        for (auto& M : cv) { fprintf(of, "MEM %d %d :", idx++, k); for (auto& m : M) { u8 b[4]; basis_of_mask(m, 4, b); for (int i = 0; i < 4; i++) fprintf(of, " %02x", b[i]); } fprintf(of, "\n"); }
        printf("Kurz set %d: %zu covers\n", k + 1, cv.size());
    }
    fclose(of); return 0;
}
static vector<Conf> read_confs(const char* file) {
    vector<Conf> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        Conf c; char ty[16];
        if (strncmp(line, "CONF", 4)) continue;
        if (sscanf(line, "CONF %d %15s hole %d aut %lld", &c.idx, ty, &c.hole, &c.aut) != 4) continue;
        c.nonext = !strcmp(ty, "nonext");
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; c.Ns.push_back((u16)v); p = e; }
        if (c.Ns.size() != 10) { printf("bad CONF line\n"); exit(1); }
        out.push_back(c);
    }
    fclose(fp); return out;
}
// Cross-run input: lines "SPREAD i ... : 14 matrices" (fourteen-check's spreads.txt), each with every choice of the
// 4th hole among its 14 -> configurations (not reduced); index = 14 * i + j.
static vector<Conf> read_spreads_as_confs(const char* file) {
    vector<Conf> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "SPREAD", 6)) continue;
        int idx = atoi(line + 7); const char* p = strchr(line, ':'); if (!p) continue; p++;
        vector<u16> Ns; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
        if (Ns.size() != 14) { printf("bad spread line\n"); exit(1); }
        for (int j = 0; j < 14; j++) { Conf c; c.idx = 14 * idx + j; c.nonext = false; c.hole = -1; c.aut = 0; c.Ns = Ns; c.Ns.erase(c.Ns.begin() + j); out.push_back(c); }
    }
    fclose(fp); return out;
}
// members of a configuration: X, Y, Z and the 10 graphs (CONF), or the 13 graphs (from spreads, holes X, Y, Z, 4th)
static vector<Mask> conf_members(const Conf& c) {
    if (!c.mem.empty()) return c.mem;
    if (c.Ns.size() == 10) return members_of(c.Ns);
    vector<Mask> M; for (u16 N : c.Ns) M.push_back(graph_mask(N)); return M;
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
static int g_input = 0; // 0 CONF, 1 SPREAD x 14, 2 MEM
static int mode_search(const char* conffile, bool fromSpreads, const char* progfile, int how, unsigned seed, int first, int last) {
    auto cf = g_input == 1 ? read_spreads_as_confs(conffile) : g_input == 2 ? read_mems(conffile) : read_confs(conffile);
    printf("search: %zu configurations%s; point order %s\n", cf.size(), fromSpreads ? " (spread file x 14)" : "", how == 0 ? "ascending" : how == 1 ? "descending" : "random");
    set<int> done; Stats tot, totExt, totNon;
    if (FILE* fp = fopen(progfile, "r")) {
        char line[1024];
        while (fgets(line, sizeof line, fp)) { int id; Stats s; if (read_stats(line, id, s) && !done.count(id)) { done.insert(id); add_stats(tot, s); } else if (!strncmp(line, "SOLUTION", 8)) printf("(previous run) %s", line); }
        fclose(fp);
    }
    FILE* pf = fopen(progfile, "a");
    auto t0 = chrono::steady_clock::now(); map<size_t, int> lineCounts;
    for (auto& c : cf) {
        if (c.idx < first || c.idx > last || done.count(c.idx)) continue;
        Problem P; build_problem(P, conf_members(c)); make_order(P, how, seed + c.idx);
        lineCounts[P.lineM.size()]++;
        if (mcount(P.D) != 60) { printf("config %d: |D| = %d\n", c.idx, mcount(P.D)); return 1; }
        Stats st = run_search(P, empty_state(), 400);
        { lock_guard<mutex> g(solmu); for (auto& kv : solutions) fprintf(pf, "SOLUTION %d %s\n", c.idx, kv.first.c_str()); }
        write_stats(pf, c.idx, st);
        add_stats(tot, st); add_stats(c.nonext ? totNon : totExt, st);
        printf("config %d (%s, %zu lines in D): nodes", c.idx, c.nonext ? "nonext" : "ext", P.lineM.size()); for (int k = 0; k <= 8; k++) printf(" %lld", st.nodes[k]);
        printf("  off-D-covered"); for (int k = 0; k <= 8; k++) printf(" %lld", st.p1[k]);
        printf("  full %lld bentext %lld  [%.0f s]\n", st.full, st.bentext, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    }
    for (auto& kv : p1spaces) { fprintf(pf, "P1 dim %zu basis", kv.second.size()); for (u64 f : kv.second) fprintf(pf, " %07x", code28(f)); fprintf(pf, "\n"); }
    fclose(pf);
    printf("lines in D:"); for (auto& kv : lineCounts) printf(" %zu:%d", kv.first, kv.second); printf("\n");
    print_stats("ext (this run)", totExt); print_stats("nonext (this run)", totNon); print_stats("search (incl. resumed)", tot);
    printf("distinct off-D-covered spans this run: %zu\n", p1spaces.size());
    printf("solutions: %zu\n", solutions.size());
    for (auto& kv : solutions) { Verdict V = verify_space(kv.second); printf("  %s APN=%d dist:", kv.first.c_str(), V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("\n"); }
    return 0;
}

// ================================================================================ positive controls
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
static bool find_clique(vector<int>& cur, const BS& cand, int need) {
    if (need == 0) return true;
    for (int wi = 0; wi < NW; wi++) { u64 w = cand.w[wi]; while (w) { int v = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
        BS nc; above(nc, cand, adj[v], v);
        cur.push_back(v); if (find_clique(cur, nc, need - 1)) return true; cur.pop_back(); } }
    return false;
}
static bool find_clique_rand(vector<int>& cur, const BS& cand, int need, mt19937& rng) {
    if (need == 0) return true;
    vector<int> vs; for (int wi = 0; wi < NW; wi++) { u64 w = cand.w[wi]; while (w) { vs.push_back(wi * 64 + (int)_tzcnt_u64(w)); w &= w - 1; } }
    if ((int)vs.size() < need) return false;
    shuffle(vs.begin(), vs.end(), rng);
    int tries = 0;
    for (int v : vs) { BS nc; band(nc, cand, adj[v]); cur.push_back(v); if (find_clique_rand(cur, nc, need - 1, rng)) return true; cur.pop_back(); if (++tries > 64) break; }
    return false;
}
struct CtlCtx { vector<Conf> cf; map<vector<u16>, int> repOf; vector<int> nonextIdx; };
static CtlCtx ctl_ctx(const char* conffile) {
    CtlCtx X; X.cf = read_confs(conffile);
    for (size_t i = 0; i < X.cf.size(); i++) {
        if (X.cf[i].idx != (int)i) { printf("config file not indexed 0..n-1\n"); exit(1); }
        vector<u16> key; for (u16 N : X.cf[i].Ns) key.push_back((u16)vid[N]); sort(key.begin(), key.end()); X.repOf[key] = (int)i;
        if (X.cf[i].nonext) X.nonextIdx.push_back((int)i);
    }
    return X;
}
struct CtlResult { int status = 1; int rep = -1; int n4 = 0, inD = 0, nsol = 0, napn = 0; bool known = false; string why; Stats st; };
// Normalise T: rank-4 radicals R1 -> X, R2 -> Y, R3 (if any) -> graph(I). Returns the profile in the new coordinates.
static bool normalise(vector<u64>& T, Profile& Q, string& why) {
    Mask Xm = XMASK(), Ym = YMASK();
    Lin Ai; basis_of_mask(Q.r4[0], 4, Ai.data()); basis_of_mask(Q.r4[1], 4, Ai.data() + 4);
    change_coords(T, lin_inv(Ai)); Q = profile(T);
    vector<int> rest; for (size_t i = 0; i < Q.r4.size(); i++) if (!meq(Q.r4[i], Xm) && !meq(Q.r4[i], Ym)) rest.push_back((int)i);
    if (!rest.empty()) {
        int M = matrix_of_graph(Q.r4[rest[0]]); if (M < 0) { why = "not a graph"; return false; }
        Lin B; for (int i = 0; i < 4; i++) { B[i] = mv4((u16)M, (u8)(1 << i)); B[4 + i] = (u8)(1 << (4 + i)); }
        change_coords(T, B); Q = profile(T);
    }
    return true;
}
// Run the search on configuration P (members in rep coordinates) for T (already in those coordinates).
static void ctl_search(const vector<Mask>& members, const vector<u64>& T, int how, CtlResult& R, bool verbose) {
    string key = space_key(T);
    Problem P; build_problem(P, members); make_order(P, how, 1);
    Profile Q = profile(T);
    for (auto& m : Q.r4) { bool in = false; for (int j = 0; j < 13; j++) if (meq(m, P.R[j])) in = true; if (!in) { R.status = 2; R.why = "rank-4 radical not a member"; return; } }
    for (size_t i = 0; i < Q.f6.size(); i++) { if (msub(Q.r6[i], P.D)) R.inD++; else { P.oracle.push_back(Q.f6[i]); P.oracleM.push_back(Q.r6[i]); } }
    map<string, vector<u64>> sol; tl_sol = &sol;
    R.st = run_search(P, empty_state(), verbose ? 400 : 1);
    tl_sol = nullptr;
    R.known = sol.count(key) > 0; R.nsol = (int)sol.size();
    for (auto& kv : sol) { Verdict V = verify_space(kv.second); R.napn += V.apn;
        if (verbose) { printf("  APN=%d dist:", V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("%s\n", kv.first == key ? "  <- known" : ""); } }
    R.status = R.known ? 0 : 3;
}
// Extendable path: a spread through T's radicals, 4 other spread elements as holes, the 13-set carried onto its
// representative by its canonical frame (this tests the classification: it must land on a listed configuration).
static CtlResult ctl_ext(const CtlCtx& X, vector<u64> T, int choice, int how, bool verbose) {
    CtlResult R; Profile Q = profile(T); R.n4 = Q.n4;
    if (Q.n4 < 2 || Q.bad) { R.why = "fewer than 2 rank-4"; return R; }
    if (!verify_space(T).apn) { R.why = "forms not an APN space (line not quadratic APN)"; return R; }
    if (!normalise(T, Q, R.why)) return R;
    Mask Xm = XMASK(), Ym = YMASK();
    vector<int> qN; bool hasI = false;
    for (auto& m : Q.r4) { if (meq(m, Xm) || meq(m, Ym)) continue; int N = matrix_of_graph(m); if (N < 0 || (N != ID4 && vid[N] < 0)) { R.why = "radical not in V5824"; return R; } if (N == ID4) hasI = true; else qN.push_back(N); }
    if (Q.n4 >= 3 && !hasI) { R.why = "normalisation to graph(I) failed"; return R; }
    BS cand; for (int i = 0; i < NW; i++) cand.w[i] = ~0ULL;
    vector<int> cur; for (int N : qN) { cur.push_back(vid[N]); BS t; band(t, cand, adj[vid[N]]); cand = t; }
    bool got = false;
    { mt19937 rng((unsigned)choice * 2654435761u + 17); vector<int> c0 = cur;
      for (int t = 0; t < 50 && !got; t++) { cur = c0; got = find_clique_rand(cur, cand, 14 - (int)cur.size(), rng); }
      if (!got) cur = c0; }
    if (!got && !find_clique(cur, cand, 14 - (int)cur.size())) { R.why = "no spread through the rank-4 radicals"; return R; }
    vector<Mask> spreadAll = { Xm, Ym, graph_mask(ID4) }; for (int v : cur) spreadAll.push_back(graph_mask(VN[v]));
    vector<int> others; for (int i = 0; i < 17; i++) { bool isQ = false; for (auto& r : Q.r4) if (meq(r, spreadAll[i])) isQ = true; if (!isQ) others.push_back(i); }
    // choice -> a 4-subset of `others`
    mt19937 rng((unsigned)choice); vector<int> o = others; shuffle(o.begin(), o.end(), rng);
    set<int> holes(o.begin(), o.begin() + 4);
    vector<Mask> members; for (int i = 0; i < 17; i++) if (!holes.count(i)) members.push_back(spreadAll[i]);
    Canon K = canon13(members);
    auto it = X.repOf.find(K.key);
    if (it == X.repOf.end()) { R.status = 2; R.why = "configuration matches no representative"; return R; }
    R.rep = it->second;
    if (X.cf[R.rep].nonext) { R.status = 2; R.why = "extendable configuration matched a non-extendable rep"; return R; }
    change_coords(T, K.frame);
    vector<Mask> repM = members_of(X.cf[R.rep].Ns);
    { set<Mask> a, b; for (auto& m : members) a.insert(mask_image(K.frame, m)); for (auto& m : repM) b.insert(m); if (a != b) { R.status = 2; R.why = "frame does not carry the members onto the rep"; return R; } }
    ctl_search(repM, T, how, R, verbose);
    if (verbose) printf("ext control -> rep %d; rank-6 radicals %d inside D, found %d spaces (%d APN), known: %s\n", R.rep, R.inD, R.nsol, R.napn, R.known ? "YES" : "NO");
    return R;
}
// Non-extendable path: T's radicals embedded into the members of a non-extendable representative.
static int g_forceRep = -1; // embed into this representative (any type) instead of cycling the non-extendable ones
static CtlResult ctl_non(const CtlCtx& X, vector<u64> T, int choice, int how, bool verbose) {
    CtlResult R; Profile Q = profile(T); R.n4 = Q.n4;
    if (Q.n4 < 2 || Q.bad) { R.why = "fewer than 2 rank-4"; return R; }
    if (!verify_space(T).apn) { R.why = "forms not an APN space (line not quadratic APN)"; return R; }
    if (X.nonextIdx.empty()) { R.why = "no non-extendable configurations"; return R; }
    if (!normalise(T, Q, R.why)) return R;
    Mask Xm = XMASK(), Ym = YMASK();
    int N4 = -1; for (auto& m : Q.r4) { if (meq(m, Xm) || meq(m, Ym)) continue; int N = matrix_of_graph(m); if (N != ID4) N4 = N; }
    if (Q.n4 > 4) { R.why = "more than 4 rank-4"; return R; }
    mt19937 rng((unsigned)choice * 7919u + 3);
    R.rep = g_forceRep >= 0 ? g_forceRep : X.nonextIdx[choice % X.nonextIdx.size()];
    vector<Mask> M = members_of(X.cf[R.rep].Ns);
    vector<array<u8, 4>> bas(13); vector<vector<u8>> pts(13); for (int i = 0; i < 13; i++) { basis_of_mask(M[i], 4, bas[i].data()); pts[i] = mpoints(M[i]); }
    vector<array<int, 3>> tr; for (int a = 0; a < 13; a++) for (int b = 0; b < 13; b++) for (int c = 0; c < 13; c++) if (a != b && a != c && b != c) tr.push_back({ a, b, c });
    shuffle(tr.begin(), tr.end(), rng);
    for (auto& t : tr) {
        Lin Cm; if (!frame_of(bas, pts, M, t[0], t[1], t[2], Cm)) { R.why = "frame"; return R; }
        Lin D;
        if (Q.n4 == 2) D = diag_ab(GL4[rng() % GL4.size()], GL4[rng() % GL4.size()]);
        else if (Q.n4 == 3) D = diag_hh(GL4[rng() % GL4.size()]);
        else {
            int tv = vid[N4]; if (tv < 0) { R.why = "4th radical not in V5824"; return R; }
            vector<int> ds; for (int d = 0; d < 13; d++) if (d != t[0] && d != t[1] && d != t[2]) { int v = vertex_in_frame(Cm, bas[d]); if (cls[v] == cls[tv]) ds.push_back(v); }
            if (ds.empty()) continue;
            int dv = ds[rng() % ds.size()];
            // h N4 h^-1 = N_d: h = h0(d)^-1 z h0(N4), z in the centraliser of the class rep
            int c = cls[tv]; u16 z = GL4[cent[c][rng() % cent[c].size()]];
            u16 h = mm4(inverse4(GL4[toRep[dv]]), mm4(z, GL4[toRep[tv]]));
            if (vid[mm4(mm4(h, (u16)N4), inverse4(h))] != dv) { R.why = "conjugator"; return R; }
            D = diag_hh(h);
        }
        Lin g = lin_comp(lin_inv(Cm), D);
        change_coords(T, g);
        ctl_search(M, T, how, R, verbose);
        if (verbose) printf("nonext control -> rep %d; rank-6 radicals %d inside D, found %d spaces (%d APN), known: %s\n", R.rep, R.inD, R.nsol, R.napn, R.known ? "YES" : "NO");
        return R;
    }
    R.status = 4; R.why = "no member 4-tuple of the 4th radical's class"; return R;
}
static int mode_ctl(const char* conffile, const char* file, int line, int choice, bool non) {
    init_conj();
    CtlCtx X = ctl_ctx(conffile);
    vector<u64> T; if (!load_line(file, line, T)) { printf("cannot load line %d\n", line); return 1; }
    g_logp1 = false;
    CtlResult R = non ? ctl_non(X, T, choice, 0, true) : ctl_ext(X, T, choice, 0, true);
    if (R.status == 1 || R.status == 2 || R.status == 4) printf("ctl: %s\n", R.why.c_str());
    print_stats("ctl", R.st);
    return R.status;
}
static int mode_ctlbatch(const char* conffile, const char* file, int minr4, long long maxn, int nch, const char* outfile) {
    init_conj();
    CtlCtx X = ctl_ctx(conffile);
    g_logp1 = false;
    FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); return 1; }
    FILE* of = fopen(outfile, "a");
    vector<pair<long long, vector<u64>>> batch; long long ln = 0, taken = 0;
    long long tot[2] = { 0 }, ok[2] = { 0 }, notfound[2] = { 0 }, nospread = 0, fail[2] = { 0 }, extra[2] = { 0 }, nonapn[2] = { 0 }, noembed = 0;
    map<int, int> repHits[2]; map<int, int> byN4[2];
    static char buf[1 << 16];
    auto flush = [&]() {
        int n = (int)batch.size(); vector<vector<CtlResult>> res(n, vector<CtlResult>(2 * nch));
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n; i++) for (int c = 0; c < 2 * nch; c++) {
            int ch = (int)(batch[i].first * 131 + c);
            res[i][c] = (c < nch) ? ctl_ext(X, batch[i].second, ch, 0, false) : ctl_non(X, batch[i].second, ch, 0, false);
        }
        for (int i = 0; i < n; i++) for (int c = 0; c < 2 * nch; c++) {
            int t = c < nch ? 0 : 1; auto& R = res[i][c];
            if (R.status == 4) { noembed++; continue; }
            if (R.why == "no spread through the rank-4 radicals") { nospread++; continue; }
            tot[t]++;
            if (R.status == 0) { ok[t]++; repHits[t][R.rep]++; byN4[t][R.n4]++; extra[t] += R.nsol - 1; nonapn[t] += R.nsol - R.napn; }
            else if (R.status == 3) { notfound[t]++; fprintf(of, "NOTFOUND %s line %lld choice %d rep %d\n", t ? "nonext" : "ext", batch[i].first, c, R.rep); }
            else { fail[t]++; fprintf(of, "FAIL %s line %lld choice %d: %s\n", t ? "nonext" : "ext", batch[i].first, c, R.why.c_str()); }
            if (R.status == 0 && R.nsol > 1) fprintf(of, "EXTRA %s line %lld choice %d rep %d: %d spaces, %d APN\n", t ? "nonext" : "ext", batch[i].first, c, R.rep, R.nsol, R.napn);
        }
        fflush(of); batch.clear();
        printf("  %lld lines read, %lld taken; ext %lld: re-found %lld, NOT %lld, fail %lld; nonext %lld: re-found %lld, NOT %lld, fail %lld; no spread %lld, no embedding %lld\n",
            ln, taken, tot[0], ok[0], notfound[0], fail[0], tot[1], ok[1], notfound[1], fail[1], nospread, noembed); fflush(stdout);
    };
    while (fgets(buf, sizeof buf, fp)) {
        ln++; u8 F[256]; if (!read_lut(buf, F)) continue;
        vector<u64> B = forms_of_lut(F);
        int n4 = 0; for (int c = 1; c < 256; c++) { u64 f = 0; for (int i = 0; i < 8; i++) if (c >> i & 1) f ^= B[i]; if (rank8(f) == 4) n4++; }
        if (n4 < minr4) continue;
        batch.push_back({ ln, B }); taken++;
        if ((int)batch.size() >= 256) flush();
        if (taken >= maxn) break;
    }
    if (!batch.empty()) flush();
    fclose(fp);
    for (int t = 0; t < 2; t++) {
        printf("ctlbatch %s %s: %lld controls; re-found %lld; NOT found %lld; setup failures %lld; other spaces found %lld (%lld not APN)\n", file, t ? "nonext" : "ext", tot[t], ok[t], notfound[t], fail[t], extra[t], nonapn[t]);
        printf("  re-found by rank-4 count:"); for (auto& kv : byN4[t]) printf(" %d:%d", kv.first, kv.second); printf("\n");
        printf("  reps hit: %zu\n", repHits[t].size());
        fprintf(of, "SUMMARY %s %s controls %lld refound %lld notfound %lld fail %lld extra %lld nonapn %lld reps %zu\n", file, t ? "nonext" : "ext", tot[t], ok[t], notfound[t], fail[t], extra[t], nonapn[t], repHits[t].size());
        fprintf(of, "REPS %s", t ? "nonext" : "ext"); for (auto& kv : repHits[t]) fprintf(of, " %d:%d", kv.first, kv.second); fprintf(of, "\n");
    }
    printf("skipped: no spread %lld, no member 4-tuple of the class (k = 4, nonext) %lld\n", nospread, noembed);
    fclose(of);
    return (notfound[0] || notfound[1] || fail[0] || fail[1]) ? 3 : 0;
}

// ================================================================================ selftest
static int mode_selftest(const char* conffile) {
    init_conj();
    auto cf = read_confs(conffile); int bad = 0;
    printf("selftest: ND4 %zu (28), ND6 %zu (13888), GL4 %zu (20160), V5824 %zu\n", ND4.size(), ND6.size(), GL4.size(), VN.size());
    bad += ND4.size() != 28 || ND6.size() != 13888 || GL4.size() != 20160 || VN.size() != 5824;
    long long checked = 0; mt19937 rng(12345); map<size_t, int> lines;
    for (auto& c : cf) {
        vector<Mask> M = conf_members(c);
        Problem P; build_problem(P, M);
        bool ok = mcount(P.D) == 60;
        for (int a = 0; a < 10; a++) { ok &= inv4(c.Ns[a]) && inv4(c.Ns[a] ^ ID4); for (int b = a + 1; b < 10; b++) ok &= inv4(c.Ns[a] ^ c.Ns[b]); }
        lines[P.lineM.size()]++;
        for (size_t li = 0; li < P.lineM.size(); li++) ok &= msub(P.lineM[li], P.D) && mcount(P.lineM[li]) == 3;
        for (int j = 0; j < 13; j++) { set<u64> s(P.f4[j].begin(), P.f4[j].end()); ok &= s.size() == 28; for (u64 f : P.f4[j]) { Mask m; ok &= radical(f, m) == 15 && meq(m, P.R[j]); checked++; } }
        if (c.idx % 16 == 0) for (size_t li = 0; li < P.lineM.size(); li++) {
            set<u64> s(P.lineF[li].begin(), P.lineF[li].end()); ok &= s.size() == 13888;
            for (u64 f : P.lineF[li]) { Mask m; ok &= rank8(f) == 6 && radical(f, m) == 3 && meq(m, P.lineM[li]); checked++; }
        }
        // canonical form: invariant under a random GL(8,2) image, and the stored aut
        { Lin g = random_gl8(rng); vector<Mask> Mg; for (auto& m : M) Mg.push_back(mask_image(g, m)); shuffle(Mg.begin(), Mg.end(), rng);
          Canon K = canon13(Mg); vector<u16> key; for (u16 N : c.Ns) key.push_back((u16)vid[N]); sort(key.begin(), key.end());
          ok &= K.key == key && K.aut == c.aut;
          set<Mask> a, b; for (auto& m : Mg) a.insert(mask_image(K.frame, m)); for (auto& m : M) b.insert(m); ok &= a == b; }
        if (!ok) { printf("config %d: FAIL\n", c.idx); bad++; }
    }
    printf("lines in D over the configurations:"); for (auto& kv : lines) printf(" %zu:%d", kv.first, kv.second); printf("\n");
    { Problem P; build_problem(P, conf_members(cf.back())); long long cnt = 0;
#pragma omp parallel for reduction(+:cnt) schedule(static, 1 << 16)
      for (int c = 0; c < (1 << 28); c++) { u64 f = form28((u32)c); if (rank8(f) == 4) { Mask m; radical(f, m); if (meq(m, P.R[5])) cnt++; } }
      printf("brute force: %lld forms with radical exactly R_5 of the last config (expect 28)\n", cnt); bad += cnt != 28; }
    printf("selftest: %zu configurations, %lld generated forms checked, %d failures\n", cf.size(), checked, bad);
    return bad;
}

int main(int argc, char** argv) {
    if (argc < 2) { printf("modes: holes <out> | configs <holes> <kurz.txt> <out> | search <confs> <prog> [asc|desc|rand <seed>] [first last] | searchsp <spreads.txt> <prog> ... | selftest <confs> | ctl <confs> <lut> <line> <choice> [non] | ctlbatch <confs> <lut> <minr4> <maxn> <nch> <out>\n"); return 1; }
    string mode = argv[1];
    ND4 = nondeg_codes(4); ND6 = nondeg_codes(6);
    init_V();
    printf("V5824: %zu vertices, %zu classes:", VN.size(), classSize.size()); for (size_t c = 0; c < classSize.size(); c++) printf(" %d (rep %04x)", classSize[c], VN[classRep[c]]); printf("\n");
    fflush(stdout);
    if (mode == "holes") return mode_holes(argv[2]);
    if (mode == "holesk") return mode_holesk(atoi(argv[2]));
    if (mode == "configs") return mode_configs(argv[2], argv[3], argv[4]);
    init_code28();
    if (mode == "covers") return mode_covers(argv[2], argv[3]);
    if (mode == "search" || mode == "searchsp" || mode == "searchmem") {
        g_input = mode == "searchsp" ? 1 : mode == "searchmem" ? 2 : 0;
        int how = 0; unsigned seed = 1; int first = 0, last = INT_MAX, ai = 4;
        if (argc > ai && !strcmp(argv[ai], "desc")) { how = 1; ai++; }
        else if (argc > ai && !strcmp(argv[ai], "rand")) { how = 2; seed = (unsigned)atoi(argv[ai + 1]); ai += 2; }
        else if (argc > ai && !strcmp(argv[ai], "asc")) ai++;
        if (argc > ai + 1) { first = atoi(argv[ai]); last = atoi(argv[ai + 1]); }
        return mode_search(argv[2], mode == "searchsp", argv[3], how, seed, first, last);
    }
    if (mode == "selftest") return mode_selftest(argv[2]);
    if (mode == "mapsp") { // spread file x 14 -> which listed configuration each is
        init_conj(); CtlCtx X = ctl_ctx(argv[2]); auto sp = read_spreads_as_confs(argv[3]);
        int n = (int)sp.size(); vector<int> hit(n, -1);
#pragma omp parallel for schedule(dynamic, 4)
        for (int i = 0; i < n; i++) { Canon K = canon13(conf_members(sp[i])); auto it = X.repOf.find(K.key); if (it != X.repOf.end()) hit[i] = it->second; }
        set<int> got; int miss = 0, nonext = 0;
        for (int i = 0; i < n; i++) { if (hit[i] < 0) miss++; else { got.insert(hit[i]); nonext += X.cf[hit[i]].nonext; } }
        int ext = 0; for (auto& c : X.cf) ext += !c.nonext;
        printf("mapsp: %d configurations; unmatched %d; matched a non-extendable rep %d; distinct reps hit %zu (extendable reps listed: %d)\n", n, miss, nonext, got.size(), ext);
        for (int i : { 737, 937 }) if (i < n) printf("  spread-file configuration %d -> rep %d\n", i, hit[i]);
        return miss || nonext;
    }
    if (mode == "ctl") return mode_ctl(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argc > 6 && !strcmp(argv[6], "non"));
    if (mode == "ctlbatch") return mode_ctlbatch(argv[2], argv[3], atoi(argv[4]), atoll(argv[5]), atoi(argv[6]), argv[7]);
    if (mode == "ctlreps") { // ctlreps <confs> <lut> <first line> <n functions with k >= 3> <rep>...: embedding controls into given reps
        init_conj(); CtlCtx X = ctl_ctx(argv[2]); g_logp1 = false;
        vector<vector<u64>> fs; vector<int> lns; FILE* fp = fopen(argv[3], "r"); static char buf[1 << 16]; int ln = 0, first = atoi(argv[4]), want = atoi(argv[5]);
        while (fgets(buf, sizeof buf, fp) && (int)fs.size() < want) { ln++; if (ln < first) continue; u8 F[256]; if (!read_lut(buf, F)) continue; auto B = forms_of_lut(F); if (profile(B).n4 >= 3) { fs.push_back(B); lns.push_back(ln); } }
        fclose(fp);
        int bad = 0;
        for (int a = 6; a < argc; a++) {
            g_forceRep = atoi(argv[a]); int ok = 0, tot = 0, extra = 0;
            vector<CtlResult> res(fs.size());
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < (int)fs.size(); i++) res[i] = ctl_non(X, fs[i], lns[i] * 31 + g_forceRep, 0, false);
            for (auto& R : res) { if (R.status == 4) continue; tot++; ok += R.status == 0; extra += R.status == 0 ? R.nsol - 1 : 0; if (R.status != 0) printf("  rep %d: status %d %s\n", g_forceRep, R.status, R.why.c_str()); }
            printf("rep %d (%s): %d of %d re-found, %d other spaces\n", g_forceRep, X.cf[g_forceRep].nonext ? "nonext" : "ext", ok, tot, extra); fflush(stdout);
            bad += tot - ok;
        }
        return bad ? 3 : 0;
    }
    return 1;
}
