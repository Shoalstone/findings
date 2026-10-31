// eight_check.cpp -- independent check of Open Problem VI.1 (arXiv 2510.12008) for 8 rank-4 components:
// [0^202, 2^45, 4^8]. Subagent (Claude Opus 5.5) for Rowen, 2026-10-30 ("search B" for i = 9). Written without reading
// pspread*.cpp, cover*.cpp, cover_prep.py, README's "Nine"/"Eight" parts, runs/i8/, runs/i9/ or notes/. Built on
// nine_check.cpp (search B for 9 rank-4): everything carried over and generalised to 8 members; new here: `count5`
// (5-cliques of V5824 counted by plain enumeration, no group action) and the configuration mass in `configs`.
//
// Coordinates: v = (x, y), x = low nibble, y = high nibble. X = {y = 0}, Y = {x = 0}, Z = graph(I) = {y = x}.
// A 4-space complementary to Y is graph(N) = {(x, Nx)}; complementary to X iff N invertible, to Z iff N + I is.
// V5824 = {N : N, N + I invertible}; two graphs are disjoint iff N + N' is invertible ("adjacent").
//
// parta   : every 5-clique of V5824 (= 8-partial spread with X, Y, Z among its members) with minimal class c, up to
//           conjugation by the minimal-class trick, each canonicalised (one class per run, or a slice of a class).
// merge   : the per-class files merged; per-orbit identity; per-orbit (not per-clique) extension invariants; the link
//           to nine-check's 5663 orbits (adding a solid), orbit by orbit.
// configs : hole sets of the orbits up to GL(8,2) (Kurz Table 9.3, r = 8: 9316 and 181), exact covers of each,
//           independent sources; writes the configuration list.
// search  : per configuration, point-branching over the forms (radicals fixed).
// ctl / ctlbatch / ctlreps : known functions through the same search path.
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
#include <memory>
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

// ================================================================================ canonical form of a partial spread
// As thirteen_check's canon13, for n >= 4 members (k = n - 3 graph vertices per frame). Frames: an ordered triple
// (A, B, C) of members sent to (X, Y, Z) by Cmap (basis a_i of A -> e_i, the B-part b_i of a_i in B + C -> e_{4+i});
// every g in GL(8,2) with g(A, B, C) = (X, Y, Z) is diag(h, h) o Cmap. Key = the least sorted vertex list over the
// frames whose ordering has the least sorted conjugacy-class vector (intrinsic, so the key is an orbit invariant; an
// image of P, so equal keys mean equivalent P). Only h moving a least-class member onto its class rep can give the
// least list (the rep is the smallest index of its class). Frames giving the key <-> Aut(P), one to one.
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
typedef array<u16, 14> VArr;
static Canon canonP(const vector<Mask>& M) {
    int n = (int)M.size(), k = n - 3; if (n < 4 || n > 17) { fprintf(stderr, "canonP: %d members\n", n); exit(1); }
    vector<array<u8, 4>> bas(n); vector<vector<u8>> pts(n);
    for (int i = 0; i < n; i++) { basis_of_mask(M[i], 4, bas[i].data()); pts[i] = mpoints(M[i]); }
    struct Ord { Lin Cm; VArr v; };
    vector<Ord> keep; array<int, 14> bestCls; bestCls.fill(INT_MAX);
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) for (int c = 0; c < n; c++) {
        if (a == b || a == c || b == c) continue;
        Ord o; if (!frame_of(bas, pts, M, a, b, c, o.Cm)) { fprintf(stderr, "canonP: members not complementary\n"); exit(1); }
        int t = 0; array<int, 14> cv; cv.fill(INT_MAX);
        for (int d = 0; d < n; d++) { if (d == a || d == b || d == c) continue; int v = vertex_in_frame(o.Cm, bas[d]); if (v < 0) { fprintf(stderr, "canonP: not in V5824\n"); exit(1); } o.v[t] = (u16)v; cv[t] = cls[v]; t++; }
        sort(cv.begin(), cv.begin() + k);
        if (cv < bestCls) { bestCls = cv; keep.clear(); }
        if (cv == bestCls) keep.push_back(o);
    }
    Canon R; R.key.assign(k, 0xffff); VArr img;
    int c0 = bestCls[0], r0 = classRep[c0];
    for (auto& o : keep) for (int e = 0; e < k; e++) {
        if (cls[o.v[e]] != c0) continue;
        u16 h0 = GL4[toRep[o.v[e]]];
        for (int z : cent[c0]) {
            int hi = gl4idx[mm4(GL4[z], h0)];
            for (int j = 0; j < k; j++) img[j] = conjTab[hi][o.v[j]];
            sort(img.begin(), img.begin() + k);
            if (img[0] != r0) { fprintf(stderr, "canonP: rep\n"); exit(1); }
            int cmp = 0; for (int j = 0; j < k && !cmp; j++) cmp = img[j] < R.key[j] ? -1 : img[j] > R.key[j] ? 1 : 0;
            if (cmp < 0) { R.key.assign(img.begin(), img.begin() + k); R.aut = 1; R.frame = lin_comp(diag_hh(GL4[hi]), o.Cm); }
            else if (cmp == 0) R.aut++;
        }
    }
    return R;
}
static vector<u16> key_matrices(const vector<u16>& key) { vector<u16> N; for (u16 v : key) N.push_back(VN[v]); return N; }

// k-cliques among a small vertex list (all pairwise adjacent), or whether one exists
static long long kcl(const vector<int>& cand, int need, bool one) {
    if (need == 0) return 1;
    if ((int)cand.size() < need) return 0;
    long long n = 0;
    for (size_t i = 0; i + need <= cand.size(); i++) {
        vector<int> nc; for (size_t j = i + 1; j < cand.size(); j++) if (bget(adj[cand[i]], cand[j])) nc.push_back(cand[j]);
        n += kcl(nc, need - 1, one); if (one && n) return n;
    }
    return n;
}
static vector<int> blist(const BS& b) { vector<int> v; for (int wi = 0; wi < NW; wi++) { u64 w = b.w[wi]; while (w) { v.push_back(wi * 64 + (int)_tzcnt_u64(w)); w &= w - 1; } } return v; }

// ================================================================================ Part A: 5-cliques, all canonicalised
// For class c: r = classRep[c], P = adj(r) & {class >= c}; 4-cliques K in P, so K + r is a 5-clique (= an 8-partial
// spread with Y, X, Z among its members) with minimal class c. Weighted count #5-cliques = sum_c |class c| *
// sum_K 1/m(K + r) (m = number of class-c vertices), kept integral by the factor 2520 (m <= 5). Per clique only cheap
// things: |Q| (Q = the common neighbourhood = every solid disjoint from all 8 members) and the canonical form.
// Everything that needs a search inside Q is done once per orbit, in `merge`. Per orbit the weighted sum must be
// 2520 * 336 * 20160 / |Aut P|. A class can be split into slices (first vertices fi with fi % nslice == slice).
static const int TYPES = 7;
static const char* TNAME[TYPES] = { "spr", "ext13", "ext12", "ext11", "ext10", "ext9", "max" }; // how far the 8 extends: 17, 13, 12, 11, 10, 9, nowhere
static const long long FR8 = 336LL * 20160, FR9 = 504LL * 20160;
static const int KN = 5; // graph members per configuration (the 8 members less X, Y, Z)
typedef array<u16, KN> Key5;
struct OrbAcc { long long w = 0, n = 0, aut = 0; int q = -1; int bad = 0; };
static int g_testFirsts = 0;
static int mode_parta(int c, const char* outfile, const char* maxfile, int slice, int nslice) {
    init_conj();
    if (c < 0 || c >= (int)classSize.size() || nslice < 1 || slice < 0 || slice >= nslice) { printf("bad class / slice\n"); return 1; }
    FILE* of = fopen(outfile, "w"); FILE* mf = fopen(maxfile, "w");
    if (!of || !mf) { printf("cannot open output\n"); return 1; }
    printf("Part A, class %d slice %d of %d: 5-cliques of V5824 (8-partial spreads with Y, X, Z among them) with minimal class %d\n", c, slice, nslice, c); fflush(stdout);
    int r = classRep[c];
    BS P; for (int i = 0; i < NW; i++) P.w[i] = 0;
    for (int v = 0; v < (int)VN.size(); v++) if (cls[v] >= c && bget(adj[r], v)) bset(P, v);
    vector<int> all = blist(P), firsts;
    for (int fi = 0; fi < (int)all.size(); fi++) if (fi % nslice == slice) firsts.push_back(all[fi]);
    if (g_testFirsts > 0 && (int)firsts.size() > g_testFirsts) firsts.resize(g_testFirsts); // quick tests only
    map<Key5, OrbAcc> orb; map<int, long long> qW; long long cw = 0, cl = 0, csq = 0, cmax = 0; CliqueStats tot; mutex mu;
    auto t0 = chrono::steady_clock::now(); int doneFirsts = 0;
#pragma omp parallel for schedule(dynamic, 1)
    for (int fi = 0; fi < (int)firsts.size(); fi++) {
        int v = firsts[fi]; BS nc; above(nc, P, adj[v], v);
        vector<int> cur = { v }; CliqueStats st;
        map<Key5, OrbAcc> lo; long long lw = 0, ll = 0, lsq = 0, lmax = 0; map<int, long long> lq; vector<string> maxLines;
        auto leaf = [&](const vector<int>& K) {
            ll++;
            int m = 1; for (int u : K) m += cls[u] == c;
            long long w = 2520 / m;
            BS Q = adj[r]; for (int u : K) { BS t; band(t, Q, adj[u]); Q = t; }
            int q = bcount(Q);
            lw += w; lsq += w * q; lq[q] += w;
            vector<u16> Ns = { VN[r] }; for (int u : K) Ns.push_back(VN[u]);
            Canon C = canonP(members_of(Ns));
            Key5 key; for (int j = 0; j < KN; j++) key[j] = C.key[j];
            OrbAcc& A = lo[key];
            if (A.n == 0) { A.aut = C.aut; A.q = q; }
            else if (A.aut != C.aut || A.q != q) A.bad++;
            A.n++; A.w += w;
            if (!q) { lmax++; char buf[128]; int o = sprintf(buf, "MAXC c%d :", c); for (u16 N : Ns) o += sprintf(buf + o, " %04x", N); maxLines.push_back(buf); }
        };
        cliques(cur, nc, KN - 2, st, leaf);
        lock_guard<mutex> g(mu);
        for (int i = 0; i < 16; i++) tot.nodes[i] += st.nodes[i];
        for (auto& kv : lo) {
            OrbAcc& A = orb[kv.first]; const OrbAcc& B = kv.second;
            if (A.n == 0) { A.aut = B.aut; A.q = B.q; }
            else if (A.aut != B.aut || A.q != B.q) A.bad++;
            A.bad += B.bad; A.n += B.n; A.w += B.w;
        }
        cw += lw; cl += ll; csq += lsq; cmax += lmax; for (auto& kv : lq) qW[kv.first] += kv.second;
        for (auto& s : maxLines) fprintf(mf, "%s\n", s.c_str());
        if (++doneFirsts % 50 == 0) { printf("  %d of %zu first vertices, %lld cliques, %zu orbits, %.0f s\n", doneFirsts, firsts.size(), cl, orb.size(), chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout); }
    }
    fclose(mf);
    double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    long long cs = classSize[c];
    printf("class %d slice %d/%d (size %d, rep %04x, |P| = %d, %zu first vertices): %lld reduced 5-cliques, %lld maximal; %zu orbits; %.0f s\n  nodes by clique size (excl. rep):", c, slice, nslice, classSize[c], VN[r], bcount(P), firsts.size(), cl, cmax, orb.size(), el);
    for (int i = 0; i < KN; i++) printf(" %lld", tot.nodes[i]); printf("\n");
    fprintf(of, "CLASS %d slice %d %d leaves %lld w %lld sumq %lld max %lld secs %.0f\n", c, slice, nslice, cl, cw * cs, csq * cs, cmax, el);
    fprintf(of, "QH %d", c); for (auto& kv : qW) fprintf(of, " %d:%lld", kv.first, kv.second * cs); fprintf(of, "\n");
    for (auto& kv : orb) {
        const OrbAcc& A = kv.second;
        fprintf(of, "ORBP aut %lld w %lld n %lld q %d bad %d :", A.aut, A.w * cs, A.n, A.q, A.bad);
        for (u16 v : kv.first) fprintf(of, " %04x", VN[v]); fprintf(of, "\n");
    }
    fprintf(of, "END\n"); fclose(of);
    return 0;
}

// ---- small clique helpers for the per-orbit invariants
static void mcrec(const vector<int>& cand, int sz, int& best, int cap) {
    if (sz > best) best = sz;
    if (best >= cap) return;
    for (size_t i = 0; i < cand.size(); i++) {
        if (sz + (int)(cand.size() - i) <= best) return;
        vector<int> nc; for (size_t j = i + 1; j < cand.size(); j++) if (bget(adj[cand[i]], cand[j])) nc.push_back(cand[j]);
        mcrec(nc, sz + 1, best, cap);
        if (best >= cap) return;
    }
}
static int max_clique(const vector<int>& cand, int cap) { int best = 0; mcrec(cand, 0, best, cap); return best; }

// earlier checks' configurations (CONF lines: n matrices, a type name), with their canonical keys
struct CPrev { int idx, type; long long aut; vector<u16> Ns; };
static const char* TNAME10[5] = { "spr", "ext13", "ext12", "ext11", "max" }; // ten-check's types
static const char* TNAME11[4] = { "spr", "ext13", "ext12", "max" }; // eleven-check's
static const char* TNAME12[3] = { "spr", "ext13", "max" };          // twelve-check's
static const char* TNAME9[6] = { "spr", "ext13", "ext12", "ext11", "ext10", "max" }; // nine-check's
static vector<CPrev> read_prev(const char* file, int nmat, const char* const* tn, int ntypes) {
    vector<CPrev> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "CONF", 4)) continue;
        CPrev c; char ty[16]; int h;
        if (sscanf(line, "CONF %d %15s hole %d aut %lld", &c.idx, ty, &h, &c.aut) != 4) { printf("bad CONF line in %s\n", file); exit(1); }
        c.type = -1; for (int t = 0; t < ntypes; t++) if (!strcmp(ty, tn[t])) c.type = t;
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; c.Ns.push_back((u16)v); p = e; }
        if ((int)c.Ns.size() != nmat || c.type < 0 || c.idx != (int)out.size()) { printf("bad CONF line in %s\n", file); exit(1); }
        out.push_back(c);
    }
    fclose(fp); return out;
}
static vector<CPrev> read_conf10(const char* file) { return read_prev(file, 7, TNAME10, 5); }
static vector<CPrev> read_conf11(const char* file) { return read_prev(file, 8, TNAME11, 4); }
static vector<CPrev> read_conf12(const char* file) { return read_prev(file, 9, TNAME12, 3); }
static vector<CPrev> read_conf9(const char* file) { return read_prev(file, 6, TNAME9, 6); }
static vector<u16> key_of_Ns(const vector<u16>& Ns) { vector<u16> k; for (u16 N : Ns) k.push_back((u16)vid[N]); sort(k.begin(), k.end()); return k; }

// ================================================================================ merge: orbits, per-orbit invariants, identities
// Identities: (1) sum over 5-cliques of |Q| = 6 * (weighted 6-cliques) = 6 * sum over nine-check's 5663 of
// 504*20160/|Aut|; (2) per orbit, weighted cliques = 336*20160/|Aut P|; (3) sum_P 336*20160/|Aut P| * (spreads through
// P) = C(14,5) * 1128864 (1128864 spreads through X, Y, Z = 14-cliques; any 5 of their 14 graph members);
// (4) per 9-orbit O: sum_P 336*20160/|Aut P| * #{q in Q(P) : P + q ~ O} = 6 * 504*20160/|Aut O|.
struct Orb8 { Key5 key; long long aut = 0, w = 0, n = 0; int q = -1, bad = 0, mc = -1, type = -1, unk = 0, tbad = 0; long long n9 = 0; map<int, int> link; };
static int mode_merge(const char* outfile, const char* conf9, const vector<const char*>& parts) {
    init_conj();
    auto c9 = read_conf9(conf9); map<vector<u16>, int> key9; long long W6 = 0;
    for (auto& c : c9) { key9[key_of_Ns(c.Ns)] = c.idx; if (FR9 % c.aut) { printf("9: aut does not divide\n"); return 1; } W6 += FR9 / c.aut; }
    printf("nine-check's configurations: %zu; sum 504*20160/|Aut| = %lld (its Part A: 47414189312)\n", c9.size(), W6);
    map<Key5, Orb8> orb; long long Wall = 0, SQ = 0, L = 0, M = 0; map<int, long long> qW; long long inconsistent = 0;
    map<int, set<int>> slicesSeen; map<int, int> nslOf; int dupSlice = 0;
    for (const char* pf : parts) {
        FILE* fp = fopen(pf, "r"); if (!fp) { printf("cannot open %s\n", pf); return 1; }
        static char line[4096]; bool end = false;
        while (fgets(line, sizeof line, fp)) {
            if (!strncmp(line, "CLASS", 5)) { int c, sl, nsl; long long l, w, sq, mx;
                if (sscanf(line, "CLASS %d slice %d %d leaves %lld w %lld sumq %lld max %lld", &c, &sl, &nsl, &l, &w, &sq, &mx) != 7) { printf("bad CLASS\n"); return 1; }
                if (nslOf.count(c) && nslOf[c] != nsl) { printf("class %d: slices of different splits\n", c); return 1; }
                nslOf[c] = nsl; if (!slicesSeen[c].insert(sl).second) dupSlice++;
                L += l; Wall += w; SQ += sq; M += mx; }
            else if (!strncmp(line, "QH", 2)) { const char* p = line + 2; char* e; strtol(p, &e, 10); p = e;
                while (true) { long q = strtol(p, &e, 10); if (e == p) break; p = e + 1; long long w = strtoll(p, &e, 10); p = e; qW[(int)q] += w; } }
            else if (!strncmp(line, "ORBP", 4)) {
                long long aut, w, n; int q, bad;
                if (sscanf(line, "ORBP aut %lld w %lld n %lld q %d bad %d", &aut, &w, &n, &q, &bad) != 5) { printf("bad ORBP\n"); return 1; }
                const char* p = strchr(line, ':') + 1; char* e; vector<u16> Ns; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
                if (Ns.size() != KN) { printf("bad ORBP key\n"); return 1; }
                Key5 k; for (int j = 0; j < KN; j++) { if (vid[Ns[j]] < 0) { printf("key not in V5824\n"); return 1; } k[j] = (u16)vid[Ns[j]]; }
                Orb8& A = orb[k];
                if (A.n == 0) { A.key = k; A.aut = aut; A.q = q; } else if (A.aut != aut || A.q != q) inconsistent++;
                A.bad += bad; A.n += n; A.w += w;
            }
            else if (!strncmp(line, "END", 3)) end = true;
        }
        fclose(fp);
        if (!end) { printf("%s: incomplete (no END)\n", pf); return 1; }
    }
    for (int c = 0; c < (int)classSize.size(); c++) if (!nslOf.count(c) || (int)slicesSeen[c].size() != nslOf[c]) { printf("class %d: slices missing\n", c); return 1; }
    if (dupSlice) { printf("%d slices given twice\n", dupSlice); return 1; }
    auto chk = [](long long x) { if (x % 2520) printf("  (not divisible by 2520!)\n"); return x / 2520; };
    long long W5 = chk(Wall);
    printf("reduced 5-cliques %lld, of them maximal (Q empty) %lld; weighted 5-cliques %lld\n", L, M, W5);
    printf("sum over 5-cliques of |Q| = %lld; 6 * (weighted 6-cliques, nine-check) = %lld  %s\n", chk(SQ), 6 * W6, chk(SQ) == 6 * W6 ? "OK" : "MISMATCH");
    printf("by |Q|:"); for (auto& kv : qW) printf(" %d:%lld", kv.first, chk(kv.second)); printf("\n");
    vector<Orb8> O; for (auto& kv : orb) O.push_back(kv.second);
    long long idbad = 0, sumAll = 0, badSum = 0;
    for (auto& A : O) { badSum += A.bad; if (FR8 % A.aut || A.w != 2520 * (FR8 / A.aut)) idbad++; sumAll += FR8 / A.aut; }
    printf("orbits %zu; clique-level inconsistencies (|Aut| or |Q|) %lld + %lld; per-orbit identity (weighted cliques = 336*20160/|Aut|) fails %lld; sum 336*20160/|Aut| = %lld (weighted 5-cliques %lld)  %s\n",
        O.size(), badSum, inconsistent, idbad, sumAll, W5, sumAll == W5 ? "OK" : "MISMATCH");
    fflush(stdout);
    // per orbit, from its key: Q, how far it extends, the number of spreads through it, and which 10-orbit P + q is
    auto t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)O.size(); i++) {
        Orb8& A = O[i];
        vector<u16> Ns; for (u16 v : A.key) Ns.push_back(VN[v]);
        BS Q; for (int j = 0; j < NW; j++) Q.w[j] = ~0ULL;
        for (u16 v : A.key) { BS t; band(t, Q, adj[v]); Q = t; }
        vector<int> qv = blist(Q);
        if ((int)qv.size() != A.q) A.tbad++;
        A.mc = max_clique(qv, 10);
        A.n9 = kcl(qv, 9, false);
        bool has[6] = { false, false, false, false, false, false };
        for (int q : qv) {
            vector<u16> N2 = Ns; N2.push_back(VN[q]);
            Canon K = canonP(members_of(N2));
            auto it = key9.find(K.key);
            if (it == key9.end()) { A.unk++; continue; }
            A.link[it->second]++; has[c9[it->second].type] = true;
        }
        // 8 -> spread (9 more), maximal 13 (5), 12 (4), 11 (3), 10 (2), 9 (1), nothing (0); nothing else is possible
        // (no maximal 14, 15, 16)
        A.type = A.mc >= 9 ? 0 : A.mc == 5 ? 1 : A.mc == 4 ? 2 : A.mc == 3 ? 3 : A.mc == 2 ? 4 : A.mc == 1 ? 5 : A.mc == 0 ? 6 : -1;
        // the same from the 9-orbits reached: an spr / ext13 / ext12 / ext11 / ext10 / max 9 makes P spr / ext13 / ext12 /
        // ext11 / ext10 / ext9
        int t2 = 6; if (!qv.empty()) { t2 = 5; for (int t = 4; t >= 0; t--) if (has[t]) t2 = t; }
        if (A.type < 0 || A.mc > 9 || t2 != A.type) A.tbad++;
        if (A.type == 0 && !A.n9) A.tbad++;
        if (A.type != 0 && A.n9) A.tbad++;
    }
    printf("per-orbit invariants: %.0f s\n", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    long long tbad = 0, unk = 0, n9sum = 0, byT[TYPES] = { 0 }, nT[TYPES] = { 0 }; map<int, int> mcH;
    vector<long long> linkSum(c9.size(), 0);
    for (auto& A : O) {
        tbad += A.tbad; unk += A.unk; n9sum += (FR8 / A.aut) * A.n9; mcH[A.mc]++;
        if (A.type >= 0) { byT[A.type] += FR8 / A.aut; nT[A.type]++; }
        for (auto& kv : A.link) linkSum[kv.first] += (FR8 / A.aut) * kv.second;
    }
    long long linkBad = 0; for (auto& c : c9) if (linkSum[c.idx] != 6 * (FR9 / c.aut)) linkBad++;
    const long long SPR = 2002LL * 1128864; // C(14,5) * spreads through X, Y, Z
    printf("largest clique in Q (how many more solids fit), orbits:"); for (auto& kv : mcH) printf(" %d:%d", kv.first, kv.second); printf("\n");
    printf("orbits by type:"); for (int t = 0; t < TYPES; t++) printf(" %s %lld", TNAME[t], nT[t]); printf("; weighted"); for (int t = 0; t < TYPES; t++) printf(" %lld", byT[t]); printf("\n");
    printf("type / Q / n9 inconsistencies %lld; P + q on no nine-check orbit %lld\n", tbad, unk);
    printf("sum 336*20160/|Aut| * (spreads through P) = %lld; C(14,5) * 1128864 = %lld  %s\n", n9sum, SPR, n9sum == SPR ? "OK" : "MISMATCH");
    printf("per 9-orbit O: sum_P 336*20160/|Aut P| * #{q : P + q in O} = 6 * 504*20160/|Aut O| fails for %lld of %zu\n", linkBad, c9.size());
    sort(O.begin(), O.end(), [](const Orb8& a, const Orb8& b) { return a.type != b.type ? a.type < b.type : a.key < b.key; });
    FILE* of = fopen(outfile, "w");
    for (size_t i = 0; i < O.size(); i++) {
        auto& A = O[i];
        fprintf(of, "ORB %zu %s aut %lld w %lld n %lld q %d mc %d n9 %lld bad %d :", i, A.type >= 0 ? TNAME[A.type] : "??", A.aut, A.w / 2520, A.n, A.q, A.mc, A.n9, A.bad + A.tbad);
        for (u16 v : A.key) fprintf(of, " %04x", VN[v]); fprintf(of, "\n");
        fprintf(of, "LINK %zu :", i); for (auto& kv : A.link) fprintf(of, " %d:%d", kv.first, kv.second); fprintf(of, "\n");
    }
    fprintf(of, "TOTAL orbits %zu weighted %lld", O.size(), W5); for (int t = 0; t < TYPES; t++) fprintf(of, " %s %lld", TNAME[t], nT[t]); fprintf(of, " reduced %lld maximal %lld\n", L, M);
    fclose(of);
    bool ok = !badSum && !inconsistent && !idbad && sumAll == W5 && chk(SQ) == 6 * W6 && !tbad && !unk && n9sum == SPR && !linkBad;

    printf("merge: %s\n", ok ? "all identities hold" : "FAILURES");
    return ok ? 0 : 2;
}

// ================================================================================ configs: hole sets, covers, sources
static const int NFAM = 6;
struct OrbRec { int idx, type; long long aut; int q; vector<u16> Ns; int hole = -1; int src[NFAM] = { 0 }; int srcCov = 0; };
static vector<OrbRec> read_orbs(const char* file) {
    vector<OrbRec> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "ORB ", 4)) continue;
        OrbRec o; char ty[16]; long long w, n;
        if (sscanf(line, "ORB %d %15s aut %lld w %lld n %lld q %d", &o.idx, ty, &o.aut, &w, &n, &o.q) != 6) { printf("bad ORB line\n"); exit(1); }
        o.type = -1; for (int t = 0; t < TYPES; t++) if (!strcmp(ty, TNAME[t])) o.type = t;
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; o.Ns.push_back((u16)v); p = e; }
        if (o.Ns.size() != KN || o.type < 0 || o.idx != (int)out.size()) { printf("bad ORB line\n"); exit(1); }
        out.push_back(o);
    }
    fclose(fp); return out;
}
static Mask holes_of(const vector<Mask>& M) { Mask cov; mclear(cov); for (auto& m : M) mor(cov, m); return mcompl(cov); }
static vector<vector<u16>> read_spread_file(const char* file) {
    vector<vector<u16>> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "SPREAD", 6)) continue;
        const char* p = strchr(line, ':'); if (!p) continue; p++;
        vector<u16> Ns; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
        if (Ns.size() != 14) { printf("bad spread line\n"); exit(1); }
        out.push_back(Ns);
    }
    fclose(fp); return out;
}
// thirteen-check's configs.txt: CONF lines, 10 matrices (members X, Y, Z + 10 graphs); returns the non-extendable ones
static vector<vector<u16>> read_conf13_nonext(const char* file) {
    vector<vector<u16>> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "CONF", 4) || !strstr(line, "nonext")) continue;
        const char* p = strchr(line, ':') + 1; vector<u16> Ns; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; }
        if (Ns.size() != 10) { printf("bad 13 CONF line\n"); exit(1); }
        out.push_back(Ns);
    }
    fclose(fp); return out;
}
// every way of removing k members from each base set (bases in order, k-subsets in lexicographic order)
static vector<vector<Mask>> drop_k(const vector<vector<Mask>>& bases, int k) {
    vector<vector<Mask>> out;
    for (auto& M : bases) {
        int n = (int)M.size(); vector<int> c(k); for (int i = 0; i < k; i++) c[i] = i;
        while (true) {
            vector<Mask> r; int j = 0; for (int i = 0; i < n; i++) { if (j < k && c[j] == i) { j++; continue; } r.push_back(M[i]); }
            out.push_back(r);
            int i = k - 1; while (i >= 0 && c[i] == n - k + i) i--;
            if (i < 0) break;
            c[i]++; for (int t = i + 1; t < k; t++) c[t] = c[t - 1] + 1;
        }
    }
    return out;
}
// the cross-run / source families (each a list of 8-member sets), independent of this file's classification:
//   0: nine-check's 5663 configurations x 9 members removed one at a time (every 8 that extends at all);
//   1: fourteen-check's 73 (spread, X Y Z) reps x C(14,6): the 14 graph members minus 6 (every 8 inside a spread: of
//      the 9 spread members an 8 misses, any 3 can be moved to X, Y, Z);
//   2: thirteen-check's 31 maximal 13s x C(13,5) (every 8 inside a maximal 13);
//   3: twelve-check's 376 configurations x C(12,4) (every 8 inside some 12);
//   4: eleven-check's 817 configurations x C(11,3) (every 8 inside some 11);
//   5: ten-check's 1777 configurations x C(10,2) (every 8 inside some 10).
static vector<vector<Mask>> familyPrev(const vector<CPrev>& c, int k) { vector<vector<Mask>> b; for (auto& x : c) b.push_back(members_of(x.Ns)); return drop_k(b, k); }
static vector<vector<Mask>> family9(const vector<CPrev>& c9) { return familyPrev(c9, 1); }
static vector<vector<Mask>> family10(const vector<CPrev>& c10) { return familyPrev(c10, 2); }
static vector<vector<Mask>> family11(const vector<CPrev>& c11) { return familyPrev(c11, 3); }
static vector<vector<Mask>> family12(const vector<CPrev>& c12) { return familyPrev(c12, 4); }
static vector<vector<Mask>> familySp(const vector<vector<u16>>& sp) {
    vector<vector<Mask>> b; for (auto& Ns : sp) { vector<Mask> M; for (u16 N : Ns) M.push_back(graph_mask(N)); b.push_back(M); }
    return drop_k(b, 6);
}
static vector<vector<Mask>> family13(const vector<vector<u16>>& c13) { vector<vector<Mask>> b; for (auto& Ns : c13) b.push_back(members_of(Ns)); return drop_k(b, 5); }
static const char* FAMNAME[NFAM] = { "nine-check's 5663 x 9", "fourteen-check's spread reps x 3003", "thirteen-check's maximal 13s x 1287", "twelve-check's 376 x 495", "eleven-check's 817 x 165", "ten-check's 1777 x 45" };
static int mode_configs(const char* orbfile, const char* conf9, const char* conf10, const char* conf11, const char* conf12, const char* spreadfile, const char* conf13, const char* outfile) {
    init_conj(); init_solids();
    printf("solids in PG(7,2): %zu (expect 200787)\n", ALLSOLIDS.size()); if (ALLSOLIDS.size() != 200787) return 1;
    auto orbs = read_orbs(orbfile); int n = (int)orbs.size();
    map<vector<u16>, int> byKey; for (int i = 0; i < n; i++) byKey[key_of_Ns(orbs[i].Ns)] = i;
    printf("orbits read: %d\n", n); fflush(stdout);
    // --- the stored key is a canonical form: recompute it (and |Aut|) from a random image
    { int bad = 0; mt19937 rng(99);
      vector<Lin> gs(n); vector<unsigned> seeds(n); for (int i = 0; i < n; i++) { gs[i] = random_gl8(rng); seeds[i] = rng(); }
#pragma omp parallel for schedule(dynamic, 8) reduction(+:bad)
      for (int i = 0; i < n; i++) {
          vector<Mask> M = members_of(orbs[i].Ns), Mg; for (auto& m : M) Mg.push_back(mask_image(gs[i], m));
          mt19937 r2(seeds[i]); shuffle(Mg.begin(), Mg.end(), r2);
          Canon K = canonP(Mg); if (K.key != key_of_Ns(orbs[i].Ns) || K.aut != orbs[i].aut) bad++;
      }
      printf("canonical form recomputed from a random GL(8,2) image of every orbit: %d mismatches\n", bad); fflush(stdout); if (bad) return 1; }
    // --- hole sets
    vector<Mask> H(n); vector<string> hk(n);
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < n; i++) { H[i] = holes_of(members_of(orbs[i].Ns)); hk[i] = hole_key(H[i]); }
    vector<Mask> holeRep; vector<string> holeKey; vector<long long> holeAut; vector<int> holeSolid, holeFirst;
    auto t0 = chrono::steady_clock::now();
    {   map<string, vector<int>> groups; for (int i = 0; i < n; i++) groups[hk[i]].push_back(i);
        printf("distinct hole invariants: %zu\n", groups.size()); fflush(stdout);
        vector<pair<string, vector<int>>> gv(groups.begin(), groups.end());
        vector<vector<pair<int, vector<int>>>> res(gv.size()); // per group: (rep orbit, members)
#pragma omp parallel for schedule(dynamic, 1)
        for (int g = 0; g < (int)gv.size(); g++) {
            vector<pair<int, vector<int>>> reps;
            for (int i : gv[g].second) {
                bool found = false;
                for (auto& rp : reps) if (iso_count(H[i], H[rp.first], true, nullptr)) { rp.second.push_back(i); found = true; break; }
                if (!found) reps.push_back({ i, { i } });
            }
            res[g] = reps;
        }
        for (size_t g = 0; g < gv.size(); g++) for (auto& rp : res[g]) {
            int h = (int)holeRep.size(); holeRep.push_back(H[rp.first]); holeKey.push_back(gv[g].first); holeFirst.push_back(rp.first);
            for (int i : rp.second) orbs[i].hole = h;
        }
        int nh = (int)holeRep.size(); holeAut.assign(nh, 0); holeSolid.assign(nh, 0);
        printf("hole classes by isomorphism: %d (%.0f s); computing |Aut H|\n", nh, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
#pragma omp parallel for schedule(dynamic, 1)
        for (int h = 0; h < nh; h++) { holeAut[h] = iso_count(holeRep[h], holeRep[h], false, nullptr); holeSolid[h] = contains_solid(holeRep[h]); }
    }
    int nh = (int)holeRep.size(), nhMax = 0;
    for (int h = 0; h < nh; h++) nhMax += !holeSolid[h];
    printf("hole sets up to GL(8,2): %d, containing no solid (maximal partial spreads): %d   [Kurz Table 9.3, r = 8: 9316 and 181] (%.0f s)\n", nh, nhMax, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    { int bad = 0; for (auto& o : orbs) if ((o.type == TYPES - 1) != !holeSolid[o.hole]) bad++; printf("orbits whose type disagrees with 'hole set contains a solid': %d\n", bad); }
    // --- exact covers of each hole set's complement by 8 solids: every cover must canonicalise onto an orbit with
    // that hole set, all of them reached, and sum over those orbits of |Aut H| / |Aut P| = number of covers
    vector<vector<int>> orbsOfHole(nh); for (int i = 0; i < n; i++) orbsOfHole[orbs[i].hole].push_back(i);
    struct CovRes { int ns = 0, m = 0; long long unk = 0, wrongHole = 0, expect = 0; vector<int> hits; bool ok = false; double secs = 0; };
    vector<CovRes> cr(nh);
    t0 = chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int h = 0; h < nh; h++) {
        auto t1 = chrono::steady_clock::now();
        CovRes& R = cr[h];
        auto cv = covers_of(holeRep[h], R.ns);
        R.m = (int)cv.size(); set<int> got;
        for (int i = 0; i < R.m; i++) {
            Canon K = canonP(cv[i]); auto it = byKey.find(K.key);
            if (it == byKey.end()) R.unk++; else if (orbs[it->second].hole != h) R.wrongHole++; else { got.insert(it->second); R.hits.push_back(it->second); }
        }
        for (int i : orbsOfHole[h]) R.expect += holeAut[h] % orbs[i].aut ? -1000000000LL : holeAut[h] / orbs[i].aut;
        R.ok = !R.unk && !R.wrongHole && got.size() == orbsOfHole[h].size() && R.expect == R.m;
        R.secs = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
    }
    long long covBad = 0, covAll = 0;
    for (int h = 0; h < nh; h++) {
        CovRes& R = cr[h]; covAll += R.m; covBad += !R.ok;
        for (int i : R.hits) orbs[i].srcCov++;
        map<int, int> types; for (int i : orbsOfHole[h]) types[orbs[i].type]++;
        printf("hole set %d: |Aut H| %lld, solid inside: %s, %d solids in the complement, %d covers, %zu orbits (", h, holeAut[h], holeSolid[h] ? "yes" : "no", R.ns, R.m, orbsOfHole[h].size());
        for (auto& kv : types) printf(" %s:%d", TNAME[kv.first], kv.second);
        printf(" ); sum |Aut H|/|Aut P| = %lld; on unknown orbits %lld, on another hole set %lld  %s [%.1f s]\n    %s\n", R.expect, R.unk, R.wrongHole, R.ok ? "OK" : "MISMATCH", R.secs, holeKey[h].c_str());
    }
    printf("hole sets with a cover mismatch: %lld; covers in all %lld (%.0f s)\n", covBad, covAll, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    { int unreached = 0; for (auto& o : orbs) unreached += !o.srcCov; printf("orbits reached by no cover: %d\n", unreached); covBad += unreached; }
    fflush(stdout);
    // --- independent sources
    auto run_family = [&](const vector<vector<Mask>>& fam, int which) {
        int m = (int)fam.size(); vector<int> hit(m, -1);
#pragma omp parallel for schedule(dynamic, 16)
        for (int i = 0; i < m; i++) { Canon K = canonP(fam[i]); auto it = byKey.find(K.key); if (it != byKey.end()) hit[i] = it->second; }
        set<int> got; int unk = 0; map<int, int> types;
        for (int i = 0; i < m; i++) { if (hit[i] < 0) { unk++; continue; } got.insert(hit[i]); types[orbs[hit[i]].type]++; orbs[hit[i]].src[which]++; }
        printf("%s: %d configurations; on no listed orbit %d; by type:", FAMNAME[which], m, unk); for (auto& kv : types) printf(" %s:%d", TNAME[kv.first], kv.second);
        int reach[TYPES] = { 0 }, all[TYPES] = { 0 }; for (int i = 0; i < n; i++) { all[orbs[i].type]++; reach[orbs[i].type] += got.count(i) > 0; }
        printf("; orbits reached:"); for (int t = 0; t < TYPES; t++) printf(" %s %d/%d", TNAME[t], reach[t], all[t]); printf("\n"); fflush(stdout);
        return unk;
    };
    int srcUnk = 0;
    srcUnk += run_family(family9(read_conf9(conf9)), 0);
    srcUnk += run_family(familySp(read_spread_file(spreadfile)), 1);
    srcUnk += run_family(family13(read_conf13_nonext(conf13)), 2);
    srcUnk += run_family(family12(read_conf12(conf12)), 3);
    srcUnk += run_family(family11(read_conf11(conf11)), 4);
    srcUnk += run_family(family10(read_conf10(conf10)), 5);
    // --- output
    FILE* of = fopen(outfile, "w"); int idx = 0;
    for (int t = 0; t < TYPES; t++) for (auto& o : orbs) if (o.type == t) {
        fprintf(of, "CONF %d %s hole %d aut %lld src", idx++, TNAME[t], o.hole, o.aut);
        for (int k = 0; k < NFAM; k++) fprintf(of, " %d", o.src[k]); fprintf(of, " %d :", o.srcCov);
        for (u16 N : o.Ns) fprintf(of, " %04x", N); fprintf(of, "\n");
    }
    // HOLE lines carry |Aut H| and the number of covers, for the exact masses (eight_check_mass.py)
    for (int h = 0; h < nh; h++) fprintf(of, "HOLE %d aut %lld solid %d covers %d solids %d | %s\n", h, holeAut[h], holeSolid[h], cr[h].m, cr[h].ns, holeKey[h].c_str());

    fclose(of);
    printf("written %d configurations to %s\n", idx, outfile);
    return (covBad || srcUnk) ? 2 : 0;
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

// The problem: 8 rank-4 radicals (members), holes D = the rest (135 points). Relaxation for controls: oracle rank-6
// forms. Forms on the lines of D are generated on first use (most configurations never branch in D).
struct Problem {
    int nr = 0; Mask R[16]; int radOf[256]; Mask D;
    vector<u64> f4[16];
    vector<Mask> lineM; vector<array<u8, 6>> lineAnn; mutable vector<vector<u64>> lineF; unique_ptr<once_flag[]> lineOnce; vector<int> linesAt[256];
    vector<int> order;
    vector<u64> oracle; vector<Mask> oracleM;
};
static void build_problem(Problem& P, const vector<Mask>& M) {
    P.nr = (int)M.size(); if (P.nr > 16) { fprintf(stderr, "too many members\n"); exit(1); }
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
    P.lineM.clear(); P.lineAnn.clear(); for (int v = 0; v < 256; v++) P.linesAt[v].clear();
    for (int a = 1; a < 256; a++) for (int b = a + 1; b < 256; b++) if (mget(P.D, a) && mget(P.D, b) && mget(P.D, a ^ b)) {
        Mask m; mclear(m); mset(m, a); mset(m, b); mset(m, a ^ b);
        if (!seen.insert(m).second) continue;
        int li = (int)P.lineM.size(); P.lineM.push_back(m);
        array<u8, 6> ann; annihilator(m, 6, ann.data()); P.lineAnn.push_back(ann);
        for (int v : { a, b, a ^ b }) P.linesAt[v].push_back(li);
    }
    P.lineF.assign(P.lineM.size(), {});
    P.lineOnce.reset(new once_flag[P.lineM.size() + 1]);
}
static const vector<u64>& line_forms(const Problem& P, int li) {
    call_once(P.lineOnce[li], [&]() { vector<u64> fs; forms_from(P.lineAnn[li].data(), 6, ND6, fs); P.lineF[li].swap(fs); });
    return P.lineF[li];
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
    else for (int li : P.linesAt[a]) if (mdisj(P.lineM[li], cov)) for (u64 f : line_forms(P, li)) cb(f);
    for (size_t i = 0; i < P.oracle.size(); i++) if (mget(P.oracleM[i], a) && !msub(P.oracleM[i], P.D) && mdisj(P.oracleM[i], cov)) cb(P.oracle[i]);
}
// p1 = nodes where every point off D is covered (all 8 rank-4 forms present); triesD = candidate forms tried at
// points of D (the rank-6 branching, which only happens at p1 nodes)
struct Stats { long long nodes[10] = { 0 }, tries[10] = { 0 }, triesD[10] = { 0 }, bentext = 0, full = 0, p1[10] = { 0 }; };
static void add_stats(Stats& a, const Stats& b) { for (int i = 0; i < 10; i++) { a.nodes[i] += b.nodes[i]; a.tries[i] += b.tries[i]; a.triesD[i] += b.triesD[i]; a.p1[i] += b.p1[i]; } a.bentext += b.bentext; a.full += b.full; }
static mutex solmu;
static map<string, vector<u64>> solutions;
static map<string, vector<u64>> p1spaces; // distinct p1 spans (capped)
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
    State T; bool inD = P.radOf[a] < 0;
    for_each_candidate(P, a, S.cov, [&](u64 f) { st.tries[S.d]++; st.triesD[S.d] += inD; if (try_add(P, S, f, T)) dfs(P, T, st); });
}
static void expand(const Problem& P, const State& S, Stats& st, vector<State>& out) {
    int a = first_uncovered(P, S);
    note_node(P, S, a, st);
    if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, st); return; }
    if (S.d == 8) return;
    State T; bool inD = P.radOf[a] < 0;
    for_each_candidate(P, a, S.cov, [&](u64 f) { st.tries[S.d]++; st.triesD[S.d] += inD; if (try_add(P, S, f, T)) out.push_back(T); });
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
    if (how == 3) { for (int v : on) P.order.push_back(v); for (int v : off) P.order.push_back(v); return; } // controls only: points of D first
    for (int v : off) P.order.push_back(v); for (int v : on) P.order.push_back(v);
}
struct Conf { int idx; int type; bool nonext; int hole; long long aut; vector<u16> Ns; vector<Mask> mem; };
static vector<Conf> read_confs(const char* file) {
    vector<Conf> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024];
    while (fgets(line, sizeof line, fp)) {
        Conf c; char ty[16];
        if (strncmp(line, "CONF", 4)) continue;
        if (sscanf(line, "CONF %d %15s hole %d aut %lld", &c.idx, ty, &c.hole, &c.aut) != 4) continue;
        c.type = -1; for (int t = 0; t < TYPES; t++) if (!strcmp(ty, TNAME[t])) c.type = t;
        c.nonext = c.type != 0;
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; c.Ns.push_back((u16)v); p = e; }
        if (c.Ns.size() != KN || c.type < 0) { printf("bad CONF line\n"); exit(1); }
        out.push_back(c);
    }
    fclose(fp); return out;
}
// Cross-run inputs that don't use this classification (families as in `configs`; type marks the family, not the orbit):
// 1: fourteen-check's spread reps x C(14,6) (index 3003 * i + 6-subset dropped): every 8 inside a spread;
// 2: thirteen-check's 31 maximal 13s x C(13,5) (1287 * i + 5-subset dropped): every 8 inside a maximal 13;
// 3: Part A's maximal 5-cliques as enumerated (only GL(4,2)-conjugation-reduced): every maximal 8;
// 4: twelve-check's 376 configurations x C(12,4) (495 * i + 4-subset dropped): every 8 inside some 12;
// 5: eleven-check's 817 configurations x C(11,3) (165 * i + triple dropped): every 8 inside some 11;
// 6: ten-check's 1777 configurations x C(10,2) (45 * i + pair dropped): every 8 inside some 10;
// 7: nine-check's 5663 configurations x 9 (9 * i + j): every 8 that extends at all.
static vector<Conf> confs_of_family(const vector<vector<Mask>>& fam, int type) {
    vector<Conf> out;
    for (size_t i = 0; i < fam.size(); i++) { Conf c; c.idx = (int)i; c.type = type; c.nonext = type != 0; c.hole = -1; c.aut = 0; c.mem = fam[i]; out.push_back(c); }
    return out;
}
static vector<Conf> read_maxc(const char* file) {
    vector<Conf> out; FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); exit(1); }
    char line[1024]; int idx = 0;
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "MAXC", 4)) continue;
        Conf c; c.idx = idx++; c.type = TYPES - 1; c.nonext = true; c.hole = -1; c.aut = 0;
        const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; c.Ns.push_back((u16)v); p = e; }
        if (c.Ns.size() != KN) { printf("bad MAXC line\n"); exit(1); }
        out.push_back(c);
    }
    fclose(fp); return out;
}
static vector<Mask> conf_members(const Conf& c) { return c.mem.empty() ? members_of(c.Ns) : c.mem; }
static void print_stats(const char* tag, const Stats& s) {
    printf("%s nodes by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.nodes[i]);
    printf("\n%s tries by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.tries[i]);
    printf("\n%s of which at points of D:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.triesD[i]);
    printf("\n%s nodes with every off-D point covered, by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.p1[i]);
    printf("\n%s bent-extension nodes %lld, full-span records %lld\n", tag, s.bentext, s.full); fflush(stdout);
}
static void write_stats(FILE* pf, int j, const Stats& s) {
    fprintf(pf, "R %d", j); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.nodes[k]); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.tries[k]);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.p1[k]); fprintf(pf, " %lld %lld", s.bentext, s.full);
    for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.triesD[k]); fprintf(pf, "\n"); fflush(pf);
}
static bool read_stats(const char* line, int& id, Stats& s) {
    const char* p = line; if (*p != 'R') return false; p++;
    char* e; id = (int)strtol(p, &e, 10); if (e == p) return false; p = e;
    long long v[38]; for (int i = 0; i < 38; i++) { v[i] = strtoll(p, &e, 10); if (e == p) return false; p = e; }
    for (int k = 0; k < 9; k++) { s.nodes[k] = v[k]; s.tries[k] = v[9 + k]; s.p1[k] = v[18 + k]; s.triesD[k] = v[29 + k]; } s.bentext = v[27]; s.full = v[28]; return true;
}
static int g_input = 0; // 0 CONF, 1 SPREAD x 3003, 2 maximal 13s x 1287, 3 MAXC, 4 twelve-check's x 495, 5 eleven-check's x 165, 6 ten-check's x 45, 7 nine-check's x 9
static int mode_search(const char* conffile, const char* progfile, int how, unsigned seed, int first, int last) {
    auto cf = g_input == 1 ? confs_of_family(familySp(read_spread_file(conffile)), 0) : g_input == 2 ? confs_of_family(family13(read_conf13_nonext(conffile)), 1)
        : g_input == 3 ? read_maxc(conffile) : g_input == 4 ? confs_of_family(family12(read_conf12(conffile)), 1)
        : g_input == 5 ? confs_of_family(family11(read_conf11(conffile)), 1) : g_input == 6 ? confs_of_family(family10(read_conf10(conffile)), 1)
        : g_input == 7 ? confs_of_family(family9(read_conf9(conffile)), 1) : read_confs(conffile);
    const char* inName[8] = { "configuration list", "spread file x 3003", "maximal 13s x 1287", "Part A maximal cliques", "twelve-check's configurations x 495", "eleven-check's configurations x 165", "ten-check's configurations x 45", "nine-check's configurations x 9" };
    printf("search: %zu configurations (%s); point order %s\n", cf.size(), inName[g_input], how == 0 ? "ascending" : how == 1 ? "descending" : "random");
    set<int> done; Stats tot, totT[TYPES];
    if (FILE* fp = fopen(progfile, "r")) {
        char line[1024];
        while (fgets(line, sizeof line, fp)) { int id; Stats s; if (read_stats(line, id, s) && !done.count(id)) { done.insert(id); add_stats(tot, s); } else if (!strncmp(line, "SOLUTION", 8)) printf("(previous run) %s", line); }
        fclose(fp);
        printf("resumed: %zu configurations already done\n", done.size());
    }
    FILE* pf = fopen(progfile, "a");
    auto t0 = chrono::steady_clock::now(); map<size_t, int> lineCounts; size_t nsolPrinted = 0;
    for (auto& c : cf) {
        if (c.idx < first || c.idx > last || done.count(c.idx)) continue;
        auto t1 = chrono::steady_clock::now();
        Problem P; build_problem(P, conf_members(c)); make_order(P, how, seed + c.idx);
        lineCounts[P.lineM.size()]++;
        if (mcount(P.D) != 135 || P.nr != 8 || P.lineM.size() != 1575) { printf("config %d: |D| = %d, %d members, %zu lines\n", c.idx, mcount(P.D), P.nr, P.lineM.size()); return 1; }
        Stats st = run_search(P, empty_state(), 400);
        { lock_guard<mutex> g(solmu); if (solutions.size() > nsolPrinted) { for (auto& kv : solutions) fprintf(pf, "SOLUTION %d %s\n", c.idx, kv.first.c_str()); nsolPrinted = solutions.size(); } }
        write_stats(pf, c.idx, st);
        add_stats(tot, st); add_stats(totT[c.type], st);
        double dt = chrono::duration<double>(chrono::steady_clock::now() - t1).count();
        printf("config %d (%s, %zu lines in D): nodes", c.idx, TNAME[c.type], P.lineM.size()); for (int k = 0; k <= 8; k++) printf(" %lld", st.nodes[k]);
        printf("  p1"); for (int k = 0; k <= 8; k++) printf(" %lld", st.p1[k]);
        long long tD = 0; for (int k = 0; k <= 8; k++) tD += st.triesD[k];
        printf("  triesD %lld full %lld bentext %lld  [%.1f s, total %.0f s]\n", tD, st.full, st.bentext, dt, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    }
    for (auto& kv : p1spaces) { fprintf(pf, "P1 dim %zu basis", kv.second.size()); for (u64 f : kv.second) fprintf(pf, " %07x", code28(f)); fprintf(pf, "\n"); }
    fclose(pf);
    printf("lines in D:"); for (auto& kv : lineCounts) printf(" %zu:%d", kv.first, kv.second); printf("\n");
    for (int t = 0; t < TYPES; t++) { string tag = string(TNAME[t]) + " (this run)"; print_stats(tag.c_str(), totT[t]); }
    print_stats("search (incl. resumed)", tot);
    printf("distinct p1 spans this run: %zu\n", p1spaces.size());
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
struct CtlResult { int status = 1; int rep = -1; int d0 = -1; int n4 = 0, inD = 0, nsol = 0, napn = 0; bool known = false; string why; Stats st; };
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
// g_startNonD (the in-D test): pick a rank-6 form g of T whose radical lies inside D, start from a random hyperplane W
// of T not containing g, with the points of g's radical first in the order. g's points are covered in T by g alone, so
// the search must find g among the main candidate generator's 13,888 forms on that line (and pass the legality and
// disjointness tests against all of W) to reach T. Status 5 if T has no rank-6 radical inside D.
static bool g_startNonD = false;
static void ctl_search(const vector<Mask>& members, const vector<u64>& T, int how, CtlResult& R, bool verbose) {
    string key = space_key(T);
    Problem P; build_problem(P, members); make_order(P, how, 1);
    Profile Q = profile(T);
    for (auto& m : Q.r4) { bool in = false; for (int j = 0; j < P.nr; j++) if (meq(m, P.R[j])) in = true; if (!in) { R.status = 2; R.why = "rank-4 radical not a member"; return; } }
    for (size_t i = 0; i < Q.f6.size(); i++) { if (msub(Q.r6[i], P.D)) R.inD++; else { P.oracle.push_back(Q.f6[i]); P.oracleM.push_back(Q.r6[i]); } }
    State S0 = empty_state();
    if (g_startNonD) {
        vector<int> inD; for (size_t i = 0; i < Q.f6.size(); i++) if (msub(Q.r6[i], P.D)) inD.push_back((int)i);
        if (inD.empty()) { R.status = 5; R.why = "no rank-6 radical inside D"; return; }
        mt19937 rng((unsigned)hash<string>()(key));
        int gi = inD[rng() % inD.size()]; u64 g = Q.f6[gi];
        int cg = -1; for (int c = 1; c < 256 && cg < 0; c++) { u64 f = 0; for (int j = 0; j < 8; j++) if (c >> j & 1) f ^= T[j]; if (f == g) cg = c; }
        if (cg < 0) { R.status = 2; R.why = "g not in T"; return; }
        int a; do a = (int)(rng() & 255); while (!par((u32)(a & cg))); // phi(x) = parity(a . coords(x)), phi(g) = 1
        for (int j = 0; j < 8; j++) {
            u64 w = T[j] ^ ((a >> j & 1) ? g : 0);
            if (find(S0.span.begin(), S0.span.end(), w) != S0.span.end()) continue;
            State T2; if (!try_add(P, S0, w, T2)) { R.status = 2; R.why = "W not legal"; return; } S0 = T2;
        }
        R.d0 = S0.d;
        if (S0.d != 7) { R.status = 2; R.why = "W not a hyperplane"; return; }
        Mask gm = Q.r6[gi]; if (!mdisj(gm, S0.cov)) { R.status = 2; R.why = "g's radical covered by W"; return; }
        vector<int> ord; for (int v = 1; v < 256; v++) if (mget(gm, v)) ord.push_back(v);
        for (int v : P.order) if (!mget(gm, v)) ord.push_back(v);
        P.order = ord;
    }
    map<string, vector<u64>> sol; tl_sol = &sol;
    R.st = run_search(P, S0, verbose ? 400 : 1);
    tl_sol = nullptr;
    R.known = sol.count(key) > 0; R.nsol = (int)sol.size();
    for (auto& kv : sol) { Verdict V = verify_space(kv.second); R.napn += V.apn;
        if (verbose) { printf("  APN=%d dist:", V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("%s\n", kv.first == key ? "  <- known" : ""); } }
    R.status = R.known ? 0 : 3;
}
// Extendable path: a spread through T's radicals, 9 other spread elements as holes, the 8-set carried onto its
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
    // choice -> a 9-subset of `others` (the holes; the other 8 spread elements are the members)
    mt19937 rng((unsigned)choice); vector<int> o = others; shuffle(o.begin(), o.end(), rng);
    if (o.size() < 9) { R.why = "too many rank-4"; return R; }
    set<int> holes(o.begin(), o.begin() + 9);
    vector<Mask> members; for (int i = 0; i < 17; i++) if (!holes.count(i)) members.push_back(spreadAll[i]);
    Canon K = canonP(members);
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
static thread_local int g_forceRep = -1; // embed into this representative (any type) instead of cycling the non-extendable ones
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
    int nm = (int)M.size();
    vector<array<u8, 4>> bas(nm); vector<vector<u8>> pts(nm); for (int i = 0; i < nm; i++) { basis_of_mask(M[i], 4, bas[i].data()); pts[i] = mpoints(M[i]); }
    vector<array<int, 3>> tr; for (int a = 0; a < nm; a++) for (int b = 0; b < nm; b++) for (int c = 0; c < nm; c++) if (a != b && a != c && b != c) tr.push_back({ a, b, c });
    shuffle(tr.begin(), tr.end(), rng);
    for (auto& t : tr) {
        Lin Cm; if (!frame_of(bas, pts, M, t[0], t[1], t[2], Cm)) { R.why = "frame"; return R; }
        Lin D;
        if (Q.n4 == 2) D = diag_ab(GL4[rng() % GL4.size()], GL4[rng() % GL4.size()]);
        else if (Q.n4 == 3) D = diag_hh(GL4[rng() % GL4.size()]);
        else {
            int tv = vid[N4]; if (tv < 0) { R.why = "4th radical not in V5824"; return R; }
            vector<int> ds; for (int d = 0; d < nm; d++) if (d != t[0] && d != t[1] && d != t[2]) { int v = vertex_in_frame(Cm, bas[d]); if (cls[v] == cls[tv]) ds.push_back(v); }
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
static int mode_ctlbatch(const char* conffile, const char* file, int minr4, long long maxn, int nch, const char* outfile, int how) {
    init_conj();
    CtlCtx X = ctl_ctx(conffile);
    g_logp1 = false;
    FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); return 1; }
    FILE* of = fopen(outfile, "a");
    vector<pair<long long, vector<u64>>> batch; long long ln = 0, taken = 0;
    Stats stT[2]; long long s0full = 0; map<int, long long> d0h[2];
    long long tot[2] = { 0 }, ok[2] = { 0 }, notfound[2] = { 0 }, nospread = 0, fail[2] = { 0 }, extra[2] = { 0 }, nonapn[2] = { 0 }, noembed = 0;
    map<int, int> repHits[2]; map<int, int> byN4[2];
    static char buf[1 << 16];
    auto flush = [&]() {
        int n = (int)batch.size(); vector<vector<CtlResult>> res(n, vector<CtlResult>(2 * nch));
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n; i++) for (int c = 0; c < 2 * nch; c++) {
            int ch = (int)(batch[i].first * 131 + c);
            res[i][c] = (c < nch) ? ctl_ext(X, batch[i].second, ch, how, false) : ctl_non(X, batch[i].second, ch, how, false);
        }
        for (int i = 0; i < n; i++) for (int c = 0; c < 2 * nch; c++) {
            int t = c < nch ? 0 : 1; auto& R = res[i][c];
            if (R.status == 4) { noembed++; continue; }
            if (R.status == 5) { s0full++; continue; }
            if (R.d0 >= 0) d0h[t][R.d0]++;
            if (R.why == "no spread through the rank-4 radicals") { nospread++; continue; }
            tot[t]++;
            add_stats(stT[t], R.st);
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
        print_stats(t ? "ctl nonext" : "ctl ext", stT[t]);
        if (g_startNonD) { printf("  started from S0 of dimension:"); for (auto& kv : d0h[t]) printf(" %d:%lld", kv.first, kv.second); printf("\n"); }
        printf("ctlbatch %s %s: %lld controls; re-found %lld; NOT found %lld; setup failures %lld; other spaces found %lld (%lld not APN)\n", file, t ? "nonext" : "ext", tot[t], ok[t], notfound[t], fail[t], extra[t], nonapn[t]);
        printf("  re-found by rank-4 count:"); for (auto& kv : byN4[t]) printf(" %d:%d", kv.first, kv.second); printf("\n");
        printf("  reps hit: %zu\n", repHits[t].size());
        fprintf(of, "SUMMARY %s %s controls %lld refound %lld notfound %lld fail %lld extra %lld nonapn %lld reps %zu\n", file, t ? "nonext" : "ext", tot[t], ok[t], notfound[t], fail[t], extra[t], nonapn[t], repHits[t].size());
        fprintf(of, "REPS %s", t ? "nonext" : "ext"); for (auto& kv : repHits[t]) fprintf(of, " %d:%d", kv.first, kv.second); fprintf(of, "\n");
    }
    printf("skipped: no spread %lld, no member 4-tuple of the class (k = 4, nonext) %lld, S0 already the whole space %lld\n", nospread, noembed, s0full);
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
        bool ok = mcount(P.D) == 135 && P.nr == 8;
        for (int a = 0; a < KN; a++) { ok &= inv4(c.Ns[a]) && inv4(c.Ns[a] ^ ID4); for (int b = a + 1; b < KN; b++) ok &= inv4(c.Ns[a] ^ c.Ns[b]); }
        lines[P.lineM.size()]++; ok &= P.lineM.size() == 1575;
        { int per[256] = { 0 }; for (auto& m : P.lineM) for (int v = 1; v < 256; v++) per[v] += mget(m, v); for (int v = 1; v < 256; v++) ok &= per[v] == (mget(P.D, v) ? 35 : 0); }
        for (size_t li = 0; li < P.lineM.size(); li++) ok &= msub(P.lineM[li], P.D) && mcount(P.lineM[li]) == 3;
        for (int j = 0; j < 8; j++) { set<u64> s(P.f4[j].begin(), P.f4[j].end()); ok &= s.size() == 28; for (u64 f : P.f4[j]) { Mask m; ok &= radical(f, m) == 15 && meq(m, P.R[j]); checked++; } }
        if (c.idx % 128 == 0) for (size_t li = 0; li < P.lineM.size(); li++) {
            const vector<u64>& lf = line_forms(P, li); set<u64> s(lf.begin(), lf.end()); ok &= s.size() == 13888;
            for (u64 f : lf) { Mask m; ok &= rank8(f) == 6 && radical(f, m) == 3 && meq(m, P.lineM[li]); checked++; }
        }
        // canonical form: invariant under a random GL(8,2) image, and the stored aut
        { Lin g = random_gl8(rng); vector<Mask> Mg; for (auto& m : M) Mg.push_back(mask_image(g, m)); shuffle(Mg.begin(), Mg.end(), rng);
          Canon K = canonP(Mg); vector<u16> key; for (u16 N : c.Ns) key.push_back((u16)vid[N]); sort(key.begin(), key.end());
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

// ================================================================================ count5: k-cliques of V5824 by plain enumeration
// No group action, no canonical forms: every clique once (increasing vertex index), 6-cliques by popcount at the last
// level. The 5-clique count must equal Part A's weighted count (merge), the 6-clique count nine-check's (47414189312).
static int mode_count5() {
    long long cnt[7] = { 0 };
    auto t0 = chrono::steady_clock::now();
    BS all; for (int i = 0; i < NW; i++) all.w[i] = ~0ULL;
    int n = (int)VN.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int v1 = 0; v1 < n; v1++) {
        long long c[7] = { 0 }; c[1] = 1;
        BS a; above(a, all, adj[v1], v1);
        for (int v2 : blist(a)) { c[2]++; BS b; above(b, a, adj[v2], v2);
            for (int v3 : blist(b)) { c[3]++; BS d; above(d, b, adj[v3], v3);
                for (int v4 : blist(d)) { c[4]++; BS e; above(e, d, adj[v4], v4);
                    for (int wi = 0; wi < NW; wi++) { u64 w = e.w[wi]; while (w) { int v5 = wi * 64 + (int)_tzcnt_u64(w); w &= w - 1;
                        c[5]++; BS f; above(f, e, adj[v5], v5); c[6] += bcount(f); } } } } }
#pragma omp critical
        for (int k = 1; k <= 6; k++) cnt[k] += c[k];
    }
    printf("count5: k-cliques of V5824 by plain enumeration:"); for (int k = 1; k <= 6; k++) printf(" %d:%lld", k, cnt[k]);
    printf("  (6-cliques: nine-check's weighted count 47414189312 %s)  [%.0f s]\n", cnt[6] == 47414189312LL ? "OK" : "MISMATCH", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    return cnt[6] == 47414189312LL ? 0 : 2;
}

int main(int argc, char** argv) {

    if (argc < 2) { printf("modes: count5 | parta <class> <orbits out> <maxcliques out> [<slice> <nslices> [<test firsts>]] | merge <orbits out> <configs.txt (check9)> <part files...> | configs <orbits> <configs.txt (check9)> <configs.txt (check10)> <configs.txt (check11)> <configs.txt (check12)> <spreads.txt (check14)> <configs.txt (check13)> <out> | search|searchsp|search13|searchmax|search12|search11|search10|search9 <input> <prog> [asc|desc|rand <seed>] [first last] | selftest <confs> | ctl <confs> <lut> <line> <choice> [non] | ctlbatch <confs> <lut> <minr4> <maxn> <nch> <out> [nond|dfirst] | ctlreps <confs> <lut> <first line> <n> <rep>...\n"); return 1; }
    string mode = argv[1];
    ND4 = nondeg_codes(4); ND6 = nondeg_codes(6);
    init_V();
    printf("V5824: %zu vertices, %zu classes:", VN.size(), classSize.size()); for (size_t c = 0; c < classSize.size(); c++) printf(" %d (rep %04x)", classSize[c], VN[classRep[c]]); printf("\n");
    fflush(stdout);
    if (mode == "parta") { int sl = 0, nsl = 1; if (argc > 6) { sl = atoi(argv[5]); nsl = atoi(argv[6]); } if (argc > 7) g_testFirsts = atoi(argv[7]); return mode_parta(atoi(argv[2]), argv[3], argv[4], sl, nsl); }
    if (mode == "count5") return mode_count5();
    if (mode == "merge") { vector<const char*> parts; for (int i = 4; i < argc; i++) parts.push_back(argv[i]); return mode_merge(argv[2], argv[3], parts); }
    if (mode == "configs") return mode_configs(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);
    init_code28();
    if (mode == "search" || mode == "searchsp" || mode == "search13" || mode == "searchmax" || mode == "search12" || mode == "search11" || mode == "search10" || mode == "search9") {
        g_input = mode == "searchsp" ? 1 : mode == "search13" ? 2 : mode == "searchmax" ? 3 : mode == "search12" ? 4 : mode == "search11" ? 5 : mode == "search10" ? 6 : mode == "search9" ? 7 : 0;
        int how = 0; unsigned seed = 1; int first = 0, last = INT_MAX, ai = 4;
        if (argc > ai && !strcmp(argv[ai], "desc")) { how = 1; ai++; }
        else if (argc > ai && !strcmp(argv[ai], "rand")) { how = 2; seed = (unsigned)atoi(argv[ai + 1]); ai += 2; }
        else if (argc > ai && !strcmp(argv[ai], "asc")) ai++;
        if (argc > ai + 1) { first = atoi(argv[ai]); last = atoi(argv[ai + 1]); }
        return mode_search(argv[2], argv[3], how, seed, first, last);
    }
    if (mode == "selftest") return mode_selftest(argv[2]);
    if (mode == "ctl") return mode_ctl(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argc > 6 && !strcmp(argv[6], "non"));
    if (mode == "ctlbatch") return mode_ctlbatch(argv[2], argv[3], atoi(argv[4]), atoll(argv[5]), atoi(argv[6]), argv[7], (g_startNonD = argc > 8 && !strcmp(argv[8], "nond"), argc > 8 && !strcmp(argv[8], "dfirst") ? 3 : 0));
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
    if (mode == "ctlrot") { // ctlrot <confs> <lut> <first line> <pool n> <rep>...: one control per rep, the functions of the
        // pool (the first n with k >= 3) in rotation, parallel over the reps
        init_conj(); CtlCtx X = ctl_ctx(argv[2]); g_logp1 = false;
        vector<vector<u64>> fs; vector<int> lns; FILE* fp = fopen(argv[3], "r"); static char buf[1 << 16]; int ln = 0, first = atoi(argv[4]), want = atoi(argv[5]);
        while (fgets(buf, sizeof buf, fp) && (int)fs.size() < want) { ln++; if (ln < first) continue; u8 F[256]; if (!read_lut(buf, F)) continue; auto B = forms_of_lut(F); if (profile(B).n4 >= 3) { fs.push_back(B); lns.push_back(ln); } }
        fclose(fp);
        vector<int> reps; for (int a = 6; a < argc; a++) reps.push_back(atoi(argv[a]));
        int n = (int)reps.size(); vector<CtlResult> res(n); vector<int> fi(n);
#pragma omp parallel for schedule(dynamic, 1)
        for (int k = 0; k < n; k++) {
            g_forceRep = reps[k];
            for (int t = 0; t < (int)fs.size(); t++) { // the next function in rotation that embeds (k = 4 needs a member of its class)
                fi[k] = (reps[k] + t) % (int)fs.size();
                res[k] = ctl_non(X, fs[fi[k]], lns[fi[k]] * 31 + reps[k], 0, false);
                if (res[k].status != 4) break;
            }
        }
        int bad = 0;
        for (int k = 0; k < n; k++) {
            auto& R = res[k];
            printf("rep %d (%s): line %d, %s, %d other spaces%s%s\n", reps[k], X.cf[reps[k]].nonext ? "nonext" : "ext", lns[fi[k]], R.status == 0 ? "re-found" : "NOT re-found",
                R.status == 0 ? R.nsol - 1 : 0, R.status ? " " : "", R.status ? R.why.c_str() : "");
            bad += R.status != 0;
        }
        fflush(stdout);
        return bad ? 3 : 0;
    }
    if (mode == "mapa") { // mapa <configs.txt> <search A's list of 9-member configurations>: which of our configurations each of A's lies on (comparison only)
        init_conj();
        map<vector<u16>, int> confOf; map<int, string> confLine; { FILE* fp = fopen(argv[2], "r"); char line[1024]; while (fgets(line, sizeof line, fp)) { int ci; if (sscanf(line, "CONF %d", &ci) != 1) continue; vector<u16> Ns; const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; } confOf[key_of_Ns(Ns)] = ci; string s(line); confLine[ci] = s.substr(0, s.find(':')); } fclose(fp); }
        vector<string> names; vector<vector<Mask>> Ms; { FILE* fp = fopen(argv[3], "r"); static char line[1 << 14]; while (fgets(line, sizeof line, fp)) {
            char nm[128]; int off; if (sscanf(line, "%127s%n", nm, &off) != 1) continue; const char* p = line + off; char* e; vector<int> v; while (true) { long x = strtol(p, &e, 16); if (e == p) break; v.push_back((int)x); p = e; }
            // A's line as at nine (name, member count "8", 4 basis bytes per member); unverified for eight: check A's file's format first
            if (v.empty() || v[0] != 0x8 || (int)v.size() != 1 + 32) { printf("unexpected line %s (%zu values)\n", nm, v.size()); return 1; }
            vector<Mask> M; for (int s = 0; s < 8; s++) {
 Mask m; mclear(m); for (int c = 1; c < 16; c++) { u8 x = 0; for (int b = 0; b < 4; b++) if (c >> b & 1) x ^= (u8)v[1 + 4 * s + b]; mset(m, x); } if (mcount(m) != 15) { printf("%s: member %d not a solid\n", nm, s); return 1; } M.push_back(m); }
            names.push_back(nm); Ms.push_back(M); } fclose(fp); }
        vector<int> hit(Ms.size(), -1); vector<long long> aut(Ms.size());
#pragma omp parallel for schedule(dynamic, 4)
        for (int i = 0; i < (int)Ms.size(); i++) { Canon K = canonP(Ms[i]); aut[i] = K.aut; auto it = confOf.find(K.key); if (it != confOf.end()) hit[i] = it->second; }
        map<int, int> mult; int unk = 0; for (size_t i = 0; i < Ms.size(); i++) { if (hit[i] < 0) unk++; else mult[hit[i]]++; printf("A %s -> %d (%s) canon aut %lld\n", names[i].c_str(), hit[i], hit[i] < 0 ? "none" : confLine[hit[i]].c_str(), aut[i]); }
        int dup = 0; for (auto& kv : mult) dup += kv.second > 1;
        printf("mapa: %zu of A's configurations, %d on none of ours, %zu of our %zu configurations hit, %d hit more than once\n", Ms.size(), unk, mult.size(), confOf.size(), dup);
        return 0;
    }
    if (mode == "orbitcov") { // orbitcov <orbits> <configs.txt> <MAXC file>: which orbit each listed clique lies on
        init_conj();
        auto orbs = read_orbs(argv[2]); map<vector<u16>, int> byKey; for (int i = 0; i < (int)orbs.size(); i++) byKey[key_of_Ns(orbs[i].Ns)] = i;
        map<vector<u16>, int> confOf; { FILE* fp = fopen(argv[3], "r"); char line[1024]; while (fgets(line, sizeof line, fp)) { int ci; if (sscanf(line, "CONF %d", &ci) != 1) continue; vector<u16> Ns; const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; } confOf[key_of_Ns(Ns)] = ci; } fclose(fp); }
        vector<vector<u16>> cl; { FILE* fp = fopen(argv[4], "r"); char line[1024]; while (fgets(line, sizeof line, fp)) { if (strncmp(line, "MAXC", 4)) continue; vector<u16> Ns; const char* p = strchr(line, ':') + 1; char* e; while (true) { long v = strtol(p, &e, 16); if (e == p) break; Ns.push_back((u16)v); p = e; } cl.push_back(Ns); } fclose(fp); }
        vector<int> hit(cl.size(), -1);
#pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < (int)cl.size(); i++) { Canon K = canonP(members_of(cl[i])); auto it = byKey.find(K.key); if (it != byKey.end()) hit[i] = it->second; }
        map<int, int> per; int unk = 0, nonmax = 0; for (int h : hit) { if (h < 0) { unk++; continue; } per[h]++; if (TNAME[orbs[h].type] != string("max")) nonmax++; }
        int maxOrbs = 0; for (auto& o : orbs) maxOrbs += TNAME[o.type] == string("max");
        printf("orbitcov: %zu cliques, %d on no listed orbit, %d on a non-max orbit; %zu orbits hit of %d max orbits\n", cl.size(), unk, nonmax, per.size(), maxOrbs);
        for (auto& o : orbs) if (TNAME[o.type] == string("max")) { auto it = confOf.find(key_of_Ns(o.Ns)); printf("  orbit %d (config %d): %d sampled\n", o.idx, it == confOf.end() ? -1 : it->second, per.count(o.idx) ? per[o.idx] : 0); }
        if (argc > 6) { // [<out> <orbit>...]: write the listed cliques lying on these orbits
            set<int> want; for (int a = 6; a < argc; a++) want.insert(atoi(argv[a]));
            FILE* of = fopen(argv[5], "w"); int nw = 0;
            for (size_t i = 0; i < cl.size(); i++) if (want.count(hit[i])) { fprintf(of, "MAXC orbit %d :", hit[i]); for (u16 N : cl[i]) fprintf(of, " %04x", N); fprintf(of, "\n"); nw++; }
            fclose(of); printf("wrote %d cliques to %s\n", nw, argv[5]);
        }
        return unk || nonmax ? 2 : 0;
    }
    return 1;
}
