// sixteen_check.cpp -- independent check of Open Problem VI.1 (arXiv 2510.12008) for [0^234, 2^5, 4^16].
// Subagent (Claude Opus 5.5) for Rowen, 2026-10-24. Shares no code with sixteen.cpp (not read).
//
// Organisation: the space S of alternating forms grows one form at a time. At each node we take the first point
// (in a fixed order) that lies in no radical of any form of the current span, and branch over every legal form
// whose radical contains that point. Every element of the span must be legal (bent, or rank 4 / rank 6 with an
// allowed radical) and radicals must be pairwise disjoint. No spread sets are enumerated: the spread appears as
// the radicals of the forms that get added. Coordinates: v = (x, y), x = low nibble, y = high nibble,
// Y = {x = 0} is the 17th spread element R17, X = {y = 0} the radical of B0 = J(y, y').
// Modes: selftest | main <progress file> [all] [K r] | known <file> | ctlB <file> <line> <s6> [hint] | ctlA <file> <line> <s6>
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
#include <intrin.h>
#include <immintrin.h>
#include <omp.h>
using namespace std;
typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16; typedef uint8_t u8;

static inline int par(u32 v) { return __popcnt(v) & 1; }
static inline u8 rowof(u64 f, int i) { return (u8)(f >> (8 * i)); }
static inline u8 applyf(u64 f, u8 v) { u8 r = 0; for (int i = 0; i < 8; i++) if (v >> i & 1) r ^= rowof(f, i); return r; }
static inline int evalf(u64 f, u8 u, u8 w) { return par(u & applyf(f, w)); }

static bool is_alt(u64 f) {
    for (int i = 0; i < 8; i++) {
        if (rowof(f, i) >> i & 1) return false;
        for (int j = 0; j < 8; j++) if ((rowof(f, i) >> j & 1) != (rowof(f, j) >> i & 1)) return false;
    }
    return true;
}

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
static inline void mor(Mask& a, const Mask& b) { for (int i = 0; i < 4; i++) a.w[i] |= b.w[i]; }
static inline int mcount(const Mask& m) { return (int)(__popcnt64(m.w[0]) + __popcnt64(m.w[1]) + __popcnt64(m.w[2]) + __popcnt64(m.w[3])); }
static inline bool meq(const Mask& a, const Mask& b) { return !memcmp(&a, &b, sizeof(Mask)); }

// nonzero points of the radical, by Gray code over all of V
static int radical(u64 f, Mask& m) {
    mclear(m); int cnt = 0; u8 cur = 0, v = 0;
    for (u32 i = 1; i < 256; i++) { int b = _tzcnt_u32(i); v ^= (u8)(1u << b); cur ^= rowof(f, b); if (!cur) { mset(m, v); cnt++; } }
    return cnt;
}

static Mask Ymask, nonYmask, Xmask, nonXmask;

// holesXY (case i = 2): rank-4 radicals must miss X too, and "in a hole" means inside Y or inside X.
// orb/minidx: every rank-4 form must be in the candidate table with orbit index >= minidx.
struct Rules { bool r4, r4compY, r6inY, r6out; bool holesXY = false; const unordered_map<u64, int>* orb = nullptr; int minidx = 0; };
static const Rules MAIN_RULES = { true, true, true, false };

static inline bool legal(u64 f, const Rules& R, Mask& m, int& cnt) {
    if (rank8(f) == 8) { cnt = 0; return true; }
    cnt = radical(f, m);
    if (cnt == 15) {
        if (!R.r4 || (R.r4compY && !mdisj(m, Ymask)) || (R.holesXY && !mdisj(m, Xmask))) return false;
        if (R.orb) { auto it = R.orb->find(f); if (it == R.orb->end() || it->second < R.minidx) return false; }
        return true;
    }
    if (cnt == 3) {
        bool inHole = mdisj(m, nonYmask) || (R.holesXY && mdisj(m, nonXmask));
        return inHole ? R.r6inY : R.r6out;
    }
    return false;
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
static inline u16 mm4(u16 A, u16 B) { // A*B
    u16 C = 0;
    for (int r = 0; r < 4; r++) { u8 row = 0; for (int k = 0; k < 4; k++) if (A >> (4 * r + k) & 1) row ^= (B >> (4 * k)) & 15; C |= (u16)row << (4 * r); }
    return C;
}
static bool inv4(u16 N) { for (u8 x = 1; x < 16; x++) if (!mv4(N, x)) return false; return true; }
static u16 fromcols4(const u8 col[4]) { u16 N = 0; for (int i = 0; i < 4; i++) for (int r = 0; r < 4; r++) if (col[i] >> r & 1) N |= 1 << (4 * r + i); return N; }

static vector<u16> GL4;
static vector<Mask> graphMask;           // per GL4 index: graph(N) nonzero points
static vector<array<u8, 4>> graphU;      // functionals (x,y) -> (Nx)_k + y_k
static vector<int> bucket[16][16];       // N with N x = y

// ---- 4x4 alternating forms by code (pairs (0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
static void code_to_rows4(u32 c, u8 R[4]) { memset(R, 0, 4); int p = 0; for (int k = 0; k < 4; k++) for (int l = k + 1; l < 4; l++) { if (c >> p & 1) { R[k] |= 1 << l; R[l] |= 1 << k; } p++; } }
static u32 rows4_to_code(const u8 R[4]) { u32 c = 0; int p = 0; for (int k = 0; k < 4; k++) for (int l = k + 1; l < 4; l++) { if (R[k] >> l & 1) c |= 1u << p; p++; } return c; }
static int eval4(const u8 R[4], u8 u, u8 w) { int s = 0; for (int k = 0; k < 4; k++) if (u >> k & 1) s ^= par(R[k] & w); return s; }
static u32 pullback4(u32 c, u16 h) { // M'(u,w) = M(hu, hw)
    u8 R[4]; code_to_rows4(c, R); u8 col[4]; for (int i = 0; i < 4; i++) col[i] = mv4(h, 1 << i);
    u8 O[4] = { 0 };
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (eval4(R, col[i], col[j])) { O[i] |= 1 << j; O[j] |= 1 << i; }
    return rows4_to_code(O);
}
static const u32 JCODE = 1 | (1 << 5); // J(e0,e1) = J(e2,e3) = 1

// ---- subspaces of V
struct Sub { Mask m; u8 basis[4]; u8 ann[6]; int dim; };
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
static void enum_subspaces(int dim, vector<Sub>& out) {
    for (u32 L = 0; L < 256; L++) if ((int)__popcnt(L) == dim) {
        int h[8]; int k = 0; for (int b = 7; b >= 0; b--) if (L >> b & 1) h[k++] = b;
        u32 freeb[8]; int nf[8]; int total = 0;
        for (int j = 0; j < dim; j++) { freeb[j] = ((1u << h[j]) - 1) & ~L; nf[j] = __popcnt(freeb[j]); total += nf[j]; }
        for (u32 t = 0; t < (1u << total); t++) {
            Sub s; s.dim = dim; u32 tt = t;
            for (int j = 0; j < dim; j++) { u32 bits = _pdep_u32(tt & ((1u << nf[j]) - 1), freeb[j]); tt >>= nf[j]; s.basis[j] = (u8)((1u << h[j]) | bits); }
            mclear(s.m);
            for (u32 c = 1; c < (1u << dim); c++) { u8 v = 0; for (int j = 0; j < dim; j++) if (c >> j & 1) v ^= s.basis[j]; mset(s.m, v); }
            vector<u8> ann;
            for (int u = 1; u < 256; u++) { bool ok = true; for (int j = 0; j < dim; j++) if (par(u & s.basis[j])) { ok = false; break; } if (ok) ann.push_back((u8)u); }
            pick_independent(ann, 8 - dim, s.ann);
            out.push_back(s);
        }
    }
}
static vector<Sub> lines, subs4;
static vector<int> lines_at[256], subs4_at[256];
static vector<char> lineInY, lineInX;

static void init_tables(bool want4any) {
    for (int v = 1; v < 256; v++) { if ((v & 15) == 0) mset(Ymask, v); else mset(nonYmask, v); }
    for (int v = 1; v < 256; v++) { if ((v >> 4) == 0) mset(Xmask, v); else mset(nonXmask, v); }
    ND4 = nondeg_codes(4); ND6 = nondeg_codes(6);
    for (u32 N = 0; N < 65536; N++) if (inv4((u16)N)) GL4.push_back((u16)N);
    graphMask.resize(GL4.size()); graphU.resize(GL4.size());
    for (size_t i = 0; i < GL4.size(); i++) {
        u16 N = GL4[i]; Mask m; mclear(m);
        for (u8 x = 1; x < 16; x++) { u8 y = mv4(N, x); mset(m, x | (y << 4)); bucket[x][y].push_back((int)i); }
        graphMask[i] = m;
        for (int k = 0; k < 4; k++) graphU[i][k] = (u8)(((N >> (4 * k)) & 15) | (1 << (4 + k)));
    }
    enum_subspaces(2, lines);
    lineInY.resize(lines.size()); lineInX.resize(lines.size());
    for (size_t i = 0; i < lines.size(); i++) {
        lineInY[i] = mdisj(lines[i].m, nonYmask);
        lineInX[i] = mdisj(lines[i].m, nonXmask);
        for (int v = 1; v < 256; v++) if (mget(lines[i].m, v)) lines_at[v].push_back((int)i);
    }
    if (want4any) {
        enum_subspaces(4, subs4);
        for (size_t i = 0; i < subs4.size(); i++) for (int v = 1; v < 256; v++) if (mget(subs4[i].m, v)) subs4_at[v].push_back((int)i);
    }
}

// ---- 28-bit codes for whole forms (pairs i<j of 0..7, lexicographic)
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

// ---- independent verification from scratch: LUT, derivative counts, Walsh transform
struct Verdict { bool apn; int dist[9]; };
static Verdict verify_space(const vector<u64>& basis) {
    Verdict V; memset(&V, 0, sizeof V);
    int n = (int)basis.size();
    u8 F[256];
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
        int k = 0; while ((1 << (8 + k)) < mx * mx) k++; // mx^2 = 2^(8+k)
        V.dist[k]++;
    }
    return V;
}

// ---- the search
struct State { vector<u64> span; Mask cov; int d; vector<u64> basis; };
static bool try_add(const State& S, u64 f, const Rules& R, State& out) {
    Mask nc = S.cov; size_t n = S.span.size();
    for (size_t i = 0; i < n; i++) {
        Mask m; int c;
        if (!legal(f ^ S.span[i], R, m, c)) return false;
        if (c) { if (!mdisj(m, nc)) return false; mor(nc, m); }
    }
    out.span.resize(2 * n);
    for (size_t i = 0; i < n; i++) { out.span[i] = S.span[i]; out.span[n + i] = S.span[i] ^ f; }
    out.cov = nc; out.d = S.d + 1; out.basis = S.basis; out.basis.push_back(f);
    return true;
}
static State empty_state() { State S; S.span = { 0 }; mclear(S.cov); S.d = 0; return S; }

struct Gen { bool g4graph, g4any, g6Y, g6any, g6oracle = false, g6X = false; };
// control C only: rank-6 forms with radical outside Y restricted to a given list (a known function's)
static vector<u64> oracle6; static vector<Mask> oracle6m;
template <class CB> static void for_each_candidate(int a, const Mask& cov, const Gen& G, CB cb) {
    int x = a & 15, y = a >> 4;
    vector<u64> buf;
    if (G.g6oracle) for (size_t i = 0; i < oracle6.size(); i++) if (mget(oracle6m[i], a) && mdisj(oracle6m[i], cov) && !mdisj(oracle6m[i], nonYmask)) cb(oracle6[i]);
    if (G.g4graph && x && y) for (int idx : bucket[x][y]) if (mdisj(graphMask[idx], cov)) {
        buf.clear(); forms_from(graphU[idx].data(), 4, ND4, buf); for (u64 f : buf) cb(f);
    }
    if (G.g4any) for (int si : subs4_at[a]) if (mdisj(subs4[si].m, cov)) {
        buf.clear(); forms_from(subs4[si].ann, 4, ND4, buf); for (u64 f : buf) cb(f);
    }
    if (G.g6Y && x == 0) for (int li : lines_at[a]) if (lineInY[li] && mdisj(lines[li].m, cov)) {
        buf.clear(); forms_from(lines[li].ann, 6, ND6, buf); for (u64 f : buf) cb(f);
    }
    if (G.g6X && y == 0) for (int li : lines_at[a]) if (lineInX[li] && mdisj(lines[li].m, cov)) {
        buf.clear(); forms_from(lines[li].ann, 6, ND6, buf); for (u64 f : buf) cb(f);
    }
    if (G.g6any) for (int li : lines_at[a]) if (!(G.g6Y && x == 0 && lineInY[li]) && !(G.g6X && y == 0 && lineInX[li]) && mdisj(lines[li].m, cov)) {
        buf.clear(); forms_from(lines[li].ann, 6, ND6, buf); for (u64 f : buf) cb(f);
    }
}

struct Config { Rules R; Gen Gnon, GY; vector<int> order; Gen GX = {}; };
static inline const Gen& gen_for(int a, const Config& C) { return (a & 15) ? ((a >> 4) ? C.Gnon : C.GX) : C.GY; }
struct Stats { long long nodes[10] = { 0 }, tries[10] = { 0 }, bentext = 0, full = 0; };
static mutex solmu;
static map<string, vector<u64>> solutions;

static void record(const State& S, Stats& st) {
    st.full++;
    string k = space_key(S.basis);
    lock_guard<mutex> g(solmu);
    if (!solutions.count(k)) {
        solutions[k] = S.basis;
        printf("SOLUTION %s\n", k.c_str()); fflush(stdout);
    }
}
static void bent_extend(const State& S, const Config& C, Stats& st) {
    st.bentext++;
    vector<u32> c; for (u64 f : S.basis) c.push_back(code28(f));
    auto B = rref28(c);
    u32 piv = 0; for (u32 b : B) { unsigned long hb; _BitScanReverse(&hb, b); piv |= 1u << hb; }
    u32 freemask = ((1u << 28) - 1) & ~piv; int nf = __popcnt(freemask);
    printf("  bent-extension node at dim %d (%d free coords)\n", S.d, nf); fflush(stdout);
    for (u64 t = 1; t < (1ULL << nf); t++) {
        u64 f = form28(_pdep_u32((u32)t, freemask));
        bool ok = true; for (u64 s : S.span) if (rank8(f ^ s) != 8) { ok = false; break; }
        if (!ok) continue;
        State T; T.span.resize(2 * S.span.size());
        for (size_t i = 0; i < S.span.size(); i++) { T.span[i] = S.span[i]; T.span[S.span.size() + i] = S.span[i] ^ f; }
        T.cov = S.cov; T.d = S.d + 1; T.basis = S.basis; T.basis.push_back(f);
        if (T.d == 8) record(T, st); else bent_extend(T, C, st);
    }
}
static int first_uncovered(const State& S, const Config& C) { for (int v : C.order) if (!mget(S.cov, v)) return v; return -1; }
// nodes where every point off the holes is covered (the order puts those points first) and holes remain
static bool g_logp1 = false; static set<string> p1keys; static long long p1nodes = 0;
static void note_p1(const State& S, int a) {
    if (!g_logp1 || a < 0 || ((a & 15) && (a >> 4))) return;
    lock_guard<mutex> g(solmu); p1nodes++; p1keys.insert(space_key(S.basis));
}
static void dfs(const State& S, const Config& C, Stats& st) {
    st.nodes[S.d]++;
    int a = first_uncovered(S, C);
    note_p1(S, a);
    if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, C, st); return; }
    if (S.d == 8) return;
    const Gen& G = gen_for(a, C);
    State T;
    for_each_candidate(a, S.cov, G, [&](u64 f) { st.tries[S.d]++; if (try_add(S, f, C.R, T)) dfs(T, C, st); });
}
// children of S (one branching step); a full span or a bent-extension node is handled in place
static void expand(const State& S, const Config& C, Stats& st, vector<State>& out) {
    st.nodes[S.d]++;
    int a = first_uncovered(S, C);
    note_p1(S, a);
    if (a < 0) { if (S.d == 8) record(S, st); else bent_extend(S, C, st); return; }
    if (S.d == 8) return;
    const Gen& G = gen_for(a, C);
    State T;
    for_each_candidate(a, S.cov, G, [&](u64 f) { st.tries[S.d]++; if (try_add(S, f, C.R, T)) out.push_back(T); });
}
static void add_stats(Stats& a, const Stats& b) { for (int i = 0; i < 10; i++) { a.nodes[i] += b.nodes[i]; a.tries[i] += b.tries[i]; } a.bentext += b.bentext; a.full += b.full; }
static void print_stats(const char* tag, const Stats& s) {
    printf("%s nodes by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.nodes[i]);
    printf("\n%s tries by dim:", tag); for (int i = 0; i <= 8; i++) printf(" %lld", s.tries[i]);
    printf("\n%s bent-extension nodes %lld, full-span records %lld\n", tag, s.bentext, s.full); fflush(stdout);
}
// expand until there are enough tasks, then DFS them in parallel
static Stats run_search(const State& S0, const Config& C, int mintasks = 2000) {
    Stats tot; vector<State> tasks = { S0 };
    while (tasks.size() && (int)tasks.size() < mintasks) {
        vector<State> nxt; for (auto& s : tasks) expand(s, C, tot, nxt);
        tasks.swap(nxt);
        if (tasks.size() && tasks[0].d >= 7) break;
    }
    printf("  %zu tasks\n", tasks.size()); fflush(stdout);
    int nt = (int)tasks.size();
    vector<Stats> per(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nt; i++) dfs(tasks[i], C, per[omp_get_thread_num()]);
    for (auto& p : per) add_stats(tot, p);
    return tot;
}

static bool g_desc = false, g_nopr = false; // alternative point order (descending) for a cross-check
static Config main_config() {
    Config C; C.R = MAIN_RULES;
    C.Gnon = { true, false, false, false }; C.GY = { false, false, true, false };
    for (int v = 1; v < 256; v++) if (v & 15) C.order.push_back(v);
    if (g_desc) reverse(C.order.begin(), C.order.end());
    size_t k = C.order.size();
    for (int v = 1; v < 256; v++) if (!(v & 15)) C.order.push_back(v);
    if (g_desc) reverse(C.order.begin() + k, C.order.end());
    return C;
}
static u64 formB0() { u8 u[4]; for (int k = 0; k < 4; k++) u[k] = (u8)(1 << (4 + k)); vector<u64> o; forms_from(u, 4, { JCODE }, o); return o[0]; }
static u64 formGraph(int idx, u32 code) { vector<u64> o; forms_from(graphU[idx].data(), 4, { code }, o); return o[0]; }
static int idxI() { for (size_t i = 0; i < GL4.size(); i++) if (GL4[i] == 0x8421) return (int)i; return -1; }

static vector<u16> symplectic() { vector<u16> S; for (u16 h : GL4) if (pullback4(JCODE, h) == JCODE) S.push_back(h); return S; }
static vector<vector<u32>> sp_orbits() {
    auto Sp = symplectic(); vector<vector<u32>> orbs; set<u32> seen;
    for (u32 c : ND4) if (!seen.count(c)) { set<u32> o; for (u16 h : Sp) o.insert(pullback4(c, h)); for (u32 x : o) seen.insert(x); orbs.push_back(vector<u32>(o.begin(), o.end())); }
    return orbs;
}

// ---------------------------------------------------------------------------------------------------- modes
static int selftest() {
    int bad = 0;
    auto chk = [&](bool ok, const char* what) { printf("%s %s\n", ok ? "ok  " : "FAIL", what); if (!ok) bad++; };
    chk(ND4.size() == 28, "28 nondegenerate alternating forms on F_2^4");
    chk(ND6.size() == 13888, "13888 nondegenerate alternating forms on F_2^6");
    chk(GL4.size() == 20160, "|GL(4,2)| = 20160");
    auto Sp = symplectic(); chk(Sp.size() == 720, "|Sp(4,2)| = 720");
    chk(lines.size() == 10795, "10795 lines (2-subspaces) of F_2^8");
    chk(subs4.size() == 200787, "200787 4-subspaces of F_2^8");
    bool okat = true; for (int v = 1; v < 256; v++) okat &= lines_at[v].size() == 127 && subs4_at[v].size() == 11811;
    chk(okat, "127 lines and 11811 4-spaces through each point");
    int ny = 0; for (char c : lineInY) ny += c; chk(ny == 35, "35 lines inside Y");
    for (int x = 1; x < 16; x++) for (int y = 1; y < 16; y++) okat &= bucket[x][y].size() == 1344;
    chk(okat, "1344 invertible N with Nx = y for each x, y != 0");
    // every generated form has exactly the intended radical, and the forms per radical are distinct
    bool okr = true; vector<u64> buf;
    for (size_t i = 0; i < GL4.size(); i++) {
        buf.clear(); forms_from(graphU[i].data(), 4, ND4, buf);
        set<u64> s(buf.begin(), buf.end()); okr &= s.size() == 28;
        for (u64 f : buf) { Mask m; okr &= is_alt(f) && radical(f, m) == 15 && meq(m, graphMask[i]); }
    }
    chk(okr, "graph(N) generator: 28 distinct forms with radical exactly graph(N), all N");
    okr = true;
    for (size_t i = 0; i < lines.size(); i += 7) {
        buf.clear(); forms_from(lines[i].ann, 6, ND6, buf);
        set<u64> s(buf.begin(), buf.end()); okr &= s.size() == 13888;
        for (u64 f : buf) { Mask m; okr &= rank8(f) == 6 && radical(f, m) == 3 && meq(m, lines[i].m); }
    }
    chk(okr, "line generator: 13888 distinct forms with radical exactly L (every 7th line)");
    okr = true;
    for (size_t i = 0; i < subs4.size(); i += 13) {
        buf.clear(); forms_from(subs4[i].ann, 4, ND4, buf);
        set<u64> s(buf.begin(), buf.end()); okr &= s.size() == 28;
        for (u64 f : buf) { Mask m; okr &= radical(f, m) == 15 && meq(m, subs4[i].m); }
    }
    chk(okr, "4-space generator: 28 distinct forms with radical exactly R (every 13th)");
    // the forms with radical exactly R are exactly 28 among all 2^28 forms? check by counting forms vanishing on R x V
    // for one R by brute force over all 2^28 codes
    {
        int idx = idxI(); long long cnt = 0, nb = 0; Mask target = graphMask[idx];
#pragma omp parallel for reduction(+:cnt, nb) schedule(static, 1 << 16)
        for (int c = 0; c < (1 << 28); c++) {
            u64 f = form28((u32)c); int rk = rank8(f);
            if (rk == 8) nb++;
            else if (rk == 4) { Mask m; radical(f, m); if (meq(m, target)) cnt++; }
        }
        chk(cnt == 28, "brute force over all 2^28 forms: exactly 28 have radical graph(I)");
        printf("     nondegenerate 8x8 alternating: %lld (expect 2^12*7*31*127 = 112881664)\n", nb);
        chk(nb == 112881664LL, "count of nondegenerate alternating forms on F_2^8");
    }
    // Desarguesian spread sets containing 0 and I = subfields F_16 of M_4(F_2): expect 336
    {
        set<vector<u16>> fields;
        for (u16 N : GL4) {
            if (N == 0x8421) continue;
            u16 I = 0x8421, N2 = mm4(N, N), N3 = mm4(N2, N);
            u16 gen[4] = { I, N, N2, N3 }; set<u16> F;
            for (int c = 0; c < 16; c++) { u16 s = 0; for (int k = 0; k < 4; k++) if (c >> k & 1) s ^= gen[k]; F.insert(s); }
            if (F.size() != 16) continue;
            bool ok = true;
            for (u16 a : F) { if (a && !inv4(a)) ok = false; for (u16 b : F) if (!F.count(mm4(a, b))) ok = false; }
            if (ok) fields.insert(vector<u16>(F.begin(), F.end()));
        }
        printf("     F_16 subfields of M_4(F_2): %zu\n", fields.size());
        chk(fields.size() == 336, "336 Desarguesian spread sets containing 0 and I");
        bool sp = true; for (auto& F : fields) for (u16 a : F) for (u16 b : F) if (a != b && !inv4(a ^ b)) sp = false;
        chk(sp, "they are spread sets (differences invertible)");
    }
    auto orbs = sp_orbits();
    printf("     Sp(J)-orbits on the 28 forms M (congruence):"); for (auto& o : orbs) printf(" %zu", o.size()); printf("\n");
    // the main problem's first two forms: B0 and graph(I) forms
    u64 B0 = formB0(); Mask m; chk(radical(B0, m) == 15 && mdisj(m, nonYmask) == false && mdisj(m, Ymask), "B0 radical is X");
    return bad;
}

static int mode_main(const char* progfile, bool all, int K, int r) {
    Config C = main_config();
    State S0 = empty_state(), S1, S2;
    if (!try_add(S0, formB0(), C.R, S1)) { printf("B0 illegal?\n"); return 1; }
    int iI = idxI();
    vector<vector<u32>> orbs = sp_orbits();
    vector<u32> Ms; vector<int> weight;
    if (all) { for (u32 c : ND4) { Ms.push_back(c); weight.push_back(1); } }
    else for (auto& o : orbs) { Ms.push_back(o[0]); weight.push_back((int)o.size()); }
    // tasks: (M, child at the second branching point), in a fixed order
    vector<State> tasks; vector<int> taskM;
    Stats pre;
    for (size_t mi = 0; mi < Ms.size(); mi++) {
        if (!try_add(S1, formGraph(iI, Ms[mi]), C.R, S2)) { printf("M %u: B0+B1 illegal\n", Ms[mi]); continue; }
        pre.nodes[S2.d]++;
        int a = first_uncovered(S2, C);
        size_t before = tasks.size();
        State T;
        for_each_candidate(a, S2.cov, C.Gnon, [&](u64 f) { pre.tries[2]++; if (try_add(S2, f, C.R, T)) { tasks.push_back(T); taskM.push_back((int)mi); } });
        printf("M %2u (orbit size %d): second point %02x, %zu dim-3 children\n", Ms[mi], weight[mi], a, tasks.size() - before);
    }
    fflush(stdout);
    // resume
    set<int> done; Stats tot = pre;
    if (FILE* fp = fopen(progfile, "r")) {
        char line[512];
        while (fgets(line, sizeof line, fp)) {
            int id; Stats s; long long n[10], t[10];
            if (sscanf(line, "T %d %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld", &id,
                &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7], &n[8], &t[0], &t[1], &t[2], &t[3], &t[4], &t[5], &t[6], &t[7], &t[8], &s.bentext, &s.full) == 21) {
                if (done.count(id)) continue;
                done.insert(id); for (int i = 0; i < 9; i++) { s.nodes[i] = n[i]; s.tries[i] = t[i]; } add_stats(tot, s);
            }
            else if (!strncmp(line, "SOLUTION", 8)) printf("(from previous run) %s", line);
        }
        fclose(fp);
    }
    printf("%zu tasks, %zu already done; this run: tasks with index %% %d == %d\n", tasks.size(), done.size(), K, r); fflush(stdout);
    FILE* pf = fopen(progfile, "a");
    mutex pmu; int nt = (int)tasks.size(); long long ndone = 0;
    auto t0 = chrono::steady_clock::now();
    vector<Stats> per(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nt; i++) {
        if (i % K != r || done.count(i)) continue;
        Stats s; dfs(tasks[i], C, s);
        lock_guard<mutex> g(pmu);
        add_stats(per[omp_get_thread_num()], s);
        fprintf(pf, "T %d", i); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.nodes[k]); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.tries[k]);
        fprintf(pf, " %lld %lld\n", s.bentext, s.full);
        if (s.full) { lock_guard<mutex> g2(solmu); for (auto& kv : solutions) fprintf(pf, "SOLUTION %s\n", kv.first.c_str()); }
        fflush(pf);
        if (++ndone % 500 == 0) { double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count(); printf("  %lld tasks done this run, %.0f s\n", ndone, el); fflush(stdout); }
    }
    fclose(pf);
    for (auto& p : per) add_stats(tot, p);
    print_stats("main", tot);
    printf("solutions (distinct spaces in these coordinates): %zu\n", solutions.size());
    for (auto& kv : solutions) {
        Verdict V = verify_space(kv.second);
        printf("  %s APN=%d dist k0..8:", kv.first.c_str(), V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("\n");
    }
    return 0;
}

// ---- known functions
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
static int mode_known(const char* file) {
    FILE* fp = fopen(file, "r"); if (!fp) { printf("cannot open %s\n", file); return 1; }
    static char buf[1 << 16]; int ln = 0; map<int, int> hist;
    while (fgets(buf, sizeof buf, fp)) {
        ln++; u8 F[256]; if (!read_lut(buf, F)) continue;
        Profile P = profile(forms_of_lut(F)); hist[P.n4]++;
        if (P.n4 >= 3) printf("line %d: rank-4 %d rank-6 %d bent %d\n", ln, P.n4, P.n6, P.n8);
    }
    for (auto& kv : hist) printf("rank-4 count %d: %d functions\n", kv.first, kv.second);
    return 0;
}
static bool load_line(const char* file, int line, vector<u64>& B) {
    FILE* fp = fopen(file, "r"); if (!fp) return false;
    static char buf[1 << 16]; int ln = 0;
    while (fgets(buf, sizeof buf, fp)) { if (++ln == line) { u8 F[256]; fclose(fp); if (!read_lut(buf, F)) return false; B = forms_of_lut(F); return true; } }
    fclose(fp); return false;
}
static u64 transform(u64 f, const u8 p[8]) { // B'(e_i, e_j) = B(p_i, p_j)
    u64 g = 0; for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (i != j && evalf(f, p[i], p[j])) g |= 1ULL << (8 * i + j); return g;
}
static void basis_of_mask(const Mask& m, int dim, u8* out) { vector<u8> c; for (int v = 1; v < 256; v++) if (mget(m, v)) c.push_back((u8)v); pick_independent(c, dim, out); }
static bool independent_of(const vector<u64>& basis, u64 f) {
    vector<u32> c; for (u64 b : basis) c.push_back(code28(b)); size_t r = rref28(c).size(); c.push_back(code28(f)); return rref28(c).size() > r;
}

// Control B: the main code path (B0 = J-form on X, graph(I), graph(N) generator, Y-line generator, main point order)
// on a relaxed problem that a known function solves: rank-6 radicals may also lie outside Y.
static int mode_ctlB(const char* file, int line, int s6, int hint) {
    vector<u64> B; if (!load_line(file, line, B)) { printf("cannot load\n"); return 1; }
    Profile P = profile(B);
    printf("known line %d: rank-4 %d rank-6 %d bent %d\n", line, P.n4, P.n6, P.n8);
    if (P.n4 < 2) { printf("need >= 2 rank-4 components\n"); return 1; }
    // Y: a 4-space complementary to every rank-4 radical, containing as many rank-6 radicals as possible
    int best = -1, bestc = -1, nY = 0;
    for (size_t i = 0; i < subs4.size(); i++) {
        bool ok = true; for (auto& r : P.r4) if (!mdisj(r, subs4[i].m)) { ok = false; break; }
        if (!ok) continue;
        nY++; int c = 0; for (auto& r : P.r6) { Mask t = r; bool in = true; for (int v = 1; v < 256; v++) if (mget(t, v) && !mget(subs4[i].m, v)) in = false; c += in; }
        if (c > bestc) { bestc = c; best = (int)i; }
    }
    printf("%d candidate Y; chosen one contains %d rank-6 radicals\n", nY, bestc);
    u8 p[8]; basis_of_mask(P.r4[0], 4, p); basis_of_mask(subs4[best].m, 4, p + 4);
    vector<u64> T; for (u64 f : B) T.push_back(transform(f, p));
    // normalise B0's form on Y to J
    Profile Q = profile(T); u64 B0 = 0;
    for (size_t i = 0; i < Q.f4.size(); i++) { bool isX = true; for (int v = 1; v < 16; v++) if (!mget(Q.r4[i], v)) isX = false; if (isX) B0 = Q.f4[i]; }
    if (!B0) { printf("no form with radical X\n"); return 1; }
    u8 Rd[4] = { 0 }; for (int k = 0; k < 4; k++) for (int l = 0; l < 4; l++) if (evalf(B0, (u8)(1 << (4 + k)), (u8)(1 << (4 + l)))) Rd[k] |= 1 << l;
    u32 dcode = rows4_to_code(Rd); u16 hh = 0; bool fh = false;
    for (u16 h : GL4) if (pullback4(dcode, h) == JCODE) { hh = h; fh = true; break; }
    if (!fh) { printf("no h\n"); return 1; }
    u8 p2[8]; for (int i = 0; i < 4; i++) { p2[i] = (u8)(1 << i); p2[4 + i] = (u8)(mv4(hh, (u8)(1 << i)) << 4); }
    for (auto& f : T) f = transform(f, p2);
    // second rank-4 radical -> graph(I)
    Q = profile(T); int i0 = -1; u64 fB0 = formB0();
    for (size_t i = 0; i < Q.f4.size(); i++) if (Q.f4[i] == fB0) i0 = (int)i;
    if (i0 < 0) { printf("B0 not normalised\n"); return 1; }
    int i1 = (i0 == 0) ? 1 : 0;
    u8 col[4]; for (int i = 0; i < 4; i++) { col[i] = 0; for (int v = 1; v < 256; v++) if (mget(Q.r4[i1], v) && (v & 15) == (1 << i)) col[i] = (u8)(v >> 4); }
    u16 N2 = fromcols4(col); if (!inv4(N2)) { printf("N2 singular\n"); return 1; }
    u16 N2inv = 0; for (u16 g : GL4) if (mm4(g, N2) == 0x8421) N2inv = g;
    u8 p3[8]; for (int i = 0; i < 4; i++) { p3[i] = mv4(N2inv, (u8)(1 << i)); p3[4 + i] = (u8)(1 << (4 + i)); }
    for (auto& f : T) f = transform(f, p3);
    Q = profile(T);
    int iI = idxI(); u32 M1 = 0xffff; u64 fB1 = 0;
    for (u32 c : ND4) { u64 f = formGraph(iI, c); if (find(Q.f4.begin(), Q.f4.end(), f) != Q.f4.end()) { M1 = c; fB1 = f; } }
    if (M1 == 0xffff) { printf("graph(I) form not found\n"); return 1; }
    string key = space_key(T);
    Verdict V0 = verify_space(T);
    printf("transformed known space: B0 = J-form, B1 = graph(I) with M = %u; APN=%d dist:", M1, V0.apn); for (int k = 0; k <= 8; k++) printf(" %d", V0.dist[k]); printf("\n  key %s\n", key.c_str());

    Config C = main_config(); C.R.r6out = true;
    C.Gnon = { true, false, false, true }; C.GY = { false, false, true, true };
    if (s6 < 0) { // control C: outside-Y rank-6 candidates from the known list, everything else in full
        C.Gnon = { true, false, false, false, true }; C.GY = { false, false, true, false, true };
        oracle6 = Q.f6; oracle6m = Q.r6;
        printf("control C: rank-6 forms outside Y restricted to the known function's; rank-4 and in-Y rank-6 generated in full\n");
    }
    State S = empty_state(), S2;
    if (!try_add(S, formB0(), C.R, S2)) { printf("seed B0 fail\n"); return 1; } S = S2;
    if (!try_add(S, fB1, C.R, S2)) { printf("seed B1 fail\n"); return 1; } S = S2;
    // seed s6 rank-6 forms of the known space (independent of the span); hint=1 also seeds the third rank-4 form
    int added = 0;
    if (hint) for (u64 f : Q.f4) if (f != fB0 && f != fB1 && independent_of(S.basis, f)) { if (try_add(S, f, C.R, S2)) { S = S2; break; } }
    for (u64 f : Q.f6) { if (added >= s6) break; if (!independent_of(S.basis, f)) continue; if (!try_add(S, f, C.R, S2)) { printf("seed fail\n"); return 1; } S = S2; added++; }
    printf("seed dimension %d (B0, B1, %s%d rank-6)\n", S.d, hint ? "a rank-4, " : "", added); fflush(stdout);
    auto t0 = chrono::steady_clock::now();
    Stats st = run_search(S, C);
    print_stats("ctlB", st);
    printf("ctlB: %zu spaces found in %.1f s; known space among them: %s\n", solutions.size(), chrono::duration<double>(chrono::steady_clock::now() - t0).count(), solutions.count(key) ? "YES" : "NO");
    int napn = 0;
    for (auto& kv : solutions) { Verdict V = verify_space(kv.second); napn += V.apn; printf("  APN=%d dist:", V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("%s\n", kv.first == key ? "  <- known" : ""); }
    printf("ctlB: %d of %zu verified APN by LUT\n", napn, solutions.size());
    return 0;
}

// Control A: fully generic rules and generators (rank 4 and rank 6 anywhere), seeded with all rank-4 forms
static int mode_ctlA(const char* file, int line, int s6) {
    vector<u64> B; if (!load_line(file, line, B)) { printf("cannot load\n"); return 1; }
    Profile P = profile(B); string key = space_key(B);
    printf("known line %d: rank-4 %d rank-6 %d bent %d\n", line, P.n4, P.n6, P.n8);
    Config C; C.R = { true, false, true, true }; C.Gnon = C.GY = C.GX = { false, true, false, true };
    for (int v = 1; v < 256; v++) C.order.push_back(v);
    State S = empty_state(), S2; int added = 0;
    for (u64 f : P.f4) if (independent_of(S.basis, f)) { if (!try_add(S, f, C.R, S2)) { printf("seed fail\n"); return 1; } S = S2; }
    for (u64 f : P.f6) { if (added >= s6) break; if (!independent_of(S.basis, f)) continue; if (!try_add(S, f, C.R, S2)) { printf("seed fail\n"); return 1; } S = S2; added++; }
    printf("seed dimension %d\n", S.d); fflush(stdout);
    auto t0 = chrono::steady_clock::now();
    Stats st = run_search(S, C);
    print_stats("ctlA", st);
    printf("ctlA: %zu spaces found in %.1f s; known space among them: %s\n", solutions.size(), chrono::duration<double>(chrono::steady_clock::now() - t0).count(), solutions.count(key) ? "YES" : "NO");
    int napn = 0;
    for (auto& kv : solutions) { Verdict V = verify_space(kv.second); napn += V.apn; printf("  APN=%d dist:", V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("%s\n", kv.first == key ? "  <- known" : ""); }
    printf("ctlA: %d of %zu verified APN by LUT\n", napn, solutions.size());
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------
// Case i = 2, [0^230, 2^10, 4^15]: holes = (X ∪ Y) \ 0, rank-6 radicals are lines of X or of Y, the 15 rank-4
// radicals miss X and Y (graphs of N with N, N + I invertible once one of them is graph(I)).
// Start: B1 = J(x+y, x'+y') on graph(I). G = {diag(h,h) : h in Sp(J)} x <swap (x,y) <-> (y,x)> (order 1440) fixes
// X ∪ Y, graph(I) and B1. B2 runs over G-orbit representatives of the other candidate rank-4 forms; from rep j, every
// rank-4 element of the span must have orbit index >= j (apply the g taking the solution's minimal-index form to its rep).
struct I2 {
    vector<u64> cand; vector<int> candN; vector<u32> candM; unordered_map<u64, int> idx;
    vector<int> orbOf, reps, orbSize; unordered_map<u64, int> orbIdx;
    vector<array<u8, 8>> group;
};
static I2 build_i2() {
    I2 T;
    for (size_t i = 0; i < GL4.size(); i++) {
        u16 N = GL4[i]; if (N == 0x8421 || !inv4(N ^ 0x8421)) continue;
        vector<u64> b; forms_from(graphU[i].data(), 4, ND4, b);
        for (size_t k = 0; k < b.size(); k++) { T.idx[b[k]] = (int)T.cand.size(); T.cand.push_back(b[k]); T.candN.push_back((int)i); T.candM.push_back(ND4[k]); }
    }
    for (u16 h : symplectic()) {
        array<u8, 8> p, q;
        for (int i = 0; i < 4; i++) { u8 he = mv4(h, (u8)(1 << i)); p[i] = he; p[4 + i] = (u8)(he << 4); q[i] = (u8)(he << 4); q[4 + i] = he; }
        T.group.push_back(p); T.group.push_back(q);
    }
    T.orbOf.assign(T.cand.size(), -1);
    for (size_t c = 0; c < T.cand.size(); c++) if (T.orbOf[c] < 0) {
        int o = (int)T.reps.size(); T.reps.push_back((int)c); int sz = 0;
        for (auto& g : T.group) {
            auto it = T.idx.find(transform(T.cand[c], g.data()));
            if (it == T.idx.end()) { fprintf(stderr, "group image is not a candidate\n"); exit(1); }
            if (T.orbOf[it->second] < 0) { T.orbOf[it->second] = o; sz++; }
            else if (T.orbOf[it->second] != o) { fprintf(stderr, "orbit clash\n"); exit(1); }
        }
        T.orbSize.push_back(sz);
    }
    for (size_t c = 0; c < T.cand.size(); c++) T.orbIdx[T.cand[c]] = T.orbOf[c];
    return T;
}
static Config i2_config() {
    Config C; C.R = { true, true, true, false }; C.R.holesXY = true;
    C.Gnon = { true, false, false, false }; C.GY = { false, false, true, false }; C.GX = { false, false, false, false, false, true };
    for (int v = 1; v < 256; v++) if ((v & 15) && (v >> 4)) C.order.push_back(v);
    if (g_desc) reverse(C.order.begin(), C.order.end());
    for (int v = 1; v < 256; v++) if (!(v & 15)) C.order.push_back(v);
    for (int v = 1; v < 256; v++) if (!(v >> 4)) C.order.push_back(v);
    return C;
}
static int i2_checks(const I2& T, u64 B1) {
    int bad = 0;
    auto chk = [&](bool ok, const char* what) { printf("%s %s\n", ok ? "ok  " : "FAIL", what); if (!ok) bad++; };
    long long nN = 0; for (u16 N : GL4) if (N != 0x8421 && inv4(N ^ 0x8421)) nN++;
    printf("     N with N, N+I invertible (N != I): %lld; candidates %zu; orbits %zu\n", nN, T.cand.size(), T.reps.size());
    chk(T.cand.size() == (size_t)nN * 28, "candidates = 28 per N");
    long long sum = 0; bool div = true; map<int, int> hist;
    for (int s : T.orbSize) { sum += s; div &= 1440 % s == 0; hist[s]++; }
    chk(sum == (long long)T.cand.size() && div, "orbit sizes sum to the candidate count and divide 1440");
    printf("     orbit sizes:"); for (auto& kv : hist) printf(" %dx%d", kv.second, kv.first); printf("\n");
    set<array<u8, 8>> G(T.group.begin(), T.group.end()); chk(G.size() == 1440, "1440 distinct group elements");
    bool fix = true;
    for (auto& g : T.group) {
        fix &= transform(B1, g.data()) == B1;
        for (int i = 0; i < 8; i++) { u8 v = g[i]; bool inX = (v >> 4) == 0, inY = (v & 15) == 0; if (!(inX || inY)) fix = false; }
    }
    chk(fix, "every group element fixes B1 and maps basis vectors into X u Y");
    return bad;
}
static void write_stats(FILE* pf, int j, const Stats& s) {
    fprintf(pf, "R %d", j); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.nodes[k]); for (int k = 0; k < 9; k++) fprintf(pf, " %lld", s.tries[k]);
    fprintf(pf, " %lld %lld\n", s.bentext, s.full); fflush(pf);
}
static int mode_i2(const char* progfile, int first, int last) {
    I2 T = build_i2();
    u64 B1 = formGraph(idxI(), JCODE);
    if (i2_checks(T, B1)) return 1;
    Config C = i2_config(); C.R.orb = &T.orbIdx;
    Rules r0 = C.R; r0.orb = nullptr;
    if (g_nopr) C.R.orb = nullptr;
    printf("i2: %s point order, orbit-index pruning %s\n", g_desc ? "descending" : "ascending", g_nopr ? "off" : "on");
    set<int> done; Stats tot;
    if (FILE* fp = fopen(progfile, "r")) {
        char line[512];
        while (fgets(line, sizeof line, fp)) {
            int id; Stats s; long long n[10], t[10];
            if (sscanf(line, "R %d %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld", &id,
                &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7], &n[8], &t[0], &t[1], &t[2], &t[3], &t[4], &t[5], &t[6], &t[7], &t[8], &s.bentext, &s.full) == 21) {
                if (done.count(id)) continue;
                done.insert(id); for (int i = 0; i < 9; i++) { s.nodes[i] = n[i]; s.tries[i] = t[i]; } add_stats(tot, s);
            }
            else if (!strncmp(line, "SOLUTION", 8)) printf("(from previous run) %s", line);
        }
        fclose(fp);
    }
    last = min(last, (int)T.reps.size() - 1);
    printf("reps %d..%d, %zu already done\n", first, last, done.size()); fflush(stdout);
    FILE* pf = fopen(progfile, "a");
    State S0 = empty_state(), S1;
    if (!try_add(S0, B1, r0, S1)) { printf("B1 illegal?\n"); return 1; }
    auto t0 = chrono::steady_clock::now();
    for (int j = first; j <= last; j++) {
        if (done.count(j)) continue;
        Stats st; State S2; C.R.minidx = j;
        int c = T.reps[j];
        if (try_add(S1, T.cand[c], C.R, S2)) { st = run_search(S2, C, 400); }
        else printf("  rep %d: B1 + B2 illegal\n", j);
        { lock_guard<mutex> g(solmu); for (auto& kv : solutions) fprintf(pf, "SOLUTION %s\n", kv.first.c_str()); }
        write_stats(pf, j, st); add_stats(tot, st);
        printf("rep %d (N %04x, M %u, orbit %d): nodes", j, GL4[T.candN[c]], T.candM[c], T.orbSize[j]); for (int k = 0; k <= 8; k++) printf(" %lld", st.nodes[k]);
        printf("  full %lld  [%.0f s]\n", st.full, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); fflush(stdout);
    }
    fclose(pf);
    print_stats("i2", tot);
    if (g_logp1) printf("this run: nodes with every non-hole point covered %lld, distinct spaces %zu\n", p1nodes, p1keys.size());
    printf("solutions: %zu\n", solutions.size());
    for (auto& kv : solutions) { Verdict V = verify_space(kv.second); printf("  %s APN=%d dist:", kv.first.c_str(), V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("\n"); }
    return 0;
}
// children at the first branching point from (B1, rep j), no orbit pruning: for the numpy recount
static int mode_i2count(int j) {
    I2 T = build_i2(); u64 B1 = formGraph(idxI(), JCODE);
    Config C = i2_config(); State S0 = empty_state(), S1, S2;
    try_add(S0, B1, C.R, S1);
    int c = T.reps[j];
    if (!try_add(S1, T.cand[c], C.R, S2)) { printf("rep %d N %04x M %u: illegal pair\n", j, GL4[T.candN[c]], T.candM[c]); return 0; }
    int a = first_uncovered(S2, C); long long n = 0; State X;
    for_each_candidate(a, S2.cov, gen_for(a, C), [&](u64 f) { if (try_add(S2, f, C.R, X)) n++; });
    printf("rep %d N %04x M %u: point %02x, %lld dim-3 children\n", j, GL4[T.candN[c]], T.candM[c], a, n);
    return 0;
}
// Control: a known function with two complementary 4-spaces, each split by 5 of its rank-6 radicals and missing
// every rank-4 radical, brought into the i = 2 coordinates, through the i = 2 code path (generators, rules, orbit
// table and pruning), relaxed only by allowing rank-6 radicals outside X ∪ Y.
static int mode_i2ctl(const char* file, int line, int s6) {
    vector<u64> B; if (!load_line(file, line, B)) { printf("cannot load\n"); return 1; }
    Profile P = profile(B);
    printf("known line %d: rank-4 %d rank-6 %d bent %d\n", line, P.n4, P.n6, P.n8);
    vector<int> hs;
    for (size_t i = 0; i < subs4.size(); i++) {
        bool ok = true; for (auto& r : P.r4) if (!mdisj(r, subs4[i].m)) { ok = false; break; }
        if (!ok) continue;
        int c = 0; for (auto& r : P.r6) { bool in = true; for (int v = 1; v < 256; v++) if (mget(r, v) && !mget(subs4[i].m, v)) in = false; c += in; }
        if (c == 5) hs.push_back((int)i);
    }
    int pi = -1, pk = -1;
    for (size_t a = 0; a < hs.size() && pi < 0; a++) for (size_t b = a + 1; b < hs.size(); b++) if (mdisj(subs4[hs[a]].m, subs4[hs[b]].m)) { pi = hs[a]; pk = hs[b]; break; }
    printf("%zu 4-spaces split by 5 rank-6 radicals and missing the rank-4 ones; complementary pair: %s\n", hs.size(), pi >= 0 ? "yes" : "no");
    if (pi < 0) return 2;
    u8 p[8]; basis_of_mask(subs4[pi].m, 4, p); basis_of_mask(subs4[pk].m, 4, p + 4);
    vector<u64> T; for (u64 f : B) T.push_back(transform(f, p));
    Profile Q = profile(T);
    u8 col[4]; for (int i = 0; i < 4; i++) { col[i] = 0; for (int v = 1; v < 256; v++) if (mget(Q.r4[0], v) && (v & 15) == (1 << i)) col[i] = (u8)(v >> 4); }
    u16 N1 = fromcols4(col); if (!inv4(N1)) { printf("N1 singular\n"); return 1; }
    u16 N1inv = 0; for (u16 g : GL4) if (mm4(g, N1) == 0x8421) N1inv = g;
    u8 p3[8]; for (int i = 0; i < 4; i++) { p3[i] = mv4(N1inv, (u8)(1 << i)); p3[4 + i] = (u8)(1 << (4 + i)); }
    for (auto& f : T) f = transform(f, p3);
    Q = profile(T); int iI = idxI(); u32 M1 = 0xffff;
    for (u32 c : ND4) if (find(Q.f4.begin(), Q.f4.end(), formGraph(iI, c)) != Q.f4.end()) M1 = c;
    if (M1 == 0xffff) { printf("graph(I) form not found\n"); return 1; }
    u16 gg = 0; for (u16 g : GL4) if (pullback4(M1, g) == JCODE) { gg = g; break; }
    u8 p4[8]; for (int i = 0; i < 4; i++) { u8 ge = mv4(gg, (u8)(1 << i)); p4[i] = ge; p4[4 + i] = (u8)(ge << 4); }
    for (auto& f : T) f = transform(f, p4);
    Q = profile(T); u64 B1 = formGraph(iI, JCODE);
    if (find(Q.f4.begin(), Q.f4.end(), B1) == Q.f4.end()) { printf("B1 not normalised\n"); return 1; }
    I2 TB = build_i2();
    int jmin = INT_MAX; u64 fmin = 0;
    for (u64 f : Q.f4) if (f != B1) { auto it = TB.orbIdx.find(f); if (it == TB.orbIdx.end()) { printf("rank-4 form not a candidate\n"); return 1; } if (it->second < jmin) { jmin = it->second; fmin = f; } }
    u64 rep = TB.cand[TB.reps[jmin]]; bool moved = false;
    for (auto& g : TB.group) if (transform(fmin, g.data()) == rep) { for (auto& f : T) f = transform(f, g.data()); moved = true; break; }
    if (!moved) { printf("no group element to the rep\n"); return 1; }
    Q = profile(T);
    string key = space_key(T); Verdict V0 = verify_space(T);
    printf("normalised: X, Y split by rank-6 radicals, B1 = J on graph(I), min orbit index %d (orbit size %d); APN=%d dist:", jmin, TB.orbSize[jmin], V0.apn);
    for (int k = 0; k <= 8; k++) printf(" %d", V0.dist[k]); printf("\n  key %s\n", key.c_str());
    Config C = i2_config(); C.R.r6out = true; C.R.orb = &TB.orbIdx; C.R.minidx = jmin;
    C.Gnon = { true, false, false, true }; C.GY = { false, false, true, true }; C.GX = { false, false, false, true, false, true };
    Rules r0 = C.R; r0.orb = nullptr;
    State S = empty_state(), S2;
    if (!try_add(S, B1, r0, S2)) { printf("seed B1 fail\n"); return 1; } S = S2;
    if (!try_add(S, rep, C.R, S2)) { printf("seed rep fail\n"); return 1; } S = S2;
    int added = 0;
    for (u64 f : Q.f6) { if (added >= s6) break; if (!independent_of(S.basis, f)) continue; if (!try_add(S, f, C.R, S2)) { printf("seed fail\n"); return 1; } S = S2; added++; }
    printf("seed dimension %d (B1, rep, %d rank-6)\n", S.d, added); fflush(stdout);
    auto t0 = chrono::steady_clock::now();
    Stats st = run_search(S, C);
    print_stats("i2ctl", st);
    printf("i2ctl: %zu spaces found in %.1f s; known space among them: %s\n", solutions.size(), chrono::duration<double>(chrono::steady_clock::now() - t0).count(), solutions.count(key) ? "YES" : "NO");
    int napn = 0;
    for (auto& kv : solutions) { Verdict V = verify_space(kv.second); napn += V.apn; printf("  APN=%d dist:", V.apn); for (int k = 0; k <= 8; k++) printf(" %d", V.dist[k]); printf("%s\n", kv.first == key ? "  <- known" : ""); }
    printf("i2ctl: %d of %zu verified APN by LUT\n", napn, solutions.size());
    return 0;
}

// Bring a space whose rank-4 radicals miss two complementary 4-spaces Xm, Ym into the i = 2 search's coordinates:
// Xm -> X, Ym -> Y, its first rank-4 radical -> graph(I) with form J, the minimal-orbit-index other rank-4 form -> its
// representative. Returns the representative's index, or -1.
static int normalise_i2(vector<u64>& T, const Mask& Xm, const Mask& Ym, const I2& TB) {
    u8 p[8]; basis_of_mask(Xm, 4, p); basis_of_mask(Ym, 4, p + 4);
    for (auto& f : T) f = transform(f, p);
    Profile Q = profile(T);
    u8 col[4]; for (int i = 0; i < 4; i++) { col[i] = 0; for (int v = 1; v < 256; v++) if (mget(Q.r4[0], v) && (v & 15) == (1 << i)) col[i] = (u8)(v >> 4); }
    u16 N1 = fromcols4(col); if (!inv4(N1)) return -1;
    u16 N1inv = 0; for (u16 g : GL4) if (mm4(g, N1) == 0x8421) N1inv = g;
    u8 p3[8]; for (int i = 0; i < 4; i++) { p3[i] = mv4(N1inv, (u8)(1 << i)); p3[4 + i] = (u8)(1 << (4 + i)); }
    for (auto& f : T) f = transform(f, p3);
    Q = profile(T); int iI = idxI(); u32 M1 = 0xffff;
    for (u32 c : ND4) if (find(Q.f4.begin(), Q.f4.end(), formGraph(iI, c)) != Q.f4.end()) M1 = c;
    if (M1 == 0xffff) return -1;
    u16 gg = 0; for (u16 g : GL4) if (pullback4(M1, g) == JCODE) { gg = g; break; }
    u8 p4[8]; for (int i = 0; i < 4; i++) { u8 ge = mv4(gg, (u8)(1 << i)); p4[i] = ge; p4[4 + i] = (u8)(ge << 4); }
    for (auto& f : T) f = transform(f, p4);
    Q = profile(T); u64 B1 = formGraph(iI, JCODE);
    if (find(Q.f4.begin(), Q.f4.end(), B1) == Q.f4.end()) return -1;
    int jmin = INT_MAX; u64 fmin = 0;
    for (u64 f : Q.f4) if (f != B1) { auto it = TB.orbIdx.find(f); if (it == TB.orbIdx.end()) return -1; if (it->second < jmin) { jmin = it->second; fmin = f; } }
    u64 rep = TB.cand[TB.reps[jmin]];
    for (auto& g : TB.group) if (transform(fmin, g.data()) == rep) { for (auto& f : T) f = transform(f, g.data()); return jmin; }
    return -1;
}
// Rowen's "SHORT" spans (runs/a/fifteen.txt): does the i = 2 search pass through each of them (normalised) at the
// moment every point off X u Y is covered?
static int mode_i2probe(const char* file) {
    FILE* fp = fopen(file, "r"); if (!fp) return 1;
    char line[4096]; vector<pair<string, int>> want; I2 TB = build_i2();
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, "SHORT", 5)) continue;
        const char* b = strstr(line, "basis"); if (!b) continue; b += 5;
        vector<u64> F; char* e; while (true) { long long c = strtoll(b, &e, 16); if (e == b) break; F.push_back(form28((u32)c)); b = e; }
        Profile P = profile(F);
        Mask cov; mclear(cov); for (auto& r : P.r4) mor(cov, r);
        vector<int> hole; for (int v = 1; v < 256; v++) if (!mget(cov, v)) hole.push_back(v);
        Mask A, Bm; mclear(A); mclear(Bm);
        int x0 = hole[0]; mset(A, x0);
        for (int y : hole) if (y != x0 && find(hole.begin(), hole.end(), x0 ^ y) != hole.end()) mset(A, y);
        for (int y : hole) if (!mget(A, y)) mset(Bm, y);
        printf("span: dim %zu, rank-4 %d, rank-6 %d, bent %d; holes %zu = %d + %d\n", F.size(), P.n4, P.n6, P.n8, hole.size(), mcount(A), mcount(Bm));
        if (P.n4 != 15 || mcount(A) != 15 || mcount(Bm) != 15) { printf("  unexpected shape\n"); continue; }
        int j = normalise_i2(F, A, Bm, TB);
        if (j < 0) { printf("  normalisation failed\n"); continue; }
        want.push_back({ space_key(F), j });
        printf("  normalised; minimal orbit index %d\n", j);
    }
    fclose(fp);
    set<int> reps; for (auto& w : want) reps.insert(w.second);
    g_logp1 = true;
    Config C = i2_config(); C.R.orb = &TB.orbIdx; Rules r0 = C.R; r0.orb = nullptr;
    State S0 = empty_state(), S1; try_add(S0, formGraph(idxI(), JCODE), r0, S1);
    for (int j : reps) { State S2; C.R.minidx = j; if (try_add(S1, TB.cand[TB.reps[j]], C.R, S2)) run_search(S2, C, 400); }
    printf("searched reps:"); for (int j : reps) printf(" %d", j);
    printf("\nnodes with every non-hole point covered: %lld (distinct spaces %zu)\n", p1nodes, p1keys.size());
    int hit = 0; for (auto& w : want) { bool h = p1keys.count(w.first) > 0; hit += h; printf("  rep %d: %s\n", w.second, h ? "visited" : "NOT VISITED"); }
    printf("%d of %zu of Rowen's spans visited by the i = 2 search\n", hit, want.size());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { printf("modes: selftest | main <progress> [all] [K r] | known <file> | ctlB <file> <line> <s6> [hint] | ctlA <file> <line> <s6> | i2 <progress> <first> <last> | i2count <j> | i2ctl <file> <line> <s6>\n"); return 1; }
    string mode = argv[1];
    init_code28();
    init_tables(mode != "main");
    if (mode == "selftest") return selftest();
    if (mode == "main") {
        bool all = false; vector<int> nums;
        for (int i = 3; i < argc; i++) {
            if (!strcmp(argv[i], "all")) all = true; else if (!strcmp(argv[i], "desc")) g_desc = true; else nums.push_back(atoi(argv[i]));
        }
        int K = nums.size() >= 2 ? nums[0] : 1, r = nums.size() >= 2 ? nums[1] : 0;
        printf("main: %s M, %s point order\n", all ? "all 28" : "Sp(J)-orbit representatives of", g_desc ? "descending" : "ascending");
        return mode_main(argv[2], all, K, r);
    }
    if (mode == "known") return mode_known(argv[2]);
    if (mode == "ctlB") return mode_ctlB(argv[2], atoi(argv[3]), atoi(argv[4]), argc > 5 ? atoi(argv[5]) : 0);
    if (mode == "ctlA") return mode_ctlA(argv[2], atoi(argv[3]), atoi(argv[4]));
    if (mode == "i2") {
        for (int i = 5; i < argc; i++) { if (!strcmp(argv[i], "desc")) g_desc = true; if (!strcmp(argv[i], "nopr")) g_nopr = true; if (!strcmp(argv[i], "p1")) g_logp1 = true; }
        return mode_i2(argv[2], argc > 3 ? atoi(argv[3]) : 0, argc > 4 ? atoi(argv[4]) : INT_MAX);
    }
    if (mode == "i2count") return mode_i2count(atoi(argv[2]));
    if (mode == "i2probe") return mode_i2probe(argv[2]);
    if (mode == "i2ctl") return mode_i2ctl(argv[2], atoi(argv[3]), atoi(argv[4]));
    return 1;
}
