// review_odd_brute: reviewer's independent check (Claude Opus 5.5, 2026-10-24). Shares no code with apn8/symp/orbit5.
// For an element A of odd prime order p in GL(8,2), enumerate EVERY A-invariant 8-dim subspace S of Alt(F_2^8)
// (no case analysis, no B1, no orbit reduction): S = sum over the isotypic parts of F_q-subspaces, one reduced-echelon
// basis per part. Each is tested for APN (rank 7 of B -> B(a,.) at every a != 0); APN ones get their rank distribution
// and how many rank-4 components are fixed by A vs moved. Pruning: only "a subspace of an APN space has rank >= dim-1
// at every a" (applied to each part's subspace alone and to the partial sums), which is necessary.
//
// usage: review_odd_brute <A> [conj seed] [count]
//   <A>: c5a, c5b, c7a, c7b, c7c (companion-matrix representatives built here), or 16 hex digits (columns, 2 each)
//   conj seed: if nonzero, replace A by P A P^-1 for a random P (counts must not change)
//   count: only print the enumeration sizes against Gaussian binomials, no APN tests
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <atomic>
#include <mutex>
#include <functional>
#include <algorithm>
#include <intrin.h>
#include <omp.h>
typedef uint8_t u8; typedef uint32_t u32; typedef uint64_t u64;

struct Mat { u8 c[8]; };
static u8 mv(const Mat& A, u8 x) { u8 r = 0; for (int i = 0; i < 8; i++) if (x >> i & 1) r ^= A.c[i]; return r; }
static Mat mmul(const Mat& A, const Mat& B) { Mat C; for (int i = 0; i < 8; i++) C.c[i] = mv(A, B.c[i]); return C; }
static Mat ident() { Mat I; for (int i = 0; i < 8; i++) I.c[i] = (u8)(1 << i); return I; }
static bool inverse(const Mat& A, Mat& R) {
    int t[256]; for (int i = 0; i < 256; i++) t[i] = -1;
    for (int x = 0; x < 256; x++) { u8 y = mv(A, (u8)x); if (t[y] >= 0) return false; t[y] = x; }
    for (int i = 0; i < 8; i++) R.c[i] = (u8)t[1 << i]; return true;
}

// forms: 28-bit coefficient vectors (bit k = pair (KI[k], KJ[k])) and 8-byte row matrices
static int KI[28], KJ[28];
static u64 C2R[4][256];
static void init_tables() {
    int k = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { KI[k] = i; KJ[k] = j; k++; }
    for (int b = 0; b < 4; b++) for (int x = 0; x < 256; x++) {
        u64 m = 0; for (int t = 0; t < 8; t++) if ((x >> t & 1) && 8 * b + t < 28) { int q = 8 * b + t; m ^= 1ull << (8 * KI[q] + KJ[q]); m ^= 1ull << (8 * KJ[q] + KI[q]); }
        C2R[b][x] = m;
    }
}
static inline u64 rows_of(u32 c) { return C2R[0][c & 255] ^ C2R[1][c >> 8 & 255] ^ C2R[2][c >> 16 & 255] ^ C2R[3][c >> 24 & 255]; }
static inline u8 rv(u64 m, unsigned a) { u8 r = 0; for (int i = 0; i < 8; i++) if (a >> i & 1) r ^= (u8)(m >> (8 * i)); return r; }
static inline int bval(u64 m, unsigned x, unsigned y) { return __popcnt(rv(m, x) & y) & 1; }
static int rank8(u64 m) {
    u8 r[8]; for (int i = 0; i < 8; i++) r[i] = (u8)(m >> (8 * i));
    int rk = 0; for (int b = 0; b < 8; b++) { int p = -1; for (int i = rk; i < 8; i++) if (r[i] >> b & 1) { p = i; break; } if (p < 0) continue; std::swap(r[p], r[rk]); for (int i = 0; i < 8; i++) if (i != rk && (r[i] >> b & 1)) r[i] ^= r[rk]; rk++; }
    return rk;
}

// linear map on coefficient space: T(B) = B o A
struct Lin { u32 col[28]; u32 tab[4][256];
    void build() { for (int b = 0; b < 4; b++) for (int x = 0; x < 256; x++) { u32 r = 0; for (int t = 0; t < 8; t++) if ((x >> t & 1) && 8 * b + t < 28) r ^= col[8 * b + t]; tab[b][x] = r; } }
    u32 operator()(u32 c) const { return tab[0][c & 255] ^ tab[1][c >> 8 & 255] ^ tab[2][c >> 16 & 255] ^ tab[3][c >> 24 & 255]; } };
static u32 coef_of_rows(u64 m) { u32 c = 0; for (int k = 0; k < 28; k++) if (m >> (8 * KI[k] + KJ[k]) & 1) c |= 1u << k; return c; }
static Lin action(const Mat& A) {
    Lin L;
    for (int k = 0; k < 28; k++) {
        u64 m = rows_of(1u << k), r = 0;
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (i != j && bval(m, A.c[i], A.c[j])) r |= 1ull << (8 * i + j);
        L.col[k] = coef_of_rows(r);
    }
    L.build(); return L;
}

// F_2 echelon on 28-bit vectors
struct Ech { u32 piv[28] = {0}; int n = 0;
    u32 reduce(u32 x) const { for (int b = 27; b >= 0 && x; b--) if ((x >> b & 1) && piv[b]) x ^= piv[b]; return x; }
    bool add(u32 x) { x = reduce(x); if (!x) return false; piv[31 - __lzcnt(x)] = x; n++; return true; } };
static std::vector<u32> kernel_of(const u32* col) {
    // rows of the matrix with tags; kernel = tags of rows reducing to 0
    u32 piv[28] = {0}, ptag[28] = {0}; std::vector<u32> ker;
    for (int k = 0; k < 28; k++) {
        u32 x = col[k], t = 1u << k;
        for (int b = 27; b >= 0 && x; b--) if ((x >> b & 1) && piv[b]) { x ^= piv[b]; t ^= ptag[b]; }
        if (!x) ker.push_back(t); else { int b = 31 - __lzcnt(x); piv[b] = x; ptag[b] = t; }
    }
    return ker;
}

struct Part { int d, m, q; u32 g; std::vector<u32> nb; std::vector<u32> tab; // tab[i*q + c] = c(T) nb[i]
};
static Lin T;
static std::vector<Part> parts;
static int P_ord;

static u32 polymod_xp(int p, u32 g, int d) { u32 r = 1; for (int i = 0; i < p; i++) { r <<= 1; if (r >> d & 1) r ^= g; } return r; }
static bool decompose(const Mat& A) {
    T = action(A);
    // order
    Mat P = A; int p = 1; while (memcmp(P.c, ident().c, 8) && p < 300) { P = mmul(A, P); p++; }
    P_ord = p;
    int d = 1; while ((1u << d) % (unsigned)p != 1) d++;
    std::vector<std::pair<u32, int>> gs; gs.push_back({3u, 1});
    for (u32 g = 1u << d; g < (2u << d); g++) { if (!(g & 1) || !(__popcnt(g) & 1)) continue; if (polymod_xp(p, g, d) == 1) gs.push_back({g, d}); }
    std::vector<u32> Tp[16]; // T^t columns
    u32 cols[17][28]; for (int k = 0; k < 28; k++) cols[0][k] = 1u << k;
    for (int t = 1; t <= 16; t++) for (int k = 0; k < 28; k++) cols[t][k] = T(cols[t - 1][k]);
    Ech all; int tot = 0;
    for (auto& gd : gs) {
        u32 g = gd.first; int dg = gd.second;
        u32 M[28]; for (int k = 0; k < 28; k++) { u32 r = 0; for (int t = 0; t <= dg; t++) if (g >> t & 1) r ^= cols[t][k]; M[k] = r; }
        std::vector<u32> ker = kernel_of(M);
        if (ker.empty()) continue;
        Part pt; pt.d = dg; pt.q = 1 << dg; pt.g = g;
        Ech sp;
        for (u32 v : ker) { if (!sp.reduce(v)) continue; pt.nb.push_back(v); u32 w = v; for (int t = 0; t < dg; t++) { if (t) w = T(w); if (!sp.add(w)) { fprintf(stderr, "F_q basis\n"); return false; } } }
        if (sp.n != (int)ker.size()) { fprintf(stderr, "part dim\n"); return false; }
        pt.m = (int)pt.nb.size();
        pt.tab.assign((size_t)pt.m * pt.q, 0);
        for (int i = 0; i < pt.m; i++) for (int c = 0; c < pt.q; c++) { u32 r = 0, w = pt.nb[i]; for (int t = 0; t < dg; t++) { if (t) w = T(w); if (c >> t & 1) r ^= w; } pt.tab[(size_t)i * pt.q + c] = r; }
        for (int i = 0; i < pt.m; i++) { u32 w = pt.nb[i]; for (int t = 0; t < dg; t++) { if (t) w = T(w); if (!all.add(w)) { fprintf(stderr, "parts overlap\n"); return false; } } }
        tot += pt.m * dg;
        parts.push_back(pt);
    }
    if (tot != 28) { fprintf(stderr, "parts sum %d\n", tot); return false; }
    // sanity: T^p = id, Fix part fixed
    for (int k = 0; k < 28; k++) { u32 x = 1u << k; for (int t = 0; t < p; t++) x = T(x); if (x != 1u << k) { fprintf(stderr, "T^p\n"); return false; } }
    printf("A order %d; parts:", p); for (auto& pt : parts) printf(" g=%x F_%d^%d", pt.g, pt.q, pt.m); printf("\n");
    return true;
}

// F_q multiplication via polynomials mod g
static u32 fmul(u32 a, u32 b, u32 g, int d) { u32 r = 0; while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a >> d & 1) a ^= g; } return r; }

// all k-dim F_q-subspaces of F_q^m (reduced echelon, leading entry = lowest index, value 1), as F_2 bases of forms
template <class F> static void enum_sub(const Part& pt, int k, F&& cb) {
    int m = pt.m, q = pt.q, d = pt.d;
    std::vector<int> piv(k);
    std::vector<u32> forms(k * d);
    std::function<void(int, int)> pick = [&](int r, int from) {
        if (r == k) {
            // free positions per row: non-pivot positions > its pivot
            std::vector<std::vector<int>> fr(k);
            for (int i = 0; i < k; i++) for (int j = piv[i] + 1; j < m; j++) if (std::find(piv.begin(), piv.end(), j) == piv.end()) fr[i].push_back(j);
            std::vector<int> nf(k); int tot = 0; for (int i = 0; i < k; i++) { nf[i] = (int)fr[i].size(); tot += nf[i]; }
            std::vector<int> val(tot, 0);
            for (;;) {
                int o = 0;
                for (int i = 0; i < k; i++) {
                    u32 f = pt.tab[(size_t)piv[i] * q + 1];
                    for (int j = 0; j < nf[i]; j++) f ^= pt.tab[(size_t)fr[i][j] * q + val[o + j]];
                    o += nf[i];
                    u32 w = f; for (int t = 0; t < d; t++) { if (t) w = T(w); forms[i * d + t] = w; }
                }
                cb(forms.data(), k * d);
                int t = 0; while (t < tot) { if (++val[t] < q) break; val[t] = 0; t++; }
                if (t == tot) break;
            }
            return;
        }
        for (int p = from; p <= m - (k - r); p++) { piv[r] = p; pick(r + 1, p + 1); }
    };
    if (k == 0) { cb(forms.data(), 0); return; }
    pick(0, 0);
}
static double gbin(int n, int k, double q) { if (k < 0 || k > n) return 0; double r = 1; for (int i = 0; i < k; i++) r = r * (pow(q, n - i) - 1) / (pow(q, i + 1) - 1); return floor(r + 0.5); }

// per-a echelon state of a set of forms
struct St { u8 piv[256][8]; u8 rk[256]; int n; };
static void st_init(St& s) { memset(&s, 0, sizeof s); }
static inline bool ins(u8* P, u8 x) { while (x) { int b = 31 - __lzcnt(x); if (!P[b]) { P[b] = x; return true; } x ^= P[b]; } return false; }
static bool st_add(St& s, const u32* f, int n) { // false if some a gets deficiency > 1
    for (int i = 0; i < n; i++) {
        u64 m = rows_of(f[i]); u8 row[8]; for (int t = 0; t < 8; t++) row[t] = (u8)(m >> (8 * t));
        u8 r[256]; r[0] = 0; for (unsigned a = 1; a < 256; a++) r[a] = r[a & (a - 1)] ^ row[_tzcnt_u32(a)];
        s.n++;
        for (unsigned a = 1; a < 256; a++) { if (ins(s.piv[a], r[a])) s.rk[a]++; if (s.rk[a] + 1 < s.n) return false; }
    }
    return true;
}
static bool self_ok(const u32* f, int n) { St s; st_init(s); return st_add(s, f, n); }
// final: does prefix + extra (total 8) have rank 7 at every a? early exit
static inline bool final_ok(const St& s, const u64* xm, int nx) {
    u8 rows[8][8]; for (int i = 0; i < nx; i++) for (int t = 0; t < 8; t++) rows[i][t] = (u8)(xm[i] >> (8 * t));
    for (unsigned a = 1; a < 256; a++) {
        u8 P[8]; memcpy(P, s.piv[a], 8); int r = s.rk[a];
        for (int i = 0; i < nx; i++) { u8 v = 0; for (int t = 0; t < 8; t++) if (a >> t & 1) v ^= rows[i][t]; if (ins(P, v)) r++; }
        if (r != 7) return false;
    }
    return true;
}

struct Result { std::mutex mu; u64 apn = 0; std::map<std::string, u64> dist; std::vector<std::vector<u32>> big; u64 bugs = 0; };
static void record(Result& R, const u32* b) {
    Ech e; for (int i = 0; i < 8; i++) e.add(b[i]);
    bool ok = e.n == 8;
    for (int i = 0; i < 8 && ok; i++) if (e.reduce(T(b[i]))) ok = false;
    u64 m[8]; for (int i = 0; i < 8; i++) m[i] = rows_of(b[i]);
    for (unsigned a = 1; a < 256 && ok; a++) { u8 P[8] = {0}; int r = 0; for (int i = 0; i < 8; i++) if (ins(P, rv(m[i], a))) r++; if (r != 7) ok = false; }
    int cnt[9] = {0}, fx4 = 0;
    for (unsigned g = 1; g < 256; g++) { u32 c = 0; for (int i = 0; i < 8; i++) if (g >> i & 1) c ^= b[i]; int rk = rank8(rows_of(c)); cnt[rk]++; if (rk == 4 && T(c) == c) fx4++; }
    char key[96]; snprintf(key, sizeof key, "8:%d 6:%d 4:%d 2:%d  (rank-4 fixed %d, moved %d)", cnt[8], cnt[6], cnt[4], cnt[2], fx4, cnt[4] - fx4);
    std::lock_guard<std::mutex> g(R.mu);
    if (!ok) { R.bugs++; return; }
    R.apn++; R.dist[key]++;
    if (cnt[4] >= 5) R.big.push_back(std::vector<u32>(b, b + 8));
}

int main(int argc, char** argv) {
    init_tables();
    if (argc < 2) { fprintf(stderr, "usage\n"); return 2; }
    std::string n = argv[1]; Mat A = ident();
    auto comp = [&](int s, u32 g, int d) { for (int i = 0; i + 1 < d; i++) A.c[s + i] = (u8)(1 << (s + i + 1)); u8 last = 0; for (int t = 0; t < d; t++) if (g >> t & 1) last |= (u8)(1 << (s + t)); A.c[s + d - 1] = last; };
    if (n == "c5a") { comp(0, 0x1F, 4); comp(4, 0x1F, 4); }
    else if (n == "c5b") comp(0, 0x1F, 4);
    else if (n == "c7a") { comp(0, 0xB, 3); comp(3, 0xB, 3); }
    else if (n == "c7b") { comp(0, 0xB, 3); comp(3, 0xD, 3); }
    else if (n == "c7c") comp(0, 0xB, 3);
    else if (n.size() == 16) { for (int i = 0; i < 8; i++) { unsigned v; sscanf(n.c_str() + 2 * i, "%2x", &v); A.c[i] = (u8)v; } }
    else { fprintf(stderr, "A?\n"); return 2; }
    unsigned long long seed = argc > 2 ? strtoull(argv[2], 0, 10) : 0;
    if (seed) {
        u64 s = seed * 0x9E3779B97F4A7C15ull; Mat Pm, Pi;
        do { for (int i = 0; i < 8; i++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; Pm.c[i] = (u8)s; } } while (!inverse(Pm, Pi));
        A = mmul(Pm, mmul(A, Pi));
    }
    printf("A ="); for (int i = 0; i < 8; i++) printf(" %02x", A.c[i]); printf("\n");
    if (!decompose(A)) return 1;
    bool countonly = argc > 3 && !strcmp(argv[3], "count");
    int np = (int)parts.size();
    // splits
    std::vector<std::vector<int>> splits; std::vector<int> k(np, 0);
    std::function<void(int, int)> sp = [&](int b, int rem) { if (b == np) { if (!rem) splits.push_back(k); return; } for (int x = 0; x <= parts[b].m && x * parts[b].d <= rem; x++) { k[b] = x; sp(b + 1, rem - x * parts[b].d); } };
    sp(0, 8);
    double total_spaces = 0;
    Result R;
    for (auto& ks : splits) {
        printf("split"); double cnt = 1; std::vector<double> cb(np);
        for (int b = 0; b < np; b++) { cb[b] = gbin(parts[b].m, ks[b], parts[b].q); cnt *= cb[b]; printf(" %d", ks[b]); }
        printf(": %.0f invariant subspaces\n", cnt); total_spaces += cnt;
        if (countonly) {
            for (int b = 0; b < np; b++) { if (cb[b] > 2e7) { printf("  part %d: %.0f (not enumerated)\n", b, cb[b]); continue; } u64 c = 0; enum_sub(parts[b], ks[b], [&](const u32*, int) { c++; }); printf("  part %d: enumerated %llu, Gaussian binomial %.0f %s\n", b, c, cb[b], (double)c == cb[b] ? "ok" : "MISMATCH"); }
            continue;
        }
        fflush(stdout);
        // stream part = the largest; others listed (filtered by their own kernel test)
        int sb = -1; for (int b = 0; b < np; b++) if (ks[b] && (sb < 0 || cb[b] > cb[sb])) sb = b;
        std::vector<std::vector<u32>> lists(np); std::vector<int> ldim(np, 0); std::vector<u64> raw(np, 0);
        for (int b = 0; b < np; b++) {
            if (!ks[b] || b == sb) continue;
            ldim[b] = ks[b] * parts[b].d;
            enum_sub(parts[b], ks[b], [&](const u32* f, int nn) { raw[b]++; if (self_ok(f, nn)) lists[b].insert(lists[b].end(), f, f + nn); });
            printf("  part %d: %llu subspaces, %zu pass their own kernel test\n", b, raw[b], lists[b].size() / ldim[b]);
        }
        // prefixes: product of the listed parts
        std::vector<int> lb; for (int b = 0; b < np; b++) if (ks[b] && b != sb) lb.push_back(b);
        std::vector<std::vector<u32>> prefixes;
        std::vector<u32> cur;
        std::function<void(int)> prod = [&](int i) {
            if (i == (int)lb.size()) { prefixes.push_back(cur); return; }
            int b = lb[i]; size_t cntb = lists[b].size() / ldim[b];
            for (size_t j = 0; j < cntb; j++) { size_t o = cur.size(); cur.insert(cur.end(), lists[b].begin() + j * ldim[b], lists[b].begin() + (j + 1) * ldim[b]); prod(i + 1); cur.resize(o); }
        };
        prod(0);
        int sdim = ks[sb] * parts[sb].d;
        std::vector<u32> slist; u64 sraw = 0; bool listed = cb[sb] <= 5e6;
        if (listed) { enum_sub(parts[sb], ks[sb], [&](const u32* f, int nn) { sraw++; if (self_ok(f, nn)) slist.insert(slist.end(), f, f + nn); }); printf("  stream part %d: %llu subspaces, %zu pass their own kernel test\n", sb, sraw, slist.size() / sdim); }
        printf("  prefixes %zu\n", prefixes.size()); fflush(stdout);
        std::atomic<u64> tested(0), pref_ok(0);
        u64 before = R.apn;
        #pragma omp parallel for schedule(dynamic, 1)
        for (long long pi = 0; pi < (long long)prefixes.size(); pi++) {
            St* s = new St; st_init(*s);
            const std::vector<u32>& pf = prefixes[pi];
            if (!st_add(*s, pf.data(), (int)pf.size())) { delete s; continue; }
            pref_ok++;
            u64 loc = 0;
            auto test = [&](const u32* f, int nn) {
                loc++;
                u64 xm[8]; for (int i = 0; i < nn; i++) xm[i] = rows_of(f[i]);
                if (final_ok(*s, xm, nn)) { u32 b8[8]; int o = 0; for (u32 x : pf) b8[o++] = x; for (int i = 0; i < nn; i++) b8[o++] = f[i]; record(R, b8); }
            };
            if (listed) { size_t c = slist.size() / sdim; for (size_t j = 0; j < c; j++) test(slist.data() + j * sdim, sdim); }
            else enum_sub(parts[sb], ks[sb], test);
            tested += loc;
            delete s;
        }
        printf("  prefixes passing %llu, final tests %llu, APN spaces in this split %llu\n", (u64)pref_ok, (u64)tested, R.apn - before); fflush(stdout);
    }
    printf("total invariant 8-dim subspaces %.0f\n", total_spaces);
    if (countonly) return 0;
    printf("APN invariant spaces: %llu (check failures %llu)\n", R.apn, R.bugs);
    for (auto& kv : R.dist) printf("  ranks %s : %llu\n", kv.first.c_str(), kv.second);
    printf("with >= 5 rank-4 components: %zu\n", R.big.size());
    for (auto& b : R.big) { printf("FOUND"); for (u32 x : b) printf(" %07x", x); printf("\n"); }
    return 0;
}
