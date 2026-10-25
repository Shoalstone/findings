// apn8: quadratic APN functions on F_2^8, searched as 8-dim subspaces of alternating forms.
// Target A: a quadratic (8,8)-APN function with >= 5 rank-4 components (open, Beneteau-Goluboff-Kolsch-Vaghasiya 2026,
// Open Problem VI.1; see notes/work.md). Rowen (Claude Opus 5.5), 2026-10-24.
//
// A homogeneous quadratic f on F_2^8 <-> its alternating form B(x,y) = f(x+y)+f(x)+f(y). Stored two ways:
//   mat: uint64, byte i = row i (bit j = B(e_i,e_j));   coef: uint32, 28 bits c_ij (i<j), index K[i][j].
// F = (f_1..f_m) is APN-so-far iff for every a != 0 the functionals r_s(a) = B_s(a,.) span a space of rank m' with
// kernel dim 8 - rank; APN (m=8) iff every kernel is {0,a} (rank 7).
//
// Search: S6 seeded with rank-4 forms (kernel dim <= 3 everywhere), exhaustive 7th coordinate over the 22-dim
// complement (Gray code, packed nibble test), then the 8th by switching (linear system: B_v(a,x_a) = 1 on each
// 2-dim kernel), then the rank distribution of the 255 components.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <atomic>
#include <omp.h>
#include <cmath>

typedef uint64_t u64; typedef uint32_t u32; typedef uint8_t u8;

static int K[8][8]; static int KI[28], KJ[28];
static void init_index() { int k = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { K[i][j] = K[j][i] = k; KI[k] = i; KJ[k] = j; k++; } }

static u64 coef_to_mat(u32 c) {
    u64 m = 0;
    for (int k = 0; k < 28; k++) if (c >> k & 1) { int i = KI[k], j = KJ[k]; m |= (u64)1 << (8 * i + j); m |= (u64)1 << (8 * j + i); }
    return m;
}
static u32 mat_to_coef(u64 m) {
    u32 c = 0;
    for (int k = 0; k < 28; k++) if (m >> (8 * KI[k] + KJ[k]) & 1) c |= 1u << k;
    return c;
}
static inline u8 row(u64 m, int i) { return (u8)(m >> (8 * i)); }
// r(a) = B(a, .) as an 8-bit mask
static inline u8 rvec(u64 m, unsigned a) { u8 r = 0; for (int i = 0; i < 8; i++) if (a >> i & 1) r ^= row(m, i); return r; }

static int rank_mat(u64 m) {
    u8 r[8]; for (int i = 0; i < 8; i++) r[i] = row(m, i);
    int rk = 0;
    for (int bit = 7; bit >= 0; bit--) {
        int p = -1; for (int i = rk; i < 8; i++) if (r[i] >> bit & 1) { p = i; break; }
        if (p < 0) continue;
        std::swap(r[p], r[rk]);
        for (int i = 0; i < 8; i++) if (i != rk && (r[i] >> bit & 1)) r[i] ^= r[rk];
        rk++;
    }
    return rk;
}

// rank of a set of 8-bit vectors
static int rank8(const u8* v, int n) {
    u8 piv[8] = {0}; int rk = 0;
    for (int t = 0; t < n; t++) {
        u8 x = v[t];
        for (int b = 7; b >= 0 && x; b--) if (x >> b & 1) { if (piv[b]) x ^= piv[b]; else { piv[b] = x; rk++; x = 0; } }
    }
    return rk;
}

// echelon basis over 28-bit coefficient space
struct Ech28 {
    u32 piv[28]; int n = 0;
    Ech28() { memset(piv, 0, sizeof piv); }
    u32 reduce(u32 x) const { for (int b = 27; b >= 0 && x; b--) if ((x >> b & 1) && piv[b]) x ^= piv[b]; return x; }
    bool add(u32 x) { x = reduce(x); if (!x) return false; int b = 31 - __lzcnt(x); piv[b] = x; n++; return true; }
};

// homogeneous quadratic f(x) from its form: f(x) = sum_{i<j} B_ij x_i x_j
static inline int fval(u64 m, unsigned x) {
    int s = 0;
    for (int i = 0; i < 8; i++) if (x >> i & 1) s ^= __popcnt((unsigned)(row(m, i) & x & (0xFEu << i) & 0xFF)) & 1;
    return s;
}

// ---------- rng ----------
struct Rng {
    u64 s[4];
    static u64 sm(u64& x) { u64 z = (x += 0x9E3779B97F4A7C15ull); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
    explicit Rng(u64 seed) { for (int i = 0; i < 4; i++) s[i] = sm(seed); }
    static u64 rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }
    u64 next() { u64 r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17; s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
    u32 below(u32 n) { return (u32)((next() >> 32) * n >> 32); }
};

static u32 wedge(u8 l, u8 m) { // coef of l^m: c_ij = l_i m_j + l_j m_i
    u32 c = 0;
    for (int k = 0; k < 28; k++) { int i = KI[k], j = KJ[k]; if (((l >> i & 1) & (m >> j & 1)) ^ ((l >> j & 1) & (m >> i & 1))) c |= 1u << k; }
    return c;
}
static u32 random_rank4(Rng& R) {
    u8 v[4];
    for (;;) { for (int i = 0; i < 4; i++) v[i] = (u8)R.next(); if (rank8(v, 4) == 4) break; }
    return wedge(v[0], v[1]) ^ wedge(v[2], v[3]);
}

// kernel dims of span(forms) at every a: returns max kernel dim
static int max_kernel(const u64* mats, int m, int* kd = nullptr) {
    int mx = 0; u8 r[8];
    for (unsigned a = 1; a < 256; a++) {
        for (int s = 0; s < m; s++) r[s] = rvec(mats[s], a);
        int k = 8 - rank8(r, m);
        if (kd) kd[a] = k;
        if (k > mx) mx = k;
    }
    return mx;
}

// full component rank distribution of an 8-form basis (255 nonzero combos): counts[rank/2]
static void rank_dist(const u64* mats, int m, int* cnt) {
    for (int i = 0; i < 5; i++) cnt[i] = 0;
    u64 cur = 0;
    for (unsigned g = 1; g < (1u << m); g++) { cur ^= mats[_tzcnt_u32(g)]; cnt[rank_mat(cur) / 2]++; }
}

// brute-force check on the look-up table: differential uniformity of F (m coordinates)
static int lut_delta(const u64* mats, int m) {
    unsigned F[256];
    for (unsigned x = 0; x < 256; x++) { unsigned y = 0; for (int s = 0; s < m; s++) y |= (unsigned)fval(mats[s], x) << s; F[x] = y; }
    int best = 0; static thread_local int cnt[256];
    for (unsigned a = 1; a < 256; a++) {
        memset(cnt, 0, sizeof(int) * (1u << m));
        for (unsigned x = 0; x < 256; x++) cnt[F[x] ^ F[x ^ a]]++;
        for (unsigned b = 0; b < (1u << m); b++) if (cnt[b] > best) best = cnt[b];
    }
    return best;
}

// ---------- the switching step: given S7 (7 forms) with all kernels <= 2, solve for 8th coordinates ----------
// returns solution space: v0 and nullspace basis (mod nothing; caller reduces mod S7). false if inconsistent.
static bool switching(const u64* mats7, u32& v0, std::vector<u32>& null) {
    // planes: kernel {0,a,x,a+x}; take each once (a minimal)
    unsigned F[256];
    for (unsigned x = 0; x < 256; x++) { unsigned y = 0; for (int s = 0; s < 7; s++) y |= (unsigned)fval(mats7[s], x) << s; F[x] = y; }
    // rows: 28-bit coef row + rhs 1, Gaussian elimination online
    u32 piv[28]; u8 rhs[28]; bool have[28]; memset(have, 0, sizeof have);
    for (unsigned a = 1; a < 256; a++) {
        unsigned xa = 0; int z = 0;
        for (unsigned x = 1; x < 256; x++) if (x != a && (F[x ^ a] ^ F[x] ^ F[a]) == 0) { z++; if (!xa) xa = x; }
        if (z == 0) continue;            // kernel {0,a}: no condition
        if (z > 2) return false;         // kernel dim > 2 (shouldn't happen)
        unsigned c = a ^ xa; if (xa < a || c < a) continue; // plane seen from its minimal element
        u32 rw = wedge((u8)a, (u8)xa); u8 r = 1;
        for (int b = 27; b >= 0 && rw; b--) if (rw >> b & 1) { if (have[b]) { rw ^= piv[b]; r ^= rhs[b]; } else { piv[b] = rw; rhs[b] = r; have[b] = true; rw = 0; r = 2; } }
        if (r == 1 && !rw) return false; // 0 = 1
    }
    // back-substitute: free vars = non-pivot bits
    // reduce to RREF
    for (int b = 0; b < 28; b++) if (have[b]) for (int c = b + 1; c < 28; c++) if (have[c] && (piv[c] >> b & 1)) { piv[c] ^= piv[b]; rhs[c] ^= rhs[b]; }
    v0 = 0; for (int b = 0; b < 28; b++) if (have[b] && rhs[b]) v0 |= 1u << b;
    null.clear();
    for (int f = 0; f < 28; f++) if (!have[f]) { u32 v = 1u << f; for (int b = 0; b < 28; b++) if (have[b] && (piv[b] >> f & 1)) v |= 1u << b; null.push_back(v); }
    return true;
}

struct Found { int cnt[5]; u32 basis[8]; };

// ---------- search from one S6 ----------
struct Stats { u64 s6 = 0, s6ok = 0, w_pass = 0, solved = 0, apn = 0; u64 hist[20] = {0}; };

static void canon8(u32* b) { // reduced echelon of 8 coef vectors, sorted: a canonical name for the subspace
    Ech28 e; for (int i = 0; i < 8; i++) e.add(b[i]);
    for (int x = 27; x >= 0; x--) if (e.piv[x]) for (int y = x + 1; y < 28; y++) if (e.piv[y] && (e.piv[y] >> x & 1)) e.piv[y] ^= e.piv[x];
    int n = 0; for (int x = 27; x >= 0; x--) if (e.piv[x]) b[n++] = e.piv[x];
}

static void search_from_s6(const u32* s6c, Stats& st, std::vector<Found>& out, int report_r4, u64 max_w_pass) {
    u64 m6[6]; for (int i = 0; i < 6; i++) m6[i] = coef_to_mat(s6c[i]);
    Ech28 e6; for (int i = 0; i < 6; i++) e6.add(s6c[i]);
    // complement basis: unit vectors not reducible
    u32 comp[22]; int nc = 0; { Ech28 e = e6; for (int b = 0; b < 28 && nc < 22; b++) if (e.add(1u << b)) comp[nc++] = 1u << b; }
    // bad a: kernel dim 3 in S6; per bad a a residual map (3 bits) whose kernel on a^perp is span_a
    int nbad = 0; unsigned bada[256]; u8 pivs[256][8];
    u8 nonpiv[256][3];
    for (unsigned a = 1; a < 256; a++) {
        u8 r[6]; for (int s = 0; s < 6; s++) r[s] = rvec(m6[s], a);
        u8 piv[8] = {0}; int rk = 0;
        for (int t = 0; t < 6; t++) { u8 x = r[t]; for (int b = 7; b >= 0 && x; b--) if (x >> b & 1) { if (piv[b]) x ^= piv[b]; else { piv[b] = x; rk++; x = 0; } } }
        if (rk < 5) return;              // kernel dim >= 4: not extendable (caller checked)
        if (rk == 5) { memcpy(pivs[nbad], piv, 8); int q = 0; for (int b = 0; b < 8; b++) if (!piv[b] && q < 3) nonpiv[nbad][q++] = (u8)b; bada[nbad++] = a; }
    }
    // packed: slot t = nibble t; value = residual bits of r_w(a)
    int nw = (nbad + 15) / 16; if (nw == 0) nw = 1;
    std::vector<u64> delta((size_t)nc * nw, 0);
    for (int k = 0; k < nc; k++) {
        u64 mk = coef_to_mat(comp[k]);
        for (int t = 0; t < nbad; t++) {
            u8 x = rvec(mk, bada[t]);
            for (int b = 7; b >= 0; b--) if ((x >> b & 1) && pivs[t][b]) x ^= pivs[t][b];
            u64 nib = 0; for (int q = 0; q < 3; q++) if (x >> nonpiv[t][q] & 1) nib |= 1ull << q;
            delta[(size_t)k * nw + t / 16] |= nib << (4 * (t % 16));
        }
    }
    u64 need[16] = {0}; for (int t = 0; t < nbad; t++) need[t / 16] |= 1ull << (4 * (t % 16));
    u64 v[16] = {0};
    const u64 ONES = 0x1111111111111111ull;
    std::unordered_set<std::string> seen;
    u32 w = 0;
    u64 passes = 0;
    for (u32 g = 1; g < (1u << nc); g++) {
        int bit = (int)_tzcnt_u32(g);
        w ^= comp[bit];
        const u64* d = &delta[(size_t)bit * nw];
        bool ok = true;
        for (int i = 0; i < nw; i++) { v[i] ^= d[i]; }
        for (int i = 0; i < nw; i++) { u64 x = v[i]; if (((x | x >> 1 | x >> 2) & ONES & need[i]) != need[i]) { ok = false; break; } }
        if (!ok) continue;
        st.w_pass++; passes++;
        if (passes > max_w_pass) break;
        u64 m7[7]; memcpy(m7, m6, sizeof m6); m7[6] = coef_to_mat(w);
        u32 v0; std::vector<u32> null;
        if (!switching(m7, v0, null)) continue;
        st.solved++;
        // reduce null and v0 modulo S7
        Ech28 e7 = e6; e7.add(w);
        std::vector<u32> nb; { Ech28 t = e7; for (u32 x : null) if (t.add(x)) nb.push_back(t.reduce(x) ? x : x); }
        // nb are directions independent of S7 (each new choice of 8th coordinate)
        int dn = (int)nb.size(); if (dn > 10) dn = 10;
        for (u32 sel = 0; sel < (1u << dn); sel++) {
            u32 v8 = v0; for (int i = 0; i < dn; i++) if (sel >> i & 1) v8 ^= nb[i];
            if (!e7.reduce(v8)) continue;
            u64 m8[8]; memcpy(m8, m7, sizeof m7); m8[7] = coef_to_mat(v8);
            int cnt[5]; rank_dist(m8, 8, cnt);
            if (cnt[0] != 0) continue;
            // APN check by kernels (fast); LUT check for reported ones
            if (max_kernel(m8, 8) != 1) continue;
            u32 b8[8]; for (int i = 0; i < 6; i++) b8[i] = s6c[i]; b8[6] = w; b8[7] = v8; canon8(b8);
            std::string key((const char*)b8, sizeof b8);
            if (!seen.insert(key).second) continue;
            st.apn++; st.hist[cnt[2] < 19 ? cnt[2] : 19]++;
            if (report_r4 < 0 || cnt[2] >= report_r4 || cnt[1] > 0) {
                Found f; memcpy(f.cnt, cnt, sizeof cnt); memcpy(f.basis, b8, sizeof b8); out.push_back(f);
            }
        }
    }
}

// random S6 seeded with k rank-4 forms; returns false if kernel condition fails
static bool make_s6(Rng& R, int k, u32* s6c) {
    Ech28 e; u64 m[6]; int n = 0;
    int tries = 0;
    while (n < 6 && tries < 200) {
        tries++;
        u32 c = n < k ? random_rank4(R) : (u32)(R.next() & 0x0FFFFFFF);
        if (!e.add(c)) continue;
        m[n] = coef_to_mat(c); s6c[n] = c; n++;
        if (max_kernel(m, n) > 9 - n) { // Lemma 3.12: kernel dim <= n - m + 1 for an APN-extendable (8,m)
            // undo: rebuild echelon without it
            n--; Ech28 e2; for (int i = 0; i < n; i++) e2.add(s6c[i]); e = e2;
        }
    }
    if (n < 6) return false;
    // no rank-2 components wanted beyond what APN allows; just require none here
    int cnt[5]; rank_dist(m, 6, cnt); if (cnt[1] > 0) return false;
    return true;
}

static void selftest() {
    // x^3 over GF(2^8) mod x^8+x^4+x^3+x+1
    auto mul = [](unsigned a, unsigned b) { unsigned r = 0; while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a & 0x100) a ^= 0x11B; } return r; };
    unsigned F[256]; for (unsigned x = 0; x < 256; x++) F[x] = mul(mul(x, x), x);
    u64 mats[8];
    for (int s = 0; s < 8; s++) {
        u64 m = 0;
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (i != j) {
            unsigned b = ((F[(1u << i) ^ (1u << j)] ^ F[1u << i] ^ F[1u << j] ^ F[0]) >> s) & 1;
            if (b) m |= (u64)1 << (8 * i + j);
        }
        mats[s] = m;
    }
    int cnt[5]; rank_dist(mats, 8, cnt);
    printf("x^3: ranks 8:%d 6:%d 4:%d 2:%d  maxker %d  lut delta %d\n", cnt[4], cnt[3], cnt[2], cnt[1], max_kernel(mats, 8), lut_delta(mats, 8));
    // fval matches LUT up to affine part: check D_a of quadratic part equals F's
    int bad = 0;
    for (unsigned a = 1; a < 256; a++) for (unsigned x = 0; x < 256; x++) {
        unsigned y = 0; for (int s = 0; s < 8; s++) y |= (unsigned)(fval(mats[s], x ^ a) ^ fval(mats[s], x) ^ fval(mats[s], a)) << s;
        if (y != (F[x ^ a] ^ F[x] ^ F[a] ^ F[0])) bad++;
    }
    printf("derivative mismatch %d (want 0)\n", bad);
    for (int t = 0; t < 1000; t++) { u32 c = (u32)(t * 2654435761u) & 0x0FFFFFFF; if (mat_to_coef(coef_to_mat(c)) != c) { printf("coef roundtrip FAIL\n"); return; } }
    Rng R(1); int r4 = 0; for (int t = 0; t < 1000; t++) if (rank_mat(coef_to_mat(random_rank4(R))) == 4) r4++;
    printf("random_rank4 rank-4: %d/1000\n", r4);
    // switching from x^3 minus one coordinate should give back an APN
    u32 v0; std::vector<u32> nl; bool ok = switching(mats, v0, nl);
    printf("switching from x^3's first 7 coords: %s, null dim %zu\n", ok ? "consistent" : "INCONSISTENT", nl.size());
}

static void lut_to_mats(const unsigned* F, u64* mats) {
    for (int s = 0; s < 8; s++) {
        u64 m = 0;
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (i != j) {
            unsigned b = ((F[(1u << i) ^ (1u << j)] ^ F[1u << i] ^ F[1u << j] ^ F[0]) >> s) & 1;
            if (b) m |= (u64)1 << (8 * i + j);
        }
        mats[s] = m;
    }
}
static void x3_basis(u32* b) {
    auto mul = [](unsigned a, unsigned c) { unsigned r = 0; while (c) { if (c & 1) r ^= a; c >>= 1; a <<= 1; if (a & 0x100) a ^= 0x11B; } return r; };
    unsigned F[256]; for (unsigned x = 0; x < 256; x++) F[x] = mul(mul(x, x), x);
    u64 m[8]; lut_to_mats(F, m); for (int i = 0; i < 8; i++) b[i] = mat_to_coef(m[i]);
    canon8(b);
}

// a random 6-dim subspace of span(b8); with prob bias, built to contain the rank-4 components first
static void pick_s6(Rng& R, const u32* b8, u32* s6, double bias) {
    std::vector<u32> r4;
    if ((R.next() >> 11) * (1.0 / 9007199254740992.0) < bias) {
        u32 cur = 0;
        for (unsigned g = 1; g < 256; g++) { cur ^= b8[_tzcnt_u32(g)]; if (rank_mat(coef_to_mat(cur)) == 4) r4.push_back(cur); }
        for (size_t i = r4.size(); i > 1; i--) std::swap(r4[i - 1], r4[R.below((u32)i)]);
    }
    Ech28 e; int n = 0;
    for (u32 c : r4) { if (n == 6) break; if (e.add(c)) s6[n++] = c; }
    while (n < 6) { u32 c = 0; u32 sel = R.below(255) + 1; for (int i = 0; i < 8; i++) if (sel >> i & 1) c ^= b8[i]; if (e.add(c)) s6[n++] = c; }
}

static void print_found(FILE* f, const Found& x);

// walk: from x^3, repeatedly drop to a 6-subspace and re-extend, moving with weight base^N4
static void walk(u64 seed, double secs, const char* outp, int th, int rep, double base, double bias) {
    omp_set_num_threads(th);
    auto t0 = std::chrono::steady_clock::now();
    std::atomic<bool> stop(false); Stats tot; int best_all = -1; double next_report = 600;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num(); Rng R(seed * 1000003ull + tid * 7919 + 1); Stats st;
        u32 cur[8]; x3_basis(cur); int curN4 = 0; u64 steps = 0;
        std::vector<Found> out;
        while (!stop.load()) {
            u32 s6[6]; pick_s6(R, cur, s6, bias);
            out.clear(); search_from_s6(s6, st, out, -1, 4000000); steps++;
            if (!out.empty()) {
                double tw = 0; std::vector<double> w(out.size());
                for (size_t i = 0; i < out.size(); i++) { w[i] = pow(base, out[i].cnt[2]); tw += w[i]; }
                double u = (R.next() >> 11) * (1.0 / 9007199254740992.0) * tw; size_t pick = 0;
                for (; pick + 1 < out.size(); pick++) { u -= w[pick]; if (u <= 0) break; }
                memcpy(cur, out[pick].basis, sizeof cur); curN4 = out[pick].cnt[2];
                std::vector<const Found*> keep;
                for (auto& x : out) if (x.cnt[2] >= rep || x.cnt[1] > 0) keep.push_back(&x);
                if (!keep.empty()) {
                    #pragma omp critical
                    { FILE* f = fopen(outp, "a"); for (auto* x : keep) print_found(f, *x); fclose(f); }
                }
                #pragma omp critical
                {
                    for (auto& x : out) {
                        if (x.cnt[2] > best_all) { best_all = x.cnt[2]; printf("[t%d step %llu] new best N4=%d: ", tid, steps, best_all); print_found(stdout, x); fflush(stdout); }
                        if (x.cnt[2] >= 5) { std::string p5 = std::string(outp) + ".FOUND5"; FILE* f = fopen(p5.c_str(), "a"); print_found(f, x); fclose(f); }
                    }
                }
            }
            double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (el > secs) stop = true;
            if (tid == 0 && el > next_report) {
                next_report += 600;
                printf("[%.0f s] thread0 steps %llu apn %llu by #rank4:", el, steps, st.apn);
                for (int i = 0; i < 20; i++) if (st.hist[i]) printf(" %d:%llu", i, st.hist[i]);
                printf("\n"); fflush(stdout);
            }
        }
        #pragma omp critical
        { tot.s6 += steps; tot.w_pass += st.w_pass; tot.solved += st.solved; tot.apn += st.apn; for (int i = 0; i < 20; i++) tot.hist[i] += st.hist[i]; }
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("walk done %.1fs steps %llu wpass %llu solved %llu apn %llu | by #rank4:", el, tot.s6, tot.w_pass, tot.solved, tot.apn);
    for (int i = 0; i < 20; i++) if (tot.hist[i]) printf(" %d:%llu", i, tot.hist[i]);
    printf("\n");
}

// walk2: like walk, but threads don't revisit their own states, move only to unseen neighbours when they can, and
// when stuck jump to a random member of the shared pool of functions with >= poolmin rank-4 components.
// Writes each globally new function with >= rep rank-4 components once.
static bool load_seeds(const char* path, std::vector<Found>& seeds) {
    FILE* f = fopen(path, "r"); if (!f) return false;
    static char line[1 << 16];
    while (fgets(line, sizeof line, f)) {
        Found x; u32 b[8];
        char* p = strstr(line, "basis");
        if (p) { if (sscanf(p + 5, "%x %x %x %x %x %x %x %x", &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7]) != 8) continue; }
        else if (line[0] == '{') {
            unsigned F[256]; char* q = line + 1; int n = 0;
            while (n < 256) { char* e; long v = strtol(q, &e, 10); if (e == q) break; F[n++] = (unsigned)v; q = e; while (*q == ',' || *q == ' ') q++; }
            if (n != 256) continue;
            u64 m[8]; lut_to_mats(F, m); for (int i = 0; i < 8; i++) b[i] = mat_to_coef(m[i]);
        } else continue;
        u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(b[i]);
        if (max_kernel(m, 8) != 1) continue;
        canon8(b); memcpy(x.basis, b, sizeof b); rank_dist(m, 8, x.cnt);
        seeds.push_back(x);
    }
    fclose(f); return true;
}

static void walk2(u64 seed, double secs, const char* outp, int th, int rep, double base, double bias, int poolmin, const char* seedfile = nullptr) {
    std::vector<Found> seeds;
    if (seedfile) { if (!load_seeds(seedfile, seeds) || seeds.empty()) { fprintf(stderr, "no seeds from %s\n", seedfile); return; } printf("%zu seeds\n", seeds.size()); fflush(stdout); }
    omp_set_num_threads(th);
    auto t0 = std::chrono::steady_clock::now();
    std::atomic<bool> stop(false);
    std::unordered_set<std::string> gseen; std::vector<Found> pool[20]; u64 gcount[20] = {0}; u64 gsteps = 0;
    double next_report = 600;
    auto key_of = [](const Found& f) { return std::string((const char*)f.basis, sizeof f.basis); };
    for (auto& x : seeds) if (gseen.insert(key_of(x)).second && x.cnt[2] >= poolmin) pool[x.cnt[2] < 19 ? x.cnt[2] : 19].push_back(x);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num(); Rng R(seed * 1000003ull + tid * 7919 + 17); Stats st;
        Found cur; x3_basis(cur.basis); rank_dist(nullptr, 0, cur.cnt);
        if (!seeds.empty()) cur = seeds[tid % seeds.size()];
        std::unordered_set<std::string> visited; visited.insert(key_of(cur));
        int stall = 0; u64 steps = 0;
        std::vector<Found> out;
        while (!stop.load()) {
            u32 s6[6]; pick_s6(R, cur.basis, s6, bias);
            out.clear(); search_from_s6(s6, st, out, -1, 4000000); steps++;
            std::vector<const Found*> cand; double tw = 0; std::vector<double> w;
            for (auto& x : out) if (!visited.count(key_of(x))) { cand.push_back(&x); w.push_back(pow(base, x.cnt[2])); tw += w.back(); }
            bool gnew = false;
            #pragma omp critical
            {
                gsteps++;
                FILE* f = nullptr;
                for (auto& x : out) {
                    if (x.cnt[2] < poolmin && x.cnt[2] < rep) continue;
                    if (!gseen.insert(key_of(x)).second) continue;
                    int n4 = x.cnt[2] < 19 ? x.cnt[2] : 19; gcount[n4]++;
                    if (x.cnt[2] >= poolmin) pool[n4].push_back(x);
                    if (x.cnt[2] >= rep) { if (!f) f = fopen(outp, "a"); print_found(f, x); gnew = true; }
                    if (x.cnt[2] >= 5) { std::string p5 = std::string(outp) + ".FOUND5"; FILE* g = fopen(p5.c_str(), "a"); print_found(g, x); fclose(g);
                                         printf("[t%d] N4=%d: ", tid, x.cnt[2]); print_found(stdout, x); fflush(stdout); }
                }
                if (f) fclose(f);
            }
            if (gnew) stall = 0; else stall++;
            if (!cand.empty() && stall < 60) {
                double u = (R.next() >> 11) * (1.0 / 9007199254740992.0) * tw; size_t pick = 0;
                for (; pick + 1 < cand.size(); pick++) { u -= w[pick]; if (u <= 0) break; }
                cur = *cand[pick]; visited.insert(key_of(cur));
            } else {
                // jump into the pool (highest classes more often)
                bool jumped = false;
                #pragma omp critical
                {
                    for (int n = 19; n >= poolmin && !jumped; n--) if (!pool[n].empty() && (R.below(10) < 7 || n == poolmin)) { cur = pool[n][R.below((u32)pool[n].size())]; jumped = true; }
                }
                if (!jumped && !out.empty()) cur = out[R.below((u32)out.size())];
                visited.insert(key_of(cur)); stall = 0;
            }
            if (visited.size() > 2000000) visited.clear();
            double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (el > secs) stop = true;
            if (tid == 0 && el > next_report) {
                next_report += 600;
                #pragma omp critical
                { printf("[%.0f s] steps %llu distinct by #rank4:", el, gsteps); for (int i = 0; i < 20; i++) if (gcount[i]) printf(" %d:%llu", i, gcount[i]); printf("\n"); fflush(stdout); }
            }
        }
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("walk2 done %.1fs steps %llu distinct by #rank4:", el, gsteps);
    for (int i = 0; i < 20; i++) if (gcount[i]) printf(" %d:%llu", i, gcount[i]);
    printf("\n");
}

static void print_found(FILE* f, const Found& x) {
    fprintf(f, "ranks 8:%d 6:%d 4:%d 2:%d basis", x.cnt[4], x.cnt[3], x.cnt[2], x.cnt[1]);
    for (int i = 0; i < 8; i++) fprintf(f, " %07x", x.basis[i]);
    fprintf(f, "\n");
}

int main(int argc, char** argv) {
    init_index();
    if (argc < 2) { fprintf(stderr, "usage: apn8 selftest | search <seed> <k_rank4> <seconds> <out> [threads] [report_r4]\n"); return 2; }
    std::string cmd = argv[1];
    if (cmd == "selftest") { selftest(); return 0; }
    if (cmd == "walk2") { // walk2 <seed> <seconds> <out> [threads] [report_r4] [base] [bias] [poolmin]
        walk2(strtoull(argv[2], 0, 10), atof(argv[3]), argv[4], argc > 5 ? atoi(argv[5]) : omp_get_max_threads(),
              argc > 6 ? atoi(argv[6]) : 4, argc > 7 ? atof(argv[7]) : 4.0, argc > 8 ? atof(argv[8]) : 0.7, argc > 9 ? atoi(argv[9]) : 3, argc > 10 ? argv[10] : nullptr);
        return 0;
    }
    if (cmd == "walk") { // walk <seed> <seconds> <out> [threads] [report_r4] [base] [bias]
        walk(strtoull(argv[2], 0, 10), atof(argv[3]), argv[4], argc > 5 ? atoi(argv[5]) : omp_get_max_threads(),
             argc > 6 ? atoi(argv[6]) : 4, argc > 7 ? atof(argv[7]) : 4.0, argc > 8 ? atof(argv[8]) : 0.7);
        return 0;
    }
    if (cmd == "search") {
        u64 seed = strtoull(argv[2], 0, 10); int k = atoi(argv[3]); double secs = atof(argv[4]); const char* outp = argv[5];
        int th = argc > 6 ? atoi(argv[6]) : omp_get_max_threads(); int rep = argc > 7 ? atoi(argv[7]) : 4;
        omp_set_num_threads(th);
        auto t0 = std::chrono::steady_clock::now();
        Stats tot; std::atomic<bool> stop(false);
        #pragma omp parallel
        {
            int tid = omp_get_thread_num(); Rng R(seed * 1000003ull + tid); Stats st; std::vector<Found> out;
            while (!stop.load()) {
                u32 s6[6]; st.s6++;
                if (!make_s6(R, k, s6)) continue;
                st.s6ok++;
                search_from_s6(s6, st, out, rep, 2000000);
                if (!out.empty()) {
                    #pragma omp critical
                    { FILE* f = fopen(outp, "a"); for (auto& x : out) { print_found(f, x); print_found(stdout, x); } fclose(f); fflush(stdout); }
                    out.clear();
                }
                double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                if (el > secs) stop = true;
            }
            #pragma omp critical
            { tot.s6 += st.s6; tot.s6ok += st.s6ok; tot.w_pass += st.w_pass; tot.solved += st.solved; tot.apn += st.apn; for (int i = 0; i < 20; i++) tot.hist[i] += st.hist[i]; }
        }
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        printf("done %.1fs s6 %llu ok %llu wpass %llu solved %llu apn %llu | by #rank4:", el, tot.s6, tot.s6ok, tot.w_pass, tot.solved, tot.apn);
        for (int i = 0; i < 20; i++) if (tot.hist[i]) printf(" %d:%llu", i, tot.hist[i]);
        printf("\n");
        return 0;
    }
    return 2;
}
