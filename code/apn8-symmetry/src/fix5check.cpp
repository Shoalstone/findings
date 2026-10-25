// fix5check: independent check of the order-5 "class 1, A fixes every rank-4 component" case.
// V = GF(256) (poly 0x11D), A = multiplication by zeta = 2^51 (order 5).
// Alt(V) = Fix (+) N. S8 = Q (+) L with Q a 4-dim subspace of Fix holding >= 5 rank-4 forms
// and L a 1-dim F16-subspace of N. Q is reduced up to GammaL(2,16) = N_GL(8,2)(<A>), and every
// orbit representative is tried with all 69,905 lines L.
// Written from scratch by Claude Opus 5.5 (2026-10-24), independent of orbit5.cpp.
//
// usage: fix5check [threads] [verify_extra_members_per_orbit] [min_rank4_in_Q, default 5; lower it for positive controls]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <array>
#include <algorithm>
#include <random>
#include <omp.h>

typedef uint64_t Mat; // alternating 8x8 matrix over F2, byte i = row i, bit j of row i = B(e_i, e_j)

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("CHECK FAILED %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

// ---------- GF(256) ----------
static uint8_t EXPT[512], LOGT[256];
static uint8_t mulSlow(uint8_t a, uint8_t b) {
    unsigned r = 0, x = a;
    for (int i = 0; i < 8; i++) { if (b >> i & 1) r ^= x; x <<= 1; if (x & 0x100) x ^= 0x11D; }
    return (uint8_t)r;
}
static void gfInit() {
    unsigned x = 1;
    for (int i = 0; i < 255; i++) { EXPT[i] = (uint8_t)x; LOGT[x] = (uint8_t)i; x = mulSlow((uint8_t)x, 2); }
    for (int i = 255; i < 512; i++) EXPT[i] = EXPT[i - 255];
}
static inline uint8_t mul(uint8_t a, uint8_t b) { return (a && b) ? EXPT[LOGT[a] + LOGT[b]] : 0; }
static uint8_t gpow(uint8_t a, unsigned e) { uint8_t r = 1; while (e--) r = mul(r, a); return r; }

// ---------- forms ----------
static inline int par8(unsigned x) { x &= 0xFF; x ^= x >> 4; x ^= x >> 2; x ^= x >> 1; return x & 1; }
static inline unsigned row(Mat m, int i) { return (unsigned)(m >> (8 * i)) & 0xFF; }
static inline int bil(Mat m, unsigned u, unsigned v) {
    unsigned acc = 0;
    for (int i = 0; i < 8; i++) if (u >> i & 1) acc ^= row(m, i) & v;
    return par8(acc);
}
static inline Mat setPair(Mat m, int i, int j) { return m | (1ULL << (8 * i + j)) | (1ULL << (8 * j + i)); }
// B o g : (x,y) -> B(g x, g y), g given by images of basis vectors
static Mat compose(Mat m, const uint8_t img[8]) {
    Mat r = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) if (bil(m, img[i], img[j])) r = setPair(r, i, j);
    return r;
}
static int rankBytes(uint8_t* r, int n) {
    int rk = 0;
    for (int c = 0; c < 8 && rk < n; c++) {
        int p = -1;
        for (int i = rk; i < n; i++) if (r[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        std::swap(r[p], r[rk]);
        for (int i = rk + 1; i < n; i++) if (r[i] >> c & 1) r[i] ^= r[rk];
        rk++;
    }
    return rk;
}
static inline int rankMat(Mat m) { uint8_t r[8]; for (int i = 0; i < 8; i++) r[i] = (uint8_t)row(m, i); return rankBytes(r, 8); }
static bool isAlt(Mat m) {
    for (int i = 0; i < 8; i++) { if (row(m, i) >> i & 1) return false; for (int j = 0; j < 8; j++) if ((row(m, i) >> j & 1) != (row(m, j) >> i & 1)) return false; }
    return true;
}

static int PIDX[8][8], PI_[28], PJ_[28];
static void pairInit() { int p = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { PIDX[i][j] = PIDX[j][i] = p; PI_[p] = i; PJ_[p] = j; p++; } }
static Mat maskToMat(uint32_t m) { Mat r = 0; for (int p = 0; p < 28; p++) if (m >> p & 1) r = setPair(r, PI_[p], PJ_[p]); return r; }
static uint32_t matToMask(Mat m) { uint32_t r = 0; for (int p = 0; p < 28; p++) if (row(m, PI_[p]) >> PJ_[p] & 1) r |= 1u << p; return r; }

// ---------- APN tests ----------
// (a) direct: for every a != 0, the 8 vectors B_k(a, .) have rank 7 (kernel of x -> (B_k(a,x))_k is {0,a})
static bool apnDirect(const Mat S[8]) {
    for (unsigned a = 1; a < 256; a++) {
        uint8_t w[8];
        for (int k = 0; k < 8; k++) { unsigned acc = 0; for (int i = 0; i < 8; i++) if (a >> i & 1) acc ^= row(S[k], i); w[k] = (uint8_t)acc; }
        if (rankBytes(w, 8) != 7) return false;
    }
    return true;
}
// rank histogram of the 255 nonzero elements of span(S)
static void rankHist(const Mat S[8], int cnt[9]) {
    for (int r = 0; r <= 8; r++) cnt[r] = 0;
    for (unsigned c = 1; c < 256; c++) { Mat m = 0; for (int k = 0; k < 8; k++) if (c >> k & 1) m ^= S[k]; cnt[rankMat(m)]++; }
}
// (b) counting: sum over nonzero B of (|rad B| - 1) = 255 (each a != 0 lies in at least one radical, so this forces exactly one)
static inline bool apnCount(const int cnt[9]) {
    return cnt[0] == 0 && (3 * cnt[6] + 15 * cnt[4] + 63 * cnt[2] == 255);
}
// (c) definition, from an explicit F
static bool apnDef(const uint8_t F[256]) {
    for (unsigned a = 1; a < 256; a++) {
        uint8_t seen[256] = {0};
        for (unsigned x = 0; x < 256; x++) if (++seen[F[x ^ a] ^ F[x]] > 2) return false;
    }
    return true;
}
static void formsOfF(const uint8_t F[256], Mat S[8]) {
    for (int k = 0; k < 8; k++) {
        Mat m = 0;
        for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) {
            unsigned ei = 1u << i, ej = 1u << j;
            if (((F[ei ^ ej] ^ F[ei] ^ F[ej] ^ F[0]) >> k) & 1) m = setPair(m, i, j);
        }
        S[k] = m;
    }
}
// quadratic F with the given forms: F_k(x) = sum_{i<j} B_k(e_i,e_j) x_i x_j
static void fOfForms(const Mat S[8], uint8_t F[256]) {
    for (unsigned x = 0; x < 256; x++) {
        unsigned y = 0;
        for (int k = 0; k < 8; k++) {
            int b = 0;
            for (int i = 0; i < 8; i++) if (x >> i & 1) b ^= par8(row(S[k], i) & x & ~((2u << i) - 1));
            y |= (unsigned)b << k;
        }
        F[x] = (uint8_t)y;
    }
}
static int dimSpan(const Mat* S, int n) {
    uint32_t v[32]; for (int i = 0; i < n; i++) v[i] = matToMask(S[i]);
    int rk = 0;
    for (int c = 0; c < 28; c++) {
        int p = -1; for (int i = rk; i < n; i++) if (v[i] >> c & 1) { p = i; break; }
        if (p < 0) continue; std::swap(v[p], v[rk]);
        for (int i = rk + 1; i < n; i++) if (v[i] >> c & 1) v[i] ^= v[rk];
        rk++;
    }
    return rk;
}

// ---------- F2 linear algebra on 28-bit vectors ----------
static uint32_t applyLin(const uint32_t col[28], uint32_t v) { uint32_t r = 0; for (int p = 0; p < 28; p++) if (v >> p & 1) r ^= col[p]; return r; }
// kernel of a linear map on F2^28 given by its column images
static std::vector<uint32_t> kernel28(const uint32_t col[28]) {
    std::vector<std::pair<uint32_t, uint32_t>> rows; // (image, preimage)
    for (int p = 0; p < 28; p++) rows.push_back({col[p], 1u << p});
    int rk = 0;
    for (int c = 0; c < 28; c++) {
        int piv = -1; for (int i = rk; i < 28; i++) if (rows[i].first >> c & 1) { piv = i; break; }
        if (piv < 0) continue; std::swap(rows[piv], rows[rk]);
        for (int i = 0; i < 28; i++) if (i != rk && (rows[i].first >> c & 1)) { rows[i].first ^= rows[rk].first; rows[i].second ^= rows[rk].second; }
        rk++;
    }
    std::vector<uint32_t> ker;
    for (int i = rk; i < 28; i++) ker.push_back(rows[i].second);
    return ker;
}
// coordinates w.r.t. a basis of F2^28: echelon table with tags
struct Coords {
    uint32_t vec[28], tag[28]; bool has[28];
    bool build(const std::vector<uint32_t>& basis) {
        memset(has, 0, sizeof has);
        for (size_t k = 0; k < basis.size(); k++) {
            uint32_t v = basis[k], t = 1u << k;
            for (int b = 27; b >= 0; b--) if (v >> b & 1) { if (has[b]) { v ^= vec[b]; t ^= tag[b]; } else { has[b] = true; vec[b] = v; tag[b] = t; break; } }
            if (!v) return false;
        }
        return true;
    }
    uint32_t of(uint32_t v) const {
        uint32_t t = 0;
        for (int b = 27; b >= 0; b--) if (v >> b & 1) { v ^= vec[b]; t ^= tag[b]; }
        return v ? 0xFFFFFFFFu : t;
    }
};

// ---------- subspaces of Fix (8-bit coordinates) ----------
static inline uint8_t applyFix(Mat G, unsigned c) { unsigned r = 0; for (int m = 0; m < 8; m++) if (c >> m & 1) r ^= row(G, m); return (uint8_t)r; }
// RREF (pivot = highest bit, rows by descending pivot), packed; 0 if rank < 4
static uint32_t rrefKey(uint8_t v[4]) {
    int cur = 0;
    for (int pos = 7; pos >= 0 && cur < 4; pos--) {
        int r = -1; for (int i = cur; i < 4; i++) if (v[i] >> pos & 1) { r = i; break; }
        if (r < 0) continue; std::swap(v[r], v[cur]);
        for (int i = 0; i < 4; i++) if (i != cur && (v[i] >> pos & 1)) v[i] ^= v[cur];
        cur++;
    }
    if (cur < 4) return 0;
    return (uint32_t)v[0] | (uint32_t)v[1] << 8 | (uint32_t)v[2] << 16 | (uint32_t)v[3] << 24;
}
static void enumQ(int k, int prevPivot, uint8_t rows[4], int piv[4], std::vector<uint32_t>& out) {
    if (k == 4) { out.push_back((uint32_t)rows[0] | (uint32_t)rows[1] << 8 | (uint32_t)rows[2] << 16 | (uint32_t)rows[3] << 24); return; }
    for (int p = prevPivot - 1; p >= 3 - k; p--) {
        bool ok = true; for (int i = 0; i < k; i++) if (rows[i] >> p & 1) ok = false; // earlier rows must be 0 at this pivot
        if (!ok) continue;
        unsigned freeMask = ((1u << p) - 1); for (int i = 0; i < k; i++) freeMask &= ~(1u << piv[i]);
        // enumerate submasks of freeMask
        unsigned s = 0;
        while (true) {
            rows[k] = (uint8_t)((1u << p) | s); piv[k] = p;
            enumQ(k + 1, p, rows, piv, out);
            if (s == freeMask) break;
            s = (s - freeMask) & freeMask;
        }
    }
}
// note: rows chosen later must be 0 at earlier pivots (guaranteed: their bits are below), and earlier rows 0 at later
// pivots is checked when the later pivot is chosen, so every output is a genuine RREF.

int main(int argc, char** argv) {
    int threads = argc > 1 ? atoi(argv[1]) : 4;
    int extra = argc > 2 ? atoi(argv[2]) : 1;
    int minR4 = argc > 3 ? atoi(argv[3]) : 5;
    omp_set_num_threads(threads);
    gfInit(); pairInit();
    double t0 = omp_get_wtime();

    // ---- field checks ----
    { bool prim = true; for (int i = 1; i < 255; i++) if (EXPT[i] == 1) prim = false; CHECK(prim, "2 not primitive"); }
    { bool ok = true; for (int a = 0; a < 256; a++) for (int b = 0; b < 256; b++) if (mul(a, b) != mulSlow(a, b)) ok = false; CHECK(ok, "mul table"); }
    const uint8_t zeta = EXPT[51], zinv = EXPT[204];
    CHECK(zeta != 1 && gpow(zeta, 5) == 1 && mul(zeta, zinv) == 1, "zeta order");
    printf("zeta = 2^51 = 0x%02X\n", zeta);

    // ---- APN test controls ----
    {
        struct { const char* name; unsigned e; bool expect; } pw[] = { {"x^3", 3, true}, {"x^9", 9, true}, {"x^5", 5, false}, {"x^17", 17, false}, {"x^33", 33, true}, {"x^129 (=x^(2^7+1), k=7)", 129, true} };
        for (auto& c : pw) {
            uint8_t F[256]; for (int x = 0; x < 256; x++) F[x] = gpow((uint8_t)x, c.e);
            Mat S[8]; formsOfF(F, S); int cnt[9]; rankHist(S, cnt);
            bool d = apnDef(F), di = apnDirect(S), co = apnCount(cnt);
            printf("control %-24s def=%d direct=%d count=%d dim=%d ranks 8/6/4/2/0 = %d/%d/%d/%d/%d\n", c.name, d, di, co, dimSpan(S, 8), cnt[8], cnt[6], cnt[4], cnt[2], cnt[0]);
            CHECK(d == c.expect && di == d && (co == d || dimSpan(S, 8) < 8), "control %s", c.name);
            // round trip forms -> F -> forms
            uint8_t G[256]; fOfForms(S, G); Mat S2[8]; formsOfF(G, S2); bool same = true; for (int k = 0; k < 8; k++) if (S2[k] != S[k]) same = false;
            CHECK(same && apnDef(G) == d, "fOfForms round trip %s", c.name);
        }
        // random 8-dim spans of alternating forms: three tests agree
        std::mt19937_64 rng(12345); int agree = 0, napn = 0, tot = 0;
        for (int t = 0; t < 3000; t++) {
            Mat S[8]; for (int k = 0; k < 8; k++) S[k] = maskToMat((uint32_t)(rng() & 0x0FFFFFFF));
            if (dimSpan(S, 8) < 8) continue;
            uint8_t F[256]; fOfForms(S, F); int cnt[9]; rankHist(S, cnt);
            bool d = apnDef(F), di = apnDirect(S), co = apnCount(cnt);
            tot++; if (d == di && di == co) agree++; if (d) napn++;
        }
        printf("random 8-dim subspaces: %d tested, %d APN, all three tests agree on %d\n", tot, napn, agree);
        CHECK(agree == tot, "random agreement");
        // perturbations of Gold x^3 (mostly non-APN, near an APN one)
        uint8_t F3[256]; for (int x = 0; x < 256; x++) F3[x] = gpow((uint8_t)x, 3);
        Mat S3[8]; formsOfF(F3, S3); agree = 0; tot = 0; napn = 0;
        for (int t = 0; t < 3000; t++) {
            Mat S[8]; memcpy(S, S3, sizeof S); S[rng() % 8] ^= maskToMat(1u << (rng() % 28));
            if (t & 1) S[rng() % 8] ^= maskToMat(1u << (rng() % 28));
            if (dimSpan(S, 8) < 8) continue;
            uint8_t F[256]; fOfForms(S, F); int cnt[9]; rankHist(S, cnt);
            bool d = apnDef(F), di = apnDirect(S), co = apnCount(cnt);
            tot++; if (d == di && di == co) agree++; if (d) napn++;
        }
        printf("perturbed Gold: %d tested, %d APN, agree %d\n", tot, napn, agree);
        CHECK(agree == tot, "perturbed agreement");
    }

    // ---- A on Alt(V): (A.B)(x,y) = B(A^-1 x, A^-1 y) ----
    uint8_t imgAinv[8]; for (int i = 0; i < 8; i++) imgAinv[i] = mul(zinv, (uint8_t)(1u << i));
    uint32_t T[28], TpI[28], Phi[28];
    for (int p = 0; p < 28; p++) T[p] = matToMask(compose(maskToMat(1u << p), imgAinv));
    for (int p = 0; p < 28; p++) {
        uint32_t v = 1u << p, s = v, w = v;
        for (int k = 1; k <= 4; k++) { w = applyLin(T, w); s ^= w; }
        Phi[p] = s; TpI[p] = T[p] ^ v;
        uint32_t w5 = applyLin(T, w); CHECK(w5 == v, "A^5 != 1 on forms");
    }
    std::vector<uint32_t> fixB = kernel28(TpI), nB = kernel28(Phi);
    printf("dim Fix = %zu, dim N = %zu\n", fixB.size(), nB.size());
    CHECK(fixB.size() == 8 && nB.size() == 20, "dims");
    for (auto v : fixB) CHECK(applyLin(T, v) == v, "Fix not fixed");
    for (auto v : nB) { CHECK(applyLin(Phi, v) == 0, "Phi5(A) != 0 on N"); CHECK(nB.size() == 20, ""); }
    // N is A-invariant: A(N) inside N (checked via coordinates below), and A has no fixed vector in N (direct sum)
    std::vector<uint32_t> allB = fixB; allB.insert(allB.end(), nB.begin(), nB.end());
    Coords co; CHECK(co.build(allB), "Fix + N is not a direct sum of Alt(V)");
    for (auto v : nB) CHECK((co.of(applyLin(T, v)) & 0xFF) == 0, "A(N) not in N");
    // N = image(A - 1)
    { std::vector<uint32_t> im; for (int p = 0; p < 28; p++) im.push_back(TpI[p]); Mat tmp[28]; for (int p = 0; p < 28; p++) tmp[p] = maskToMat(im[p]); int d = dimSpan(tmp, 28);
      bool inN = true; for (auto v : im) if (co.of(v) & 0xFF) inN = false; CHECK(d == 20 && inN, "image(A-1) != N"); }

    Mat fixMat[8], nMat[20];
    for (int k = 0; k < 8; k++) fixMat[k] = maskToMat(fixB[k]);
    for (int k = 0; k < 20; k++) nMat[k] = maskToMat(nB[k]);
    Mat fixElem[256]; int fixRank[256]; int fixHist[9] = {0};
    for (unsigned c = 0; c < 256; c++) { Mat m = 0; for (int k = 0; k < 8; k++) if (c >> k & 1) m ^= fixMat[k]; fixElem[c] = m; fixRank[c] = rankMat(m); if (c) fixHist[fixRank[c]]++; CHECK(isAlt(m), "not alternating"); }
    printf("Fix rank distribution 8/6/4/2: %d/%d/%d/%d\n", fixHist[8], fixHist[6], fixHist[4], fixHist[2]);
    {   Mat S[8]; memcpy(S, fixMat, sizeof S); int cnt[9]; rankHist(S, cnt); uint8_t F[256]; fOfForms(S, F);
        bool d = apnDef(F), di = apnDirect(S), c2 = apnCount(cnt);
        printf("S8 = Fix: APN def=%d direct=%d count=%d\n", d, di, c2); CHECK(d == di && di == c2, "Fix tests disagree"); }
    Mat nLo[1024], nHi[1024];
    for (unsigned c = 0; c < 1024; c++) { Mat a = 0, b = 0; for (int k = 0; k < 10; k++) if (c >> k & 1) { a ^= nMat[k]; b ^= nMat[10 + k]; } nLo[c] = a; nHi[c] = b; }
    auto nToMat = [&](uint32_t x) { return nLo[x & 1023] ^ nHi[x >> 10]; };
    uint32_t TN[20]; for (int k = 0; k < 20; k++) TN[k] = co.of(applyLin(T, nB[k])) >> 8;
    auto applyTN = [&](uint32_t x) { uint32_t r = 0; for (int k = 0; k < 20; k++) if (x >> k & 1) r ^= TN[k]; return r; };

    // Gold's S8 is invariant under multiplication by zeta: sanity that it lies in N (it has no fixed forms)
    {   uint8_t F[256]; for (int x = 0; x < 256; x++) F[x] = gpow((uint8_t)x, 3); Mat S[8]; formsOfF(F, S);
        int inN = 0; for (int k = 0; k < 8; k++) { uint32_t t = co.of(matToMask(S[k])); if ((t & 0xFF) == 0) inN++; }
        printf("Gold x^3 forms lying in N: %d of 8\n", inN); CHECK(inN == 8, "Gold not in N"); }

    // ---- the group GammaL(2,16), acting on Fix ----
    // F16 inside GF(256): 0 and 2^(17j). theta = 2 (not in F16). v = a + b*theta.
    uint8_t F16[16]; F16[0] = 0; for (int j = 0; j < 15; j++) F16[j + 1] = EXPT[17 * j];
    uint8_t decA[256], decB[256]; bool decSeen[256] = {0};
    for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) { uint8_t v = F16[a] ^ mul(F16[b], 2); CHECK(!decSeen[v], "dec collision"); decSeen[v] = true; decA[v] = F16[a]; decB[v] = F16[b]; }
    for (int i = 0; i < 16; i++) CHECK(gpow(F16[i], 16) == F16[i], "F16 elt");
    CHECK(gpow(zeta, 16) == zeta, "zeta not in F16");
    std::vector<std::array<uint8_t, 3>> gens; // (k, u, w): g(v) = a u + b w where v^(2^k) = a + b theta
    for (int k = 0; k < 4; k++) for (int u = 1; u < 256; u++) for (int w = 0; w < 256; w++) {
        bool dep = false; for (int i = 0; i < 16; i++) if (mul(F16[i], (uint8_t)u) == w) dep = true;
        if (!dep) gens.push_back({(uint8_t)k, (uint8_t)u, (uint8_t)w});
    }
    const int G = (int)gens.size();
    printf("|GammaL(2,16)| = %d (expect %d), |GL(2,16)| = %d\n", G, 255 * 240 * 4, G / 4);
    CHECK(G == 244800, "group order");
    std::vector<Mat> gfix(G);
    int badNorm = 0, badBij = 0, badFix = 0;
    std::map<int, int> conjCount; // which power of A each g conjugates A to
    #pragma omp parallel for reduction(+:badNorm,badBij,badFix) schedule(static)
    for (int gi = 0; gi < G; gi++) {
        int k = gens[gi][0]; uint8_t u = gens[gi][1], w = gens[gi][2];
        uint8_t gv[256];
        for (int v = 0; v < 256; v++) { uint8_t x = gpow((uint8_t)v, 1u << k); gv[v] = mul(decA[x], u) ^ mul(decB[x], w); }
        bool seen[256] = {0}; for (int v = 0; v < 256; v++) { if (seen[gv[v]]) badBij++; seen[gv[v]] = true; }
        // additivity and normalizing: g(zeta v) = zeta^(2^k) g(v)
        uint8_t zk = gpow(zeta, 1u << k);
        for (int v = 0; v < 256; v++) { if (gv[mul(zeta, (uint8_t)v)] != mul(zk, gv[v])) badNorm++; if (gv[v ^ 0x5A] != (gv[v] ^ gv[0x5A])) badNorm++; }
        uint8_t img[8]; for (int i = 0; i < 8; i++) img[i] = gv[1u << i];
        Mat Gm = 0;
        for (int m = 0; m < 8; m++) {
            uint32_t t = co.of(matToMask(compose(fixMat[m], img)));
            if (t >> 8) badFix++;
            Gm |= (Mat)(t & 0xFF) << (8 * m);
        }
        gfix[gi] = Gm;
    }
    printf("group checks: non-bijective %d, not normalizing/additive %d, Fix not preserved %d\n", badBij, badNorm, badFix);
    CHECK(badBij == 0 && badNorm == 0 && badFix == 0, "group");
    { std::vector<std::array<uint8_t,3>> s = gens; std::sort(s.begin(), s.end()); CHECK(std::unique(s.begin(), s.end()) == s.end(), "dup gens"); }
    { int bad = 0; for (int gi = 0; gi < G; gi++) for (unsigned c = 1; c < 256; c++) if (fixRank[applyFix(gfix[gi], c)] != fixRank[c]) bad++; CHECK(bad == 0, "rank not preserved on Fix"); }
    { int kern = 0; for (int gi = 0; gi < G; gi++) { bool id = true; for (int m = 0; m < 8; m++) if (row(gfix[gi], m) != (1u << m)) id = false; if (id) kern++; }
      printf("elements acting trivially on Fix: %d\n", kern); }

    // ---- 4-dim subspaces of Fix ----
    std::vector<uint32_t> allQ; { uint8_t r[4]; int pv[4]; enumQ(0, 8, r, pv, allQ); }
    printf("4-dim subspaces of Fix: %zu (expect 200787)\n", allQ.size());
    CHECK(allQ.size() == 200787, "Gaussian binomial");
    { std::vector<uint32_t> s = allQ; std::sort(s.begin(), s.end()); CHECK(std::unique(s.begin(), s.end()) == s.end(), "dup Q");
      for (uint32_t key : allQ) { uint8_t v[4] = {(uint8_t)key, (uint8_t)(key >> 8), (uint8_t)(key >> 16), (uint8_t)(key >> 24)}; CHECK(rrefKey(v) == key, "rref mismatch"); } }
    auto qRows = [](uint32_t key, uint8_t v[4]) { for (int i = 0; i < 4; i++) v[i] = (uint8_t)(key >> (8 * i)); };
    auto r4count = [&](uint32_t key) { uint8_t v[4]; qRows(key, v); int n = 0; for (unsigned s = 1; s < 16; s++) { unsigned c = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) c ^= v[i]; if (fixRank[c] == 4) n++; } return n; };
    std::vector<uint32_t> qual; std::map<int, int> r4dist, r4distAll;
    for (uint32_t key : allQ) { int n = r4count(key); r4distAll[n]++; if (n >= minR4) { qual.push_back(key); r4dist[n]++; } }
    std::sort(qual.begin(), qual.end());
    printf("rank-4 counts over all 4-dim Q:"); for (auto& p : r4distAll) printf(" %d:%d", p.first, p.second); printf("\n");
    printf("qualifying Q (>= %d rank-4 forms): %zu; by rank-4 count:", minR4, qual.size()); for (auto& p : r4dist) printf(" %d:%d", p.first, p.second); printf("\n");
    {   // radicals of the rank-4 forms in Fix, and whether a qualifying Q has two rank-4 forms with the same radical
        // (that alone rules out APN: every nonzero a in the shared radical would lie in two radicals)
        std::map<std::vector<uint8_t>, int> radCount; int badRad = 0;
        std::vector<std::vector<uint8_t>> radOf(256);
        for (unsigned c = 1; c < 256; c++) if (fixRank[c] == 4) {
            std::vector<uint8_t> r; for (unsigned x = 1; x < 256; x++) { unsigned acc = 0; for (int i = 0; i < 8; i++) if (x >> i & 1) acc ^= row(fixElem[c], i); if (!acc) r.push_back((uint8_t)x); }
            bool finv = true; for (uint8_t x : r) if (!std::binary_search(r.begin(), r.end(), mul(zeta, x))) finv = false;
            if (r.size() != 15 || !finv) badRad++;
            radOf[c] = r; radCount[r]++;
        }
        printf("rank-4 forms of Fix: %d distinct radicals (F16-lines of V), forms per radical:", (int)radCount.size());
        std::map<int, int> per; for (auto& p : radCount) per[p.second]++; for (auto& p : per) printf(" %d:%d", p.first, p.second); printf("\n");
        CHECK(badRad == 0, "rank-4 radical not a 4-dim A-invariant subspace");
        int shared = 0, distinct = 0;
        for (uint32_t key : qual) {
            uint8_t v[4]; for (int i = 0; i < 4; i++) v[i] = (uint8_t)(key >> (8 * i));
            std::map<std::vector<uint8_t>, int> seen; bool dup = false;
            for (unsigned s = 1; s < 16; s++) { unsigned c = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) c ^= v[i]; if (fixRank[c] == 4 && ++seen[radOf[c]] > 1) dup = true; }
            if (dup) shared++; else distinct++;
        }
        printf("qualifying Q with two rank-4 forms sharing a radical: %d; with all rank-4 radicals distinct: %d\n", shared, distinct);
    }
    auto qIndex = [&](uint32_t key) { auto it = std::lower_bound(qual.begin(), qual.end(), key); return (it != qual.end() && *it == key) ? (int)(it - qual.begin()) : -1; };

    // ---- orbits ----
    struct Orbit { uint32_t rep; int size, stab, r4; std::vector<int> members; };
    auto orbits = [&](int gcount, const char* name) {
        std::vector<char> vis(qual.size(), 0); std::vector<Orbit> out; long long total = 0; int bad = 0;
        for (size_t qi = 0; qi < qual.size(); qi++) {
            if (vis[qi]) continue;
            std::vector<int> img(gcount);
            #pragma omp parallel for schedule(static)
            for (int gi = 0; gi < gcount; gi++) {
                uint8_t v[4]; qRows(qual[qi], v); for (int i = 0; i < 4; i++) v[i] = applyFix(gfix[gi], v[i]);
                img[gi] = qIndex(rrefKey(v));
            }
            Orbit o; o.rep = qual[qi]; o.stab = 0; o.r4 = r4count(qual[qi]);
            for (int gi = 0; gi < gcount; gi++) { if (img[gi] < 0) { bad++; continue; } if (img[gi] == (int)qi) o.stab++; }
            std::sort(img.begin(), img.end()); img.erase(std::unique(img.begin(), img.end()), img.end());
            for (int m : img) if (m >= 0) { CHECK(!vis[m], "orbits overlap"); vis[m] = 1; o.members.push_back(m); }
            o.size = (int)o.members.size();
            CHECK((long long)o.size * o.stab == gcount, "orbit-stabilizer %s: %d * %d != %d", name, o.size, o.stab, gcount);
            total += o.size; out.push_back(o);
        }
        CHECK(bad == 0, "image of qualifying Q not qualifying");
        CHECK(total == (long long)qual.size(), "orbit sizes don't sum");
        printf("%s: %zu orbits, sizes sum to %lld\n", name, out.size(), total);
        for (auto& o : out) printf("  rep %08X  rank4=%d  size=%d  stab=%d\n", o.rep, o.r4, o.size, o.stab);
        return out;
    };
    std::vector<Orbit> orbGL = orbits(G / 4, "GL(2,16) orbits"); // gens with k = 0 come first
    for (int gi = 0; gi < G / 4; gi++) CHECK(gens[gi][0] == 0, "GL ordering");
    std::vector<Orbit> orbGamma = orbits(G, "GammaL(2,16) orbits");
    printf("[%.1f s]\n", omp_get_wtime() - t0);

    // ---- lines L in N ----
    std::vector<std::array<uint32_t, 4>> lines;
    for (uint32_t b = 1; b < (1u << 20); b++) {
        uint32_t v[5]; v[0] = b; for (int i = 1; i < 5; i++) v[i] = applyTN(v[i - 1]);
        if ((v[4] ^ v[3] ^ v[2] ^ v[1] ^ v[0]) != 0) { CHECK(false, "Phi5 on N coords"); }
        uint32_t mn = 0xFFFFFFFF; int nz = 0;
        for (unsigned s = 1; s < 16; s++) { uint32_t x = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) x ^= v[i]; if (x) nz++; mn = std::min(mn, x); }
        CHECK(nz == 15, "line not 4-dim");
        if (mn == b) lines.push_back({v[0], v[1], v[2], v[3]});
    }
    printf("lines L (1-dim F16-subspaces of N): %zu (expect 69905)\n", lines.size());
    CHECK(lines.size() == 69905, "line count");
    const int NL = (int)lines.size();

    // ---- search ----
    struct Res { long long apn = 0, napn = 0, mismatch = 0, defChecked = 0, defBad = 0; std::map<std::array<int, 5>, long long> hist; std::map<int, long long> r4all; int maxR4 = 0; };
    auto search = [&](uint32_t qkey) {
        uint8_t qv[4]; qRows(qkey, qv);
        Mat qm[16]; for (unsigned s = 0; s < 16; s++) { unsigned c = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) c ^= qv[i]; qm[s] = fixElem[c]; }
        int nT = omp_get_max_threads(); std::vector<Res> part(nT);
        #pragma omp parallel
        {
            Res& R = part[omp_get_thread_num()];
            #pragma omp for schedule(dynamic, 256)
            for (int li = 0; li < NL; li++) {
                Mat lm[16]; for (unsigned s = 0; s < 16; s++) { uint32_t x = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) x ^= lines[li][i]; lm[s] = nToMat(x); }
                int cnt[9] = {0};
                for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) if (a | b) cnt[rankMat(qm[a] ^ lm[b])]++;
                bool c1 = apnCount(cnt);
                Mat S[8] = {qm[1], qm[2], qm[4], qm[8], lm[1], lm[2], lm[4], lm[8]};
                bool c2 = apnDirect(S);
                if (c1 != c2) R.mismatch++;
                R.r4all[cnt[4]]++; R.maxR4 = std::max(R.maxR4, cnt[4]);
                if (c1) {
                    R.apn++;
                    int r4Q = 0; for (int a = 1; a < 16; a++) if (rankMat(qm[a]) == 4) r4Q++;
                    R.hist[{cnt[8], cnt[6], cnt[4], cnt[2], r4Q}]++;
                    if (R.defChecked < 5) { uint8_t F[256]; fOfForms(S, F); R.defChecked++; if (!apnDef(F)) R.defBad++; }
                } else R.napn++;
            }
        }
        Res T_; for (auto& R : part) { T_.apn += R.apn; T_.napn += R.napn; T_.mismatch += R.mismatch; T_.defChecked += R.defChecked; T_.defBad += R.defBad; T_.maxR4 = std::max(T_.maxR4, R.maxR4);
            for (auto& h : R.hist) T_.hist[h.first] += h.second; for (auto& h : R.r4all) T_.r4all[h.first] += h.second; }
        return T_;
    };
    auto sameRes = [](const Res& a, const Res& b) { return a.apn == b.apn && a.napn == b.napn && a.hist == b.hist && a.r4all == b.r4all; };
    {   // positive control through the search's own plumbing: Gold x^3's S8 is a 2-dim F16-subspace of N.
        // Rebuild it from N-coordinates via nToMat, check it is F16-closed (A maps it into itself), and count its ranks.
        uint8_t F[256]; for (int x = 0; x < 256; x++) F[x] = gpow((uint8_t)x, 3); Mat S[8]; formsOfF(F, S);
        Mat S2[8]; uint32_t c[8]; for (int k = 0; k < 8; k++) { c[k] = co.of(matToMask(S[k])) >> 8; S2[k] = nToMat(c[k]); }
        bool same = true; for (int k = 0; k < 8; k++) if (S2[k] != S[k]) same = false;
        Mat S3[16]; for (int k = 0; k < 8; k++) { S3[k] = S2[k]; S3[8 + k] = nToMat(applyTN(c[k])); }
        int cnt[9]; rankHist(S2, cnt);
        printf("Gold via N-coords: rebuilt=%d, dim(S8 + A S8)=%d, ranks 8/6/4/2 = %d/%d/%d/%d, count-test APN=%d, direct=%d\n",
               same, dimSpan(S3, 16), cnt[8], cnt[6], cnt[4], cnt[2], apnCount(cnt), apnDirect(S2));
        CHECK(same && dimSpan(S3, 16) == 8 && apnCount(cnt) && apnDirect(S2), "Gold via N coords");
    }
    long long totalApnPairs = 0, totalPairs = 0;
    std::mt19937 rng(777);
    for (size_t oi = 0; oi < orbGamma.size(); oi++) {
        Orbit& o = orbGamma[oi];
        Res R = search(o.rep);
        printf("orbit %zu rep %08X rank4(Q)=%d size=%d: APN %lld of %d, test mismatches %lld, def-checked %lld (bad %lld), max rank-4 count over all L = %d\n",
               oi, o.rep, o.r4, o.size, R.apn, NL, R.mismatch, R.defChecked, R.defBad, R.maxR4);
        printf("    rank-4 count of S8 over all L:"); for (auto& h : R.r4all) printf(" %d:%lld", h.first, h.second); printf("\n");
        for (auto& h : R.hist) printf("    APN ranks 8/6/4/2 = %d/%d/%d/%d (rank-4 in Q %d): %lld\n", h.first[0], h.first[1], h.first[2], h.first[3], h.first[4], h.second);
        CHECK(R.mismatch == 0 && R.defBad == 0, "APN test disagreement");
        {   // the definition on a sample of this rep's S8 (built the same way as in the search)
            uint8_t qv[4]; qRows(o.rep, qv); int agreeDef = 0, n = 500;
            for (int t = 0; t < n; t++) {
                const auto& L = lines[rng() % NL];
                Mat S[8]; for (int i = 0; i < 4; i++) { S[i] = fixElem[qv[i]]; S[4 + i] = nToMat(L[i]); }
                uint8_t F[256]; fOfForms(S, F); if (apnDef(F) == apnDirect(S)) agreeDef++;
            }
            printf("    definition agrees with the form test on %d of %d sampled L\n", agreeDef, n);
            CHECK(agreeDef == n, "definition sample");
        }
        totalApnPairs += R.apn * o.size; totalPairs += (long long)NL * o.size;
        for (int e = 0; e < extra && o.size > 1; e++) {
            int m; do { m = o.members[rng() % o.size]; } while (qual[m] == o.rep);
            Res R2 = search(qual[m]);
            printf("    member %08X gives the same result: %s\n", qual[m], sameRes(R, R2) ? "yes" : "NO");
            CHECK(sameRes(R, R2), "orbit member result differs");
        }
        fflush(stdout);
    }
    printf("total (Q, L) pairs %lld, APN pairs %lld\n", totalPairs, totalApnPairs);
    printf("[%.1f s] failures: %d\n", omp_get_wtime() - t0, failures);
    return failures ? 1 : 0;
}
