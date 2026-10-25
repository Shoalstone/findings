// symp: every quadratic (8,8)-APN component space S8 invariant under a given A in GL(8,2) of odd prime order p that
// contains a rank-4 form B1 whose <A>-orbit has p elements with pairwise disjoint radicals (the "p-cycle" case of
// Open Problem VI.1's symmetric search). Written independently of orbit5's search (only apn8's helpers are shared).
// Claude Opus 5.5, for Rowen, 2026-10-24.
//
// Math. p odd, so F_2[A] is semisimple on Alt(V) (28-dim; A acts by (A.B)(x,y) = B(A^-1 x, A^-1 y)):
//   Alt(V) = Fix + I_1 + .. + I_r,  I_j = ker g_j(A) for the irreducible factors g_j != x+1 of x^p - 1.
// I_j is an F_q-space (q = 2^deg g_j, x acting as A). An A-invariant S8 is the sum of its intersections with these
// parts (the projections are polynomials in A), each an F_q-subspace; it contains F_2[A]B1 = <pi_0 B1> + sum_j F_q pi_j B1.
// So: for each split (k_0 = dim S8 & Fix, k_j = F_q-dim of S8 & I_j, sum d_j k_j = 8), enumerate the subspaces of
// each part containing pi_j B1, each once (reduced echelon over F_q in a complement of pi_j B1), pruning on kernel
// defects: an n-dim subspace of an APN S8 has rank{B(a,.) : B in it} >= n - 1 at every a != 0 (apn8's
// max_kernel <= 9 - n). At n = 8 that is APN. A form is kept in the blocks' coordinates, so "S8 & I_j" is a set of
// coordinate vectors over F_q.
//
// usage:
//   symp info <A>                       decomposition of Alt(V) under A, possible splits
//   symp enumtest <A>                   the subspace enumerator against Gaussian binomials (each subspace once)
//   symp run <A> <out> [threads] [M r]  the p-cycle search (B1 over all rank-4 forms; 1/M sample, residue r)
//   symp ext <A> <space> <g> <out>      extension step alone from B1 = component g (1..255, first one not in Fix
//                                       from g on) of a known space; reports whether the space itself comes back.
//                                       <space>: x3 (Gold x^3, GF(256) mod 0x11B) or <seedfile>:<line> (apn8 format)
//   symp autfind <seedfile> <line> <p>  order-p elements of Aut(S8) for a seed with >= 3 rank-4 components
// <A>: c5a (Phi5^2), c5b (Phi5 (x+1)^4), c7a (p1^2 (x+1)^2), c7b (p1 p2 (x+1)^2), c7c (p1 (x+1)^5),
//      gf<p> (multiplication by an order-p element of GF(256) mod 0x11B, p | 255), or 16 hex digits (columns A e_0..A e_7,
//      two digits each, bit r = row r). p1 = x^3+x+1, p2 = x^3+x^2+1.
#define main apn8_main
#include "apn8.cpp"
#undef main
#include <mutex>
#include <set>
#include <array>
#include <map>
#include <functional>
typedef uint16_t u16;

// ---------- 8x8 matrices over F_2 as column bytes ----------
struct M8 { u8 c[8]; };
static inline u8 mv(const M8& A, u8 x) { u8 r = 0; for (int i = 0; i < 8; i++) if (x >> i & 1) r ^= A.c[i]; return r; }
static M8 mm(const M8& A, const M8& B) { M8 C; for (int i = 0; i < 8; i++) C.c[i] = mv(A, B.c[i]); return C; }
static M8 ident8() { M8 I; for (int i = 0; i < 8; i++) I.c[i] = (u8)(1 << i); return I; }
static bool eq8(const M8& A, const M8& B) { return memcmp(A.c, B.c, 8) == 0; }
static bool inv8(const M8& A, M8& R) {
    int t[256]; for (int i = 0; i < 256; i++) t[i] = -1;
    for (unsigned x = 0; x < 256; x++) { u8 y = mv(A, (u8)x); if (t[y] >= 0) return false; t[y] = (int)x; }
    for (int i = 0; i < 8; i++) R.c[i] = (u8)t[1 << i];
    return true;
}
static int order8(const M8& A) { M8 P = A; for (int k = 1; k <= 4096; k++) { if (eq8(P, ident8())) return k; P = mm(A, P); } return -1; }
static void print_m8(FILE* f, const M8& A) { for (int i = 0; i < 8; i++) fprintf(f, "%02x", A.c[i]); }

// ---------- linear maps on the 28-bit coefficient space ----------
struct L28 {
    u32 img[28]; u32 tab[4][256];
    void build() {
        for (int b = 0; b < 4; b++) for (int x = 0; x < 256; x++) {
            u32 r = 0; for (int k = 0; k < 8; k++) if ((x >> k & 1) && 8 * b + k < 28) r ^= img[8 * b + k];
            tab[b][x] = r;
        }
    }
    u32 operator()(u32 c) const { return tab[0][c & 255] ^ tab[1][c >> 8 & 255] ^ tab[2][c >> 16 & 255] ^ tab[3][c >> 24 & 255]; }
};
static L28 l28_id() { L28 L; for (int k = 0; k < 28; k++) L.img[k] = 1u << k; L.build(); return L; }
static L28 l28_comp(const L28& X, const L28& Y) { L28 L; for (int k = 0; k < 28; k++) L.img[k] = X(Y.img[k]); L.build(); return L; }

static u64 CM[4][256];
static void init_cm() {
    for (int b = 0; b < 4; b++) for (int x = 0; x < 256; x++) {
        u64 m = 0; for (int k = 0; k < 8; k++) if ((x >> k & 1) && 8 * b + k < 28) m ^= coef_to_mat(1u << (8 * b + k));
        CM[b][x] = m;
    }
}
static inline u64 c2m(u32 c) { return CM[0][c & 255] ^ CM[1][c >> 8 & 255] ^ CM[2][c >> 16 & 255] ^ CM[3][c >> 24 & 255]; }

static inline int bform(u64 m, u8 u, u8 v) { return __popcnt((unsigned)(rvec(m, u) & v)) & 1; }
// the form B(Ai x, Ai y) as a coefficient vector
static u32 act_form(u64 m, const M8& Ai) {
    u32 c = 0; for (int k = 0; k < 28; k++) if (bform(m, Ai.c[KI[k]], Ai.c[KJ[k]])) c |= 1u << k;
    return c;
}
static L28 alt_action(const M8& A) { // B -> A.B, (A.B)(x,y) = B(A^-1 x, A^-1 y)
    M8 Ai; inv8(A, Ai); L28 L;
    for (int k = 0; k < 28; k++) L.img[k] = act_form(coef_to_mat(1u << k), Ai);
    L.build(); return L;
}

// echelon with tags (which combination of inputs each pivot row is)
struct EchT {
    u32 piv[28], tag[28];
    EchT() { memset(piv, 0, sizeof piv); memset(tag, 0, sizeof tag); }
    void reduce(u32& x, u32& t) const { for (int b = 27; b >= 0 && x; b--) if ((x >> b & 1) && piv[b]) { x ^= piv[b]; t ^= tag[b]; } }
    bool add(u32 x, u32 t) { reduce(x, t); if (!x) return false; int b = 31 - (int)__lzcnt(x); piv[b] = x; tag[b] = t; return true; }
};
static std::vector<u32> kernel28(const L28& L) {
    EchT e; std::vector<u32> ker;
    for (int k = 0; k < 28; k++) { u32 x = L.img[k], t = 1u << k; e.reduce(x, t); if (!x) ker.push_back(t); else e.add(x, t); }
    return ker;
}

// ---------- polynomials over F_2 (bit i = x^i) ----------
static u32 pmulmod(u32 a, u32 b, u32 g, int d) {
    u32 r = 0; while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a >> d & 1) a ^= g; }
    return r;
}

// ---------- the decomposition ----------
struct Block {
    int d = 1, m = 0, off = 0, q = 2; u32 g = 3;
    std::vector<u32> elem;          // elem[i*d + t] = A^t v_i: F_2 basis of the part, coordinate bit off + i*d + t
    std::vector<u16> mul; std::vector<u8> mulx;
    std::vector<u32> tab;           // tab[i*q + c] = the form of c in F_q-coordinate i
    u32 form(const u8* w) const { u32 r = 0; for (int i = 0; i < m; i++) r ^= tab[(size_t)i * q + w[i]]; return r; }
};
struct Decomp {
    int p = 0; M8 A{}; std::vector<L28> pw;   // pw[i] = (B -> A^i.B)
    std::vector<Block> B;                      // B[0] = Fix (d = 1)
    L28 coord;                                 // form -> coordinates in the blocks' basis
    std::vector<M8> Apow;                      // A^i on V
    void coords(u32 c, int b, u8* v) const {
        u32 t = coord(c); const Block& bl = B[b];
        for (int i = 0; i < bl.m; i++) v[i] = (u8)((t >> (bl.off + i * bl.d)) & ((1u << bl.d) - 1));
    }
};

static bool is_odd_prime(int p) { if (p < 3 || p % 2 == 0) return false; for (int i = 3; i * i <= p; i += 2) if (p % i == 0) return false; return true; }

static bool build_decomp(const M8& A, Decomp& D, bool verbose) {
    int p = order8(A);
    if (!is_odd_prime(p)) { fprintf(stderr, "A has order %d, not an odd prime\n", p); return false; }
    D.p = p; D.A = A;
    D.Apow.resize(p); D.Apow[0] = ident8(); for (int i = 1; i < p; i++) D.Apow[i] = mm(A, D.Apow[i - 1]);
    L28 a1 = alt_action(A);
    D.pw.resize(p); D.pw[0] = l28_id(); for (int i = 1; i < p; i++) D.pw[i] = l28_comp(a1, D.pw[i - 1]);
    { L28 t = l28_comp(a1, D.pw[p - 1]); for (int k = 0; k < 28; k++) if (t.img[k] != 1u << k) { fprintf(stderr, "A_alt^p != 1\n"); return false; } }
    int d = 1; while ((1u << d) % (unsigned)p != 1) d++;
    std::vector<u32> gs; gs.push_back(3);
    for (u32 g = 1u << d; g < (2u << d); g++) {
        if (!(g & 1) || !(__popcnt(g) & 1)) continue;          // x | g or (x+1) | g
        u32 r = 1; for (int i = 0; i < p; i++) r = pmulmod(r, 2, g, d);  // x^p mod g
        if (r == 1) gs.push_back(g);
    }
    D.B.clear(); int tot = 0;
    for (u32 g : gs) {
        int dg = g == 3 ? 1 : d;
        L28 L; for (int k = 0; k < 28; k++) { u32 r = 0; for (int t = 0; t <= dg; t++) if (g >> t & 1) r ^= D.pw[t].img[k]; L.img[k] = r; }
        L.build();
        std::vector<u32> ker = kernel28(L);
        if (ker.empty() && g != 3) continue;
        Block bl; bl.d = dg; bl.g = g; bl.q = 1 << dg;
        Ech28 span;
        for (u32 v : ker) {
            if (!span.reduce(v)) continue;
            u32 w = v;
            for (int t = 0; t < dg; t++) { if (t) w = D.pw[1](w); bl.elem.push_back(w); if (!span.add(w)) { fprintf(stderr, "F_q basis failed\n"); return false; } }
        }
        if (bl.elem.size() != ker.size()) { fprintf(stderr, "part for g=%x: dim %zu not a multiple of %d\n", g, ker.size(), dg); return false; }
        bl.m = (int)ker.size() / dg; bl.off = tot; tot += (int)ker.size();
        bl.mul.resize((size_t)bl.q * bl.q); bl.mulx.resize(bl.q);
        for (int a = 0; a < bl.q; a++) for (int b = 0; b < bl.q; b++) bl.mul[(size_t)a * bl.q + b] = (u16)(dg == 1 ? (a & b) : pmulmod(a, b, g, dg));
        for (int a = 0; a < bl.q; a++) bl.mulx[a] = (u8)(dg == 1 ? a : pmulmod(a, 2, g, dg));
        bl.tab.assign((size_t)bl.m * bl.q, 0);
        for (int i = 0; i < bl.m; i++) for (int c = 0; c < bl.q; c++) { u32 r = 0; for (int t = 0; t < dg; t++) if (c >> t & 1) r ^= bl.elem[i * dg + t]; bl.tab[(size_t)i * bl.q + c] = r; }
        D.B.push_back(bl);
    }
    if (tot != 28) { fprintf(stderr, "parts sum to %d, not 28\n", tot); return false; }
    EchT e; int idx = 0;
    for (auto& bl : D.B) for (u32 w : bl.elem) { if (!e.add(w, 1u << idx)) { fprintf(stderr, "parts not independent\n"); return false; } idx++; }
    for (int k = 0; k < 28; k++) { u32 x = 1u << k, t = 0; e.reduce(x, t); D.coord.img[k] = t; }
    D.coord.build();
    // self-checks: coordinates round-trip; A acts as x on each part; Fix really fixed
    Rng R(12345);
    for (int it = 0; it < 200; it++) {
        u32 c = (u32)R.next() & 0x0FFFFFFF, back = 0;
        for (int b = 0; b < (int)D.B.size(); b++) { u8 v[28]; D.coords(c, b, v); back ^= D.B[b].form(v); }
        if (back != c) { fprintf(stderr, "coordinate round trip failed\n"); return false; }
        for (int b = 0; b < (int)D.B.size(); b++) {
            const Block& bl = D.B[b]; u8 v[28], w[28];
            for (int i = 0; i < bl.m; i++) { v[i] = (u8)R.below(bl.q); w[i] = bl.mulx[v[i]]; }
            if (D.pw[1](bl.form(v)) != bl.form(w)) { fprintf(stderr, "A != x on part %d\n", b); return false; }
        }
    }
    if (verbose) {
        printf("A = "); print_m8(stdout, A); printf("  order %d, deg of factors %d\n", p, d);
        for (auto& bl : D.B) printf("  part g=%x: F_%d-dim %d (F_2-dim %d)\n", bl.g, bl.q, bl.m, bl.m * bl.d);
    }
    return true;
}

// ---------- named classes ----------
static void companion(M8& A, int s, u32 g, int d) {
    for (int i = 0; i + 1 < d; i++) A.c[s + i] = (u8)(1 << (s + i + 1));
    u8 last = 0; for (int t = 0; t < d; t++) if (g >> t & 1) last |= (u8)(1 << (s + t));
    A.c[s + d - 1] = last;
}
static unsigned gmul11b(unsigned a, unsigned b) { unsigned r = 0; while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a & 0x100) a ^= 0x11B; } return r; }
static bool parse_A(const char* s, M8& A) {
    A = ident8(); std::string n = s;
    const u32 PHI5 = 0x1F, P1 = 0xB, P2 = 0xD;
    if (n == "c5a") { companion(A, 0, PHI5, 4); companion(A, 4, PHI5, 4); return true; }
    if (n == "c5b") { companion(A, 0, PHI5, 4); return true; }
    if (n == "c7a") { companion(A, 0, P1, 3); companion(A, 3, P1, 3); return true; }
    if (n == "c7b") { companion(A, 0, P1, 3); companion(A, 3, P2, 3); return true; }
    if (n == "c7c") { companion(A, 0, P1, 3); return true; }
    if (n.size() == 4 && n.compare(0, 3, "c3k") == 0 && n[3] >= '1' && n[3] <= '4') { // (x^2+x+1)^k (x+1)^(8-2k)
        for (int i = 0; i < n[3] - '0'; i++) companion(A, 2 * i, 7, 2);
        return true;
    }
    if (n.size() > 2 && n.compare(0, 2, "gf") == 0) {
        int p = atoi(n.c_str() + 2);
        if (p < 2 || 255 % p) { fprintf(stderr, "gf<p>: p must divide 255\n"); return false; }
        for (unsigned l = 2; l < 256; l++) {
            unsigned x = 1; for (int i = 0; i < p; i++) x = gmul11b(x, l);
            if (x == 1) { for (int i = 0; i < 8; i++) A.c[i] = (u8)gmul11b(l, 1u << i); return true; }
        }
        return false;
    }
    if (n.size() == 16) {
        for (int i = 0; i < 8; i++) { unsigned v; if (sscanf(n.c_str() + 2 * i, "%2x", &v) != 1) return false; A.c[i] = (u8)v; }
        M8 Ai; if (!inv8(A, Ai)) { fprintf(stderr, "A not invertible\n"); return false; }
        return true;
    }
    fprintf(stderr, "unknown A '%s'\n", s); return false;
}

// ---------- incremental kernel-defect state ----------
// for every a != 0: an echelon of {B(a,.)} over the forms so far, and whether one form already failed to raise its rank
struct KS {
    u8 piv[256][8]; u8 def[256]; int n; u32 basis[8];
    KS() { memset(this, 0, sizeof *this); }
    bool add(u32 c) {
        if (n >= 8) return false;
        u64 m = c2m(c); u8 rows[8]; for (int i = 0; i < 8; i++) rows[i] = row(m, i);
        u8 rv[256]; rv[0] = 0;
        for (unsigned a = 1; a < 256; a++) rv[a] = rv[a & (a - 1)] ^ rows[_tzcnt_u32(a)];
        for (unsigned a = 1; a < 256; a++) {
            u8 x = rv[a]; u8* P = piv[a];
            while (x) { int b = 31 - (int)__lzcnt(x); if (!P[b]) { P[b] = x; break; } x ^= P[b]; }
            if (!x) { if (def[a]) return false; def[a] = 1; }
        }
        basis[n++] = c; return true;
    }
};
struct FormList { // no pruning: for the enumerator test
    u32 basis[28]; int n = 0;
    bool add(u32 c) { basis[n++] = c; return true; }
};

// ---------- the subspace enumerator ----------
// Rows of each block in reduced echelon form over F_q on the allowed positions P (a complement of pi_b B1):
// pivot = highest allowed position with a nonzero entry, normalised to 1; rows in increasing pivot order; a new row
// is zero above its pivot and at earlier pivots (earlier rows are zero at its pivot automatically). Each row adds the
// d forms x^s * row (s < d), i.e. its F_q-line.
struct Search {
    const Decomp* D = nullptr;
    std::vector<int> pos[16]; int extra[16] = {0}; int ord[16] = {0}; int nord = 0;
    int pivsel[16][28];
    u64 adds[29] = {0}; u64 tries = 0;
};
template <class St, class Leaf>
static void dfs(Search& S, int oi, int r, int lastpi, const St& st, Leaf& leaf) {
    if (oi == S.nord) { leaf(st); return; }
    int b = S.ord[oi]; const Block& bl = S.D->B[b]; const std::vector<int>& P = S.pos[b];
    int n = (int)P.size(), k = S.extra[b];
    if (r == k) { dfs(S, oi + 1, 0, -1, st, leaf); return; }
    for (int pi = lastpi + 1; pi <= n - (k - r); pi++) {
        int fr[32], nf = 0, ip = 0;
        for (int j = 0; j < pi; j++) { if (ip < r && S.pivsel[b][ip] == j) { ip++; continue; } fr[nf++] = P[j]; }
        int v[32] = {0}; v[P[pi]] = 1;
        S.pivsel[b][r] = pi;
        for (;;) {
            S.tries++;
            St c = st; bool ok = true; u8 w[32];
            for (int i = 0; i < bl.m; i++) w[i] = (u8)v[i];
            for (int s = 0; s < bl.d; s++) {
                if (s) for (int i = 0; i < bl.m; i++) w[i] = bl.mulx[w[i]];
                if (!c.add(bl.form(w))) { ok = false; break; }
                S.adds[c.n]++;
            }
            if (ok) dfs(S, oi, r + 1, pi, c, leaf);
            int t = 0; while (t < nf) { if (++v[fr[t]] < bl.q) break; v[fr[t]] = 0; t++; }
            if (t == nf) break;
        }
    }
}

static std::vector<int> g_only; // if set (--split a,b,..), only this split is searched
static void all_splits(const Decomp& D, const int* mink, int b, int rem, std::vector<int>& k, std::vector<std::vector<int>>& out) {
    if (b == (int)D.B.size()) { if (rem == 0 && (g_only.empty() || g_only == k)) out.push_back(k); return; }
    const Block& bl = D.B[b];
    for (int kb = mink[b]; kb <= bl.m && kb * bl.d <= rem; kb++) { k[b] = kb; all_splits(D, mink, b + 1, rem - kb * bl.d, k, out); }
}
static std::string split_name(const Decomp& D, const std::vector<int>& k, const bool* pres) {
    std::string s; char buf[64];
    for (int b = 0; b < (int)D.B.size(); b++) {
        if (b == 0) snprintf(buf, sizeof buf, "Fix:%d%s", k[b], pres[b] ? "(f)" : "");
        else snprintf(buf, sizeof buf, " F%d[%x]:%d%s", D.B[b].q, D.B[b].g, k[b], pres[b] ? "(n)" : "");
        s += buf;
    }
    return s;
}

// a search's per-split counters
struct SplitStat { u64 b1 = 0, tries = 0, leaves = 0; u64 adds[29] = {0}; };

// prepare: coordinates of B1 in each part, which parts it meets, W = F_2[A]B1 in a KS; false if W already fails
// The start of a search: an A-invariant subspace W0 (F_2[A]B1, or a given one), its F_q-dim sdim[b] in each part and
// the pivot positions excl[b] of its echelon form there (the search's rows live on the other positions, a complement).
struct Prep { u8 co[16][28]; bool pres[16]; int h[16]; int sdim[16]; std::vector<int> excl[16]; KS w; int wdim = 0; };
static bool prepare(const Decomp& D, u32 B1, Prep& P) {
    P.w = KS();
    for (int b = 0; b < (int)D.B.size(); b++) {
        const Block& bl = D.B[b]; D.coords(B1, b, P.co[b]);
        P.pres[b] = false; P.h[b] = -1; P.excl[b].clear();
        for (int i = 0; i < bl.m; i++) if (P.co[b][i]) { P.pres[b] = true; P.h[b] = i; }
        P.sdim[b] = P.pres[b] ? 1 : 0;
        if (!P.pres[b]) continue;
        P.excl[b].push_back(P.h[b]);
        u8 w[28]; memcpy(w, P.co[b], bl.m);
        for (int s = 0; s < bl.d; s++) { if (s) for (int i = 0; i < bl.m; i++) w[i] = bl.mulx[w[i]]; if (!P.w.add(bl.form(w))) return false; }
    }
    P.wdim = P.w.n;
    return true;
}
// start from a given A-invariant subspace (basis of n forms): its projections onto the parts are its intersections
// with them; echelon form over F_q in each. Returns false if it isn't invariant or fails the kernel test.
static bool prepare_sub(const Decomp& D, const u32* basis, int n, Prep& P) {
    P.w = KS(); int tot = 0;
    for (int b = 0; b < (int)D.B.size(); b++) {
        const Block& bl = D.B[b]; P.excl[b].clear(); P.sdim[b] = 0; P.pres[b] = false; P.h[b] = -1;
        std::vector<std::vector<u8>> rows; // echelon over F_q, pivot = highest nonzero position, normalised to 1
        std::vector<int> piv;
        for (int j = 0; j < n; j++) {
            std::vector<u8> v(bl.m); D.coords(basis[j], b, v.data());
            for (size_t r = 0; r < rows.size(); r++) { u8 c = v[piv[r]]; if (c) for (int i = 0; i < bl.m; i++) v[i] ^= (u8)bl.mul[(size_t)c * bl.q + rows[r][i]]; }
            int h = -1; for (int i = 0; i < bl.m; i++) if (v[i]) h = i;
            if (h < 0) continue;
            u8 inv = 0; for (int c = 1; c < bl.q; c++) if (bl.mul[(size_t)v[h] * bl.q + c] == 1) inv = (u8)c;
            for (int i = 0; i < bl.m; i++) v[i] = (u8)bl.mul[(size_t)inv * bl.q + v[i]];
            for (size_t r = 0; r < rows.size(); r++) { u8 c = rows[r][h]; if (c) for (int i = 0; i < bl.m; i++) rows[r][i] ^= (u8)bl.mul[(size_t)c * bl.q + v[i]]; }
            rows.push_back(v); piv.push_back(h);
        }
        P.sdim[b] = (int)rows.size(); P.pres[b] = !rows.empty(); P.excl[b] = piv;
        for (auto& v : rows) {
            u8 w[28]; for (int i = 0; i < bl.m; i++) w[i] = v[i];
            for (int s = 0; s < bl.d; s++) { if (s) for (int i = 0; i < bl.m; i++) w[i] = bl.mulx[w[i]]; if (!P.w.add(bl.form(w))) return false; tot++; }
        }
    }
    P.wdim = P.w.n;
    if (tot != n) { fprintf(stderr, "prepare_sub: start space not A-invariant (%d vs %d)\n", tot, n); return false; }
    return true;
}
static void search_setup(const Decomp& D, const Prep& P, const std::vector<int>& ks, Search& S) {
    S.D = &D; S.nord = 0;
    for (int b = (int)D.B.size() - 1; b >= 0; b--) { // big-field parts first
        S.pos[b].clear();
        for (int i = 0; i < D.B[b].m; i++) if (std::find(P.excl[b].begin(), P.excl[b].end(), i) == P.excl[b].end()) S.pos[b].push_back(i);
        S.extra[b] = ks[b] - P.sdim[b];
        if (S.extra[b] > 0) S.ord[S.nord++] = b;
    }
}
// run every split from a prepared start; leaf(KS) gets each APN S8
template <class Leaf>
static void extend(const Decomp& D, const Prep& P, Leaf& leaf, std::map<std::string, SplitStat>& stats) {
    std::vector<std::vector<int>> sp; std::vector<int> k(D.B.size(), 0);
    all_splits(D, P.sdim, 0, 8, k, sp);
    for (auto& ks : sp) {
        Search S; search_setup(D, P, ks, S);
        SplitStat& st = stats[split_name(D, ks, P.pres)];
        u64 leaves = 0;
        auto lf = [&](const KS& s) { leaves++; leaf(s); };
        dfs(S, 0, 0, -1, P.w, lf);
        st.b1++; st.tries += S.tries; st.leaves += leaves;
        for (int i = 0; i < 29; i++) st.adds[i] += S.adds[i];
    }
}

// ---------- results ----------
struct Global {
    std::mutex mu; std::set<std::array<u32, 8>> seen; std::vector<Found> found;
    std::map<std::string, u64> dist; FILE* out = nullptr; u64 leaves = 0, bugs = 0;
};
static std::string dist_name(const int* cnt) { char b[64]; snprintf(b, sizeof b, "8:%d 6:%d 4:%d 2:%d", cnt[4], cnt[3], cnt[2], cnt[1]); return b; }
static bool invariant(const Decomp& D, const u32* basis) {
    Ech28 e; for (int i = 0; i < 8; i++) e.add(basis[i]);
    if (e.n != 8) return false;
    for (int i = 0; i < 8; i++) if (e.reduce(D.pw[1](basis[i]))) return false;
    return true;
}
static void record(const Decomp& D, Global& G, const KS& s, u32 B1) {
    u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(s.basis[i]);
    bool ok = s.n == 8 && max_kernel(m, 8) == 1 && invariant(D, s.basis);
    Found f; rank_dist(m, 8, f.cnt);
    memcpy(f.basis, s.basis, sizeof f.basis); canon8(f.basis);
    { Ech28 e; for (int i = 0; i < 8; i++) e.add(f.basis[i]); if (e.reduce(B1)) ok = false; }
    std::array<u32, 8> key; memcpy(key.data(), f.basis, sizeof f.basis);
    std::lock_guard<std::mutex> g(G.mu);
    G.leaves++;
    if (!ok) { G.bugs++; printf("BUG: leaf fails the APN/invariance/B1 check\n"); return; }
    if (!G.seen.insert(key).second) return;
    G.found.push_back(f); G.dist[dist_name(f.cnt)]++;
    if (G.out) { print_found(G.out, f); fflush(G.out); }
    if (f.cnt[2] >= 5) { printf("FOUND "); print_found(stdout, f); fflush(stdout); }
}

// ---------- radicals ----------
static int radical(u64 m, u64* mask, u8* basis) { // mask: 256-bit membership
    u8 rows[8]; for (int i = 0; i < 8; i++) rows[i] = row(m, i);
    u8 rv[256]; rv[0] = 0; for (unsigned a = 1; a < 256; a++) rv[a] = rv[a & (a - 1)] ^ rows[_tzcnt_u32(a)];
    mask[0] = mask[1] = mask[2] = mask[3] = 0; int d = 0; u8 piv[8] = {0};
    for (unsigned a = 0; a < 256; a++) if (!rv[a]) {
        mask[a >> 6] |= 1ull << (a & 63);
        u8 x = (u8)a; for (int b = 7; b >= 0 && x; b--) if (x >> b & 1) { if (piv[b]) x ^= piv[b]; else { piv[b] = x; basis[d++] = (u8)a; x = 0; } }
    }
    return d;
}
static inline bool inmask(const u64* mask, u8 x) { return mask[x >> 6] >> (x & 63) & 1; }

// ---------- modes ----------
static double gauss_binom(int n, int k, double q) {
    if (k < 0 || k > n) return 0;
    double r = 1; for (int i = 0; i < k; i++) r = r * (pow(q, n - i) - 1) / (pow(q, i + 1) - 1);
    return std::floor(r + 0.5);
}
static std::vector<u32> canon_span(const u32* b, int n) {
    Ech28 e; for (int i = 0; i < n; i++) e.add(b[i]);
    for (int x = 27; x >= 0; x--) if (e.piv[x]) for (int y = x + 1; y < 28; y++) if (e.piv[y] && (e.piv[y] >> x & 1)) e.piv[y] ^= e.piv[x];
    std::vector<u32> v; for (int x = 27; x >= 0; x--) if (e.piv[x]) v.push_back(e.piv[x]);
    return v;
}

static int mode_enumtest(const Decomp& D) {
    int fails = 0; Rng R(7);
    // one block at a time: k-dim F_q-subspaces, free or containing a random vector u
    for (int b = 0; b < (int)D.B.size(); b++) {
        const Block& bl = D.B[b];
        for (int with = 0; with < 2; with++) for (int k = with; k <= bl.m; k++) {
            double expect = with ? gauss_binom(bl.m - 1, k - 1, bl.q) : gauss_binom(bl.m, k, bl.q);
            if (expect > 3e5 || expect < 1) continue;
            Search S; S.D = &D; S.nord = 1; S.ord[0] = b;
            FormList st0; u8 u[28] = {0}; int h = -1;
            if (with) { while (h < 0) for (int i = 0; i < bl.m; i++) { u[i] = (u8)R.below(bl.q); if (u[i]) h = i; }
                        for (int s = 0; s < bl.d; s++) { if (s) for (int i = 0; i < bl.m; i++) u[i] = bl.mulx[u[i]]; st0.add(bl.form(u)); } }
            for (int i = 0; i < bl.m; i++) if (i != h) S.pos[b].push_back(i);
            S.extra[b] = k - with;
            std::set<std::vector<u32>> seen; u64 leaves = 0; bool badinv = false;
            auto lf = [&](const FormList& s) {
                leaves++; seen.insert(canon_span(s.basis, s.n));
                if (s.n != k * bl.d) badinv = true;
                for (int i = 0; i < s.n && !badinv; i++) { Ech28 e; for (int j = 0; j < s.n; j++) e.add(s.basis[j]); if (e.n != s.n || e.reduce(D.pw[1](s.basis[i]))) badinv = true; }
            };
            if (S.extra[b] == 0) lf(st0); else dfs(S, 0, 0, -1, st0, lf);
            bool ok = leaves == (u64)expect && seen.size() == leaves && !badinv;
            fails += !ok;
            printf("%s part %d (F_%d^%d) k=%d%s: leaves %llu distinct %zu expected %.0f\n", ok ? "OK  " : "FAIL", b, bl.q, bl.m, k,
                   with ? " containing u" : "", leaves, seen.size(), expect);
        }
    }
    // two blocks at once (product), free
    if (D.B.size() >= 2) {
        int b0 = 0, b1 = 1; const Block& x = D.B[b0]; const Block& y = D.B[b1];
        int k0 = std::min(2, x.m), k1 = 1;
        double expect = gauss_binom(x.m, k0, x.q) * gauss_binom(y.m, k1, y.q);
        if (expect <= 3e6) {
            Search S; S.D = &D; S.nord = 2; S.ord[0] = b1; S.ord[1] = b0;
            for (int i = 0; i < x.m; i++) S.pos[b0].push_back(i);
            for (int i = 0; i < y.m; i++) S.pos[b1].push_back(i);
            S.extra[b0] = k0; S.extra[b1] = k1;
            std::set<std::vector<u32>> seen; u64 leaves = 0;
            auto lf = [&](const FormList& s) { leaves++; seen.insert(canon_span(s.basis, s.n)); };
            FormList st0; dfs(S, 0, 0, -1, st0, lf);
            bool ok = leaves == (u64)expect && seen.size() == leaves;
            fails += !ok;
            printf("%s parts 0+1, k=(%d,%d): leaves %llu distinct %zu expected %.0f\n", ok ? "OK  " : "FAIL", k0, k1, leaves, seen.size(), expect);
        }
    }
    printf("enumtest: %d failures\n", fails);
    return fails ? 1 : 0;
}

static void print_stats(FILE* f, const std::map<std::string, SplitStat>& stats) {
    for (auto& kv : stats) {
        fprintf(f, "  split %-40s B1s %llu tries %llu leaves %llu  adds by dim:", kv.first.c_str(), kv.second.b1, kv.second.tries, kv.second.leaves);
        for (int i = 1; i <= 8; i++) fprintf(f, " %llu", kv.second.adds[i]);
        fprintf(f, "\n");
    }
}
static void print_summary(FILE* f, Global& G) {
    fprintf(f, "APN leaves %llu, distinct APN spaces %zu, bugs %llu\n", G.leaves, G.found.size(), G.bugs);
    for (auto& kv : G.dist) fprintf(f, "  ranks %s : %llu\n", kv.first.c_str(), kv.second);
    int n5 = 0; for (auto& x : G.found) if (x.cnt[2] >= 5) n5++;
    fprintf(f, "with >= 5 rank-4 components: %d\n", n5);
}

static bool load_space(const char* spec, u32* basis) {
    std::string s = spec;
    if (s == "x3") { x3_basis(basis); return true; }
    size_t c = s.rfind(':'); if (c == std::string::npos) return false;
    std::vector<Found> v; if (!load_seeds(s.substr(0, c).c_str(), v)) return false;
    int line = atoi(s.c_str() + c + 1); if (line < 0 || line >= (int)v.size()) return false;
    memcpy(basis, v[line].basis, 8 * sizeof(u32)); return true;
}

static int mode_ext(const Decomp& D, const char* spec, int g0, const char* outp) {
    u32 basis[8]; if (!load_space(spec, basis)) { fprintf(stderr, "can't load %s\n", spec); return 2; }
    canon8(basis);
    printf("space %s invariant under A: %s\n", spec, invariant(D, basis) ? "yes" : "NO");
    u32 B1 = 0; int g = g0;
    for (; g < 256; g++) { u32 c = 0; for (int i = 0; i < 8; i++) if (g >> i & 1) c ^= basis[i]; if (D.pw[1](c) != c) { B1 = c; break; } }
    if (!B1) { fprintf(stderr, "no component outside Fix from g=%d on\n", g0); return 2; }
    printf("B1 = component %d = %07x, rank %d\n", g, B1, rank_mat(coef_to_mat(B1)));
    { Ech28 e; for (int i = 0; i < D.p; i++) e.add(D.pw[i](B1)); printf("orbit span dim %d\n", e.n); }
    Prep P; bool okw = prepare(D, B1, P);
    printf("W = F_2[A]B1: dim %d, kernel condition %s; B1 meets parts:", P.wdim, okw ? "passes" : "FAILS");
    for (int b = 0; b < (int)D.B.size(); b++) printf(" %d", (int)P.pres[b]);
    printf("\n");
    if (!okw) return 1;
    Global G; G.out = fopen(outp, "w");
    std::map<std::string, SplitStat> stats;
    auto leaf = [&](const KS& s) { record(D, G, s, B1); };
    auto t0 = std::chrono::steady_clock::now();
    extend(D, P, leaf, stats);
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    print_stats(stdout, stats); print_summary(stdout, G);
    std::array<u32, 8> key; memcpy(key.data(), basis, sizeof basis);
    printf("the space itself found: %s  (%.2f s)\n", G.seen.count(key) ? "YES" : "NO", el);
    for (auto& x : G.found) print_found(stdout, x);
    if (G.out) fclose(G.out);
    return G.seen.count(key) ? 0 : 1;
}

// stage 1: B1 over all rank-4 forms, orbit-minimal, radicals of the orbit pairwise disjoint, W = F_2[A]B1 passing
// the kernel test. The radical test is done twice: as R & A^k R = 0, and on the radicals of the orbit forms
// themselves (a check of rad(A.B) = A rad(B) and of the action); disagreements are counted.
struct Stage1 { u64 r4 = 0, fixed = 0, min = 0, rad = 0, w = 0, radmismatch = 0; std::vector<u32> cand; };
static void stage1(const Decomp& D, Stage1& S1) {
    const int p = D.p;
    std::atomic<u64> n_r4(0), n_fixed(0), n_min(0), n_rad(0), n_w(0), n_mis(0);
    std::vector<u32>& cand = S1.cand; std::mutex cm;
    const u32 CH = 1u << 16;
    #pragma omp parallel
    {
        std::vector<u32> loc; u64 r4 = 0, fx = 0, mn = 0, rd = 0, wp = 0, mis = 0;
        #pragma omp for schedule(dynamic, 16)
        for (int ch = 0; ch < (int)((1u << 28) / CH); ch++) {
            for (u32 c = (u32)ch * CH; c < (u32)(ch + 1) * CH; c++) {
                if (!c) continue;
                u64 m = c2m(c);
                if (rank_mat(m) != 4) continue;
                r4++;
                if (D.pw[1](c) == c) { fx++; continue; }
                bool isMin = true; for (int i = 1; i < p && isMin; i++) if (D.pw[i](c) < c) isMin = false;
                if (!isMin) continue;
                mn++;
                u64 mask[4]; u8 rb[8]; radical(m, mask, rb);
                bool disj = true;
                for (int k = 1; k <= (p - 1) / 2 && disj; k++) {
                    // R & A^k R = 0 iff no nonzero combination of the images lies in R
                    u8 img[4]; for (int i = 0; i < 4; i++) img[i] = mv(D.Apow[k], rb[i]);
                    for (unsigned s = 1; s < 16 && disj; s++) { u8 x = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) x ^= img[i]; if (inmask(mask, x)) disj = false; }
                }
                {
                    u64 rm[16][4]; bool d2 = true;
                    for (int i = 0; i < p && i < 16; i++) { u8 bb[8]; radical(c2m(D.pw[i](c)), rm[i], bb); }
                    for (int i = 0; i < p && i < 16 && d2; i++) for (int j = i + 1; j < p && j < 16 && d2; j++)
                        if ((rm[i][0] & rm[j][0] & ~1ull) | (rm[i][1] & rm[j][1]) | (rm[i][2] & rm[j][2]) | (rm[i][3] & rm[j][3])) d2 = false;
                    if (p <= 16 && d2 != disj) mis++;
                }
                if (!disj) continue;
                rd++;
                Prep P; bool okw = prepare(D, c, P);
                { // the same test on the plain orbit span, with apn8's max_kernel
                    Ech28 e; u64 ms[8]; int n = 0; bool fit = true;
                    for (int i = 0; i < p; i++) { u32 x = D.pw[i](c); if (e.add(x)) { if (n == 8) { fit = false; break; } ms[n++] = c2m(x); } }
                    bool ok2 = fit && max_kernel(ms, n) <= 9 - n;
                    if (ok2 != okw || (okw && n != P.wdim)) mis++;
                }
                if (!okw) continue;
                wp++; loc.push_back(c);
            }
        }
        n_r4 += r4; n_fixed += fx; n_min += mn; n_rad += rd; n_w += wp; n_mis += mis;
        std::lock_guard<std::mutex> g(cm); cand.insert(cand.end(), loc.begin(), loc.end());
    }
    std::sort(cand.begin(), cand.end());
    S1.r4 = n_r4; S1.fixed = n_fixed; S1.min = n_min; S1.rad = n_rad; S1.w = n_w; S1.radmismatch = n_mis;
}

// the pruned search against no pruning: every leaf of every split for every M-th B1, APN tested by apn8's max_kernel
static int mode_brute(const Decomp& D, int th, int M) {
    omp_set_num_threads(th);
    Stage1 S1; stage1(D, S1);
    printf("stage 1: %zu B1s (radical-test mismatches %llu); taking every %d-th\n", S1.cand.size(), S1.radmismatch, M);
    std::vector<u32> todo; for (size_t i = 0; i < S1.cand.size(); i += M) todo.push_back(S1.cand[i]);
    std::atomic<u64> leaves(0), apn_brute(0), apn_pruned(0), mism(0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)todo.size(); i++) {
        Prep P; if (!prepare(D, todo[i], P)) continue;
        std::vector<std::vector<int>> sp; std::vector<int> k(D.B.size(), 0);
        all_splits(D, P.sdim, 0, 8, k, sp);
        u64 lv = 0, ab = 0, ap = 0;
        for (auto& ks : sp) {
            Search S; search_setup(D, P, ks, S);
            FormList f0; for (int j = 0; j < P.w.n; j++) f0.add(P.w.basis[j]);
            auto lf = [&](const FormList& s) { lv++; u64 m[8]; for (int j = 0; j < 8; j++) m[j] = coef_to_mat(s.basis[j]); if (s.n == 8 && max_kernel(m, 8) == 1) ab++; };
            Search S2 = S;
            dfs(S, 0, 0, -1, f0, lf);
            auto lp = [&](const KS&) { ap++; };
            dfs(S2, 0, 0, -1, P.w, lp);
        }
        leaves += lv; apn_brute += ab; apn_pruned += ap; if (ab != ap) mism++;
    }
    printf("B1s %zu, unpruned leaves %llu, APN by brute force %llu, by the pruned search %llu, B1s disagreeing %llu\n",
           todo.size(), (u64)leaves, (u64)apn_brute, (u64)apn_pruned, (u64)mism);
    return mism ? 1 : 0;
}

static int mode_run(const Decomp& D, const char* outp, int th, int M, int rr) {
    omp_set_num_threads(th);
    auto t0 = std::chrono::steady_clock::now();
    auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    const int p = D.p;
    Stage1 S1; stage1(D, S1);
    const std::vector<u32>& cand = S1.cand;
    std::vector<u32> todo; for (size_t i = 0; i < cand.size(); i++) if ((int)(i % M) == rr) todo.push_back(cand[i]);
    char hdr[600];
    snprintf(hdr, sizeof hdr, "# symp p=%d A=%02x%02x%02x%02x%02x%02x%02x%02x\n# stage 1 (%.1f s): rank-4 %llu, fixed by A %llu, orbit-min %llu, radicals pairwise disjoint %llu (second test disagrees %llu), W passes kernel test %llu; this run takes %zu (1/%d, residue %d)\n",
             p, D.A.c[0], D.A.c[1], D.A.c[2], D.A.c[3], D.A.c[4], D.A.c[5], D.A.c[6], D.A.c[7], el(),
             S1.r4, S1.fixed, S1.min, S1.rad, S1.radmismatch, S1.w, todo.size(), M, rr);
    printf("%s", hdr); fflush(stdout);
    Global G; G.out = fopen(outp, "w"); if (G.out) { fputs(hdr, G.out); fflush(G.out); }
    // stage 2: extensions
    std::map<std::string, SplitStat> stats; std::mutex sm;
    std::map<int, u64> wdims; // W dim and whether f != 0: key = 10*dim + (f != 0)
    std::atomic<u64> done(0); double next = 60;
    #pragma omp parallel
    {
        std::map<std::string, SplitStat> loc; std::map<int, u64> wl;
        #pragma omp for schedule(dynamic, 4)
        for (int i = 0; i < (int)todo.size(); i++) {
            u32 B1 = todo[i]; Prep P;
            if (!prepare(D, B1, P)) continue;
            wl[10 * P.wdim + (P.pres[0] ? 1 : 0)]++;
            auto leaf = [&](const KS& s) { record(D, G, s, B1); };
            extend(D, P, leaf, loc);
            u64 dn = ++done;
            if (omp_get_thread_num() == 0 && el() > next) {
                next = el() + 60;
                std::lock_guard<std::mutex> g(G.mu);
                printf("[%.0f s] %llu / %zu B1s, APN spaces %zu\n", el(), dn, todo.size(), G.found.size()); fflush(stdout);
            }
        }
        std::lock_guard<std::mutex> g(sm);
        for (auto& kv : loc) { SplitStat& t = stats[kv.first]; t.b1 += kv.second.b1; t.tries += kv.second.tries; t.leaves += kv.second.leaves; for (int j = 0; j < 29; j++) t.adds[j] += kv.second.adds[j]; }
        for (auto& kv : wl) wdims[kv.first] += kv.second;
    }
    FILE* fs[2] = {stdout, G.out};
    for (FILE* f : fs) {
        if (!f) continue;
        fprintf(f, "# stage 2 (%.1f s total): W types (dim, f!=0):", el());
        for (auto& kv : wdims) fprintf(f, " (%d,%d):%llu", kv.first / 10, kv.first % 10, kv.second);
        fprintf(f, "\n");
        print_stats(f, stats); print_summary(f, G);
    }
    if (G.out) fclose(G.out);
    return 0;
}

// order-p elements of Aut(S8): A is fixed by where it sends three rank-4 radicals R1, R2, R3 (R3 = graph of
// T: R1 -> R2) and by A|R1; then A|R2 = T' A|R1 T^-1.
static int mode_autfind(const char* file, int line, int p) {
    std::vector<Found> v; if (!load_seeds(file, v) || line < 0 || line >= (int)v.size()) { fprintf(stderr, "can't load\n"); return 2; }
    const u32* S = v[line].basis; u64 Sm[8]; for (int i = 0; i < 8; i++) Sm[i] = coef_to_mat(S[i]);
    Ech28 E; for (int i = 0; i < 8; i++) E.add(S[i]);
    std::vector<std::array<u8, 4>> rads; std::vector<std::array<u64, 4>> masks;
    u32 cur = 0;
    for (unsigned g = 1; g < 256; g++) {
        cur ^= S[_tzcnt_u32(g)];
        u64 m = coef_to_mat(cur); if (rank_mat(m) != 4) continue;
        std::array<u8, 4> b; std::array<u64, 4> mk; u8 bb[8]; radical(m, mk.data(), bb); for (int i = 0; i < 4; i++) b[i] = bb[i];
        rads.push_back(b); masks.push_back(mk);
    }
    int k = (int)rads.size(); printf("ranks %s, %d rank-4 radicals\n", dist_name(v[line].cnt).c_str(), k);
    if (k < 3) return 2;
    // y in R_j with x + y in R_l
    auto graph = [&](int i, int j, int l, u8 x) -> u8 {
        for (unsigned s = 0; s < 16; s++) { u8 y = 0; for (int t = 0; t < 4; t++) if (s >> t & 1) y ^= rads[j][t]; if (inmask(masks[l].data(), (u8)(x ^ y))) return y; }
        return 0;
    };
    u8 r1[4], T[4]; for (int i = 0; i < 4; i++) { r1[i] = rads[0][i]; T[i] = graph(0, 1, 2, r1[i]); }
    M8 Src; for (int i = 0; i < 4; i++) { Src.c[i] = r1[i]; Src.c[4 + i] = T[i]; }
    M8 Si; if (!inv8(Src, Si)) { fprintf(stderr, "R1 + R2 != V\n"); return 2; }
    std::vector<u16> gl4; // GL(4,2) as 4 nibble columns
    for (unsigned G = 0; G < 65536; G++) { u8 c[4]; for (int i = 0; i < 4; i++) c[i] = (u8)(G >> (4 * i) & 15); if (rank8(c, 4) == 4) gl4.push_back((u16)G); }
    std::set<std::array<u8, 8>> outs;
    for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) for (int c = 0; c < k; c++) {
        if (a == b || a == c || b == c) continue;
        u8 s1[4], Tp[4]; for (int i = 0; i < 4; i++) { s1[i] = rads[a][i]; Tp[i] = graph(a, b, c, s1[i]); }
        for (u16 G : gl4) {
            M8 Img;
            for (int i = 0; i < 4; i++) {
                u8 u = 0, tu = 0; for (int j = 0; j < 4; j++) if (G >> (4 * i + j) & 1) { u ^= s1[j]; tu ^= Tp[j]; }
                Img.c[i] = u; Img.c[4 + i] = tu;
            }
            M8 A; for (int i = 0; i < 8; i++) A.c[i] = mv(Img, Si.c[i]);
            if (order8(A) != p) continue;
            M8 Ai; inv8(A, Ai);
            bool inv = true; for (int i = 0; i < 8 && inv; i++) if (E.reduce(act_form(Sm[i], Ai))) inv = false;
            if (!inv) continue;
            std::array<u8, 8> key; memcpy(key.data(), A.c, 8); outs.insert(key);
        }
    }
    printf("%zu elements of order %d\n", outs.size(), p);
    for (auto& x : outs) { M8 A; memcpy(A.c, x.data(), 8); int fixdim = 0; for (unsigned y = 1; y < 256; y++) if (mv(A, (u8)y) == y) fixdim++; print_m8(stdout, A); printf("  fixed vectors %d\n", fixdim); }
    return 0;
}

// ---------- order 3, class c3k<k>: V = W + T, W = F_4^k on bits 0..2k-1 (A = omega), T = F_2^(8-2k) fixed ----------
static u8 f4mul(u8 a, u8 b) { return (u8)pmulmod(a, b, 7, 2); }

// lemma checks for c3k<k>: fixed forms vanish on T x W; the other part vanishes on T x T and on (x, omega x), x in W;
// the smallest subspace of Herm(W) = Fix & Alt(W) whose quadrics x -> B(x, omega x) have no common zero on W \ 0.
static int mode_lemmas3(const Decomp& D, int k) {
    const int nw = 2 * k; const u8 Wm = (u8)((1u << nw) - 1);
    const Block& F = D.B[0];
    int bad = 0;
    for (u32 c : F.elem) { u64 m = c2m(c); for (int w = 0; w < nw; w++) for (int t = nw; t < 8; t++) if (bform(m, (u8)(1 << t), (u8)(1 << w))) bad++; }
    printf("fixed forms nonzero on T x W: %d (want 0)\n", bad);
    bad = 0;
    for (size_t b = 1; b < D.B.size(); b++) for (u32 c : D.B[b].elem) {
        u64 m = c2m(c);
        for (int s = nw; s < 8; s++) for (int t = nw; t < 8; t++) if (bform(m, (u8)(1 << s), (u8)(1 << t))) bad++;
        for (unsigned x = 1; x <= Wm; x++) if (bform(m, (u8)x, mv(D.A, (u8)x))) bad++;
    }
    printf("non-fixed forms nonzero on T x T or on some (x, omega x): %d (want 0)\n", bad);
    // Herm(W): fixed forms supported on W x W
    std::vector<u32> H; { u32 Tmask = 0; for (int q = 0; q < 28; q++) if (KI[q] >= nw || KJ[q] >= nw) Tmask |= 1u << q;
        Ech28 e; for (u32 g = 1; g < (1u << F.m); g++) { u32 c = 0; for (int i = 0; i < F.m; i++) if (g >> i & 1) c ^= F.elem[i]; if (!(c & Tmask) && e.add(c)) H.push_back(c); } }
    int nt = 8 - nw;
    printf("dim Fix %d = dim Herm(W) %zu + dim Alt(T) %d; other part F_4-dim %d\n", F.m, H.size(), nt * (nt - 1) / 2, D.B.size() > 1 ? D.B[1].m : 0);
    // zero sets Z(B) over x in W \ 0, as bitmasks
    int nh = (int)H.size(); std::vector<std::array<u64, 4>> Z(1u << nh);
    for (u32 g = 0; g < (1u << nh); g++) {
        u32 c = 0; for (int i = 0; i < nh; i++) if (g >> i & 1) c ^= H[i];
        u64 m = c2m(c); std::array<u64, 4> z = {0, 0, 0, 0};
        for (unsigned x = 1; x <= Wm; x++) if (!bform(m, (u8)x, mv(D.A, (u8)x))) z[x >> 6] |= 1ull << (x & 63);
        Z[g] = z;
    }
    auto nonempty = [](const std::array<u64, 4>& z) { return (z[0] | z[1] | z[2] | z[3]) != 0; };
    int jmin = -1;
    for (int j = 1; j <= nh && jmin < 0; j++) {
        // j-dim subspaces of Herm as j independent elements (by index tuples, increasing): common zero = AND of zero sets
        std::vector<u32> idx(j); bool found = false;
        std::function<void(int, u32, std::array<u64, 4>, int)> rec = [&](int depth, u32 start, std::array<u64, 4> z, int dummy) {
            if (found) return;
            if (depth == j) { if (!nonempty(z)) found = true; return; }
            for (u32 g = start; g < (1u << nh) && !found; g++) {
                if (!g) continue;
                // independence: g not in the span of the chosen ones
                bool dep = false; for (u32 s = 1; s < (1u << depth) && !dep; s++) { u32 c = 0; for (int i = 0; i < depth; i++) if (s >> i & 1) c ^= idx[i]; if (c == g) dep = true; }
                if (dep) continue;
                idx[depth] = g; std::array<u64, 4> z2; for (int i = 0; i < 4; i++) z2[i] = z[i] & Z[g][i];
                if (depth + 1 < j || !nonempty(z2)) rec(depth + 1, g + 1, z2, dummy);
                if (depth + 1 == j && !nonempty(z2)) found = true;
            }
        };
        std::array<u64, 4> all; for (int i = 0; i < 4; i++) all[i] = ~0ull; all[0] &= ~1ull;
        rec(0, 1, all, 0);
        printf("  a %d-dim subspace of Herm(W) with no common zero: %s\n", j, found ? "exists" : "none");
        if (found) jmin = j;
    }
    int bmin = nt > 0 ? k : 0, amin = std::max(jmin, nt > 0 ? nt - 1 : 0);
    printf("so a >= %d (Herm projection >= %d, Alt(T) projection >= %d), b >= %d; splits a + 2b = 8 left:", amin, jmin, nt > 0 ? nt - 1 : 0, bmin);
    int any = 0; for (int b = bmin; 2 * b <= 8; b++) { int a = 8 - 2 * b; if (a >= amin && a <= F.m && (D.B.size() < 2 || b <= D.B[1].m)) { printf(" (%d,%d)", a, b); any++; } }
    printf("%s\n", any ? "" : " none");
    return 0;
}

// C(A) for c3k2, with the Frobenius of W (which conjugates A to A^2, so it preserves <A>-invariance)
static std::vector<M8> c3k2_gens() {
    std::vector<M8> G;
    auto f4lin = [&](u8 al, u8 be, u8 ga, u8 de) {
        M8 M = ident8();
        for (int i = 0; i < 4; i++) {
            u8 x = (u8)(1 << i), z1 = x & 3, z2 = x >> 2 & 3;
            u8 y1 = f4mul(al, z1) ^ f4mul(be, z2), y2 = f4mul(ga, z1) ^ f4mul(de, z2);
            M.c[i] = (u8)(y1 | y2 << 2);
        }
        return M;
    };
    G.push_back(f4lin(2, 0, 0, 1)); G.push_back(f4lin(1, 1, 0, 1)); G.push_back(f4lin(1, 2, 0, 1)); G.push_back(f4lin(0, 1, 1, 0));
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) if (i != j) { M8 M = ident8(); M.c[4 + i] = (u8)((1 << (4 + i)) | (1 << (4 + j))); G.push_back(M); }
    M8 Fr = ident8(); Fr.c[1] = 3; Fr.c[3] = 12; G.push_back(Fr); // z -> z^2 on each F_4 coordinate: omega -> 1 + omega
    return G;
}

// 3-cycle case for c3k2: B1 over all rank-4 forms with R & AR = 0 and W = F_2[A]B1 passing, taken up to the group
// generated by c3k2_gens (orbits by BFS), then the extension with --split (default 4,2).
static int mode_cycle3(const Decomp& D, const char* outp, int th, int maxreps) {
    omp_set_num_threads(th);
    auto t0 = std::chrono::steady_clock::now();
    auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    std::vector<M8> gens = c3k2_gens();
    M8 A2 = mm(D.A, D.A);
    for (auto& g : gens) { M8 gi; inv8(g, gi); M8 c = mm(g, mm(D.A, gi)); if (!eq8(c, D.A) && !eq8(c, A2)) { fprintf(stderr, "generator doesn't normalise <A>\n"); return 2; } }
    std::vector<L28> gl; for (auto& g : gens) gl.push_back(alt_action(g));
    std::vector<u32> X; std::mutex xm; std::atomic<u64> n_r4(0), n_rad(0);
    const u32 CH = 1u << 16;
    #pragma omp parallel
    {
        std::vector<u32> loc; u64 r4 = 0, rd = 0;
        #pragma omp for schedule(dynamic, 16)
        for (int ch = 0; ch < (int)((1u << 28) / CH); ch++) for (u32 c = (u32)ch * CH; c < (u32)(ch + 1) * CH; c++) {
            if (!c) continue; u64 m = c2m(c); if (rank_mat(m) != 4) continue;
            r4++;
            u64 mask[4]; u8 rb[8]; radical(m, mask, rb);
            bool disj = true; u8 img[4]; for (int i = 0; i < 4; i++) img[i] = mv(D.A, rb[i]);
            for (unsigned s = 1; s < 16 && disj; s++) { u8 x = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) x ^= img[i]; if (inmask(mask, x)) disj = false; }
            if (!disj) continue;
            rd++;
            Prep P; if (!prepare(D, c, P)) continue;
            loc.push_back(c);
        }
        n_r4 += r4; n_rad += rd;
        std::lock_guard<std::mutex> g(xm); X.insert(X.end(), loc.begin(), loc.end());
    }
    std::sort(X.begin(), X.end());
    std::vector<u64> inX(1u << 22, 0), seen(1u << 22, 0);
    for (u32 c : X) inX[c >> 6] |= 1ull << (c & 63);
    std::vector<u32> reps; std::vector<u64> sizes; std::map<u64, u64> sizehist; u64 escaped = 0;
    std::vector<u32> q;
    for (u32 c : X) {
        if (seen[c >> 6] >> (c & 63) & 1) continue;
        reps.push_back(c); q.clear(); q.push_back(c); seen[c >> 6] |= 1ull << (c & 63);
        for (size_t h = 0; h < q.size(); h++) for (auto& L : gl) {
            u32 y = L(q[h]);
            if (seen[y >> 6] >> (y & 63) & 1) continue;
            if (!(inX[y >> 6] >> (y & 63) & 1)) { escaped++; continue; }
            seen[y >> 6] |= 1ull << (y & 63); q.push_back(y);
        }
        sizes.push_back(q.size()); sizehist[q.size()]++;
    }
    FILE* out = fopen(outp, "w");
    FILE* fs[2] = {stdout, out};
    for (FILE* f : fs) {
        if (!f) continue;
        fprintf(f, "# cycle3 c3k2 (%.1f s): rank-4 %llu, R & AR = 0: %llu, W passes: %zu, orbits under C(A)xFrob: %zu (images leaving the set: %llu, want 0)\n#  orbit sizes:",
                el(), (u64)n_r4, (u64)n_rad, X.size(), reps.size(), escaped);
        for (auto& kv : sizehist) fprintf(f, " %llu^%llu", kv.first, kv.second);
        fprintf(f, "\n"); fflush(f);
    }
    if (maxreps >= 0 && (int)reps.size() > maxreps) reps.resize(maxreps);
    Global G; G.out = out;
    std::map<std::string, SplitStat> stats; std::mutex sm; std::atomic<u64> done(0); double next = 60;
    #pragma omp parallel
    {
        std::map<std::string, SplitStat> loc;
        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < (int)reps.size(); i++) {
            u32 B1 = reps[i]; Prep P; if (!prepare(D, B1, P)) continue;
            auto leaf = [&](const KS& s) { record(D, G, s, B1); };
            extend(D, P, leaf, loc);
            u64 dn = ++done;
            if (omp_get_thread_num() == 0 && el() > next) { next = el() + 60; std::lock_guard<std::mutex> g(G.mu); printf("[%.0f s] %llu / %zu reps, APN %zu\n", el(), dn, reps.size(), G.found.size()); fflush(stdout); }
        }
        std::lock_guard<std::mutex> g(sm);
        for (auto& kv : loc) { SplitStat& t = stats[kv.first]; t.b1 += kv.second.b1; t.tries += kv.second.tries; t.leaves += kv.second.leaves; for (int j = 0; j < 29; j++) t.adds[j] += kv.second.adds[j]; }
    }
    for (FILE* f : fs) { if (!f) continue; fprintf(f, "# extension (%.1f s total), %zu reps\n", el(), reps.size()); print_stats(f, stats); print_summary(f, G); }
    if (out) fclose(out);
    return 0;
}

// the no-3-cycle case for c3k2: five fixed rank-4 forms with radicals L_i + S_i (L_i the five F_4-lines of W, S_i a
// spread of T, matched by a bijection), Q = their span (or any 4-dim Q in Fix containing it), then P.
static int mode_fixed5(const Decomp& D, const char* outp, int th) {
    omp_set_num_threads(th);
    const Block& F = D.B[0];
    // fixed rank-4 forms by radical type
    std::map<std::string, int> types; std::vector<std::pair<u32, std::array<u64, 4>>> r4;
    for (u32 g = 1; g < (1u << F.m); g++) {
        u32 c = 0; for (int i = 0; i < F.m; i++) if (g >> i & 1) c ^= F.elem[i];
        u64 m = c2m(c); if (rank_mat(m) != 4) continue;
        std::array<u64, 4> mk; u8 rb[8]; radical(m, mk.data(), rb);
        int dw = 0, dt = 0; for (unsigned x = 1; x < 16; x++) if (inmask(mk.data(), (u8)x)) dw++; for (unsigned t = 1; t < 16; t++) if (inmask(mk.data(), (u8)(t << 4))) dt++;
        dw = 31 - (int)__lzcnt((unsigned)dw + 1); dt = 31 - (int)__lzcnt((unsigned)dt + 1); // 2^dim - 1 nonzero elements
        char b[32]; snprintf(b, sizeof b, "W%d+T%d", dw, dt); types[b]++;
        r4.push_back({c, mk});
    }
    printf("fixed rank-4 forms by radical type (dim R&W + dim R&T):"); for (auto& kv : types) printf(" %s:%d", kv.first.c_str(), kv.second); printf("\n");
    // largest pairwise complementary family of their radicals (distinct radicals), and largest containing a non-(2,2) one
    std::vector<std::array<u64, 4>> rads; std::vector<bool> is22;
    { std::set<std::array<u64, 4>> s; for (auto& x : r4) if (s.insert(x.second).second) rads.push_back(x.second); }
    for (auto& mk : rads) { int dw = 0, dt = 0; for (unsigned x = 1; x < 16; x++) dw += inmask(mk.data(), (u8)x); for (unsigned t = 1; t < 16; t++) dt += inmask(mk.data(), (u8)(t << 4)); is22.push_back(dw == 3 && dt == 3); }
    auto disj = [](const std::array<u64, 4>& a, const std::array<u64, 4>& b) { return !((a[0] & b[0] & ~1ull) | (a[1] & b[1]) | (a[2] & b[2]) | (a[3] & b[3])); };
    int nr = (int)rads.size(), best = 0, bestNon = 0; std::vector<int> cl;
    std::function<void(int)> clique = [&](int from) {
        int sz = (int)cl.size(); bool non = false; for (int i : cl) if (!is22[i]) non = true;
        best = std::max(best, sz); if (non) bestNon = std::max(bestNon, sz);
        for (int j = from; j < nr; j++) { bool ok = true; for (int i : cl) if (!disj(rads[i], rads[j])) { ok = false; break; } if (ok) { cl.push_back(j); clique(j + 1); cl.pop_back(); } }
    };
    clique(0);
    printf("distinct radicals %d; largest pairwise complementary family %d, largest with a radical not of type (F_4-line + 2-dim of T) %d\n", nr, best, bestNon);
    // lines of W, 2-dim subspaces of T, spreads
    std::vector<u16> lines; for (unsigned v = 1; v < 16; v++) { u8 w = mv(D.A, (u8)v); u16 L = (u16)(1 | 1u << v | 1u << w | 1u << (v ^ w)); if (std::find(lines.begin(), lines.end(), L) == lines.end()) lines.push_back(L); }
    std::vector<u16> planes; for (unsigned a = 1; a < 16; a++) for (unsigned b = a + 1; b < 16; b++) { u16 S = (u16)(1 | 1u << a | 1u << b | 1u << (a ^ b)); if (std::find(planes.begin(), planes.end(), S) == planes.end()) planes.push_back(S); }
    std::vector<std::array<int, 5>> spreads; std::array<int, 5> cur;
    std::function<void(int, int)> sp = [&](int d, int from) {
        if (d == 5) { spreads.push_back(cur); return; }
        for (int j = from; j < (int)planes.size(); j++) { bool ok = true; for (int i = 0; i < d; i++) if ((planes[cur[i]] & planes[j]) != 1) ok = false; if (ok) { cur[d] = j; sp(d + 1, j + 1); } }
    };
    sp(0, 0);
    printf("F_4-lines %zu, 2-dim subspaces of T %zu, spreads %zu\n", lines.size(), planes.size(), spreads.size());
    auto radmask = [&](u16 L, u16 S) { std::array<u64, 4> mk = {0, 0, 0, 0}; for (unsigned w = 0; w < 16; w++) if (L >> w & 1) for (unsigned t = 0; t < 16; t++) if (S >> t & 1) { u8 x = (u8)(w | t << 4); mk[x >> 6] |= 1ull << (x & 63); } return mk; };
    std::vector<std::vector<std::vector<u32>>> form(lines.size(), std::vector<std::vector<u32>>(planes.size()));
    int multi = 0;
    for (size_t l = 0; l < lines.size(); l++) for (size_t s = 0; s < planes.size(); s++) {
        auto mk = radmask(lines[l], planes[s]);
        for (auto& x : r4) if (x.second == mk) form[l][s].push_back(x.first);
        if (form[l][s].size() != 1) multi++;
    }
    printf("(line, plane) pairs whose fixed rank-4 form isn't unique: %d (want 0)\n", multi);
    // configurations
    std::set<std::vector<u32>> starts; u64 configs = 0, pass = 0; std::map<int, u64> spandim;
    int perm[5] = {0, 1, 2, 3, 4};
    for (auto& S : spreads) {
        int p[5]; memcpy(p, perm, sizeof p);
        do {
            configs++;
            KS st; Ech28 e; bool ok = true; u32 bs[5]; int n = 0;
            for (int i = 0; i < 5 && ok; i++) { u32 c = form[i][S[p[i]]][0]; if (e.add(c)) { bs[n++] = c; if (!st.add(c)) ok = false; } }
            spandim[e.n]++;
            if (!ok || e.n > 4) continue;
            pass++; starts.insert(canon_span(bs, n));
        } while (std::next_permutation(p, p + 5));
    }
    printf("configurations %llu, span dims:", configs); for (auto& kv : spandim) printf(" %d:%llu", kv.first, kv.second);
    printf("; passing the kernel test with dim <= 4: %llu, distinct starting spaces %zu\n", pass, starts.size()); fflush(stdout);
    Global G; G.out = fopen(outp, "w");
    std::map<std::string, SplitStat> stats; std::vector<std::vector<u32>> sv(starts.begin(), starts.end()); std::mutex sm;
    auto t0 = std::chrono::steady_clock::now();
    #pragma omp parallel
    {
        std::map<std::string, SplitStat> loc;
        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < (int)sv.size(); i++) {
            Prep P; if (!prepare_sub(D, sv[i].data(), (int)sv[i].size(), P)) continue;
            u32 b0 = sv[i][0];
            auto leaf = [&](const KS& s) { record(D, G, s, b0); };
            extend(D, P, leaf, loc);
        }
        std::lock_guard<std::mutex> g(sm);
        for (auto& kv : loc) { SplitStat& t = stats[kv.first]; t.b1 += kv.second.b1; t.tries += kv.second.tries; t.leaves += kv.second.leaves; for (int j = 0; j < 29; j++) t.adds[j] += kv.second.adds[j]; }
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    FILE* fs[2] = {stdout, G.out};
    for (FILE* f : fs) { if (!f) continue; fprintf(f, "# fixed5 extension %.1f s\n", el); print_stats(f, stats); print_summary(f, G); }
    if (G.out) fclose(G.out);
    return 0;
}

// start from a known invariant space's intersection with the parts in <mask> (bit b = part b) and search the rest
// with the space's own split; reports whether the space comes back.
static int mode_extq(const Decomp& D, const char* spec, int mask, const char* outp) {
    u32 basis[8]; if (!load_space(spec, basis)) { fprintf(stderr, "can't load %s\n", spec); return 2; }
    canon8(basis);
    if (!invariant(D, basis)) { printf("space not invariant under A\n"); return 1; }
    std::vector<u32> start; std::vector<int> ks(D.B.size());
    for (int b = 0; b < (int)D.B.size(); b++) {
        Ech28 e; std::vector<u32> part;
        for (int j = 0; j < 8; j++) { u8 v[28]; D.coords(basis[j], b, v); u32 c = D.B[b].form(v); if (e.add(c)) part.push_back(c); }
        ks[b] = (int)part.size() / D.B[b].d;
        if (mask >> b & 1) start.insert(start.end(), part.begin(), part.end());
    }
    if (g_only.empty()) g_only = ks;
    printf("space's split:"); for (int x : ks) printf(" %d", x); printf("; start dim %zu\n", start.size());
    Prep P; if (!prepare_sub(D, start.data(), (int)start.size(), P)) { printf("start fails\n"); return 1; }
    Global G; G.out = fopen(outp, "w"); std::map<std::string, SplitStat> stats;
    u32 b0 = start.empty() ? 0 : start[0];
    auto leaf = [&](const KS& s) { record(D, G, s, b0); };
    auto t0 = std::chrono::steady_clock::now();
    extend(D, P, leaf, stats);
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    print_stats(stdout, stats); print_summary(stdout, G);
    std::array<u32, 8> key; memcpy(key.data(), basis, sizeof basis);
    printf("the space itself found: %s  (%.2f s)\n", G.seen.count(key) ? "YES" : "NO", el);
    if (G.out) fclose(G.out);
    return G.seen.count(key) ? 0 : 1;
}

// random starts: <count> random <dim>-dim F_q-subspaces of part <b> (passing the kernel test), each extended
static int mode_randsub(const Decomp& D, int b, int dim, int count, u64 seed, const char* outp, int th) {
    omp_set_num_threads(th);
    Global G; G.out = fopen(outp, "w"); std::map<std::string, SplitStat> stats; std::mutex sm;
    std::atomic<int> tried(0), started(0);
    #pragma omp parallel
    {
        Rng R(seed * 7919 + omp_get_thread_num() + 1); std::map<std::string, SplitStat> loc;
        #pragma omp for schedule(dynamic, 1)
        for (int it = 0; it < count; it++) {
            const Block& bl = D.B[b];
            for (int att = 0; att < 1000; att++) {
                tried++;
                std::vector<u32> st; Ech28 e; bool ok = true;
                for (int r = 0; r < dim && ok; r++) {
                    u8 v[28]; for (int i = 0; i < bl.m; i++) v[i] = (u8)R.below(bl.q);
                    for (int s = 0; s < bl.d; s++) { if (s) for (int i = 0; i < bl.m; i++) v[i] = bl.mulx[v[i]]; u32 c = bl.form(v); if (!e.add(c)) { ok = false; break; } st.push_back(c); }
                }
                if (!ok) continue;
                Prep P; if (!prepare_sub(D, st.data(), (int)st.size(), P)) continue;
                started++;
                auto leaf = [&](const KS& s) { record(D, G, s, st[0]); };
                extend(D, P, leaf, loc);
                break;
            }
        }
        std::lock_guard<std::mutex> g(sm);
        for (auto& kv : loc) { SplitStat& t = stats[kv.first]; t.b1 += kv.second.b1; t.tries += kv.second.tries; t.leaves += kv.second.leaves; for (int j = 0; j < 29; j++) t.adds[j] += kv.second.adds[j]; }
    }
    printf("random starts: %d tried, %d passed the kernel test\n", (int)tried, (int)started);
    print_stats(stdout, stats); print_summary(stdout, G);
    for (auto& x : G.found) print_found(stdout, x);
    if (G.out) fclose(G.out);
    return 0;
}

int main(int argc, char** argv) {
    init_index(); init_cm();
    for (int i = 1; i + 1 < argc; i++) if (std::string(argv[i]) == "--split") { // strip "--split a,b,.."
        g_only.clear(); for (const char* s = argv[i + 1]; *s; ) { g_only.push_back(atoi(s)); while (*s && *s != ',') s++; if (*s) s++; }
        for (int j = i; j + 2 < argc; j++) argv[j] = argv[j + 2];
        argc -= 2; break;
    }
    if (argc < 2) { fprintf(stderr, "usage: see the header of symp.cpp\n"); return 2; }
    std::string cmd = argv[1];
    if (cmd == "autfind" && argc >= 5) return mode_autfind(argv[2], atoi(argv[3]), atoi(argv[4]));
    if (argc < 3) return 2;
    M8 A; if (!parse_A(argv[2], A)) return 2;
    Decomp D; if (!build_decomp(A, D, true)) return 1;
    if (cmd == "info") {
        bool pres[16] = {false}; int mink[16] = {0}; std::vector<std::vector<int>> sp; std::vector<int> k(D.B.size(), 0);
        all_splits(D, mink, 0, 8, k, sp);
        printf("splits of 8 (any B1):\n"); for (auto& s : sp) printf("  %s\n", split_name(D, s, pres).c_str());
        int fv = 0; for (unsigned y = 1; y < 256; y++) if (mv(A, (u8)y) == y) fv++;
        printf("fixed nonzero vectors of A on V: %d\n", fv);
        return 0;
    }
    if (cmd == "enumtest") return mode_enumtest(D);
    if (cmd == "run" && argc >= 4) {
        int th = argc > 4 ? atoi(argv[4]) : omp_get_max_threads();
        int M = argc > 6 ? atoi(argv[5]) : 1, rr = argc > 6 ? atoi(argv[6]) : 0;
        return mode_run(D, argv[3], th, M, rr);
    }
    if (cmd == "lemmas3") return mode_lemmas3(D, argv[2][3] - '0');
    if (cmd == "cycle3" && argc >= 4) { if (g_only.empty()) g_only = {4, 2}; return mode_cycle3(D, argv[3], argc > 4 ? atoi(argv[4]) : omp_get_max_threads(), argc > 5 ? atoi(argv[5]) : -1); }
    if (cmd == "fixed5" && argc >= 4) { if (g_only.empty()) g_only = {4, 2}; return mode_fixed5(D, argv[3], argc > 4 ? atoi(argv[4]) : omp_get_max_threads()); }
    if (cmd == "extq" && argc >= 6) return mode_extq(D, argv[3], atoi(argv[4]), argv[5]);
    if (cmd == "randsub" && argc >= 8) return mode_randsub(D, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), strtoull(argv[6], 0, 10), argv[7], argc > 8 ? atoi(argv[8]) : omp_get_max_threads());
    if (cmd == "brute" && argc >= 4) return mode_brute(D, argc > 4 ? atoi(argv[4]) : omp_get_max_threads(), atoi(argv[3]));
    if (cmd == "ext" && argc >= 6) return mode_ext(D, argv[3], atoi(argv[4]), argv[5]);
    fprintf(stderr, "bad arguments\n"); return 2;
}
