// sixteen: is there a quadratic APN F on F_2^8 with 16 rank-4 components (amplitude distribution [0^234, 2^5, 4^16],
// the case Bénéteau et al. single out in Open Problem VI.1)? Rowen (Claude Opus 5.5), 2026-10-24.
//
// Reduction (README, "Sixteen"): the 16 rank-4 radicals are pairwise complementary 4-spaces; the 15 points they miss
// form a 4-space R17 (hyperplane count), so with R17 they make a spread, and the 5 rank-6 radicals are lines in R17.
// GL(8,2) is transitive on ordered triples of pairwise complementary 4-spaces, so take R17 = Y = {(0, y)}, one
// rank-4 radical X = graph(0), another graph(I). The spread is then {Y} + {graph(N) : N in Sigma}, Sigma a spread set
// of 16 matrices containing 0 and I (pairwise differences invertible). Enumerate every such Sigma (cliques), and for
// each, search the forms: the component with radical graph(N) is B(x, y) = d(N x1 + x2, N y1 + y2), d one of the 28
// nondegenerate alternating forms on F_2^4. The form on X is normalised to d = J by diag(g, g) (which fixes X, Y,
// graph(I) and permutes the Sigmas, all of which are enumerated).
// Every form in the span W must have rank 8, or rank 4 with radical graph(N) for an N in Sigma (each N once), or
// rank 6 with radical inside Y; all radicals pairwise disjoint (that is the kernel test). APN <=> radicals partition
// V \ 0. W reaches dim 8 with all 16 covered => solution; all 16 covered at dim < 8 => "short" (needs extension;
// reported).
//
// usage: sixteen [maxcliques] [i2 | i3]   prints clique counts (spread sets; fields among them: must be 336), DFS stats,
//                                 and any solution / short case as 'basis' lines (28-bit coefficient vectors).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <atomic>
#include <mutex>
#include <intrin.h>
#include <omp.h>

typedef uint16_t M4;  // row i in bits 4i..4i+3; (N x)_i = parity(row_i & x)
typedef uint64_t u64;

static uint8_t g_inv[65536];
static inline int mv(M4 m, int x) {
    int y = 0;
    for (int i = 0; i < 4; i++) y |= (__popcnt((m >> (4 * i)) & 15 & x) & 1) << i;
    return y;
}
static bool invertible(M4 m) {
    int r[4];
    for (int i = 0; i < 4; i++) r[i] = (m >> (4 * i)) & 15;
    int rank = 0;
    for (int c = 0; c < 4; c++) {
        int p = -1;
        for (int i = rank; i < 4; i++) if (r[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        std::swap(r[p], r[rank]);
        for (int i = 0; i < 4; i++) if (i != rank && (r[i] >> c & 1)) r[i] ^= r[rank];
        rank++;
    }
    return rank == 4;
}
static M4 mmul(M4 a, M4 b) {  // (a b) x = a (b x)
    M4 c = 0;
    for (int j = 0; j < 4; j++) {  // column j of c = a (b e_j)
        int col = mv(a, mv(b, 1 << j));
        for (int i = 0; i < 4; i++) if (col >> i & 1) c |= 1 << (4 * i + j);
    }
    return c;
}
static const M4 ID = 0x8421;

// Alternating forms on F_2^4: 6 bits, pairs (0,1),(0,2),(0,3),(1,2),(1,3),(2,3). As 4 row nibbles.
static int P4[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
static std::vector<M4> g_D28;  // nondegenerate ones, as row nibbles (symmetric with zero diagonal)
// holes h = i of [0^(238-4i), 2^(5i), 4^(17-i)]: 17 - h rank-4 radicals. For h <= 3 every partial spread of 17 - h
// solids in PG(7,2) extends to a spread (README "Fifteen and fourteen"), so the holes are h spread elements: Y; X
// (graph(0)) from h = 2; for h = 3 also graph(N_m), each m = 2..15 tried. The normalised form sits on X (h = 1) or
// graph(I). Rank-6 radicals must lie in the union of the holes (for h = 3 a line can cross all three).
static int g_holes = 1;
static M4 alt4(int bits) {
    M4 m = 0;
    for (int k = 0; k < 6; k++) if (bits >> k & 1) { int i = P4[k][0], j = P4[k][1]; m |= 1 << (4 * i + j); m |= 1 << (4 * j + i); }
    return m;
}

// Forms on V = F_2^8 as 8 row bytes packed in a u64 (row a in bits 8a..8a+7). Point x = x1 | x2 << 4 (x1: X coords).
static inline int bval(u64 B, int x, int y) {
    int v = 0;
    for (int a = 0; a < 8; a++) if (x >> a & 1) v ^= (B >> (8 * a)) & 255;
    return __popcnt(v & y) & 1;
}
static u64 form_of(M4 N, M4 d) {  // B(x, y) = d(N x1 + x2, N y1 + y2)
    int L[8];
    for (int a = 0; a < 8; a++) { int x = 1 << a; L[a] = mv(N, x & 15) ^ (x >> 4); }
    u64 B = 0;
    for (int a = 0; a < 8; a++) for (int b = 0; b < 8; b++) if (a != b && (__popcnt(mv(d, L[a]) & L[b]) & 1)) B |= (u64)1 << (8 * a + b);
    return B;
}
struct Rad { u64 m[4]; int size; };
static inline void radical(u64 B, Rad& r) {  // all x with B(x, .) = 0, by Gray code over the row XORs
    r.m[0] = r.m[1] = r.m[2] = r.m[3] = 0;
    r.m[0] = 1;  // x = 0
    r.size = 1;
    int v = 0, x = 0;
    for (unsigned g = 1; g < 256; g++) {
        int a = _tzcnt_u32(g);
        x ^= 1 << a;
        v ^= (B >> (8 * a)) & 255;
        if (!v) { r.m[x >> 6] |= (u64)1 << (x & 63); r.size++; }
    }
}
static u64 to_coef(u64 B) {  // 28-bit vector, bit k = pair (i, j), i < j in lexicographic order
    u64 c = 0; int k = 0;
    for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++, k++) if (B >> (8 * i + j) & 1) c |= (u64)1 << k;
    return c;
}

struct Sigma { M4 N[16]; };

struct Search {
    const Sigma* S;
    u64 elems[256]; int nel;   // all elements of W (elems[0] = 0)
    u64 cov[4];                // union of radicals (as point sets), including 0
    int covered;               // bitmask over Sigma indices with a form in W
    int dim;
    long long nodes = 0, sols = 0, shorts = 0;
    std::vector<std::vector<u64>>* out;
    std::mutex* mu;
    int holek = 0;             // Sigma indices whose graph is a hole (not a rank-4 radical)
    u64 holept[4] = {0, 0, 0, 0};  // all hole points: Y, plus the graphs in holek (rank-6 radicals lie in here)

    int find_N(M4 N) const { for (int k = 0; k < 16; k++) if (S->N[k] == N) return k; return -1; }

    void set_holes(int hk) {
        holek = hk;
        holept[0] = holept[1] = holept[2] = holept[3] = 0;
        for (int y = 1; y < 16; y++) { int x = y << 4; holept[x >> 6] |= (u64)1 << (x & 63); }
        for (int k = 0; k < 16; k++) if (hk >> k & 1)
            for (int x1 = 1; x1 < 16; x1++) { int x = x1 | mv(S->N[k], x1) << 4; holept[x >> 6] |= (u64)1 << (x & 63); }
    }

    // classify C's radical: returns -2 reject, -1 ok (rank 8, or 6 inside the holes), k >= 0 rank 4 with radical
    // graph(N_k), N_k not a hole
    int classify(u64 C, Rad& r) const {
        radical(C, r);
        if (r.size == 1) return -1;
        if (r.size == 4) {  // rank 6: radical must lie in the union of the holes
            u64 rm0 = r.m[0] & ~(u64)1;
            if ((rm0 & ~holept[0]) | (r.m[1] & ~holept[1]) | (r.m[2] & ~holept[2]) | (r.m[3] & ~holept[3])) return -2;
            return -1;
        }
        if (r.size == 16) {  // rank 4: radical must be graph(N) for N in Sigma
            int col[4] = {-1, -1, -1, -1};
            for (int w = 0; w < 4; w++) {
                u64 m = r.m[w];
                while (m) {
                    int x = w * 64 + (int)_tzcnt_u64(m); m &= m - 1;
                    if (!x) continue;
                    if (!(x & 15)) return -2;  // meets Y
                    int x1 = x & 15;
                    if (__popcnt(x1) == 1) col[_tzcnt_u32(x1)] = x >> 4;
                }
            }
            M4 N = 0;
            for (int j = 0; j < 4; j++) { if (col[j] < 0) return -2; for (int i = 0; i < 4; i++) if (col[j] >> i & 1) N |= 1 << (4 * i + j); }
            int k = find_N(N);
            if (k < 0 || (holek >> k & 1)) return -2;
            return k;
        }
        return -2;
    }

    void report(const char* what) {
        std::vector<u64> basis;  // greedy basis of W
        u64 sp[8]; int n = 0;
        for (int i = 1; i < nel && n < 8; i++) {
            u64 v = elems[i];
            for (int j = 0; j < n; j++) if ((v ^ sp[j]) < v) v ^= sp[j];
            if (v) { sp[n++] = v; basis.push_back(elems[i]); }
        }
        std::lock_guard<std::mutex> g(*mu);
        printf("%s dim %d covered %04x basis", what, dim, covered);
        for (u64 b : basis) printf(" %07llx", (unsigned long long)to_coef(b));
        printf("\n"); fflush(stdout);
        out->push_back(basis);
    }

    void dfs() {
        nodes++;
        const int target = 0xffff & ~holek;
        if (covered == target) {
            if (dim == 8) { sols++; report("SOLUTION"); }
            else { shorts++; report("SHORT"); }
            return;
        }
        if (dim == 8) return;
        int j = _tzcnt_u32(~covered & target);
        for (M4 d : g_D28) {
            u64 B = form_of(S->N[j], d);
            u64 ncov[4] = {cov[0], cov[1], cov[2], cov[3]};
            int ncovered = covered;
            bool ok = true;
            Rad r;
            for (int i = 0; i < nel && ok; i++) {
                u64 C = B ^ elems[i];
                int k = classify(C, r);
                if (k == -2) { ok = false; break; }
                if (k >= 0) { if (ncovered >> k & 1) { ok = false; break; } ncovered |= 1 << k; }
                // radical minus 0 must avoid everything covered so far
                u64 rm0 = r.m[0] & ~(u64)1;
                if ((rm0 & ncov[0]) | (r.m[1] & ncov[1]) | (r.m[2] & ncov[2]) | (r.m[3] & ncov[3])) { ok = false; break; }
                ncov[0] |= rm0; ncov[1] |= r.m[1]; ncov[2] |= r.m[2]; ncov[3] |= r.m[3];
            }
            if (!ok) continue;
            int n0 = nel;
            for (int i = 0; i < n0; i++) elems[n0 + i] = B ^ elems[i];
            u64 sv[4] = {cov[0], cov[1], cov[2], cov[3]}; int sc = covered;
            memcpy(cov, ncov, sizeof cov); covered = ncovered; nel = 2 * n0; dim++;
            dfs();
            dim--; nel = n0; covered = sc; memcpy(cov, sv, sizeof cov);
        }
    }
};

struct Clique {
    const std::vector<M4>& cand; const std::vector<u64>& adj; int W; long maxcl;
    std::atomic<long long>& ncl; std::atomic<long long>& nfield; std::atomic<long long>& nsemi;
    std::atomic<long long>& nodes; std::atomic<long long>& sols; std::atomic<long long>& shorts;
    std::vector<std::vector<u64>>& found; std::mutex& mu; M4 J; double t0;
    int pick[14];

    // out = in & (indices > v)
    void after(std::vector<u64>& out, const u64* in, int v) const {
        for (int w = 0; w < W; w++) {
            int lo = w * 64;
            u64 m = in[w];
            if (lo + 64 <= v + 1) m = 0;
            else if (lo <= v) m &= ~0ull << ((v + 1) - lo);
            out[w] = m;
        }
    }
    void rec(int depth, const std::vector<u64>& P) {
        if (maxcl >= 0 && ncl.load() >= maxcl) return;
        int rem = 0;
        for (int w = 0; w < W; w++) rem += (int)__popcnt64(P[w]);
        if (rem < 14 - depth) return;
        if (depth == 14) { leaf(); return; }
        std::vector<u64> Q(W), T(W);
        for (int w = 0; w < W; w++) {
            u64 m = P[w];
            while (m) {
                int v = w * 64 + (int)_tzcnt_u64(m); m &= m - 1;
                pick[depth] = v;
                for (int k = 0; k < W; k++) T[k] = P[k] & adj[(size_t)v * W + k];
                after(Q, T.data(), v);
                rec(depth + 1, Q);
            }
        }
    }
    void leaf() {
        long long id = ncl++;
        Sigma S;
        S.N[0] = 0; S.N[1] = ID;
        for (int i = 0; i < 14; i++) S.N[2 + i] = cand[pick[i]];
        auto in = [&](M4 x) { for (int c = 0; c < 16; c++) if (S.N[c] == x) return true; return false; };
        bool add = true, mul = true;
        for (int a = 0; a < 16 && add; a++) for (int b = a + 1; b < 16 && add; b++) add = in(S.N[a] ^ S.N[b]);
        if (add) {
            nsemi++;
            for (int a = 2; a < 16 && mul; a++) for (int b = 2; b < 16 && mul; b++) mul = in(mmul(S.N[a], S.N[b]));
            if (mul) nfield++;
        }
        // holes: Y always; X (index 0) from two holes on; with three, one more graph(N_m), m = 2..15, tried in turn.
        // The normalised form sits on X, or on graph(I) (index 1) once X is a hole.
        int k0 = g_holes >= 2 ? 1 : 0;
        std::vector<int> hks;
        if (g_holes == 1) hks.push_back(0);
        else if (g_holes == 2) hks.push_back(1);
        else for (int m = 2; m < 16; m++) hks.push_back(1 | 1 << m);
        for (int hk : hks) {
            Search se;
            se.S = &S; se.out = &found; se.mu = &mu;
            se.set_holes(hk);
            u64 B0 = form_of(S.N[k0], J);
            Rad r; radical(B0, r);
            se.elems[0] = 0; se.elems[1] = B0; se.nel = 2; se.dim = 1; se.covered = 1 << k0;
            se.cov[0] = r.m[0] | 1; se.cov[1] = r.m[1]; se.cov[2] = r.m[2]; se.cov[3] = r.m[3];
            se.dfs();
            nodes += se.nodes; sols += se.sols; shorts += se.shorts;
        }
        if (id % 20000 == 0) {
            std::lock_guard<std::mutex> g(mu);
            printf("cliques %lld  %.0f s  nodes %lld\n", id, omp_get_wtime() - t0, nodes.load()); fflush(stdout);
        }
    }
};

int main(int argc, char** argv) {
    long maxcl = -1;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "i2")) g_holes = 2;
        else if (!strcmp(argv[i], "i3")) g_holes = 3;
        else maxcl = atol(argv[i]);
    }
    printf("holes %d (%d rank-4 components)\n", g_holes, 17 - g_holes);
    for (int m = 0; m < 65536; m++) g_inv[m] = invertible((M4)m);
    for (int b = 0; b < 64; b++) { M4 a = alt4(b); if (g_inv[a]) g_D28.push_back(a); }
    printf("nondegenerate alternating 4x4: %zu (expect 28)\n", g_D28.size());
    // J = the first one is fine as the normalised form on X (any nondegenerate one; all are GL(4,2)-equivalent)
    M4 J = g_D28[0];

    std::vector<M4> cand;
    for (int m = 0; m < 65536; m++) if (g_inv[m] && g_inv[m ^ ID] && m != ID) cand.push_back((M4)m);
    int nc = (int)cand.size();
    int W = (nc + 63) / 64;
    printf("candidates (N, N + I invertible): %d\n", nc); fflush(stdout);
    std::vector<u64> adj((size_t)nc * W, 0);
    for (int a = 0; a < nc; a++) for (int b = a + 1; b < nc; b++) if (g_inv[cand[a] ^ cand[b]]) {
        adj[(size_t)a * W + b / 64] |= (u64)1 << (b % 64);
        adj[(size_t)b * W + a / 64] |= (u64)1 << (a % 64);
    }

    std::atomic<long long> ncl(0), nfield(0), nsemi(0), nodes(0), sols(0), shorts(0);
    std::vector<std::vector<u64>> found;
    std::mutex mu;
    double t0 = omp_get_wtime();

    // all 14-cliques (ordered: increasing candidate index), parallel over the first vertex
    #pragma omp parallel for schedule(dynamic, 1)
    for (int v0 = 0; v0 < nc; v0++) {
        Clique cq{cand, adj, W, maxcl, ncl, nfield, nsemi, nodes, sols, shorts, found, mu, J, t0};
        std::vector<u64> P0(W);
        cq.after(P0, &adj[(size_t)v0 * W], v0);
        cq.pick[0] = v0;
        cq.rec(1, P0);
    }
    printf("spread sets (cliques) %lld, semifield %lld, field %lld (expect 336)\n", ncl.load(), nsemi.load(), nfield.load());
    printf("form search: nodes %lld, solutions %lld, short %lld, %.0f s\n", nodes.load(), sols.load(), shorts.load(), omp_get_wtime() - t0);
    return 0;
}
