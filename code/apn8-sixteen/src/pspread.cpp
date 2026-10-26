// pspread: partial spreads of K solids (pairwise complementary 4-spaces) in F_2^8, up to GL(8,2).
// Rowen (Claude Opus 5.5), 2026-10-25. Written for i = 5 (K = 12 rank-4 radicals); README "Twelve".
// Written without reading thirteen_check.cpp / fourteen_check.cpp (only their write-ups' numbers, used as checks).
//
// 1. Three members go to Y = {(0,y)}, X = {(x,0)}, Z = graph(I); the other K - 3 are graph(N), N in V5824 = {N : N, N + I
//    invertible}, pairwise N + N' invertible: a (K-3)-clique. Conjugation by GL(4,2) (diag(h,h)) fixes X, Y, Z.
//    Reduction: c = the least conjugacy class present; one class-c member is moved to the class rep r_c; the rest are
//    a clique in adj(r_c) with classes >= c. A reduced clique K stands for |class c| / m(K) cliques (m = number of
//    class-c members), kept exact as multiples of 360360.
// 2. Each clique is classified by how far it extends: w = the largest clique in its common neighbourhood in V5824.
// 3. Canonical form over GL(8,2): over every ordered triple (A,B,C) of members, the basis change g with A -> X, B -> Y,
//    C -> Z (unique up to diag(h,h)) turns the other members into graph(N)'s. Frames with the least sorted class vector
//    are kept; over those, h runs through the h with h N h^-1 = r_c for a least-class member N; key = least sorted
//    index list of the conjugated members. The (frame, h) giving the key number |Aut P| (every g in GL(8,2) carrying P
//    onto the key's configuration is among them).
// 4. Check: sum over orbits of K(K-1)(K-2) * 20160 / |Aut P| = the weighted clique count, per extension class.
//
// usage: pspread enum K out      configurations (cover.cpp input lines: tag, K, 4K hex bytes) + summary on stdout
//        pspread selftest K file  re-canonicalise each line's configuration after a random GL(8,2) map and shuffle
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <mutex>
#include <random>
#include <functional>
#include <intrin.h>
#include <omp.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;

static inline int par(u32 x) { return __popcnt(x) & 1; }

// ---- 4x4 matrices: 16 bits, row i = bits 4i..4i+3; (m x)_i = parity(row_i & x)
static inline int mv4(u32 m, int x) {
    int y = 0;
    for (int i = 0; i < 4; i++) y |= par((m >> (4 * i)) & 15 & x) << i;
    return y;
}
static inline u32 from_cols4(const int* c) {
    u32 m = 0;
    for (int j = 0; j < 4; j++) for (int i = 0; i < 4; i++) if (c[j] >> i & 1) m |= 1u << (4 * i + j);
    return m;
}
static inline u32 mul4(u32 a, u32 b) {
    int c[4];
    for (int j = 0; j < 4; j++) c[j] = mv4(a, mv4(b, 1 << j));
    return from_cols4(c);
}
static const u32 I4 = 0x8421;

// ---- 8x8 matrices as 8 column bytes
struct M8 { u8 c[8]; };
static inline int ap8(const M8& m, int x) { int y = 0; for (int j = 0; j < 8; j++) if (x >> j & 1) y ^= m.c[j]; return y; }
static bool inv8(const M8& m, M8& r) {
    // rows of [m | I] as 16-bit: row i = (bits of row i of m) | (row i of I) << 8
    u32 row[8];
    for (int i = 0; i < 8; i++) {
        u32 v = 0;
        for (int j = 0; j < 8; j++) if (m.c[j] >> i & 1) v |= 1u << j;
        row[i] = v | (1u << (8 + i));
    }
    for (int col = 0; col < 8; col++) {
        int p = -1;
        for (int i = col; i < 8; i++) if (row[i] >> col & 1) { p = i; break; }
        if (p < 0) return false;
        std::swap(row[p], row[col]);
        for (int i = 0; i < 8; i++) if (i != col && (row[i] >> col & 1)) row[i] ^= row[col];
    }
    // row i of inverse = row[i] >> 8; columns:
    for (int j = 0; j < 8; j++) {
        int v = 0;
        for (int i = 0; i < 8; i++) if (row[i] >> (8 + j) & 1) v |= 1 << i;
        r.c[j] = (u8)v;
    }
    return true;
}

// ---- globals
static std::vector<u32> GL;                 // GL(4,2)
static std::vector<u32> INV(65536, 0);      // inverse of invertible matrices
static std::vector<int> IDX(65536, -1);     // V5824 index
static std::vector<u32> V;                  // V5824 members
static std::vector<int> CLS;                // class of V[i] (0..4, ordered by rep index)
static std::vector<int> REP;                // rep index per class
static std::vector<int> CSIZE;
static std::vector<u32> TOREP;              // per V index: h with h N h^-1 = rep
static std::vector<std::vector<u32>> CENT;  // centraliser of each rep
static int NV, NW;
static std::vector<u64> ADJ;                // NV x NW bitset
static inline const u64* adj(int i) { return &ADJ[(size_t)i * NW]; }

static bool invertible4(u32 m) {
    int r[4];
    for (int i = 0; i < 4; i++) r[i] = (m >> (4 * i)) & 15;
    int rank = 0;
    for (int col = 0; col < 4; col++) {
        int p = -1;
        for (int i = rank; i < 4; i++) if (r[i] >> col & 1) { p = i; break; }
        if (p < 0) continue;
        std::swap(r[p], r[rank]);
        for (int i = 0; i < 4; i++) if (i != rank && (r[i] >> col & 1)) r[i] ^= r[rank];
        rank++;
    }
    return rank == 4;
}
static inline u32 conj(u32 h, u32 n) { return mul4(mul4(h, n), INV[h]); }

static void setup() {
    for (u32 m = 0; m < 65536; m++) if (invertible4(m)) GL.push_back(m);
    for (u32 a : GL) for (u32 b : GL) if (mul4(a, b) == I4) { INV[a] = b; break; }
    for (u32 m = 0; m < 65536; m++) if (invertible4(m) && invertible4(m ^ I4)) { IDX[m] = (int)V.size(); V.push_back(m); }
    NV = (int)V.size(); NW = (NV + 63) / 64;
    CLS.assign(NV, -1); TOREP.assign(NV, 0);
    for (int i = 0; i < NV; i++) {
        if (CLS[i] >= 0) continue;
        int c = (int)REP.size(); REP.push_back(i); CSIZE.push_back(0);
        std::vector<u32> cent;
        for (u32 h : GL) {
            u32 n2 = conj(h, V[i]);
            int j = IDX[n2];
            if (j == i) cent.push_back(h);
            if (CLS[j] < 0) { CLS[j] = c; TOREP[j] = INV[h]; CSIZE[c]++; }  // h V[i] h^-1 = V[j]  =>  h^-1 V[j] h = V[i]
        }
        CENT.push_back(cent);
    }
    ADJ.assign((size_t)NV * NW, 0);
    for (int i = 0; i < NV; i++) for (int j = 0; j < NV; j++)
        if (i != j && invertible4(V[i] ^ V[j])) ADJ[(size_t)i * NW + (j >> 6)] |= 1ull << (j & 63);
}

// ---- solids as 4 basis bytes; point x = x1 | x2 << 4
struct Sol { u8 b[4]; };
static Sol graphN(u32 n) { Sol s; for (int i = 0; i < 4; i++) s.b[i] = (u8)((1 << i) | (mv4(n, 1 << i) << 4)); return s; }
static Sol solY() { Sol s; for (int i = 0; i < 4; i++) s.b[i] = (u8)(1 << (4 + i)); return s; }

// N with g(S) = graph(N), or -1 if g(S) isn't complementary to Y
static int graph_of(const M8& g, const Sol& s) {
    int lo[4], hi[4];
    for (int j = 0; j < 4; j++) { int p = ap8(g, s.b[j]); lo[j] = p & 15; hi[j] = p >> 4; }
    // L = cols lo, H = cols hi; N = H L^-1
    u32 L = from_cols4(lo), H = from_cols4(hi);
    if (!invertible4(L)) return -1;
    return (int)mul4(H, INV[L]);
}

struct Canon { std::vector<u16> key; int aut; };

static Canon canon(const std::vector<Sol>& P) {
    int K = (int)P.size(), R = K - 3;
    std::vector<int> bestcv;
    std::vector<std::vector<int>> frames;  // V indices of the other members, per minimal frame
    std::vector<int> cv(R), ns(R);
    for (int a = 0; a < K; a++) for (int b = 0; b < K; b++) {
        if (b == a) continue;
        M8 AB;
        for (int j = 0; j < 4; j++) { AB.c[j] = P[a].b[j]; AB.c[4 + j] = P[b].b[j]; }
        M8 ABi;
        if (!inv8(AB, ABi)) { fprintf(stderr, "not complementary\n"); exit(1); }
        for (int c = 0; c < K; c++) {
            if (c == a || c == b) continue;
            // phi(a_i): the B-part of the point of C whose A-part is a_i. C's points: span of its basis.
            int phi[16];  // indexed by A-coefficient vector
            for (int t = 0; t < 16; t++) phi[t] = -1;
            for (int t = 1; t < 16; t++) {
                int p = 0;
                for (int j = 0; j < 4; j++) if (t >> j & 1) p ^= P[c].b[j];
                int co = ap8(ABi, p);       // coordinates in (a1..a4, b1..b4)
                int ac = co & 15, bc = co >> 4;
                int bp = 0;
                for (int j = 0; j < 4; j++) if (bc >> j & 1) bp ^= P[b].b[j];
                phi[ac] = bp;
            }
            M8 Gi;  // e_i -> a_i, e_{4+i} -> phi(a_i)
            bool ok = true;
            for (int j = 0; j < 4; j++) { Gi.c[j] = P[a].b[j]; if (phi[1 << j] < 0) ok = false; Gi.c[4 + j] = (u8)phi[1 << j]; }
            if (!ok) { fprintf(stderr, "C not a graph over A+B\n"); exit(1); }
            M8 G;
            if (!inv8(Gi, G)) { fprintf(stderr, "frame singular\n"); exit(1); }
            int r = 0;
            for (int m = 0; m < K; m++) {
                if (m == a || m == b || m == c) continue;
                int n = graph_of(G, P[m]);
                if (n < 0 || IDX[n] < 0) { fprintf(stderr, "member not a graph in V5824\n"); exit(1); }
                ns[r] = IDX[n]; cv[r] = CLS[IDX[n]]; r++;
            }
            std::vector<int> s = cv;
            std::sort(s.begin(), s.end());
            if (bestcv.empty() || s < bestcv) { bestcv = s; frames.clear(); }
            if (s == bestcv) frames.push_back(ns);
        }
    }
    int c0 = bestcv[0];
    Canon out; out.aut = 0;
    std::vector<u16> cur(R);
    for (auto& f : frames) {
        for (int j = 0; j < R; j++) {
            if (CLS[f[j]] != c0) continue;
            for (u32 z : CENT[c0]) {
                u32 h = mul4(z, TOREP[f[j]]);
                for (int l = 0; l < R; l++) cur[l] = (u16)IDX[conj(h, V[f[l]])];
                std::sort(cur.begin(), cur.end());
                if (out.key.empty() || cur < out.key) { out.key = cur; out.aut = 1; }
                else if (cur == out.key) out.aut++;
            }
        }
    }
    return out;
}

static std::vector<Sol> config_of(const std::vector<u16>& key) {
    std::vector<Sol> P;
    P.push_back(solY()); P.push_back(graphN(0)); P.push_back(graphN(I4));
    for (u16 k : key) P.push_back(graphN(V[k]));
    return P;
}

// ---- largest clique (capped) inside a candidate bitset
static int maxclique(std::vector<u64> cand, int cap) {
    int best = 0;
    for (int w = 0; w < NW; w++) while (cand[w]) {
        int i = w * 64 + (int)_tzcnt_u64(cand[w]);
        cand[w] &= cand[w] - 1;
        std::vector<u64> nc(NW);
        bool any = false;
        for (int u = 0; u < NW; u++) { nc[u] = cand[u] & adj(i)[u]; any |= nc[u] != 0; }
        int s = 1 + (any && cap > 1 ? maxclique(nc, cap - 1) : 0);
        if (s > best) best = s;
        if (best >= cap) return best;
    }
    return best;
}

struct Orbit { int aut; int ext; u64 reduced; };

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: pspread enum K out | selftest K file\n"); return 1; }
    setup();
    printf("GL(4,2) %zu, V5824 %d, classes", GL.size(), NV);
    for (size_t c = 0; c < REP.size(); c++) printf(" %d(cent %zu)", CSIZE[c], CENT[c].size());
    printf("\n"); fflush(stdout);
    std::string mode = argv[1];
    int K = atoi(argv[2]), R = K - 3;
    const u64 SCALE = 360360;

    if (mode == "selftest") {
        FILE* f = fopen(argv[3], "r");
        char line[4096]; int n = 0, bad = 0;
        std::mt19937_64 rng(12345);
        while (fgets(line, sizeof line, f)) {
            char tag[256]; int k, off;
            if (sscanf(line, "%255s %d%n", tag, &k, &off) != 2 || k != K) continue;
            std::vector<Sol> P(K);
            char* p = line + off;
            for (int s = 0; s < K; s++) for (int j = 0; j < 4; j++) { char* e; P[s].b[j] = (u8)strtoul(p, &e, 16); p = e; }
            Canon c1 = canon(P);
            M8 g, gi;
            do { for (int j = 0; j < 8; j++) g.c[j] = (u8)rng(); } while (!inv8(g, gi));
            std::vector<Sol> Q = P;
            for (auto& s : Q) {
                u32 C = GL[rng() % GL.size()];  // random basis of the image: rows of C as coefficients
                Sol t;
                for (int j = 0; j < 4; j++) {
                    int v = 0;
                    for (int b = 0; b < 4; b++) if (C >> (4 * j + b) & 1) v ^= s.b[b];
                    t.b[j] = (u8)ap8(g, v);
                }
                s = t;
            }
            std::shuffle(Q.begin(), Q.end(), rng);
            Canon c2 = canon(Q);
            n++;
            if (c1.key != c2.key || c1.aut != c2.aut) { bad++; printf("MISMATCH %s aut %d %d\n", tag, c1.aut, c2.aut); }
        }
        printf("selftest: %d configurations, %d mismatches\n", n, bad);
        return bad != 0;
    }

    // ---- enumeration
    std::vector<u64> wsum(18, 0), rsum(18, 0);
    std::map<std::vector<u16>, Orbit> orbits;
    std::mutex mu;
    double t0 = omp_get_wtime();
    for (int c = 0; c < (int)REP.size(); c++) {
        int r0 = REP[c];
        // candidates: adj(r0), class >= c
        std::vector<u64> cand(NW, 0);
        for (int j = 0; j < NV; j++) if (CLS[j] >= c && (adj(r0)[j >> 6] >> (j & 63) & 1)) cand[j >> 6] |= 1ull << (j & 63);
        std::vector<int> first;
        for (int j = 0; j < NV; j++) if (cand[j >> 6] >> (j & 63) & 1) first.push_back(j);
        u64 wc = 0;
        #pragma omp parallel
        {
            std::vector<u64> lw(18, 0), lr(18, 0);
            std::map<std::vector<u16>, Orbit> lo;
            std::vector<int> stack;
            #pragma omp for schedule(dynamic, 1)
            for (int fi = 0; fi < (int)first.size(); fi++) {
                // cliques of size R - 1 in cand with first element first[fi] (members in increasing index)
                int f0 = first[fi];
                std::vector<std::vector<u64>> lev(R + 1, std::vector<u64>(NW));
                std::vector<int> mem(R);
                mem[0] = r0; mem[1] = f0;
                for (int u = 0; u < NW; u++) lev[2][u] = cand[u] & adj(f0)[u];
                // only later indices
                for (int j = 0; j <= f0; j++) lev[2][j >> 6] &= ~(1ull << (j & 63));
                // recursive DFS via explicit lambda
                std::function<void(int)> dfs = [&](int d) {
                    if (d == R) {
                        // clique complete: mem[0..R-1]
                        int m = 0;
                        for (int t = 0; t < R; t++) m += CLS[mem[t]] == c;
                        u64 w = (u64)CSIZE[c] * SCALE / m;
                        std::vector<u64> q(NW);
                        for (int u = 0; u < NW; u++) q[u] = ~0ull;
                        for (int t = 0; t < R; t++) for (int u = 0; u < NW; u++) q[u] &= adj(mem[t])[u];
                        if (NV % 64) q[NW - 1] &= (1ull << (NV % 64)) - 1;
                        int e = maxclique(q, 17 - K);
                        lw[e] += w; lr[e]++;
                        std::vector<Sol> P;
                        P.push_back(solY()); P.push_back(graphN(0)); P.push_back(graphN(I4));
                        for (int t = 0; t < R; t++) P.push_back(graphN(V[mem[t]]));
                        Canon cn = canon(P);
                        auto it = lo.find(cn.key);
                        if (it == lo.end()) lo[cn.key] = Orbit{cn.aut, e, 1};
                        else {
                            it->second.reduced++;
                            if (it->second.aut != cn.aut || it->second.ext != e) {
                                std::lock_guard<std::mutex> g(mu);
                                printf("INCONSISTENT orbit (aut %d vs %d, ext %d vs %d)\n", it->second.aut, cn.aut, it->second.ext, e);
                            }
                        }
                        return;
                    }
                    std::vector<u64>& cs = lev[d];
                    for (int w = 0; w < NW; w++) {
                        u64 bits = cs[w];
                        while (bits) {
                            int i = w * 64 + (int)_tzcnt_u64(bits);
                            bits &= bits - 1;
                            mem[d] = i;
                            if (d + 1 < R) {
                                bool any = false;
                                for (int u = 0; u < NW; u++) {
                                    u64 v = cs[u] & adj(i)[u];
                                    if (u < (i >> 6)) v = 0;
                                    else if (u == (i >> 6)) v &= (i & 63) == 63 ? 0 : ~((2ull << (i & 63)) - 1);
                                    lev[d + 1][u] = v; any |= v != 0;
                                }
                                if (!any) continue;
                            }
                            dfs(d + 1);
                        }
                    }
                };
                if (R == 1) { /* not used */ }
                else if (R == 2) { for (int u = 0; u < NW; u++) lev[2][u] = 0; mem[1] = f0; dfs(2); }
                else dfs(2);
            }
            std::lock_guard<std::mutex> g(mu);
            for (int e = 0; e < 18; e++) { wsum[e] += lw[e]; rsum[e] += lr[e]; }
            for (auto& kv : lo) {
                auto it = orbits.find(kv.first);
                if (it == orbits.end()) orbits[kv.first] = kv.second;
                else {
                    it->second.reduced += kv.second.reduced;
                    if (it->second.aut != kv.second.aut || it->second.ext != kv.second.ext)
                        printf("INCONSISTENT orbit at merge\n");
                }
            }
        }
        printf("class %d done, %.0f s, orbits so far %zu\n", c, omp_get_wtime() - t0, orbits.size()); fflush(stdout);
    }
    // ---- summary + check
    u64 frames = (u64)K * (K - 1) * (K - 2) * 20160;
    printf("reduced cliques by extension w (largest clique in the common neighbourhood):\n");
    for (int e = 0; e < 18; e++) if (rsum[e]) {
        u64 osum = 0; int no = 0;
        for (auto& kv : orbits) if (kv.second.ext == e) { no++; osum += frames / kv.second.aut; if (frames % kv.second.aut) printf("aut does not divide frames\n"); }
        u64 wtot = wsum[e] / SCALE;
        printf("  w=%d: reduced %llu, weighted %llu%s, orbits %d, orbit sum %llu %s\n", e, (unsigned long long)rsum[e],
               (unsigned long long)wtot, wsum[e] % SCALE ? " (non-integer!)" : "", no, (unsigned long long)osum,
               osum == wtot ? "MATCH" : "MISMATCH");
    }
    printf("orbits %zu\n", orbits.size());
    FILE* out = fopen(argv[3], "w");
    int n = 0;
    for (auto& kv : orbits) {
        std::vector<Sol> P = config_of(kv.first);
        fprintf(out, "p%d-%d-a%d-w%d %d", K, n++, kv.second.aut, kv.second.ext, K);
        for (auto& s : P) for (int j = 0; j < 4; j++) fprintf(out, " %02x", s.b[j]);
        fprintf(out, "\n");
    }
    fclose(out);
    printf("wrote %d configurations, %.0f s\n", n, omp_get_wtime() - t0);
    return 0;
}
