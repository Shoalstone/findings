// aut8: the linear automorphism group Aut(S8) = {A in GL(8,2) : S8 o A = S8} of quadratic APN functions with >= 3
// rank-4 components (for quadratic F, this is the input side of its EA self-equivalences). Any such A permutes the
// rank-4 radicals; with R1, R2, R3 three of them (pairwise complementary, R3 = graph of T : R1 -> R2), A is fixed by
// where it sends R1, R2, R3 and by A|R1, since A|R2 = T' A|R1 T^-1. So: k(k-1)(k-2) x 20160 candidates, each checked.
// aut8 <seed file>...  (lines as load_seeds: {LUT} or "... basis b0..b7"); prints one line per function.
// AUT8_PRINT=<order> in the environment also prints each automorphism of that order as "A <col0..col7>" (A e_i, hex).
// Rowen (Claude Opus 5.5), 2026-10-24.
#define main apn8_main
#include "apn8.cpp"
#undef main
#include <map>
#include <string>

static inline int bform(u64 m, unsigned x, unsigned y) { return __popcnt(rvec(m, x) & y) & 1; }

struct Rad { std::vector<u8> el; u8 mem[256]; };
static Rad radical_of(u64 m) { Rad r; memset(r.mem, 0, 256); r.mem[0] = 1; for (unsigned x = 1; x < 256; x++) if (rvec(m, x) == 0) { r.el.push_back((u8)x); r.mem[x] = 1; } return r; }

// T : R1 -> R2 with R3 = {x + T x}; as a 256-table on R1's elements
static bool graph_map(const Rad& R1, const Rad& R2, const Rad& R3, u8* T) {
    u8 p1[256], p2[256]; memset(p1, 0, 256); memset(p2, 0, 256);
    std::vector<u8> a1(R1.el), a2(R2.el); a1.push_back(0); a2.push_back(0);
    for (u8 a : a1) for (u8 b : a2) { p1[a ^ b] = a; p2[a ^ b] = b; }
    memset(T, 0, 256);
    for (u8 r : R3.el) { if (!p1[r] || !p2[r]) return false; T[p1[r]] = p2[r]; }
    return true;
}

int main(int argc, char** argv) {
    init_index();
    const char* po = getenv("AUT8_PRINT"); int print_order = po ? atoi(po) : 0;
    for (int fi = 1; fi < argc; fi++) {
        std::vector<Found> fs; if (!load_seeds(argv[fi], fs)) { fprintf(stderr, "can't read %s\n", argv[fi]); continue; }
        for (size_t idx = 0; idx < fs.size(); idx++) {
            Found& F = fs[idx];
            if (F.cnt[2] < 3 || F.cnt[1]) continue;
            u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(F.basis[i]);
            Ech28 S; for (int i = 0; i < 8; i++) S.add(F.basis[i]);
            std::vector<Rad> R; u64 cur = 0;
            for (unsigned g = 1; g < 256; g++) { cur ^= m[_tzcnt_u32(g)]; if (rank_mat(cur) == 4) R.push_back(radical_of(cur)); }
            int k = (int)R.size();
            u8 e[4]; { u8 piv[8] = {0}; int n = 0; for (u8 x : R[0].el) { u8 y = x; for (int b = 7; b >= 0 && y; b--) if (y >> b & 1) { if (piv[b]) y ^= piv[b]; else { piv[b] = y; y = 0; e[n++] = x; } } if (n == 4) break; } }
            u8 T[256]; if (!graph_map(R[0], R[1], R[2], T)) { fprintf(stderr, "graph\n"); return 1; }
            std::map<int, long> ord, fixv, fixs; std::map<std::string, long> perm; long total = 0;
            for (int s1 = 0; s1 < k; s1++) for (int s2 = 0; s2 < k; s2++) for (int s3 = 0; s3 < k; s3++) {
                if (s1 == s2 || s1 == s3 || s2 == s3) continue;
                u8 Tp[256]; if (!graph_map(R[s1], R[s2], R[s3], Tp)) continue;
                const auto& E = R[s1].el;
                for (u8 u0 : E) for (u8 u1 : E) for (u8 u2 : E) for (u8 u3 : E) {
                    u8 u[4] = {u0, u1, u2, u3}; if (rank8(u, 4) != 4) continue;
                    u8 b[8], im[8]; for (int i = 0; i < 4; i++) { b[i] = e[i]; im[i] = u[i]; b[4 + i] = T[e[i]]; im[4 + i] = Tp[u[i]]; }
                    // LUT of A from the basis b -> im
                    u8 A[256]; { u8 inv[256]; u8 x = 0, y = 0; A[0] = 0; for (unsigned g = 1; g < 256; g++) { int t = _tzcnt_u32(g); x ^= b[t]; y ^= im[t]; A[x] = y; } (void)inv; }
                    bool ok = true;
                    for (int t = 0; t < 8 && ok; t++) {
                        u32 c = 0; for (int q = 0; q < 28; q++) if (bform(m[t], A[1u << KI[q]], A[1u << KJ[q]])) c |= 1u << q;
                        if (S.reduce(c)) ok = false;
                    }
                    if (!ok) continue;
                    total++;
                    int o = 1; { u8 P[256]; memcpy(P, A, 256); for (;;) { bool id = true; for (int x = 1; x < 256; x <<= 1) if (P[x] != x) { id = false; break; } if (id) break; u8 Q[256]; for (int x = 0; x < 256; x++) Q[x] = A[P[x]]; memcpy(P, Q, 256); o++; } }
                    ord[o]++;
                    if (print_order == o) { printf("A"); for (int i = 0; i < 8; i++) printf(" %02x", A[1u << i]); printf("\n"); }
                    { int fx = 0; for (int x = 1; x < 256; x++) if (A[x] == x) fx++; int d = 0; while ((1 << d) - 1 < fx) d++; fixv[o * 100 + d]++; }
                    { // dim of S8 & Fix(A) (forms with B o A = B)
                        int fs = 0; u64 cm = 0;
                        for (unsigned g = 1; g < 256; g++) { cm ^= m[_tzcnt_u32(g)]; u32 c = mat_to_coef(cm), c2 = 0;
                            for (int q = 0; q < 28; q++) if (bform(cm, A[1u << KI[q]], A[1u << KJ[q]])) c2 |= 1u << q;
                            if (c2 == c) fs++; }
                        int d = 0; while ((1 << d) - 1 < fs) d++; fixs[o * 100 + d]++;
                    }
                    // permutation of the rank-4 radicals: R_i -> A R_i
                    std::vector<int> pi(k); for (int i = 0; i < k; i++) { pi[i] = -1; for (int j = 0; j < k; j++) if (R[j].mem[A[R[i].el[0]]]) pi[i] = j; }
                    std::vector<int> cyc; std::vector<int> vis(k, 0);
                    for (int i = 0; i < k; i++) if (!vis[i]) { int l = 0, j = i; while (!vis[j]) { vis[j] = 1; j = pi[j]; l++; } cyc.push_back(l); }
                    std::sort(cyc.rbegin(), cyc.rend()); std::string cs; for (int l : cyc) cs += std::to_string(l);
                    perm[cs]++;
                }
            }
            printf("%s:%zu ranks 8:%d 6:%d 4:%d |Aut| %ld orders", argv[fi], idx, F.cnt[4], F.cnt[3], F.cnt[2], total);
            for (auto& kv : ord) printf(" %d^%ld", kv.first, kv.second);
            printf(" fixV(order:dim)"); for (auto& kv : fixv) printf(" %d:%d^%ld", kv.first / 100, kv.first % 100, kv.second);
            printf(" fixS8(order:dim)"); for (auto& kv : fixs) printf(" %d:%d^%ld", kv.first / 100, kv.first % 100, kv.second);
            printf(" perms"); for (auto& kv : perm) printf(" %s:%ld", kv.first.c_str(), kv.second);
            printf("\n"); fflush(stdout);
        }
    }
    return 0;
}
