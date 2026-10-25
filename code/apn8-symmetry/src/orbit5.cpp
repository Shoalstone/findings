// orbit5: every quadratic (8,8)-APN function whose component space S8 is invariant under an order-5 linear map A
// (a linear self-equivalence, input side) and contains a rank-4 form B1 whose radical has an A-orbit of 5 pairwise
// complementary radicals, i.e. 5 rank-4 components permuted cyclically. A case of Open Problem VI.1.
//
// GL(8,2) has two classes of order-5 elements: char poly Phi5^2 (V = F_16^2, A = zeta scalar) and Phi5 (x+1)^4
// (A = zeta on F_16 + identity on F_2^4). A acts on forms by (A.B)(x, y) = B(A^-1 x, A^-1 y). Alt(V) = Fix + N
// (trivial and nontrivial isotypic parts; projection onto Fix is 1 + A + .. + A^4). An A-invariant S8 is
// (S8 & Fix) + (S8 & N), and S8 & N is an F_16-subspace (F_16 = F_2[A] acting on N), of F_2-dim 4k.
// B1 = f + n (f in Fix, n in N, n != 0 since B1's radical isn't fixed), and S8 contains f and F_16 n. So:
//   k = 1: S8 = F_16 n + Q, Q a 4-dim subspace of Fix containing f;
//   k = 2: f = 0, S8 = F_16 n + F_16 n'.
// Output: each APN S8 found (deduplicated by its canonical basis) with its rank distribution.
// Rowen (Claude Opus 5.5), 2026-10-24.
#define main apn8_main
#include "apn8.cpp"
#undef main
#include <set>
#include <mutex>
#include <map>
#include <functional>
#include <array>

struct Set256 { u64 w[4] = {0, 0, 0, 0}; void add(unsigned x) { w[x >> 6] |= 1ull << (x & 63); } bool meets(const Set256& o) const { return (w[0] & o.w[0] & ~1ull) | (w[1] & o.w[1]) | (w[2] & o.w[2]) | (w[3] & o.w[3]); } };
static int radical(u64 m, Set256& S) {
    u8 r[8], c[8]; for (int i = 0; i < 8; i++) { r[i] = row(m, i); c[i] = (u8)(1 << i); }
    int rk = 0;
    for (int bit = 7; bit >= 0; bit--) {
        int p = -1; for (int i = rk; i < 8; i++) if (r[i] >> bit & 1) { p = i; break; }
        if (p < 0) continue;
        std::swap(r[p], r[rk]); std::swap(c[p], c[rk]);
        for (int i = 0; i < 8; i++) if (i != rk && (r[i] >> bit & 1)) { r[i] ^= r[rk]; c[i] ^= c[rk]; }
        rk++;
    }
    int d = 8 - rk; S = Set256();
    for (unsigned g = 0; g < (1u << d); g++) { u8 x = 0; for (int k = 0; k < d; k++) if (g >> k & 1) x ^= c[rk + k]; S.add(x); }
    return d;
}

// GF(16) with x^4 + x + 1, GF(256) with x^8 + x^4 + x^3 + x^2 + 1 (0x11D)
static unsigned gmul(unsigned a, unsigned b, unsigned poly, int deg) {
    unsigned r = 0; while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a >> deg & 1) a ^= poly; } return r;
}

static u8 Av[256], Ai[256];
static u32 actcol[28]; // A. on coef basis vectors
static u32 act(u32 c) { u32 r = 0; while (c) { int k = (int)_tzcnt_u32(c); r ^= actcol[k]; c &= c - 1; } return r; }
static inline int bform(u64 m, unsigned x, unsigned y) { return __popcnt(rvec(m, x) & y) & 1; }

static void setup_A(int cls) {
    if (cls == 1) { // V = GF(256), A = multiplication by an element of order 5
        unsigned z = 1; for (int i = 0; i < 51; i++) z = gmul(z, 2, 0x11D, 8); // 2 is primitive for 0x11D
        for (unsigned x = 0; x < 256; x++) Av[x] = (u8)gmul(x, z, 0x11D, 8);
    } else {        // low nibble GF(16) times zeta (order 5: 2^3), high nibble fixed
        unsigned z = 8;
        for (unsigned x = 0; x < 256; x++) Av[x] = (u8)(gmul(x & 15, z, 0x13, 4) | (x & 0xF0));
    }
    for (unsigned x = 0; x < 256; x++) Ai[Av[x]] = (u8)x;
    // order check
    for (unsigned x = 1; x < 256; x++) { unsigned y = x; for (int i = 0; i < 5; i++) y = Av[y]; if (y != x) { fprintf(stderr, "A order\n"); exit(1); } }
    for (int k = 0; k < 28; k++) {
        u64 m = coef_to_mat(1u << k); u32 r = 0;
        for (int q = 0; q < 28; q++) if (bform(m, Ai[1u << KI[q]], Ai[1u << KJ[q]])) r |= 1u << q;
        actcol[k] = r;
    }
}

static bool g_enumtest = false; // orbit5 <class> <out> <B1> enumtest: k1 enumerator check, no pruning, count subspaces
static std::set<std::vector<u32>> g_enum; static long g_enumcount = 0;
struct Elem { u64 m; Set256 rad; int rk; };
static bool partition_ok(const u32* fs, int n) { // all nonzero combos: rank >= 2, radicals pairwise disjoint
    std::vector<Elem> el; Ech28 e;
    for (int i = 0; i < n; i++) {
        if (!e.add(fs[i])) continue;
        u64 fm = coef_to_mat(fs[i]); size_t old = el.size();
        std::vector<u64> news; news.push_back(fm); for (size_t j = 0; j < old; j++) news.push_back(fm ^ el[j].m);
        for (u64 m : news) {
            Elem x; x.m = m; x.rk = rank_mat(m);
            if (x.rk == 0) return false;
            if (x.rk < 8) { radical(m, x.rad); for (auto& y : el) if (y.rk < 8 && x.rad.meets(y.rad)) return false; }
            el.push_back(x);
        }
    }
    return true;
}

int main(int argc, char** argv) {
    init_index();
    if (argc < 3) { fprintf(stderr, "usage: orbit5 <class 1|2> <out>\n"); return 2; }
    int cls = atoi(argv[1]); setup_A(cls);
    if (argc > 4 && !strcmp(argv[4], "enumtest")) g_enumtest = true;
    // Fix basis: kernel of act - I
    std::vector<u32> fix;
    {
        u32 rows[28], comb[28]; for (int k = 0; k < 28; k++) { rows[k] = actcol[k] ^ (1u << k); comb[k] = 1u << k; }
        int rk = 0;
        for (int bit = 27; bit >= 0; bit--) {
            int p = -1; for (int i = rk; i < 28; i++) if (rows[i] >> bit & 1) { p = i; break; }
            if (p < 0) continue;
            std::swap(rows[p], rows[rk]); std::swap(comb[p], comb[rk]);
            for (int i = 0; i < 28; i++) if (i != rk && (rows[i] >> bit & 1)) { rows[i] ^= rows[rk]; comb[i] ^= comb[rk]; }
            rk++;
        }
        for (int i = rk; i < 28; i++) fix.push_back(comb[i]);
    }
    int df = (int)fix.size();
    std::vector<u32> fixel(1u << df); for (u32 g = 0; g < (1u << df); g++) { u32 v = 0; for (int i = 0; i < df; i++) if (g >> i & 1) v ^= fix[i]; fixel[g] = v; }
    for (u32 v : fixel) if (act(v) != v) { fprintf(stderr, "fix wrong\n"); return 1; }
    auto pfix = [&](u32 v) { u32 a1 = act(v), a2 = act(a1), a3 = act(a2), a4 = act(a3); return v ^ a1 ^ a2 ^ a3 ^ a4; };
    std::vector<u32> nbas; { Ech28 e; for (int k = 0; k < 28; k++) { u32 v = (1u << k) ^ pfix(1u << k); if (e.add(v)) nbas.push_back(v); } }
    printf("class %d: dim Fix = %d, dim N = %zu\n", cls, df, nbas.size()); fflush(stdout);
    if (df + (int)nbas.size() != 28) { fprintf(stderr, "dims\n"); return 1; }

    std::mutex mu; std::set<std::vector<u32>> seen; FILE* out = fopen(argv[2], "w");
    std::atomic<long> nB1(0), nOrbit(0), nPart(0), nApn(0);
    std::map<std::string, long> types;
    auto report = [&](u32* b8, const char* how) {
        u64 m[8]; for (int i = 0; i < 8; i++) m[i] = coef_to_mat(b8[i]);
        if (max_kernel(m, 8) != 1) return;
        int cnt[5]; rank_dist(m, 8, cnt); if (cnt[0]) return;
        canon8(b8); std::vector<u32> key(b8, b8 + 8);
        std::lock_guard<std::mutex> g(mu);
        if (!seen.insert(key).second) return;
        nApn++;
        char t[64]; snprintf(t, sizeof t, "ranks 8:%d 6:%d 4:%d 2:%d", cnt[4], cnt[3], cnt[2], cnt[1]); types[t]++;
        fprintf(out, "%s %s basis", t, how); for (int i = 0; i < 8; i++) fprintf(out, " %07x", b8[i]); fprintf(out, "\n"); fflush(out);
        if (cnt[2] >= 5) { printf("FOUND %s\n", t); fflush(stdout); }
    };
    auto kernel_ok = [&](const u32* fs, int n) { u64 m[8]; for (int i = 0; i < n; i++) m[i] = coef_to_mat(fs[i]); return max_kernel(m, n) <= 9 - n; };

    // orbit5 <class> <out> fixed: the other order-5 case, A fixing every rank-4 component (>= 5 of them). They lie in
    // Q = S8 & Fix, of dim 0, 4 or 8, so either S8 = Fix, or Q is 4-dim with >= 5 rank-4 forms and S8 = Q + F16 n.
    if (argc > 3 && !strcmp(argv[3], "fixed")) {
        u32 fb[8]; for (int i = 0; i < 8; i++) fb[i] = fix[i];
        report(fb, "fix8");
        std::vector<std::array<u32, 4>> Qs; long nQ = 0;
        for (int p0 = 0; p0 < 8; p0++) for (int p1 = p0 + 1; p1 < 8; p1++) for (int p2 = p1 + 1; p2 < 8; p2++) for (int p3 = p2 + 1; p3 < 8; p3++) {
            int P[4] = {p0, p1, p2, p3}; unsigned piv = (1u << p0) | (1u << p1) | (1u << p2) | (1u << p3);
            unsigned M[4]; for (int i = 0; i < 4; i++) M[i] = ((1u << P[i]) - 1) & ~piv;
            int bits = 0; for (int i = 0; i < 4; i++) bits += __popcnt(M[i]);
            for (u32 g = 0; g < (1u << bits); g++) { // spread g's bits over the free positions of the 4 rows
                unsigned r[4]; u32 gg = g;
                for (int i = 0; i < 4; i++) { r[i] = 1u << P[i]; for (int b = 0; b < 8; b++) if (M[i] >> b & 1) { if (gg & 1) r[i] |= 1u << b; gg >>= 1; } }
                std::array<u32, 4> q; for (int i = 0; i < 4; i++) { u32 v = 0; for (int b = 0; b < 8; b++) if (r[i] >> b & 1) v ^= fix[b]; q[i] = v; }
                nQ++;
                int r4 = 0; u32 cur = 0; for (unsigned h = 1; h < 16; h++) { cur ^= q[_tzcnt_u32(h)]; if (rank_mat(coef_to_mat(cur)) == 4) r4++; }
                if (r4 < 5 || !partition_ok(q.data(), 4)) continue;
                Qs.push_back(q);
            }
        }
        printf("fixed: %ld subspaces Q of Fix, %zu with >= 5 rank-4 forms and partition ok\n", nQ, Qs.size()); fflush(stdout);
        u32 nb[5][4]; int m = 0; { Ech28 e; for (u32 v : nbas) { if (m == 5) break; if (!e.reduce(v)) continue; nb[m][0] = v; for (int i = 1; i < 4; i++) nb[m][i] = act(nb[m][i - 1]); for (int i = 0; i < 4; i++) e.add(nb[m][i]); m++; } }
        if (m != 5) { fprintf(stderr, "F16 basis of N\n"); return 1; }
        auto smul = [&](unsigned c, int j) { u32 v = 0; for (int i = 0; i < 4; i++) if (c >> i & 1) v ^= nb[j][i]; return v; };
        // optional: orbit5 <c> <out> fixed <qfirst> <qlast>  (inclusive range of Q indices, for chunked/resumable runs)
        long qfirst = argc > 4 ? atol(argv[4]) : 0, qlast = argc > 5 ? atol(argv[5]) : (long)Qs.size() - 1;
        if (qlast >= (long)Qs.size()) qlast = (long)Qs.size() - 1;
        std::atomic<long> qdone(0); double tq0 = omp_get_wtime();
        #pragma omp parallel for schedule(dynamic, 1)
        for (long qi = qfirst; qi <= qlast; qi++) {
            struct Done { std::atomic<long>& d; long tot; double t0; ~Done() { long k = ++d; if (k % 100 == 0 || k == tot) { printf("Q done %ld/%ld %.0f s\n", k, tot, omp_get_wtime() - t0); fflush(stdout); } } } done_{qdone, qlast - qfirst + 1, tq0};
            for (int lead = 0; lead < 5; lead++) for (u32 t = 0; t < (1u << (4 * (4 - lead))); t++) {
                u32 v = nb[lead][0]; for (int j = 0; j < 4 - lead; j++) v ^= smul(t >> (4 * j) & 15, lead + 1 + j);
                u32 b8[8] = {Qs[qi][0], Qs[qi][1], Qs[qi][2], Qs[qi][3], v, act(v), act(act(v)), act(act(act(v)))};
                report(b8, "fix4");
            }
        }
        printf("class %d fixed (Q %ld..%ld): APN spaces %ld\n", cls, qfirst, qlast, nApn.load());
        for (auto& kv : types) printf("  %s: %ld\n", kv.first.c_str(), kv.second);
        fclose(out); return 0;
    }

    auto process = [&](u32 b1, bool validate) {
        u32 orb[5]; orb[0] = b1; for (int i = 1; i < 5; i++) orb[i] = act(orb[i - 1]);
        if (act(orb[4]) != b1) { fprintf(stderr, "orbit\n"); exit(1); }
        if (!validate) {
            if (*std::min_element(orb, orb + 5) != b1) return; // one B1 per orbit
            nB1++;
            Set256 R[5]; for (int i = 0; i < 5; i++) radical(coef_to_mat(orb[i]), R[i]);
            for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) if (R[i].meets(R[j])) return;
            nOrbit++;
        }
        if (!g_enumtest && !partition_ok(orb, 5)) return;
        nPart++;
        u32 f = b1 ^ orb[1] ^ orb[2] ^ orb[3] ^ orb[4]; u32 n = b1 ^ f; // projection onto Fix is the orbit sum
        u32 Fn[4] = {n, act(n), act(act(n)), act(act(act(n)))};
        // k = 1: S8 = F16 n + Q, Q 4-dim in Fix containing f
        {
            // Fix basis fb with f first; Q = <f> + (Q & H), H = span(fb[lo..7]), so each Q is one subspace of H,
            // enumerated once in reduced echelon form (pivot = top bit, rows in increasing pivot order), pruned per row.
            std::vector<u32> fb; { Ech28 t; if (f) { t.add(f); fb.push_back(f); } for (u32 v : fix) if (t.add(v)) fb.push_back(v); }
            int lo = f ? 1 : 0;
            auto comb = [&](unsigned r) { u32 v = 0; for (int i = 0; i < 8; i++) if (r >> i & 1) v ^= fb[i]; return v; };
            std::vector<u32> cur(Fn, Fn + 4); if (f) cur.push_back(f);
            std::function<void(int, unsigned)> rec = [&](int pmin, unsigned piv) {
                if ((int)cur.size() == 8) {
                    u32 b8[8]; for (int i = 0; i < 8; i++) b8[i] = cur[i];
                    if (g_enumtest) { canon8(b8); std::lock_guard<std::mutex> g(mu); g_enum.insert(std::vector<u32>(b8, b8 + 8)); g_enumcount++; return; }
                    report(b8, "k1"); return;
                }
                int left = 8 - (int)cur.size();
                for (int p = pmin; p <= 8 - left; p++) {
                    unsigned M = ((1u << p) - 1) & ~((1u << lo) - 1) & ~piv;
                    for (unsigned s = M;; s = (s - 1) & M) {
                        cur.push_back(comb((1u << p) | s));
                        if (g_enumtest || kernel_ok(cur.data(), (int)cur.size())) rec(p + 1, piv | 1u << p);
                        cur.pop_back();
                        if (!s) break;
                    }
                }
            };
            rec(lo, 0);
        }
        // k = 2: f = 0, S8 = F16 n + F16 n', n' in N
        // n' runs over the projective points of N / F16 n (F16-dim 4: 4369 points), so each space comes once.
        if (!f) {
            Ech28 e; for (int i = 0; i < 4; i++) e.add(Fn[i]);
            u32 nb[4][4]; int m = 0; // F16-basis n2..n5 of a complement, with its A-powers
            for (u32 v : nbas) { if (m == 4) break; if (!e.reduce(v)) continue; nb[m][0] = v; for (int i = 1; i < 4; i++) nb[m][i] = act(nb[m][i - 1]); for (int i = 0; i < 4; i++) e.add(nb[m][i]); m++; }
            if (m != 4) { fprintf(stderr, "F16 basis\n"); exit(1); }
            auto smul = [&](unsigned c, int j) { u32 v = 0; for (int i = 0; i < 4; i++) if (c >> i & 1) v ^= nb[j][i]; return v; };
            for (int lead = 0; lead < 4; lead++) {
                int tail = 3 - lead;
                for (unsigned t = 0; t < (1u << (4 * tail)); t++) {
                    u32 v = nb[lead][0]; for (int j = 0; j < tail; j++) v ^= smul(t >> (4 * j) & 15, lead + 1 + j);
                    u32 b8[8] = {Fn[0], Fn[1], Fn[2], Fn[3], v, act(v), act(act(v)), act(act(act(v)))};
                    report(b8, "k2");
                }
            }
        }
    };
    if (argc > 5 && argv[3][0] == 's') { // sample: c = r mod M only, timed, to estimate a full run
        long M = atol(argv[4]), r = atol(argv[5]); double t0 = omp_get_wtime();
        #pragma omp parallel for schedule(dynamic, 64)
        for (long c = r ? r : M; c < (1L << 28); c += M) { if (rank_mat(coef_to_mat((u32)c)) == 4) process((u32)c, false); }
        printf("sample 1/%ld (r=%ld): %.1f s wall\n", M, r, omp_get_wtime() - t0);
    } else if (argc > 3) { // validation: one given B1, no radical-orbit condition
        u32 b1 = (u32)strtoul(argv[3], 0, 16); process(b1, true);
        if (g_enumtest) printf("enumtest: k1 leaves %ld, distinct %zu\n", g_enumcount, g_enum.size());
    } else {
        #pragma omp parallel for schedule(dynamic, 1 << 14)
        for (long c = 1; c < (1L << 28); c++) { if (rank_mat(coef_to_mat((u32)c)) == 4) process((u32)c, false); }
    }
    printf("class %d: B1 orbits %ld, radical orbits complementary %ld, partition ok %ld, APN spaces %ld\n", cls, nB1.load(), nOrbit.load(), nPart.load(), nApn.load());
    for (auto& kv : types) printf("  %s: %ld\n", kv.first.c_str(), kv.second);
    fclose(out);
    return 0;
}
