// review_odd_engine: reviewer's positive control (Claude Opus 5.5, 2026-10-24) for c3k4.cpp's extension Engine at
// every r = 1..4 (search2 uses r = 1, 2; searchf 2 uses r = 2; the others 3, 4; only r = 3, 4 had positive controls).
// The Engine is c3k4.cpp's own code, copied verbatim (lines 14-205) into c3k4_engine.inc by the build step:
//   sed -n 14,205p src/c3k4.cpp > src/c3k4_engine.inc   (or the Python line in README.md)
// For known (6,1) APN spaces S8 = P + Q from results/c3k4-ff.txt: seed = P + (6 - r) random vectors of Q, then the
// Engine must give S8 back. usage: review_odd_engine <c3k4-ff.txt> [max spaces] [spaces also tried at r = 4, default 2]
#define main apn8_main
#include "apn8.cpp"
#undef main
#include "c3k4_engine.inc"

int main(int argc, char** argv) {
    init_index(); setup();
    std::vector<Found> v; if (argc < 2 || !load_seeds(argv[1], v)) { fprintf(stderr, "can't load\n"); return 2; }
    int maxn = argc > 2 ? atoi(argv[2]) : 1 << 30;
    int max4 = argc > 3 ? atoi(argv[3]) : 2; // r = 4 costs minutes per space
    setvbuf(stdout, nullptr, _IONBF, 0);
    Rng R(777);
    long tried[5] = {0}, hit[5] = {0}, bad = 0;
    // order: all spaces with 3 rank-4 first, then the rest
    std::vector<int> order; for (int i = 0; i < (int)v.size(); i++) if (v[i].cnt[2] == 3) order.push_back(i);
    for (int i = 0; i < (int)v.size(); i++) if (v[i].cnt[2] != 3) order.push_back(i);
    int done = 0;
    for (int idx : order) {
        if (done++ >= maxn) break;
        const u32* b = v[idx].basis;
        // P = S8 & N, Q = S8 & Fix (projections of the basis; both lie in S8 since S8 is A-invariant)
        Ech28 eq, ep; std::vector<u32> Q, P;
        for (int i = 0; i < 8; i++) { u32 f = fixpart(b[i]), n = b[i] ^ f; if (f && eq.add(f)) Q.push_back(f); if (n && ep.add(n)) P.push_back(n); }
        if (Q.size() != 6 || P.size() != 2) { bad++; continue; }
        u32 target[8]; memcpy(target, b, sizeof target); canon8(target);
        for (int r = 1; r <= 4; r++) {
            if (r == 4 && done > max4) continue;
            // random basis of Q: random invertible combination
            std::vector<u32> qb;
            for (;;) { qb.clear(); Ech28 e; for (int t = 0; t < 6; t++) { u32 c = 0; unsigned g = (unsigned)R.below(63) + 1; for (int i = 0; i < 6; i++) if (g >> i & 1) c ^= Q[i]; if (e.add(c)) qb.push_back(c); } if (qb.size() == 6) break; }
            std::vector<u32> seed = P; for (int t = 0; t < 6 - r; t++) seed.push_back(qb[t]);
            Engine E; tried[r]++;
            if (!E.init(seed, r)) { printf("space %d r=%d: init failed\n", idx, r); continue; }
            bool found = false;
            E.run([&](u32* b8) { u32 c[8]; memcpy(c, b8, sizeof c); canon8(c); if (!memcmp(c, target, sizeof c)) found = true; });
            if (found) hit[r]++; else printf("space %d r=%d: NOT re-found\n", idx, r);
        }
    }
    printf("spaces skipped (not (6,1)): %ld\n", bad);
    for (int r = 1; r <= 4; r++) printf("r=%d: %ld tried, %ld re-found\n", r, tried[r], hit[r]);
    return 0;
}
