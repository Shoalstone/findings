// classify8: EA-invariants of quadratic (8,8)-APN functions from apn8 output lines, to count inequivalent ones.
// For quadratic APN F, EA-equivalence = CCZ-equivalence, and EA-equivalent F give EA-equivalent ortho-derivatives
// (Canteaut-Couvreur-Perrin 2022, Prop. 36). So the ortho-derivative's differential spectrum and extended Walsh
// spectrum are CCZ-invariants of F: distinct invariant -> inequivalent. Equal invariant proves nothing.
// Usage: classify8 <in> <out>   out lines: '<invariant hash> <rank type> <basis>' for the first function of each
// distinct invariant, then a summary per rank type on stdout. Rowen (Claude Opus 5.5), 2026-10-24.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <omp.h>
typedef uint64_t u64; typedef uint32_t u32;

static int KI[28], KJ[28];

static void lut(const u32* basis, unsigned* F) {
    for (unsigned x = 0; x < 256; x++) {
        unsigned y = 0;
        for (int s = 0; s < 8; s++) {
            int v = 0;
            for (int k = 0; k < 28; k++) if (basis[s] >> k & 1) v ^= (x >> KI[k] & 1) & (x >> KJ[k] & 1);
            y |= (unsigned)v << s;
        }
        F[x] = y;
    }
}

static bool ortho(const unsigned* F, unsigned* pi) {
    pi[0] = 0;
    for (unsigned a = 1; a < 256; a++) {
        unsigned d[8]; for (int j = 0; j < 8; j++) { unsigned x = 1u << j; d[j] = F[x ^ a] ^ F[x] ^ F[a] ^ F[0]; }
        unsigned found = 0; int n = 0;
        for (unsigned c = 1; c < 256; c++) { bool ok = true; for (int j = 0; j < 8 && ok; j++) if (__popcnt(c & d[j]) & 1) ok = false; if (ok) { found = c; n++; } }
        if (n != 1) return false;
        pi[a] = found;
    }
    return true;
}

static u64 fnv(u64 h, u64 v) { for (int i = 0; i < 8; i++) { h ^= (v >> (8 * i)) & 0xFF; h *= 0x100000001B3ull; } return h; }

static std::string invariant(const unsigned* pi) {
    // differential spectrum of pi: histogram of N(a,b), a != 0
    std::map<int, int> ds; static thread_local int cnt[256];
    for (unsigned a = 1; a < 256; a++) {
        memset(cnt, 0, sizeof cnt);
        for (unsigned x = 0; x < 256; x++) cnt[pi[x] ^ pi[x ^ a]]++;
        for (int b = 0; b < 256; b++) ds[cnt[b]]++;
    }
    // extended Walsh spectrum: histogram of |W(u,v)|, v != 0
    std::map<int, int> ws; int w[256];
    for (unsigned v = 1; v < 256; v++) {
        for (unsigned x = 0; x < 256; x++) w[x] = (__popcnt(v & pi[x]) & 1) ? -1 : 1;
        for (int h = 1; h < 256; h <<= 1) for (int i = 0; i < 256; i += 2 * h) for (int j = i; j < i + h; j++) { int a = w[j], b = w[j + h]; w[j] = a + b; w[j + h] = a - b; }
        for (int u = 0; u < 256; u++) ws[abs(w[u])]++;
    }
    std::string s = "D";
    for (auto& kv : ds) s += " " + std::to_string(kv.first) + ":" + std::to_string(kv.second);
    s += " W";
    for (auto& kv : ws) s += " " + std::to_string(kv.first) + ":" + std::to_string(kv.second);
    return s;
}

static int walsh_rank_dist(const unsigned* F, int* cnt) { // cnt[r/2]; returns 0 if some component isn't quadratic-shaped
    for (int i = 0; i < 5; i++) cnt[i] = 0;
    int w[256];
    for (unsigned b = 1; b < 256; b++) {
        for (unsigned x = 0; x < 256; x++) w[x] = (__popcnt(b & F[x]) & 1) ? -1 : 1;
        for (int h = 1; h < 256; h <<= 1) for (int i = 0; i < 256; i += 2 * h) for (int j = i; j < i + h; j++) { int a = w[j], c = w[j + h]; w[j] = a + c; w[j + h] = a - c; }
        int m = 0; for (int u = 0; u < 256; u++) if (abs(w[u]) > m) m = abs(w[u]);
        int lg = 0; while ((1 << (lg + 1)) <= m) lg++;
        int r = 2 * (8 - lg); if (r < 0 || r > 8 || (1 << lg) != m) return 0;
        cnt[r / 2]++;
    }
    return 1;
}

struct Item { int fi; std::string src; unsigned F[256]; };
struct Seen { std::string src; unsigned mask; };

// Streams the inputs in batches (the 3.8M file is 3.4 GB). --min4 k keeps only types with >= k rank-4 components
// or a rank-2 one; the totals per type are counted for every type regardless.
int main(int argc, char** argv) {
    { int k = 0; for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) { KI[k] = i; KJ[k] = j; k++; } }
    int min4 = 0; int a0 = 1;
    if (argc > 2 && !strcmp(argv[1], "--min4")) { min4 = atoi(argv[2]); a0 = 3; }
    if (argc < a0 + 2) { fprintf(stderr, "usage: classify8 [--min4 k] <out> <in>...   (apn8 lines or {256 LUT values} lines)\n"); return 2; }
    const char* outp = argv[a0]; int f0 = a0 + 1, nf = argc - f0;
    std::unordered_set<std::string> seenb;
    std::map<std::string, std::unordered_map<std::string, Seen>> per; // type -> invariant -> (first src, file mask)
    std::map<std::string, long> total; long nread = 0, bad = 0;
    static char line[1 << 16];
    const size_t BATCH = 1 << 16;
    std::vector<Item> items; items.reserve(BATCH);
    auto flush = [&]() {
        size_t n = items.size(); std::vector<std::string> type(n), inv(n); std::vector<char> okv(n);
        #pragma omp parallel for schedule(dynamic, 16)
        for (long i = 0; i < (long)n; i++) {
            int cnt[5]; unsigned pi[256];
            okv[i] = walsh_rank_dist(items[i].F, cnt) && ortho(items[i].F, pi);
            if (!okv[i]) continue;
            char t[96]; snprintf(t, sizeof t, "ranks 8:%d 6:%d 4:%d 2:%d", cnt[4], cnt[3], cnt[2], cnt[1]); type[i] = t;
            if (cnt[2] >= min4 || cnt[1] > 0) inv[i] = invariant(pi);
        }
        for (size_t i = 0; i < n; i++) {
            nread++;
            if (!okv[i]) { bad++; continue; }
            total[type[i]]++;
            if (inv[i].empty()) continue;
            auto& m = per[type[i]];
            auto it = m.find(inv[i]);
            if (it == m.end()) m.emplace(inv[i], Seen{items[i].src, 1u << items[i].fi}); else it->second.mask |= 1u << items[i].fi;
        }
        items.clear();
    };
    for (int fi = 0; fi < nf; fi++) {
        FILE* f = fopen(argv[f0 + fi], "r"); if (!f) { fprintf(stderr, "can't open %s\n", argv[f0 + fi]); return 1; }
        long ln = 0;
        while (fgets(line, sizeof line, f)) {
            ln++;
            Item it; it.fi = fi; it.src = std::string(argv[f0 + fi]) + ":" + std::to_string(ln);
            char* p = strstr(line, "basis");
            if (p) {
                u32 b[8];
                if (sscanf(p + 5, "%x %x %x %x %x %x %x %x", &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7]) != 8) continue;
                if (!seenb.insert(std::string(p + 5)).second) continue;
                lut(b, it.F);
                it.src += " basis" + std::string(p + 5); while (!it.src.empty() && (it.src.back() == '\n' || it.src.back() == '\r')) it.src.pop_back();
            } else if (line[0] == '{' || line[0] == '[') {
                char* q = line + 1; int n = 0;
                while (n < 256) { char* e; long v = strtol(q, &e, 10); if (e == q) break; it.F[n++] = (unsigned)v; q = e; while (*q == ',' || *q == ' ') q++; }
                if (n != 256) continue;
            } else continue;
            items.push_back(std::move(it));
            if (items.size() == BATCH) { flush(); if (nread % (1 << 20) == 0) fprintf(stderr, "%ld read\n", nread); }
        }
        fclose(f);
    }
    flush();
    FILE* o = fopen(outp, "w");
    printf("%ld functions read, %ld not quadratic APN-shaped\n", nread, bad);
    for (auto& kv : total) {
        auto pit = per.find(kv.first);
        printf("  %s: %ld functions", kv.first.c_str(), kv.second);
        if (pit == per.end()) { printf("\n"); continue; }
        std::vector<std::pair<std::string, Seen>> sorted(pit->second.begin(), pit->second.end());
        std::sort(sorted.begin(), sorted.end(), [](auto& x, auto& y) { return x.first < y.first; });
        std::map<unsigned, int> bymask;
        for (auto& iv : sorted) {
            bymask[iv.second.mask]++;
            u64 h = 0xCBF29CE484222325ull; for (char c : iv.first) h = fnv(h, (unsigned char)c);
            fprintf(o, "%016llx files=%x %s src %s\n", h, iv.second.mask, kv.first.c_str(), iv.second.src.c_str());
        }
        printf(", %zu distinct invariants;", sorted.size());
        for (auto& bm : bymask) { printf(" in{"); for (int k = 0; k < nf; k++) if (bm.first >> k & 1) printf("%d", k); printf("}:%d", bm.second); }
        printf("\n");
    }
    fclose(o);
    for (int k = 0; k < nf; k++) printf("  file %d = %s\n", k, argv[f0 + k]);
    return 0;
}
