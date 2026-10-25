// extend6: every APN 8-space containing a given 6-space, with no caps (apn8's search_from_s6: Gray code over the
// 22-dim complement for the 7th coordinate, the 8th by the switching linear system, all of its solutions).
// Used for sixteen.cpp's SHORT cases (a span of rank-4 forms covering all their radicals in fewer than 8 dims).
// usage: extend6 <file with 'SHORT dim 6 ... basis h1..h6' lines> [min rank-4 to print, default 0]
// Rowen (Claude Opus 5.5), 2026-10-24.
#define main apn8_main
#include "apn8.cpp"
#undef main

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: extend6 <file> [min4]\n"); return 2; }
    int min4 = argc > 2 ? atoi(argv[2]) : 0;
    init_index();  // apn8's pair tables (its own main does this)
    g_dn_cap = 64;
    FILE* f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }
    char line[8192];
    long ncase = 0, nskip = 0;
    unsigned long long tot_apn = 0;
    while (fgets(line, sizeof line, f)) {
        if (strncmp(line, "SHORT", 5)) continue;
        int dim = 0;
        const char* p = strstr(line, "dim ");
        if (p) dim = atoi(p + 4);
        const char* b = strstr(line, "basis");
        if (!b || dim != 6) { printf("skip (dim %d): %s", dim, line); nskip++; continue; }
        b += 5;
        u32 s6[6]; int n = 0;
        while (n < 6) { char* end; unsigned long v = strtoul(b, &end, 16); if (end == b) break; s6[n++] = (u32)v; b = end; }
        if (n != 6) { printf("bad line: %s", line); nskip++; continue; }
        Stats st; std::vector<Found> out;
        search_from_s6(s6, st, out, -1, ~0ull);
        ncase++;
        tot_apn += st.apn;
        int best = 0; for (auto& fd : out) best = std::max(best, fd.cnt[2]);
        printf("case %ld: w_pass %llu solved %llu apn %llu (max rank-4 %d)\n", ncase, (unsigned long long)st.w_pass,
               (unsigned long long)st.solved, (unsigned long long)st.apn, best);
        for (auto& fd : out) if (fd.cnt[2] >= min4) {
            printf("  ranks 8:%d 6:%d 4:%d 2:%d basis", fd.cnt[4], fd.cnt[3], fd.cnt[2], fd.cnt[1]);
            for (int i = 0; i < 8; i++) printf(" %07x", fd.basis[i]);
            printf("\n");
        }
        fflush(stdout);
    }
    printf("cases %ld, skipped %ld, APN extensions %llu\n", ncase, nskip, tot_apn);
    return 0;
}
