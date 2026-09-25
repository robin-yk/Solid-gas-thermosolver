/* Transition-matrix Monte Carlo for oxygen vacancies on the rutile O sublattice.
 *
 * Energy: pairwise additive, E = sum_{i<j} J_ij s_i s_j (ZHA2017 pair energy
 * within the cutoff), site energy 0 (bulk reference).
 * Constraint: every vacancy holds two Ti3+ among its three Ti neighbours and a
 * Ti holds at most one extra electron. A configuration is allowed only if such
 * an assignment exists (bipartite b-matching, kept by augmenting paths).
 * Moves: grand-canonical flip of a random site, and vacancy relocation.
 * Collected: unbiased acceptance for N -> N +/- 1, which gives the macrostate
 * distribution P(N) at mu0, hence ln Q(N) for every N.
 */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int nO, nTi;
    const int *nb_start, *nb_idx;
    const double *nb_J;
    const int *oti;
    int *occ, *owner, *slot, *vis, *vac, *vpos;
    double *h;
    int stamp, N;
    int *log_t, *log_prev, nlog;
    int *slog_v, *slog_k, *slog_prev, nslog;
    uint64_t rng;
} Sys;

static double urand(Sys *s) {
    s->rng ^= s->rng >> 12; s->rng ^= s->rng << 25; s->rng ^= s->rng >> 27;
    return ((s->rng * 2685821657736338717ULL) >> 11) * (1.0 / 9007199254740992.0);
}

static void set_owner(Sys *s, int t, int v) {
    s->log_t[s->nlog] = t; s->log_prev[s->nlog] = s->owner[t]; s->nlog++;
    s->owner[t] = v;
}

static void set_slot(Sys *s, int v, int k, int t) {
    s->slog_v[s->nslog] = v; s->slog_k[s->nslog] = k; s->slog_prev[s->nslog] = s->slot[2 * v + k]; s->nslog++;
    s->slot[2 * v + k] = t;
}

static void rollback(Sys *s) {
    while (s->nslog) { s->nslog--; s->slot[2 * s->slog_v[s->nslog] + s->slog_k[s->nslog]] = s->slog_prev[s->nslog]; }
    while (s->nlog) { s->nlog--; s->owner[s->log_t[s->nlog]] = s->log_prev[s->nlog]; }
}

static int owns(Sys *s, int v, int t) { return s->slot[2 * v] == t || s->slot[2 * v + 1] == t; }

/* u gives up t_old and takes another Ti, possibly displacing others. */
static int augment(Sys *s, int u, int t_old) {
    for (int k = 0; k < 3; k++) {
        int t = s->oti[3 * u + k];
        if (s->vis[t] == s->stamp || owns(s, u, t)) continue;
        s->vis[t] = s->stamp;
        int w = s->owner[t];
        if (w < 0 || augment(s, w, t)) {
            set_owner(s, t, u);
            set_slot(s, u, s->slot[2 * u] == t_old ? 0 : 1, t);
            return 1;
        }
    }
    return 0;
}

/* x takes one more Ti into an empty slot. */
static int grab(Sys *s, int x) {
    int k_free = s->slot[2 * x] < 0 ? 0 : 1;
    for (int k = 0; k < 3; k++) {
        int t = s->oti[3 * x + k];
        if (s->vis[t] == s->stamp || owns(s, x, t)) continue;
        s->vis[t] = s->stamp;
        int w = s->owner[t];
        if (w < 0 || augment(s, w, t)) {
            set_owner(s, t, x);
            set_slot(s, x, k_free, t);
            return 1;
        }
    }
    return 0;
}

static int try_insert(Sys *s, int i) {
    s->nlog = s->nslog = 0;
    for (int r = 0; r < 2; r++) {
        s->stamp++;
        if (!grab(s, i)) { rollback(s); return 0; }
    }
    return 1;
}

static void flip_field(Sys *s, int i, double sign) {
    for (int p = s->nb_start[i]; p < s->nb_start[i + 1]; p++) s->h[s->nb_idx[p]] += sign * s->nb_J[p];
}

static void add_vac(Sys *s, int i) {
    s->occ[i] = 1; flip_field(s, i, 1.0);
    s->vpos[i] = s->N; s->vac[s->N++] = i;
}

static void remove_vac(Sys *s, int i) {
    s->occ[i] = 0; flip_field(s, i, -1.0);
    int last = s->vac[--s->N];
    s->vac[s->vpos[i]] = last; s->vpos[last] = s->vpos[i]; s->vpos[i] = -1;
    for (int k = 0; k < 2; k++) { int t = s->slot[2 * i + k]; if (t >= 0) s->owner[t] = -1; s->slot[2 * i + k] = -1; }
}

static double pair_J(Sys *s, int i, int j) {
    for (int p = s->nb_start[i]; p < s->nb_start[i + 1]; p++) if (s->nb_idx[p] == j) return s->nb_J[p];
    return 0.0;
}

/* Widom sum over every empty site that could take a vacancy:
 * W = sum_i exp(-beta h_i), so Q(N+1)/Q(N) = <W>_N / (N+1). */
static double widom_log(Sys *s, double beta, double *wsite) {
    double hmin = 1e300;
    for (int i = 0; i < s->nO; i++) if (!s->occ[i] && s->h[i] < hmin) hmin = s->h[i];
    double sum = 0.0;
    for (int i = 0; i < s->nO; i++) {
        double w = 0.0;
        if (!s->occ[i]) {
            int nfree = 0;
            for (int k = 0; k < 3; k++) nfree += s->owner[s->oti[3 * i + k]] < 0;
            int ok = nfree >= 2;
            if (!ok) { ok = try_insert(s, i); if (ok) rollback(s); }
            if (ok) w = exp(-beta * (s->h[i] - hmin));
        }
        if (wsite) wsite[i] = w;
        sum += w;
    }
    return log(sum) - beta * hmin;
}

static void setup(Sys *s, int nO, int nTi, const int *nb_start, const int *nb_idx, const double *nb_J,
                  const int *oti, int *occ, int *owner, int *slot, double *h, uint64_t seed) {
    memset(s, 0, sizeof *s);
    s->nO = nO; s->nTi = nTi; s->nb_start = nb_start; s->nb_idx = nb_idx; s->nb_J = nb_J; s->oti = oti;
    s->occ = occ; s->owner = owner; s->slot = slot; s->h = h; s->rng = seed ? seed : 88172645463325252ULL;
    s->vis = calloc(nTi, sizeof(int)); s->vac = malloc(nO * sizeof(int)); s->vpos = malloc(nO * sizeof(int));
    int cap = 4 * nO + 16;
    s->log_t = malloc(cap * sizeof(int)); s->log_prev = malloc(cap * sizeof(int));
    s->slog_v = malloc(cap * sizeof(int)); s->slog_k = malloc(cap * sizeof(int)); s->slog_prev = malloc(cap * sizeof(int));
    for (int i = 0; i < nO; i++) { s->vpos[i] = -1; if (occ[i]) { s->vpos[i] = s->N; s->vac[s->N++] = i; } }
}

static void teardown(Sys *s) {
    free(s->vis); free(s->vac); free(s->vpos); free(s->log_t); free(s->log_prev);
    free(s->slog_v); free(s->slog_k); free(s->slog_prev);
}

/* Canonical sampling at fixed N: relocation of a vacancy to a random site or
 * to a random site within the cutoff. Every sample_every moves, the log Widom
 * sum is added (log-sum-exp) into *lnW_acc and *nsamp is incremented. */
int canon_run(int nO, int nTi, const int *nb_start, const int *nb_idx, const double *nb_J,
              const int *oti, int *occ, int *owner, int *slot, double *h,
              double beta, long nsteps, int sample_every, uint64_t seed,
              double *lnW_acc, long *nsamp) {
    Sys s;
    setup(&s, nO, nTi, nb_start, nb_idx, nb_J, oti, occ, owner, slot, h, seed);
    for (long step = 0; step < nsteps; step++) {
        if (s.N > 0) {
            int i = s.vac[(int)(urand(&s) * s.N)];
            int j;
            if (urand(&s) < 0.5) j = (int)(urand(&s) * nO);
            else {
                int a = nb_start[i], b = nb_start[i + 1];
                j = nb_idx[a + (int)(urand(&s) * (b - a))];
            }
            if (!occ[j]) {
                double dE = -h[i] + h[j] - pair_J(&s, i, j);
                int t0 = slot[2 * i], t1 = slot[2 * i + 1];
                remove_vac(&s, i);
                if (try_insert(&s, j) && (dE <= 0 || urand(&s) < exp(-beta * dE))) {
                    add_vac(&s, j);
                } else {
                    rollback(&s);
                    owner[t0] = i; owner[t1] = i; slot[2 * i] = t0; slot[2 * i + 1] = t1;
                    add_vac(&s, i);
                }
            }
        }
        if ((step + 1) % sample_every == 0) {
            double lw = widom_log(&s, beta, NULL);
            double m = lw > *lnW_acc ? lw : *lnW_acc;
            *lnW_acc = (*nsamp == 0) ? lw : m + log(exp(*lnW_acc - m) + exp(lw - m));
            (*nsamp)++;
        }
    }
    teardown(&s);
    return s.N;
}

/* Add one vacancy by heat bath over all feasible empty sites. */
int grow(int nO, int nTi, const int *nb_start, const int *nb_idx, const double *nb_J,
         const int *oti, int *occ, int *owner, int *slot, double *h, double beta, uint64_t seed) {
    Sys s;
    setup(&s, nO, nTi, nb_start, nb_idx, nb_J, oti, occ, owner, slot, h, seed);
    double *w = malloc(nO * sizeof(double));
    widom_log(&s, beta, w);
    double tot = 0.0;
    for (int i = 0; i < nO; i++) tot += w[i];
    int pick = -1;
    if (tot > 0) {
        double r = urand(&s) * tot, c = 0.0;
        for (int i = 0; i < nO; i++) { c += w[i]; if (w[i] > 0 && c >= r) { pick = i; break; } }
        if (pick >= 0 && try_insert(&s, pick)) add_vac(&s, pick); else pick = -1;
    }
    free(w);
    int N = s.N;
    teardown(&s);
    return pick < 0 ? -1 : N;
}
