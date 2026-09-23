#!/usr/bin/env python3
"""Rigorous LOWER bound L(t) on the fractional packing number nu_f(t) by a direct packing LP
with column generation (poses) and constraint generation (points), certified EXACTLY at the
vertices of the arrangement of the support squares.

    nu_f(t) = sup mu(Pi)   s.t.  mu({Q : p in Q}) <= 1 for every p in C=[0,t]^2,   mu >= 0

Pi = placements of a CLOSED unit square inside C (the convention of certificates/FORMAT.md and
verify/).  The measure is D4-symmetrised: a pose r = (cx, cy, theta) with mass mu_r puts mu_r/8 on
each of its 8 dihedral images, so the LP's coverage at a point is the ORBIT-AVERAGED coverage and
the constraint set can live in the fundamental domain {0 <= x <= y <= t/2}.

Certification (the only rigorous step).  For a finite set of closed squares the coverage function
cov(p) = sum_r mu_r [p in Q_r] is piecewise constant on the arrangement of the squares' edges and,
because the squares are CLOSED, cov is upper semicontinuous: cov at a vertex of the arrangement is
>= cov on every open face/edge whose closure contains it.  Hence

    max_{p in C} cov(p) = max over V,   V = {square corners} u {pairwise edge intersections} u
                                            {container corners}

(a face of the arrangement inside C is a convex polygon whose vertices are in V; a square inside C
meets the container boundary only at its own corners or along an edge, whose endpoints are corners).
We compute the candidates of V in floating point and evaluate cov at each with a tolerance
tol = 1e-8 in the square's own frame ( |R(p)-u|_inf <= 1/2 + tol ), which can only OVER-count:
every square containing the true vertex v contains the tol-dilation of the computed v' (|v'-v| is
below 1e-9: pose angles are snapped to multiples of 1e-5 rad, so any two non-parallel edge lines
meet at an angle >= 1e-5 rad and the 2x2 solve is well conditioned; pairs with |det| < 1e-6 are
exactly parallel and contribute no vertex).  So M := max_{v' in V'} cov_tol(v') >= max_C cov and
mu/M is a feasible packing measure:   L = sum(mu)/M <= nu_f(t).   Every support pose is checked
to be admissible (closed unit square inside C with margin 1e-9).

Everything else (LP, pricing, which points are rows) is heuristic and only decides WHICH measure
gets certified.  Pricing: the dual y of the point constraints is a D4-symmetric weighted point set
(a certificate file), and an improving column is a closed unit square capturing < 1 of it; we use
the exact Rust verifier (topk worst placements per angle bin) when the atom count allows, and a
fast raster pricer (2-D prefix sums over a rotated raster, then exact local refinement) otherwise.

Usage
    python3 search/packing_dual.py T TAG [--rounds R] [--time SEC] [--seed-pitch H] [--seed-dth DEG]
                                          [--warm runs/dual_OTHER_support.txt] [--threads K] ...
    -> runs/dual_TAG.log, runs/dual_TAG.json, runs/dual_TAG_support.txt
"""
import sys, os, math, time, json, subprocess, ctypes, hashlib, argparse
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
VERIFY = os.path.join(REPO, 'verify', 'target', 'release', 'verify')
RUNS = os.path.join(REPO, 'runs')
TOL = 1e-8            # containment tolerance (square frame), used identically in LP rows and certification
ADM = 1e-9            # admissibility margin: cx in [wid/2 + ADM, t - wid/2 - ADM]
ANG_UNIT = 1e-5       # pose angles are multiples of this (radians)

# ============================================================================== C kernels
CSRC = r'''
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
typedef struct { double cx, cy, ct, st, u0, u1, w; int col; } Sq;
typedef struct { int nb; double cell; int *start; int *idx; } Grid;

static void build_grid(const Sq *sq, int n, double s, double cell, Grid *g) {
    int nb = (int)(s / cell) + 3; g->nb = nb; g->cell = cell;
    g->start = (int*)calloc((size_t)nb * nb + 1, sizeof(int));
    g->idx = (int*)malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    int *cnt = (int*)calloc((size_t)nb * nb, sizeof(int));
    for (int i = 0; i < n; i++) { int bx = (int)(sq[i].cx / cell), by = (int)(sq[i].cy / cell);
        if (bx < 0) bx = 0; if (by < 0) by = 0; if (bx >= nb) bx = nb - 1; if (by >= nb) by = nb - 1; cnt[bx * nb + by]++; }
    for (int b = 0; b < nb * nb; b++) g->start[b + 1] = g->start[b] + cnt[b];
    memset(cnt, 0, sizeof(int) * (size_t)nb * nb);
    for (int i = 0; i < n; i++) { int bx = (int)(sq[i].cx / cell), by = (int)(sq[i].cy / cell);
        if (bx < 0) bx = 0; if (by < 0) by = 0; if (bx >= nb) bx = nb - 1; if (by >= nb) by = nb - 1;
        int b = bx * nb + by; g->idx[g->start[b] + cnt[b]++] = i; }
    free(cnt);
}
static void free_grid(Grid *g) { free(g->start); free(g->idx); }

static inline void fill_sq(Sq *q, const double *a, int i, int has_w) {
    const double *p = a + (size_t)i * (has_w ? 6 : 5);
    q->cx = p[0]; q->cy = p[1]; q->ct = p[2]; q->st = p[3]; q->col = (int)p[4]; q->w = has_w ? p[5] : 0.0;
    q->u0 = q->cx * q->ct + q->cy * q->st; q->u1 = -q->cx * q->st + q->cy * q->ct;
}

/* coverage (tolerant) of the weighted squares at point (x,y) */
static inline double cov_at(const Sq *sq, const Grid *g, double x, double y, double tol) {
    double r = 0.7072 + tol; int nb = g->nb; double c = g->cell;
    int bx0 = (int)floor((x - r) / c), bx1 = (int)floor((x + r) / c), by0 = (int)floor((y - r) / c), by1 = (int)floor((y + r) / c);
    if (bx0 < 0) bx0 = 0; if (by0 < 0) by0 = 0; if (bx1 >= nb) bx1 = nb - 1; if (by1 >= nb) by1 = nb - 1;
    double tot = 0.0;
    for (int bx = bx0; bx <= bx1; bx++) for (int by = by0; by <= by1; by++) {
        int b = bx * nb + by;
        for (int k = g->start[b]; k < g->start[b + 1]; k++) { const Sq *q = &sq[g->idx[k]];
            double q0 = x * q->ct + y * q->st, q1 = -x * q->st + y * q->ct;
            if (fabs(q0 - q->u0) <= 0.5 + tol && fabs(q1 - q->u1) <= 0.5 + tol) tot += q->w; }
    }
    return tot;
}

/* ---- certification: max of cov over all arrangement vertices.  sq: n x 6 (cx,cy,ct,st,col,w).
   Returns number of vertices evaluated; *maxcov = M; bad vertices with cov > thresh written to out
   (x,y,cov) up to cap (nout = count written; *overflow=1 if some were dropped). */
typedef struct { double *buf; long n; long total; double M; } TB;
static void tb_flush(TB *tb, double *out, long cap, long *cnt, int *ovf) {
    #pragma omp critical
    { long room = cap - *cnt; long take = tb->n < room ? tb->n : room;
      if (take > 0) memcpy(out + 3 * (*cnt), tb->buf, sizeof(double) * 3 * (size_t)take);
      *cnt += take; if (take < tb->n) *ovf = 1; }
    tb->n = 0;
}
static inline void emit(const Sq *sq, const Grid *g, double s, double tol, int fund_only, double thresh,
                        double x, double y, TB *tb, double *out, long cap, long *cnt, int *ovf) {
    if (x < -1e-7 || y < -1e-7 || x > s + 1e-7 || y > s + 1e-7) return;
    if (fund_only && !(x <= y + 1e-7 && y <= 0.5 * s + 1e-7)) return;
    double cv = cov_at(sq, g, x, y, tol); tb->total++; if (cv > tb->M) tb->M = cv;
    if (cv > thresh) { tb->buf[3 * tb->n] = x; tb->buf[3 * tb->n + 1] = y; tb->buf[3 * tb->n + 2] = cv; tb->n++;
        if (tb->n == 4096) tb_flush(tb, out, cap, cnt, ovf); }
}
long vertex_cov(int n, const double *a, double s, double tol, int fund_only, double thresh,
                double *out, long cap, long *nout, double *maxcov, int *overflow, int nthreads) {
    Sq *sq = (Sq*)malloc(sizeof(Sq) * (size_t)(n > 0 ? n : 1));
    for (int i = 0; i < n; i++) fill_sq(&sq[i], a, i, 1);
    Grid g; build_grid(sq, n, s, 0.25, &g);
    long total = 0; double M = 0.0; long cnt = 0; int ovf = 0;
    omp_set_num_threads(nthreads);
    #pragma omp parallel reduction(+:total) reduction(max:M)
    {
        TB tb; tb.buf = (double*)malloc(sizeof(double) * 3 * 4096); tb.n = 0; tb.total = 0; tb.M = 0.0;
        #pragma omp for schedule(dynamic, 16)
        for (int A = 0; A < n; A++) {
            const Sq *qa = &sq[A];
            for (int e1 = -1; e1 <= 1; e1 += 2) for (int e2 = -1; e2 <= 1; e2 += 2) {
                double v0 = qa->u0 + 0.5 * e1, v1 = qa->u1 + 0.5 * e2;
                emit(sq, &g, s, tol, fund_only, thresh, qa->ct * v0 - qa->st * v1, qa->st * v0 + qa->ct * v1, &tb, out, cap, &cnt, &ovf); }
            double la[4][3] = {{qa->ct, qa->st, qa->u0 - 0.5}, {qa->ct, qa->st, qa->u0 + 0.5}, {-qa->st, qa->ct, qa->u1 - 0.5}, {-qa->st, qa->ct, qa->u1 + 0.5}};
            double r = 1.4143; int nbk = g.nb; double c = g.cell;
            int bx0 = (int)floor((qa->cx - r) / c), bx1 = (int)floor((qa->cx + r) / c), by0 = (int)floor((qa->cy - r) / c), by1 = (int)floor((qa->cy + r) / c);
            if (bx0 < 0) bx0 = 0; if (by0 < 0) by0 = 0; if (bx1 >= nbk) bx1 = nbk - 1; if (by1 >= nbk) by1 = nbk - 1;
            for (int bx = bx0; bx <= bx1; bx++) for (int by = by0; by <= by1; by++) {
                int b = bx * nbk + by;
                for (int k = g.start[b]; k < g.start[b + 1]; k++) { int B = g.idx[k]; if (B <= A) continue;
                    const Sq *qb = &sq[B];
                    if (fabs(qb->cx - qa->cx) > 1.4143 || fabs(qb->cy - qa->cy) > 1.4143) continue;
                    double lb[4][3] = {{qb->ct, qb->st, qb->u0 - 0.5}, {qb->ct, qb->st, qb->u0 + 0.5}, {-qb->st, qb->ct, qb->u1 - 0.5}, {-qb->st, qb->ct, qb->u1 + 0.5}};
                    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
                        double a1 = la[i][0], b1 = la[i][1], c1 = la[i][2], a2 = lb[j][0], b2 = lb[j][1], c2 = lb[j][2];
                        double det = a1 * b2 - a2 * b1; if (fabs(det) < 1e-6) continue;
                        double x = (c1 * b2 - c2 * b1) / det, y = (a1 * c2 - a2 * c1) / det;
                        double oa = (i < 2) ? (-x * qa->st + y * qa->ct - qa->u1) : (x * qa->ct + y * qa->st - qa->u0);
                        if (fabs(oa) > 0.5 + 1e-7) continue;
                        double ob = (j < 2) ? (-x * qb->st + y * qb->ct - qb->u1) : (x * qb->ct + y * qb->st - qb->u0);
                        if (fabs(ob) > 0.5 + 1e-7) continue;
                        emit(sq, &g, s, tol, fund_only, thresh, x, y, &tb, out, cap, &cnt, &ovf);
                    }
                }
            }
        }
        if (tb.n > 0) tb_flush(&tb, out, cap, &cnt, &ovf);
        total += tb.total; if (tb.M > M) M = tb.M;
        free(tb.buf);
    }
    { double cc[4][2] = {{0, 0}, {s, 0}, {0, s}, {s, s}}; int nn = fund_only ? 1 : 4;
      for (int i = 0; i < nn; i++) { double cv = cov_at(sq, &g, cc[i][0], cc[i][1], tol); total++; if (cv > M) M = cv;
        if (cv > thresh && cnt < cap) { out[3 * cnt] = cc[i][0]; out[3 * cnt + 1] = cc[i][1]; out[3 * cnt + 2] = cv; cnt++; } } }
    *nout = cnt; *maxcov = M; *overflow = ovf;
    free_grid(&g); free(sq); return total;
}

/* ---- incidence: rows = points (npts x 2), cols = squares' col ids (sq: n x 5, cx,cy,ct,st,col).
   pass 0: count per point -> cnt[npts]; pass 1: fill (rowidx, colidx) using offsets. */
void incidence(int npts, const double *pts, int n, const double *a, double s, double tol,
               int pass, long *cnt, const long *off, int *ri, int *ci, int nthreads) {
    Sq *sq = (Sq*)malloc(sizeof(Sq) * (size_t)(n > 0 ? n : 1));
    for (int i = 0; i < n; i++) fill_sq(&sq[i], a, i, 0);
    Grid g; build_grid(sq, n, s, 0.25, &g);
    omp_set_num_threads(nthreads);
    #pragma omp parallel for schedule(dynamic, 64)
    for (int p = 0; p < npts; p++) {
        double x = pts[2 * p], y = pts[2 * p + 1]; double r = 0.7072 + tol; int nb = g.nb; double c = g.cell;
        int bx0 = (int)floor((x - r) / c), bx1 = (int)floor((x + r) / c), by0 = (int)floor((y - r) / c), by1 = (int)floor((y + r) / c);
        if (bx0 < 0) bx0 = 0; if (by0 < 0) by0 = 0; if (bx1 >= nb) bx1 = nb - 1; if (by1 >= nb) by1 = nb - 1;
        long k = pass ? off[p] : 0; long m = 0;
        for (int bx = bx0; bx <= bx1; bx++) for (int by = by0; by <= by1; by++) {
            int b = bx * nb + by;
            for (int t = g.start[b]; t < g.start[b + 1]; t++) { const Sq *q = &sq[g.idx[t]];
                double q0 = x * q->ct + y * q->st, q1 = -x * q->st + y * q->ct;
                if (fabs(q0 - q->u0) <= 0.5 + tol && fabs(q1 - q->u1) <= 0.5 + tol) {
                    if (pass) { ri[k] = p; ci[k] = q->col; k++; } m++; } }
        }
        if (!pass) cnt[p] = m;
    }
    free_grid(&g); free(sq);
}

/* ---- exact capture of a weighted atom set by closed unit squares at poses (m x 3: cx,cy,th) */
void capture(int na, const double *ax, const double *ay, const double *aw, int m, const double *poses,
             double tol, double *out, int nthreads) {
    /* atoms sorted by x are assumed (ax ascending) for pruning */
    omp_set_num_threads(nthreads);
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < m; i++) {
        double cx = poses[3 * i], cy = poses[3 * i + 1], th = poses[3 * i + 2];
        double ct = cos(th), st = sin(th); double u0 = cx * ct + cy * st, u1 = -cx * st + cy * ct;
        double lo = cx - 0.7072 - tol, hi = cx + 0.7072 + tol;
        int a = 0, b = na;                          /* binary search on ax */
        { int l = 0, r = na; while (l < r) { int mid = (l + r) / 2; if (ax[mid] < lo) l = mid + 1; else r = mid; } a = l; }
        { int l = a, r = na; while (l < r) { int mid = (l + r) / 2; if (ax[mid] <= hi) l = mid + 1; else r = mid; } b = l; }
        double tot = 0.0;
        for (int k = a; k < b; k++) { double x = ax[k], y = ay[k]; if (fabs(y - cy) > 0.7072 + tol) continue;
            double q0 = x * ct + y * st, q1 = -x * st + y * ct;
            if (fabs(q0 - u0) <= 0.5 + tol && fabs(q1 - u1) <= 0.5 + tol) tot += aw[k]; }
        out[i] = tot;
    }
}

/* ---- exact closed-square sweep pricer.  For each angle: rotate the atoms, sweep the arrangement
   of breakpoints q0 +- (h+tol) (strips) and, inside each strip, q1 +- (h+tol) (cells); the captured
   weight of the CLOSED unit square is constant on each open cell of centres; report the K best cells
   (admissible centre, separation >= sep) per angle.  out: nang*K x 4 (cap, cx, cy, th), -1 if unused. */
static int cmp_d(const void *a, const void *b) { double x = *(const double*)a, y = *(const double*)b; return (x < y) ? -1 : (x > y); }
typedef struct { double q0, q1, w; } At;
static int cmp_at0(const void *a, const void *b) { double x = ((const At*)a)->q0, y = ((const At*)b)->q0; return (x < y) ? -1 : (x > y); }
static int cmp_at1(const void *a, const void *b) { double x = ((const At*)a)->q1, y = ((const At*)b)->q1; return (x < y) ? -1 : (x > y); }
void sweep_price(int na, const double *ax, const double *ay, const double *aw, double s, int nang, const double *angs,
                 double tol, int K, double sep, double *out, int nthreads) {
    omp_set_num_threads(nthreads);
    #pragma omp parallel
    {
        At *at = (At*)malloc(sizeof(At) * (na + 1)); At *act = (At*)malloc(sizeof(At) * (na + 1));
        double *bx = (double*)malloc(sizeof(double) * (2 * na + 4)); double *by = (double*)malloc(sizeof(double) * (2 * na + 4));
        double *pre = (double*)malloc(sizeof(double) * (na + 2));
        double *best = (double*)malloc(sizeof(double) * 4 * (K + 1));
        #pragma omp for schedule(dynamic, 1)
        for (int ia = 0; ia < nang; ia++) {
            double th = angs[ia], ct = cos(th), st = sin(th), h = 0.5 + tol;
            double wid = fabs(ct) + fabs(st), lo = wid / 2 + 1e-9, hi = s - wid / 2 - 1e-9;
            int nbest = 0;
            for (int k = 0; k < K; k++) out[(ia * K + k) * 4] = -1.0;
            if (hi <= lo) continue;
            for (int k = 0; k < na; k++) { at[k].q0 = ax[k] * ct + ay[k] * st; at[k].q1 = -ax[k] * st + ay[k] * ct; at[k].w = aw[k]; }
            qsort(at, na, sizeof(At), cmp_at0);
            /* u0 range of admissible centres: rotated box corners */
            double u0min = 1e9, u0max = -1e9; double cc[4][2] = {{lo, lo}, {hi, lo}, {hi, hi}, {lo, hi}};
            for (int i = 0; i < 4; i++) { double u = cc[i][0] * ct + cc[i][1] * st; if (u < u0min) u0min = u; if (u > u0max) u0max = u; }
            int nb = 0; for (int k = 0; k < na; k++) { bx[nb++] = at[k].q0 - h; bx[nb++] = at[k].q0 + h; } bx[nb++] = u0min; bx[nb++] = u0max;
            qsort(bx, nb, sizeof(double), cmp_d);
            for (int w = 0; w + 1 < nb; w++) {
                double a = bx[w], b = bx[w + 1]; if (b - a < 1e-12) continue; if (b <= u0min || a >= u0max) continue;
                double m0 = 0.5 * (a + b);
                /* active atoms: q0 in [b-h, a+h] */
                int i0 = 0, i1 = na; { int l = 0, r = na; while (l < r) { int mid = (l + r) / 2; if (at[mid].q0 < b - h) l = mid + 1; else r = mid; } i0 = l; }
                { int l = i0, r = na; while (l < r) { int mid = (l + r) / 2; if (at[mid].q0 <= a + h) l = mid + 1; else r = mid; } i1 = l; }
                int m = i1 - i0;
                /* u1 range: rotated box restricted to this strip (widened: conservative superset) */
                double u1min = 1e9, u1max = -1e9;
                for (int i = 0; i < 4; i++) { double p0 = cc[i][0] * ct + cc[i][1] * st, p1 = -cc[i][0] * st + cc[i][1] * ct;
                    double r0 = cc[(i + 1) % 4][0] * ct + cc[(i + 1) % 4][1] * st, r1 = -cc[(i + 1) % 4][0] * st + cc[(i + 1) % 4][1] * ct;
                    if (p0 >= a - 1e-9 && p0 <= b + 1e-9) { if (p1 < u1min) u1min = p1; if (p1 > u1max) u1max = p1; }
                    for (int e = 0; e < 2; e++) { double xc = e ? b : a;
                        if ((p0 <= xc && xc <= r0) || (r0 <= xc && xc <= p0)) { double yv = (fabs(r0 - p0) > 1e-12) ? p1 + (r1 - p1) * (xc - p0) / (r0 - p0) : p1;
                            if (yv < u1min) u1min = yv; if (yv > u1max) u1max = yv; if (fabs(r0 - p0) <= 1e-12) { if (r1 < u1min) u1min = r1; if (r1 > u1max) u1max = r1; } } } }
                if (!(u1min <= u1max)) continue;
                u1min -= 1e-7; u1max += 1e-7;
                if (m == 0) {   /* whole strip empty: value 0 at any admissible centre */
                    double u1 = 0.5 * (u1min + u1max); double cx = m0 * ct - u1 * st, cy = m0 * st + u1 * ct;
                    if (cx < lo) cx = lo; if (cx > hi) cx = hi; if (cy < lo) cy = lo; if (cy > hi) cy = hi;
                    double v = 0.0; int ok = 1;
                    for (int k = 0; k < nbest; k++) if (fabs(best[4 * k + 1] - cx) < sep && fabs(best[4 * k + 2] - cy) < sep) { if (v < best[4 * k]) { best[4 * k] = v; best[4 * k + 1] = cx; best[4 * k + 2] = cy; } ok = 0; break; }
                    if (ok) { if (nbest < K) { best[4 * nbest] = v; best[4 * nbest + 1] = cx; best[4 * nbest + 2] = cy; nbest++; }
                              else { int worst = 0; for (int k = 1; k < K; k++) if (best[4 * k] > best[4 * worst]) worst = k; if (v < best[4 * worst]) { best[4 * worst] = v; best[4 * worst + 1] = cx; best[4 * worst + 2] = cy; } } }
                    continue; }
                for (int k = 0; k < m; k++) act[k] = at[i0 + k];
                qsort(act, m, sizeof(At), cmp_at1);
                pre[0] = 0; for (int k = 0; k < m; k++) pre[k + 1] = pre[k] + act[k].w;
                int nby = 0; for (int k = 0; k < m; k++) { by[nby++] = act[k].q1 - h; by[nby++] = act[k].q1 + h; } by[nby++] = u1min; by[nby++] = u1max;
                qsort(by, nby, sizeof(double), cmp_d);
                for (int v_ = 0; v_ + 1 < nby; v_++) {
                    double c = by[v_], d = by[v_ + 1]; if (d - c < 1e-12) continue; if (d <= u1min || c >= u1max) continue;
                    int j0, j1; { int l = 0, r = m; while (l < r) { int mid = (l + r) / 2; if (act[mid].q1 < d - h) l = mid + 1; else r = mid; } j0 = l; }
                    { int l = j0, r = m; while (l < r) { int mid = (l + r) / 2; if (act[mid].q1 <= c + h) l = mid + 1; else r = mid; } j1 = l; }
                    double v = (j1 > j0) ? pre[j1] - pre[j0] : 0.0;
                    if (nbest == K) { int worst = 0; for (int k = 1; k < K; k++) if (best[4 * k] > best[4 * worst]) worst = k; if (v >= best[4 * worst]) continue; }
                    double u1 = 0.5 * (c + d); double cx = m0 * ct - u1 * st, cy = m0 * st + u1 * ct;
                    if (cx < lo - 1e-6 || cx > hi + 1e-6 || cy < lo - 1e-6 || cy > hi + 1e-6) continue;
                    if (cx < lo) cx = lo; if (cx > hi) cx = hi; if (cy < lo) cy = lo; if (cy > hi) cy = hi;
                    int ok = 1;
                    for (int k = 0; k < nbest; k++) if (fabs(best[4 * k + 1] - cx) < sep && fabs(best[4 * k + 2] - cy) < sep) { if (v < best[4 * k]) { best[4 * k] = v; best[4 * k + 1] = cx; best[4 * k + 2] = cy; } ok = 0; break; }
                    if (!ok) continue;
                    if (nbest < K) { best[4 * nbest] = v; best[4 * nbest + 1] = cx; best[4 * nbest + 2] = cy; nbest++; }
                    else { int worst = 0; for (int k = 1; k < K; k++) if (best[4 * k] > best[4 * worst]) worst = k; best[4 * worst] = v; best[4 * worst + 1] = cx; best[4 * worst + 2] = cy; }
                }
            }
            for (int k = 0; k < nbest; k++) { out[(ia * K + k) * 4] = best[4 * k]; out[(ia * K + k) * 4 + 1] = best[4 * k + 1]; out[(ia * K + k) * 4 + 2] = best[4 * k + 2]; out[(ia * K + k) * 4 + 3] = th; }
        }
        free(at); free(act); free(bx); free(by); free(pre); free(best);
    }
}

/* ---- raster pricer: for each angle, rotate atoms, bin on a raster of pitch p (rotated frame),
   2-D prefix sums, evaluate the capture of the unit square at every raster centre (cells whose
   centres lie in the square), block minima (B x B blocks), return up to K block-minima per angle
   with capture < thr and admissible centre.  out: nang*K x 4 (cap, cx, cy, th); -1 if unused. */
void raster_price(int na, const double *ax, const double *ay, const double *aw, double s,
                  int nang, const double *angs, double pitch, int B, int K, double thr, double *out, int nthreads) {
    omp_set_num_threads(nthreads);
    #pragma omp parallel
    {
        double *S = NULL; long Scap = 0;
        #pragma omp for schedule(dynamic, 1)
        for (int ia = 0; ia < nang; ia++) {
            double th = angs[ia], ct = cos(th), st = sin(th);
            for (int k = 0; k < K; k++) out[(ia * K + k) * 4] = -1.0;
            /* rotated container corners */
            double q0min = 1e9, q0max = -1e9, q1min = 1e9, q1max = -1e9;
            double cc[4][2] = {{0, 0}, {s, 0}, {0, s}, {s, s}};
            for (int i = 0; i < 4; i++) { double q0 = cc[i][0] * ct + cc[i][1] * st, q1 = -cc[i][0] * st + cc[i][1] * ct;
                if (q0 < q0min) q0min = q0; if (q0 > q0max) q0max = q0; if (q1 < q1min) q1min = q1; if (q1 > q1max) q1max = q1; }
            int n0 = (int)((q0max - q0min) / pitch) + 2, n1 = (int)((q1max - q1min) / pitch) + 2;
            long need = (long)(n0 + 1) * (n1 + 1);
            if (need > Scap) { free(S); S = (double*)malloc(sizeof(double) * need); Scap = need; }
            memset(S, 0, sizeof(double) * need);
            #define SS(i, j) S[(long)(i) * (n1 + 1) + (j)]
            for (int k = 0; k < na; k++) { double q0 = ax[k] * ct + ay[k] * st, q1 = -ax[k] * st + ay[k] * ct;
                int i = (int)((q0 - q0min) / pitch), j = (int)((q1 - q1min) / pitch);
                if (i < 0) i = 0; if (j < 0) j = 0; if (i >= n0) i = n0 - 1; if (j >= n1) j = n1 - 1; SS(i + 1, j + 1) += aw[k]; }
            for (int i = 1; i <= n0; i++) for (int j = 1; j <= n1; j++) SS(i, j) += SS(i - 1, j) + SS(i, j - 1) - SS(i - 1, j - 1);
            /* square of side 1 covers hw = round(0.5/pitch) cells on each side of the centre cell */
            int hw = (int)floor(0.5 / pitch + 1e-9);
            double wid = fabs(ct) + fabs(st), lo = wid / 2, hi = s - wid / 2;
            int nb0 = (n0 + B - 1) / B, nb1 = (n1 + B - 1) / B;
            long nblk = (long)nb0 * nb1;
            double *bmin = (double*)malloc(sizeof(double) * nblk); int *barg = (int*)malloc(sizeof(int) * 2 * nblk);
            for (long b = 0; b < nblk; b++) bmin[b] = 1e9;
            for (int i = 0; i < n0; i++) { int i0 = i - hw, i1 = i + hw; if (i0 < 0 || i1 >= n0) continue;
                for (int j = 0; j < n1; j++) { int j0 = j - hw, j1 = j + hw; if (j0 < 0 || j1 >= n1) continue;
                    double u0 = q0min + (i + 0.5) * pitch, u1 = q1min + (j + 0.5) * pitch;
                    double cx = u0 * ct - u1 * st, cy = u0 * st + u1 * ct;
                    if (cx < lo || cx > hi || cy < lo || cy > hi) continue;
                    double v = SS(i1 + 1, j1 + 1) - SS(i0, j1 + 1) - SS(i1 + 1, j0) + SS(i0, j0);
                    long b = (long)(i / B) * nb1 + (j / B);
                    if (v < bmin[b]) { bmin[b] = v; barg[2 * b] = i; barg[2 * b + 1] = j; } } }
            /* K best blocks with min < thr */
            for (int k = 0; k < K; k++) { long best = -1; double bv = thr;
                for (long b = 0; b < nblk; b++) if (bmin[b] < bv) { bv = bmin[b]; best = b; }
                if (best < 0) break;
                int i = barg[2 * best], j = barg[2 * best + 1];
                double u0 = q0min + (i + 0.5) * pitch, u1 = q1min + (j + 0.5) * pitch;
                out[(ia * K + k) * 4] = bv; out[(ia * K + k) * 4 + 1] = u0 * ct - u1 * st; out[(ia * K + k) * 4 + 2] = u0 * st + u1 * ct; out[(ia * K + k) * 4 + 3] = th;
                bmin[best] = 1e9; }
            free(bmin); free(barg);
            #undef SS
        }
        free(S);
    }
}
'''


def build_lib():
    os.makedirs(RUNS, exist_ok=True)
    h = hashlib.sha256(CSRC.encode()).hexdigest()[:12]
    so = os.path.join(RUNS, f'_pdual_{h}.so'); src = os.path.join(RUNS, f'_pdual_{h}.c')
    if not os.path.exists(so):
        open(src, 'w').write(CSRC)
        subprocess.run(['gcc', '-O2', '-march=native', '-fopenmp', '-shared', '-fPIC', '-o', so, src, '-lm'], check=True)
    lib = ctypes.CDLL(so)
    dp = ctypes.POINTER(ctypes.c_double); ip = ctypes.POINTER(ctypes.c_int); lp = ctypes.POINTER(ctypes.c_long)
    lib.vertex_cov.restype = ctypes.c_long
    lib.vertex_cov.argtypes = [ctypes.c_int, dp, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double, dp, ctypes.c_long, lp, dp, ip, ctypes.c_int]
    lib.incidence.restype = None
    lib.incidence.argtypes = [ctypes.c_int, dp, ctypes.c_int, dp, ctypes.c_double, ctypes.c_double, ctypes.c_int, lp, lp, ip, ip, ctypes.c_int]
    lib.capture.restype = None
    lib.capture.argtypes = [ctypes.c_int, dp, dp, dp, ctypes.c_int, dp, ctypes.c_double, dp, ctypes.c_int]
    lib.sweep_price.restype = None
    lib.sweep_price.argtypes = [ctypes.c_int, dp, dp, dp, ctypes.c_double, ctypes.c_int, dp, ctypes.c_double, ctypes.c_int, ctypes.c_double, dp, ctypes.c_int]
    lib.raster_price.restype = None
    lib.raster_price.argtypes = [ctypes.c_int, dp, dp, dp, ctypes.c_double, ctypes.c_int, dp, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_double, dp, ctypes.c_int]
    return lib


def _dp(a): return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


# ============================================================================== poses
def snap_angle(th):
    return round(th / ANG_UNIT) * ANG_UNIT


def wid_of(th):
    return abs(math.cos(th)) + abs(math.sin(th))


def admissible(cx, cy, th, t, margin=ADM):
    w2 = wid_of(th) / 2
    return (w2 + margin <= cx <= t - w2 - margin) and (w2 + margin <= cy <= t - w2 - margin)


def clamp_pose(cx, cy, th, t):
    w2 = wid_of(th) / 2; lo = w2 + ADM; hi = t - w2 - ADM
    if hi < lo: return None
    return (min(max(cx, lo), hi), min(max(cy, lo), hi), th)


def pose_images(cx, cy, th, t):
    """the 8 dihedral images of a pose (angle sign as in nu_f.max_coverage)"""
    return [(cx, cy, th), (t - cx, cy, -th), (cx, t - cy, -th), (t - cx, t - cy, th),
            (cy, cx, -th), (t - cy, cx, th), (cy, t - cx, th), (t - cy, t - cx, -th)]


def canon_pose(cx, cy, th, t):
    """canonical representative of the orbit: angle |th| in [0, pi/4] (angles are mod pi/2),
    lexicographically smallest centre among the images at that angle"""
    th = snap_angle(th)
    a = abs(th)
    if a > math.pi / 4 + 1e-12:           # fold into (-pi/4, pi/4]
        th = th - math.copysign(math.pi / 2, th); th = snap_angle(th); a = abs(th)
    best = None
    for (x, y, u) in pose_images(cx, cy, th, t):
        uu = u if u >= 0 else u + math.pi / 2
        if abs(uu - a) < 1e-9:
            key = (round(x, 9), round(y, 9))
            if best is None or key < best: best = key
    return (best[0], best[1], a)


def pose_key(p):
    return (round(p[0] * 1e7), round(p[1] * 1e7), round(p[2] / ANG_UNIT))


def images_array(poses, t, weights=None):
    """n x 5 (cx,cy,ct,st,col) or n x 6 (+w) array of all images"""
    rows = []
    for k, (cx, cy, th) in enumerate(poses):
        for (x, y, u) in pose_images(cx, cy, th, t):
            r = [x, y, math.cos(u), math.sin(u), float(k)]
            if weights is not None: r.append(weights[k] / 8.0)
            rows.append(r)
    return np.ascontiguousarray(np.array(rows, dtype=float)) if rows else np.zeros((0, 6 if weights is not None else 5))


def fund(pts, t):
    """map points to the fundamental domain 0 <= x <= y <= t/2"""
    x = np.minimum(pts[:, 0], t - pts[:, 0]); y = np.minimum(pts[:, 1], t - pts[:, 1])
    return np.c_[np.minimum(x, y), np.maximum(x, y)]


# ============================================================================== model
class Model:
    def __init__(self, t, lib, threads, log):
        self.t = t; self.lib = lib; self.threads = threads; self.log = log
        self.poses = []; self.pkey = {}; self.col_zero = []       # columns
        self.pts = np.zeros((0, 2)); self.pkeys = set(); self.row_zero = np.zeros(0, dtype=int)   # rows
        self.A = None                                              # csr rows x cols (values 1/8 per image hit)

    # ---- columns
    def add_poses(self, cand):
        new = []
        for (cx, cy, th) in cand:
            p = clamp_pose(cx, cy, th, self.t)
            if p is None: continue
            p = canon_pose(*p, self.t)
            if not admissible(*p, self.t): continue
            k = pose_key(p)
            if k in self.pkey: continue
            self.pkey[k] = len(self.poses); self.poses.append(p); self.col_zero.append(0); new.append(p)
        if new and len(self.pts):
            imgs = images_array(new, self.t)
            Anew = self._incidence(self.pts, imgs, ncols=len(new))
            self.A = Anew if self.A is None else sp.hstack([self.A, Anew], format='csr')
        return len(new)

    def add_points(self, pts):
        pts = np.asarray(pts, dtype=float).reshape(-1, 2)
        if len(pts) == 0: return 0
        pts = fund(pts, self.t)
        keep = []
        for (x, y) in pts:
            k = (round(x * 1e9), round(y * 1e9))
            if k in self.pkeys: continue
            self.pkeys.add(k); keep.append((x, y))
        if not keep: return 0
        P = np.array(keep)
        if self.poses:
            imgs = images_array(self.poses, self.t)
            Anew = self._incidence(P, imgs, ncols=len(self.poses))
            self.A = Anew if self.A is None else sp.vstack([self.A, Anew], format='csr')
        self.pts = np.vstack([self.pts, P]); self.row_zero = np.concatenate([self.row_zero, np.zeros(len(P), dtype=int)])
        return len(P)

    def _incidence(self, P, imgs, ncols):
        P = np.ascontiguousarray(P, dtype=float); imgs = np.ascontiguousarray(imgs, dtype=float)
        n = len(P); cnt = np.zeros(n, dtype=np.int64); off = np.zeros(n, dtype=np.int64)
        self.lib.incidence(n, _dp(P), len(imgs), _dp(imgs), self.t, TOL, 0, cnt.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
                           off.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), None, None, self.threads)
        off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64); nnz = int(off[-1])
        ri = np.zeros(max(nnz, 1), dtype=np.int32); ci = np.zeros(max(nnz, 1), dtype=np.int32)
        self.lib.incidence(n, _dp(P), len(imgs), _dp(imgs), self.t, TOL, 1, cnt.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
                           off.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), ri.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                           ci.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), self.threads)
        return sp.coo_matrix((np.full(nnz, 0.125), (ri[:nnz], ci[:nnz])), shape=(n, ncols)).tocsr()

    def drop_columns(self, keep):
        idx = np.nonzero(keep)[0]
        self.poses = [self.poses[i] for i in idx]; self.col_zero = [self.col_zero[i] for i in idx]
        self.pkey = {pose_key(p): i for i, p in enumerate(self.poses)}
        self.A = self.A[:, idx]

    def drop_rows(self, keep):
        idx = np.nonzero(keep)[0]
        self.pts = self.pts[idx]; self.row_zero = self.row_zero[idx]; self.A = self.A[idx]
        self.pkeys = set((round(x * 1e9), round(y * 1e9)) for x, y in self.pts)

    # ---- LP
    def solve(self):
        n = len(self.poses); m = len(self.pts)
        res = linprog(c=-np.ones(n), A_ub=self.A, b_ub=np.ones(m), bounds=(0, None), method='highs')
        if not res.success: return None
        mu = np.maximum(res.x, 0.0); y = np.maximum(-res.ineqlin.marginals, 0.0)
        return mu, y, -res.fun

    # ---- certification (rigorous)
    def certify(self, mu, thresh=1.0 + 1e-9, cap=400_000, fund_only=True):
        sup = np.nonzero(mu > 1e-12)[0]
        for i in sup: assert admissible(*self.poses[i], self.t, margin=0.0), self.poses[i]   # poses are clamped with margin ADM at creation
        imgs = images_array([self.poses[i] for i in sup], self.t, weights=mu[sup])
        out = np.zeros((cap, 3)); nout = ctypes.c_long(0); M = ctypes.c_double(0.0); ovf = ctypes.c_int(0)
        nv = self.lib.vertex_cov(len(imgs), _dp(imgs), self.t, TOL, 1 if fund_only else 0, thresh, _dp(out), cap,
                                 ctypes.byref(nout), ctypes.byref(M), ctypes.byref(ovf), self.threads)
        bad = out[:nout.value]
        return M.value, bad, nv, bool(ovf.value)


# ============================================================================== pricing
def atoms_from_dual(model, y, ymin=0.0, maxrows=None):
    """D4-symmetrised weighted atoms (x, y, w) of the dual point weights y (mass y_p/8 per image);
    with maxrows, only the heaviest rows are kept (whole orbits, so the set stays D4-symmetric)"""
    sel = np.nonzero(y > ymin)[0]
    if maxrows is not None and len(sel) > maxrows: sel = sel[np.argsort(y[sel])[::-1][:maxrows]]
    if len(sel) == 0: return np.zeros(0), np.zeros(0), np.zeros(0)
    P = model.pts[sel]; w = y[sel]; t = model.t
    X = []; Y = []; W = []
    for (a, b) in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        px = t - P[:, 0] if a else P[:, 0]; py = t - P[:, 1] if b else P[:, 1]
        X += [px, py]; Y += [py, px]; W += [w / 8, w / 8]
    X = np.concatenate(X); Y = np.concatenate(Y); W = np.concatenate(W)
    o = np.argsort(X, kind='stable')
    return np.ascontiguousarray(X[o]), np.ascontiguousarray(Y[o]), np.ascontiguousarray(W[o])


def capture(lib, ax, ay, aw, poses, threads, tol=TOL):
    poses = np.ascontiguousarray(np.asarray(poses, dtype=float).reshape(-1, 3)); out = np.zeros(len(poses))
    if len(poses): lib.capture(len(ax), _dp(ax), _dp(ay), _dp(aw), len(poses), _dp(poses), tol, _dp(out), threads)
    return out


def refine(lib, ax, ay, aw, poses, t, threads, steps=(4e-3, 1e-3, 2.5e-4, 6e-5), dth=(2e-3, 5e-4, 1.2e-4, 3e-5)):
    """coordinate-descent on the exact capture (heuristic), keeping poses admissible"""
    cur = [clamp_pose(*p, t) for p in poses]; cur = [p for p in cur if p is not None]
    if not cur: return []
    cur = np.array(cur); val = capture(lib, ax, ay, aw, cur, threads)
    for h, dh in zip(steps, dth):
        for _ in range(2):
            moves = np.array([[h, 0, 0], [-h, 0, 0], [0, h, 0], [0, -h, 0], [0, 0, dh], [0, 0, -dh], [h, h, 0], [-h, -h, 0], [h, -h, 0], [-h, h, 0]])
            cand = cur[:, None, :] + moves[None, :, :]; cand = cand.reshape(-1, 3)
            cl = [clamp_pose(*c, t) for c in cand]
            cand = np.array([c if c is not None else (np.nan, np.nan, np.nan) for c in cl])
            v = capture(lib, ax, ay, aw, np.nan_to_num(cand), threads); v[np.isnan(cand[:, 0])] = 9e9
            v = v.reshape(len(cur), len(moves)); j = np.argmin(v, axis=1); better = v[np.arange(len(cur)), j] < val - 1e-12
            cur[better] = cand.reshape(len(cur), len(moves), 3)[np.arange(len(cur)), j][better]; val = np.minimum(val, v[np.arange(len(cur)), j])
    return [(float(a), float(b), float(c), float(v)) for (a, b, c), v in zip(cur, val)]


def write_cert(ax, ay, aw, t, path, D=100000, WD=10 ** 8):
    """the symmetrised atom set as an exact certificate (weights rounded UP)"""
    tf = Fraction(t).limit_denominator(1000); assert (tf * D).denominator == 1
    SD = int(tf * D)
    X = np.clip(np.round(ax * D).astype(np.int64), 0, SD); Y = np.clip(np.round(ay * D).astype(np.int64), 0, SD)
    W = np.ceil(aw * WD).astype(np.int64)
    keep = W > 0
    with open(path, 'w') as f:
        f.write(f"{tf.numerator} {tf.denominator}\n{D}\n{WD}\n{int(keep.sum())}\n")
        for a, b, c in zip(X[keep], Y[keep], W[keep]): f.write(f"{a} {b} {c}\n")
    return W[keep].sum() / WD


def run_verifier(path, N, threads, topk, sep):
    args = [VERIFY, path, '12', str(N), str(threads), str(topk), sep]
    r = subprocess.run(args, capture_output=True, text=True, timeout=7200)
    ln = [l for l in r.stdout.split('\n') if l.startswith('min covered')]
    if not ln: raise RuntimeError(r.stdout[-500:] + r.stderr[-500:])
    frac = ln[0].split('=')[1].strip().split()[0]; a, b = frac.split('/')
    sym = 'D4-symmetric atom set: true' in r.stdout
    return Fraction(int(a), int(b)), sym


def price_verifier(model, ax, ay, aw, N, topk, tag, maxatoms, threads, log):
    """exact pricing: the verifier's worst placements (shrunk sigma-square at bin angles) of the
    dual certificate; atoms are pruned to the heaviest `maxatoms` (only ADDS witnesses)."""
    t = model.t
    path = os.path.join(RUNS, f'dual_{tag}_price.txt'); sep = os.path.join(RUNS, f'dual_{tag}_sep.txt')
    tot = write_cert(ax, ay, aw, t, path)
    t0 = time.time(); mv, sym = run_verifier(path, N, threads, topk, sep)
    wit = []
    if os.path.exists(sep):
        for line in open(sep):
            q = line.split()
            if len(q) == 4: wit.append((float(q[2]), float(q[3]), float(q[1]), float(q[0])))
        os.remove(sep)
    log(f"   verifier N={N} atoms={len(ax)} total={tot:.5f} min={float(mv):.6f} sym={sym} witnesses={len(wit)} ({time.time()-t0:.0f}s)")
    return float(mv), tot, wit


def price_sweep(model, lib, ax, ay, aw, angs, K, sep, threads, tol=TOL):
    """exact closed-square minimum capture per angle (K best well-separated centres per angle)"""
    out = np.zeros((len(angs) * K, 4)); angs = np.ascontiguousarray(angs, dtype=float)
    lib.sweep_price(len(ax), _dp(ax), _dp(ay), _dp(aw), model.t, len(angs), _dp(angs), tol, K, sep, _dp(out), threads)
    out = out[out[:, 0] >= 0]
    return [(float(cx), float(cy), float(th), float(v)) for v, cx, cy, th in out]


def neighbours(poses, t):
    """perturbations of poses (for the LP to fine-tune the support)"""
    out = []
    for (cx, cy, th) in poses:
        for d, dth in ((0.02, math.radians(1.0)), (0.004, math.radians(0.2)), (0.0008, math.radians(0.04))):
            for (a, b, c) in ((d, 0, 0), (-d, 0, 0), (0, d, 0), (0, -d, 0), (0, 0, dth), (0, 0, -dth), (d, d, 0), (-d, -d, 0), (d, -d, 0), (-d, d, 0)):
                out.append((cx + a, cy + b, th + c))
    return out


def price_raster(model, lib, ax, ay, aw, angs, pitch, B, K, thr, threads):
    out = np.zeros((len(angs) * K, 4)); angs = np.ascontiguousarray(angs, dtype=float)
    lib.raster_price(len(ax), _dp(ax), _dp(ay), _dp(aw), model.t, len(angs), _dp(angs), pitch, B, K, thr, _dp(out), threads)
    out = out[out[:, 0] >= 0]
    return [(float(cx), float(cy), float(th), float(v)) for v, cx, cy, th in out]


# ============================================================================== main loop
def seed_poses(t, pitch, dth_deg):
    out = []
    for th in np.arange(0.0, math.pi / 4 + 1e-9, math.radians(dth_deg)):
        th = snap_angle(min(th, math.pi / 4)); w2 = wid_of(th) / 2
        k = int(math.floor((t / 2 - w2 - ADM) / pitch))            # symmetric about t/2 so that orbits dedupe
        g = t / 2 + pitch * np.arange(-k, k + 1)
        g = np.unique(np.concatenate([g, [w2 + ADM, t - w2 - ADM]]))
        for cx in g:
            for cy in g: out.append((float(cx), float(cy), th))
    return out


def seed_points(t, pitch):
    n0 = int(math.ceil(t / pitch)); g = (np.arange(n0) + 0.5) * (t / n0)
    G0, G1 = np.meshgrid(g, g, indexing='ij'); G = np.c_[G0.ravel(), G1.ravel()]
    return G[(G[:, 0] <= G[:, 1] + 1e-12) & (G[:, 1] <= t / 2 + 1e-12)]


def corners(poses, t):
    imgs = images_array(poses, t)
    out = []
    for e1 in (-0.5, 0.5):
        for e2 in (-0.5, 0.5):
            u0 = imgs[:, 0] * imgs[:, 2] + imgs[:, 1] * imgs[:, 3] + e1; u1 = -imgs[:, 0] * imgs[:, 3] + imgs[:, 1] * imgs[:, 2] + e2
            out.append(np.c_[imgs[:, 2] * u0 - imgs[:, 3] * u1, imgs[:, 3] * u0 + imgs[:, 2] * u1])
    return np.concatenate(out)


def read_support(path):
    poses = []
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        q = line.split()
        if q[0] == 'pose': poses.append((float(q[1]), float(q[2]), math.radians(float(q[3]))))
    return poses


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('T'); ap.add_argument('TAG')
    ap.add_argument('--rounds', type=int, default=60); ap.add_argument('--time', type=float, default=5400)
    ap.add_argument('--seed-pitch', type=float, default=0.05); ap.add_argument('--seed-dth', type=float, default=5.0)
    ap.add_argument('--row-pitch', type=float, default=0.02)
    ap.add_argument('--seed-file', default=None, help='extra seed poses: lines "cx cy theta_deg"')
    ap.add_argument('--warm', default=None, help='support file of another run (poses rescaled to T)')
    ap.add_argument('--threads', type=int, default=32)
    ap.add_argument('--N', type=int, default=1000, help='verifier angle parameter for pricing')
    ap.add_argument('--topk', type=int, default=8)
    ap.add_argument('--maxatoms', type=int, default=6000, help='prune the pricing certificate to this many atoms')
    ap.add_argument('--verifier-every', type=int, default=1)
    ap.add_argument('--cg-want', type=int, default=1500, help='max new columns per round')
    ap.add_argument('--row-cap', type=int, default=30000, help='max new rows (bad vertices) per round')
    ap.add_argument('--raster-dth', type=float, default=0.1); ap.add_argument('--raster-pitch', type=float, default=0.004)
    ap.add_argument('--col-age', type=int, default=4, help='drop columns at mu=0 for this many rounds')
    ap.add_argument('--row-age', type=int, default=6, help='drop rows with y=0 and slack>0.02 for this many rounds')
    ap.add_argument('--smooth', type=float, default=0.5, help='dual smoothing factor for pricing (0 = off)')
    ap.add_argument('--sweep-rows', type=int, default=350, help='heaviest dual rows used by the sweep pricer')
    ap.add_argument('--sweep-k', type=int, default=6, help='candidates per angle from the exact sweep pricer')
    ap.add_argument('--polish-after', type=int, default=10**9, help='round after which no more columns are added')
    ap.add_argument('--polish-time', type=float, default=900, help='extra seconds for the rows-only polish phase')
    ap.add_argument('--final-full', action='store_true', help='certify the best measure over the whole container too')
    a = ap.parse_args()
    t = float(Fraction(a.T)); tag = a.TAG
    os.makedirs(RUNS, exist_ok=True)
    lf = open(os.path.join(RUNS, f'dual_{tag}.log'), 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"packing_dual t={t} tag={tag} args={vars(a)}")
    lib = build_lib(); T0 = time.time()
    m = Model(t, lib, a.threads, log)
    # seed columns
    seeds = seed_poses(t, a.seed_pitch, a.seed_dth)
    if a.warm:
        lam = None
        for line in open(a.warm):
            if line.startswith('# t='): lam = t / float(line.split('=')[1].split()[0])
        wp = read_support(a.warm); lam = lam or 1.0
        seeds += [(cx * lam, cy * lam, th) for cx, cy, th in wp]
        log(f"warm start: {len(wp)} poses from {a.warm} (centres scaled by {lam:.6f})")
    if a.seed_file:
        ext = []
        for l in open(a.seed_file):
            q = l.split()
            if not q or l.startswith('#'): continue
            if q[0] == 'pose': q = q[1:]                       # a support file of another run
            elif q[0] == 'point': continue
            if len(q) >= 3: ext.append((float(q[0]), float(q[1]), math.radians(float(q[2]))))
        seeds += ext; log(f"seed file: {len(ext)} poses from {a.seed_file}")
    m.add_poses(seeds)
    # seed rows: grid + corners of the seed squares
    m.add_points(seed_points(t, a.row_pitch))
    m.add_points(corners(m.poses, t))
    log(f"seed: {len(m.poses)} poses, {len(m.pts)} rows, nnz={m.A.nnz} ({time.time()-T0:.0f}s)")
    hist = []; best = dict(L=0.0); angs = np.arange(0.0, math.pi / 4 + 1e-9, math.radians(a.raster_dth)); ysm = {}
    for rd in range(a.rounds):
        t0 = time.time(); out = m.solve()
        if out is None: log("LP failed"); break
        mu, y, mass = out; tlp = time.time() - t0
        # ---- rigorous certification of the LP measure
        t0 = time.time(); M, bad, nv, ovf = m.certify(mu, cap=2_000_000)
        L = mass / M; tc = time.time() - t0
        nsup = int((mu > 1e-12).sum()); ysup = int((y > 1e-12).sum())
        rec = dict(round=rd, t=time.time() - T0, mass=mass, M=M, L=L, support=nsup, dual_support=ysup, cols=len(m.poses), rows=len(m.pts),
                   nnz=int(m.A.nnz), vertices=int(nv), bad=int(len(bad)), lp_s=tlp, cert_s=tc)
        log(f"round {rd}: LP mass={mass:.5f} M={M:.6f} -> L={L:.5f} | support={nsup} dual={ysup} cols={len(m.poses)} rows={len(m.pts)} nnz={m.A.nnz} "
            f"vertices={nv} bad={len(bad)}{' (capped)' if ovf else ''} | lp {tlp:.0f}s cert {tc:.0f}s total {time.time()-T0:.0f}s")
        write_support(m, mu, y, t, tag + '_last', L, M, mass, rd)          # every round (provenance / warm starts)
        if L > best['L']:
            best = dict(L=L, M=M, mass=mass, round=rd, support=nsup)
            write_support(m, mu, y, t, tag, L, M, mass, rd)
        # ---- rows: the violated vertices (worst first, capped)
        t0 = time.time(); n_old_cols = len(mu); n_old_rows = len(y)
        slack = 1.0 - np.asarray(m.A @ mu).ravel()
        if len(bad):
            o = np.argsort(bad[:, 2])[::-1][:a.row_cap]; nrow = m.add_points(bad[o, :2])
        else: nrow = 0
        # ---- pricing (skipped in the polish phase)
        ax, ay, aw = atoms_from_dual(m, y, ymin=1e-12)
        # dual smoothing (Wentges): price also against a running average of the duals
        keys = [(round(x * 1e9), round(yy * 1e9)) for x, yy in m.pts[:n_old_rows]]
        ysm_new = {}
        for k_, yv in zip(keys, y): ysm_new[k_] = a.smooth * ysm.get(k_, 0.0) + (1 - a.smooth) * yv
        ysm = ysm_new; ys_vec = np.array([ysm[k_] for k_ in keys] + [0.0] * (len(m.pts) - n_old_rows))
        cand = []; vmin = None; U = None; rmin = None; ncol = 0
        polish = (rd >= a.polish_after) or (time.time() - T0 > a.time)
        if not polish:
            if a.verifier_every and rd % a.verifier_every == 0 and len(ax):
                pruned = (y > 1e-12).sum() > a.maxatoms // 8
                try:
                    px, py, pw = atoms_from_dual(m, y, ymin=1e-12, maxrows=a.maxatoms // 8) if pruned else (ax, ay, aw)
                    vmin, tot, wit = price_verifier(m, px, py, pw, a.N, a.topk, tag, a.maxatoms, a.threads, log)
                    if not pruned and vmin > 0: U = tot / vmin
                    cand += [(cx, cy, th) for cx, cy, th, v in wit]
                except Exception as e:
                    log(f"   [verifier failed: {e}]")
            t1 = time.time()
            sx, sy, sww = atoms_from_dual(m, y, ymin=1e-12, maxrows=a.sweep_rows)
            sw = price_sweep(m, lib, sx, sy, sww, angs, a.sweep_k, 0.05, a.threads)
            swmin = min([v for cx, cy, th, v in sw], default=None)
            cand += [(cx, cy, th) for cx, cy, th, v in sw]
            if a.smooth > 0:
                sx, sy, sww = atoms_from_dual(m, ys_vec, ymin=1e-12, maxrows=a.sweep_rows)
                if len(sx):
                    sw2 = price_sweep(m, lib, sx, sy, sww, angs, a.sweep_k, 0.05, a.threads); cand += [(cx, cy, th) for cx, cy, th, v in sw2]
            sup_poses = [m.poses[i] for i in np.nonzero(mu > 1e-9)[0]]
            nb = neighbours(sup_poses, t); tsw = time.time() - t1
            raw = [(cx, cy, th, v) for (cx, cy, th), v in zip(cand, capture(lib, ax, ay, aw, cand, a.threads))] if cand else []
            ref = refine(lib, ax, ay, aw, cand, t, a.threads) if cand else []
            nbv = [(cx, cy, th, v) for (cx, cy, th), v in zip(nb, capture(lib, ax, ay, aw, nb, a.threads))] if nb else []
            allc = sorted(raw + ref + nbv, key=lambda r: r[3])
            rmin = allc[0][3] if allc else None
            newc = [(cx, cy, th) for cx, cy, th, v in allc if v < 1.0 - 1e-7][:a.cg_want]
            ncol = m.add_poses(newc)
            log(f"   sweep min {swmin} ({tsw:.0f}s) | candidates {len(cand)} raw+refined+neighbours {len(allc)} improving {sum(1 for r in allc if r[3] < 1 - 1e-7)}")
        # ---- housekeeping: age and drop stale columns / rows (never the ones just added)
        for i in range(n_old_cols): m.col_zero[i] = m.col_zero[i] + 1 if mu[i] <= 1e-12 else 0
        keepc = np.array([cz < a.col_age for cz in m.col_zero], dtype=bool)
        if (~keepc).sum() > 0: m.drop_columns(keepc)
        stale = (y <= 1e-12) & (slack > 0.02)
        m.row_zero[:n_old_rows] = np.where(stale, m.row_zero[:n_old_rows] + 1, 0)
        keepr = m.row_zero < a.row_age
        if (~keepr).sum() > 0: m.drop_rows(keepr)
        rec.update(new_rows=nrow, new_cols=ncol, verifier_min=vmin, U=U, refined_min=rmin, cols_after=len(m.poses), rows_after=len(m.pts))
        hist.append(rec)
        log(f"   +rows {nrow} +cols {ncol} (verifier min {vmin}, refined min {rmin}, U={U}) -> cols={len(m.poses)} rows={len(m.pts)} ({time.time()-t0:.0f}s)")
        json.dump(dict(t=t, tag=tag, args=vars(a), best=best, hist=hist), open(os.path.join(RUNS, f'dual_{tag}.json'), 'w'), indent=1)
        if time.time() - T0 > a.time + a.polish_time: log("time limit"); break
        if nrow == 0 and ncol == 0:
            log("converged: no violated vertex" + ("" if polish else ", no improving column")); break
    log(f"BEST t={t} L={best['L']:.5f} (round {best.get('round')}, mass {best.get('mass')}, M {best.get('M')}) wall {time.time()-T0:.0f}s")
    if a.final_full and best['L'] > 0:
        # re-certify the saved best measure over the WHOLE container (no symmetry reduction)
        poses = []; mus = []
        for line in open(os.path.join(RUNS, f'dual_{tag}_support.txt')):
            q = line.split()
            if q and q[0] == 'pose': poses.append((float(q[1]), float(q[2]), math.radians(float(q[3])))); mus.append(float(q[4]))
        mm = Model(t, lib, a.threads, log); mm.poses = poses
        M2, bad2, nv2, _ = mm.certify(np.array(mus), cap=10, fund_only=False)
        log(f"full-container certification of the best measure: M={M2:.9f} over {nv2} vertices -> L={sum(mus)/M2:.5f}")
        best['M_full'] = M2; best['L_full'] = sum(mus) / M2
        json.dump(dict(t=t, tag=tag, args=vars(a), best=best, hist=hist), open(os.path.join(RUNS, f'dual_{tag}.json'), 'w'), indent=1)


def write_support(m, mu, y, t, tag, L, M, mass, rd):
    with open(os.path.join(RUNS, f'dual_{tag}_support.txt'), 'w') as f:
        f.write(f"# t={t} round={rd} mass={mass:.9f} M={M:.9f} L={L:.9f}\n")
        f.write("# D4-symmetrised measure: pose cx cy theta_deg mu  (mass mu/8 on each of the 8 dihedral images)\n")
        for i in np.nonzero(mu > 1e-6)[0]:
            cx, cy, th = m.poses[i]; f.write(f"pose {cx:.13f} {cy:.13f} {math.degrees(th):.11f} {mu[i]:.12f}\n")
        f.write("# constraint points (fundamental domain 0<=x<=y<=t/2) with dual weight y\n")
        for i in range(len(m.pts)):
            f.write(f"point {m.pts[i,0]:.10f} {m.pts[i,1]:.10f} {y[i]:.10f}\n")


def analyse(path):
    """describe a support file: mass by angle, where the centres sit, concentration"""
    poses = []; mus = []; t = None
    for line in open(path):
        q = line.split()
        if line.startswith('# t='): t = float(line.split('=')[1].split()[0])
        if q and q[0] == 'pose': poses.append((float(q[1]), float(q[2]), float(q[3]))); mus.append(float(q[4]))
    mus = np.array(mus); P = np.array(poses); tot = mus.sum(); o = np.argsort(mus)[::-1]
    print(f"{path}: t={t} support={len(mus)} mass={tot:.5f}")
    cum = np.cumsum(mus[o]) / tot
    for f in (0.25, 0.5, 0.75, 0.9, 0.99): print(f"  {f*100:.0f}% of the mass on the {int(np.searchsorted(cum, f)) + 1} heaviest poses")
    print("  angle histogram (mass per 5 deg bin):")
    for lo in range(0, 46, 5):
        sel = (P[:, 2] >= lo) & (P[:, 2] < lo + 5) if lo < 45 else (P[:, 2] >= 45)
        if sel.any(): print(f"    [{lo:2d},{lo+5:2d}) deg: mass {mus[sel].sum():7.4f} ({100*mus[sel].sum()/tot:5.1f}%)  poses {sel.sum()}")
    near0 = P[:, 2] < 1.0; print(f"  axis-aligned (|theta|<1 deg): mass {mus[near0].sum():.4f} ({100*mus[near0].sum()/tot:.1f}%)")
    # where the centres sit (canonical image, so cx <= cy, cy <= t/2 roughly): distance to the nearest wall
    d = np.minimum(np.minimum(P[:, 0], t - P[:, 0]), np.minimum(P[:, 1], t - P[:, 1]))
    for lo, hi in ((0, 0.52), (0.52, 0.8), (0.8, 1.2), (1.2, 1.6), (1.6, 2.0)):
        sel = (d >= lo) & (d < hi)
        if sel.any(): print(f"  centre distance to nearest wall in [{lo},{hi}): mass {mus[sel].sum():7.4f} ({100*mus[sel].sum()/tot:5.1f}%) poses {sel.sum()}")
    print("  heaviest poses (cx cy theta_deg mu):")
    for i in o[:20]: print(f"    {P[i,0]:.5f} {P[i,1]:.5f} {P[i,2]:9.4f} {mus[i]:.5f}")


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'analyse': analyse(sys.argv[2])
    else: main()
