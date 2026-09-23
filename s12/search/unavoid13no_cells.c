/* unavoid13-no: maximal cells of an arrangement of closed unit squares (float, lenient).
 *
 * Input: n squares (cx, cy, cos, sin), container [0,m]^2, tolerance tol.
 * Streams every arrangement vertex v (square corners and intersections of edges of two squares),
 * computes S(v) = {Q : dist(v, Q) <= tol} (lenient incidence, a superset of the truth), the
 * polygon P(v) = intersection of the squares of S(v) each enlarged by tol, and keeps v iff no
 * square outside S(v) (enlarged by tol) meets P(v).
 *
 * Why this is the dominance filter [proved for tol = 0]: v is dominated by an arrangement vertex
 * v' (S(v') a strict superset of S(v)) iff some square Q not in S(v) meets P(v): if Q meets P(v),
 * Q ∩ P(v) is a convex polygon inside every square of S(v) and inside Q, so any of its vertices
 * (arrangement vertices) dominates v; conversely v' in Q ∩ P(v).  With tol > 0 the test is
 * lenient in the safe direction (fewer drops).  Survivors are deduped on S(v).
 *
 * Build: gcc -O2 -fopenmp -shared -fPIC -o unavoid13no_cells.so unavoid13no_cells.c
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct { double x, y; } pt;

static int n_sq; static const double *CX, *CY, *CC, *SS; static double M, TOL;
static int ncell; static double hcell;
static int *cell_start, *cell_list;   /* squares by grid cell (of centre) */
static int NW;                        /* words per bitset */

static inline int cell_of(double x) { int c = (int)floor(x / hcell); if (c < 0) c = 0; if (c >= ncell) c = ncell - 1; return c; }

/* lenient containment: dist(v, Q) <= tol  <=>  |x'| <= 1/2 + tol and |y'| <= 1/2 + tol */
static inline int inside(int q, double px, double py, double tol) {
    double dx = px - CX[q], dy = py - CY[q];
    double xp = dx * CC[q] + dy * SS[q], yp = -dx * SS[q] + dy * CC[q];
    return fabs(xp) <= 0.5 + tol && fabs(yp) <= 0.5 + tol;
}

/* clip convex polygon (in/out arrays, n points) by half-plane a.x + b.y <= c ; returns new count */
static int clip(pt *in, int n, double a, double b, double c, pt *out) {
    int m = 0;
    for (int i = 0; i < n; i++) {
        pt P = in[i], Q = in[(i + 1) % n];
        double dp = a * P.x + b * P.y - c, dq = a * Q.x + b * Q.y - c;
        if (dp <= 0) out[m++] = P;
        if ((dp < 0 && dq > 0) || (dp > 0 && dq < 0)) {
            double t = dp / (dp - dq);
            pt R = { P.x + t * (Q.x - P.x), P.y + t * (Q.y - P.y) };
            out[m++] = R;
        }
    }
    return m;
}

/* clip polygon by square q enlarged by tol; poly buffers must hold >= 64 points */
static int clip_square(pt *poly, int n, int q, double tol, pt *tmp) {
    double c = CC[q], s = SS[q], cx = CX[q], cy = CY[q];
    double off1 = cx * c + cy * s, off2 = -cx * s + cy * c;
    n = clip(poly, n, c, s, off1 + 0.5 + tol, tmp);   if (n == 0) return 0;
    n = clip(tmp, n, -c, -s, -off1 + 0.5 + tol, poly); if (n == 0) return 0;
    n = clip(poly, n, -s, c, off2 + 0.5 + tol, tmp);  if (n == 0) return 0;
    n = clip(tmp, n, s, -c, -off2 + 0.5 + tol, poly);
    return n;
}

typedef struct {
    double *coords; uint64_t *bits; int cap, cnt;
} outbuf;

static void push(outbuf *ob, double x, double y, const uint64_t *S) {
    if (ob->cnt >= ob->cap) {
        ob->cap = ob->cap ? ob->cap * 2 : 4096;
        ob->coords = realloc(ob->coords, sizeof(double) * 2 * ob->cap);
        ob->bits = realloc(ob->bits, sizeof(uint64_t) * NW * ob->cap);
    }
    ob->coords[2 * ob->cnt] = x; ob->coords[2 * ob->cnt + 1] = y;
    memcpy(ob->bits + (size_t)NW * ob->cnt, S, sizeof(uint64_t) * NW);
    ob->cnt++;
}

/* process one vertex: returns 1 if kept (and pushes) */
static int process(double px, double py, outbuf *ob, uint64_t *S, int *members, int *others,
                   pt *poly, pt *tmp) {
    if (px < -TOL || px > M + TOL || py < -TOL || py > M + TOL) return 0;
    if (px < 0) px = 0; if (px > M) px = M; if (py < 0) py = 0; if (py > M) py = M;
    memset(S, 0, sizeof(uint64_t) * NW);
    int nm = 0, no = 0;
    int cx0 = cell_of(px), cy0 = cell_of(py);
    for (int a = cx0 - 1; a <= cx0 + 1; a++) for (int b = cy0 - 1; b <= cy0 + 1; b++) {
        if (a < 0 || b < 0 || a >= ncell || b >= ncell) continue;
        int cid = a * ncell + b;
        for (int k = cell_start[cid]; k < cell_start[cid + 1]; k++) {
            int q = cell_list[k];
            if (inside(q, px, py, TOL)) { members[nm++] = q; S[q >> 6] |= (1ULL << (q & 63)); }
            else others[no++] = q;
        }
    }
    if (nm == 0) return 0;
    /* P(v): start from the first member square enlarged by tol */
    int q0 = members[0];
    double c = CC[q0], s = SS[q0], h = 0.5 + 2 * TOL;
    poly[0].x = CX[q0] + (-h) * c - (-h) * s; poly[0].y = CY[q0] + (-h) * s + (-h) * c;
    poly[1].x = CX[q0] + ( h) * c - (-h) * s; poly[1].y = CY[q0] + ( h) * s + (-h) * c;
    poly[2].x = CX[q0] + ( h) * c - ( h) * s; poly[2].y = CY[q0] + ( h) * s + ( h) * c;
    poly[3].x = CX[q0] + (-h) * c - ( h) * s; poly[3].y = CY[q0] + (-h) * s + ( h) * c;
    int np = 4;
    for (int i = 1; i < nm; i++) {
        np = clip_square(poly, np, members[i], 2 * TOL, tmp);
        if (np == 0) { /* numerically empty: keep v conservatively */
            push(ob, px, py, S); return 1;
        }
        if (np > 60) { push(ob, px, py, S); return 1; }
    }
    /* bounding box of P(v) */
    double bx0 = 1e9, bx1 = -1e9, by0 = 1e9, by1 = -1e9;
    for (int i = 0; i < np; i++) {
        if (poly[i].x < bx0) bx0 = poly[i].x; if (poly[i].x > bx1) bx1 = poly[i].x;
        if (poly[i].y < by0) by0 = poly[i].y; if (poly[i].y > by1) by1 = poly[i].y;
    }
    /* any non-member square (enlarged by tol) meeting P(v)?  bbox prefilter (SAT on 4 axes), then clip */
    pt p2[64], t2[64];
    for (int i = 0; i < no; i++) {
        int q = others[i];
        double w2 = (fabs(CC[q]) + fabs(SS[q])) * 0.5 + TOL;
        double mx = (bx0 + bx1) * 0.5, my = (by0 + by1) * 0.5, hx = (bx1 - bx0) * 0.5, hy = (by1 - by0) * 0.5;
        if (fabs(mx - CX[q]) > w2 + hx) continue;
        if (fabs(my - CY[q]) > w2 + hy) continue;
        double c1 = CC[q], s1 = SS[q];
        if (fabs(mx * c1 + my * s1 - (CX[q] * c1 + CY[q] * s1)) > 0.5 + TOL + hx * fabs(c1) + hy * fabs(s1)) continue;
        if (fabs(-mx * s1 + my * c1 - (-CX[q] * s1 + CY[q] * c1)) > 0.5 + TOL + hx * fabs(s1) + hy * fabs(c1)) continue;
        memcpy(p2, poly, sizeof(pt) * np);
        int r = clip_square(p2, np, q, TOL, t2);
        if (r > 0) return 0;      /* dominated */
    }
    push(ob, px, py, S);
    return 1;
}

/* segment intersection (closed, lenient): returns 1 and the point */
static int seg_isect(pt A, pt B, pt C, pt D, double tol, pt *out) {
    double rx = B.x - A.x, ry = B.y - A.y, sx = D.x - C.x, sy = D.y - C.y;
    double den = rx * sy - ry * sx;
    if (fabs(den) < 1e-12) return 0;
    double qx = C.x - A.x, qy = C.y - A.y;
    double t = (qx * sy - qy * sx) / den, u = (qx * ry - qy * rx) / den;
    if (t < -tol || t > 1 + tol || u < -tol || u > 1 + tol) return 0;
    out->x = A.x + t * rx; out->y = A.y + t * ry;
    return 1;
}

/* main entry.  Returns number of survivors (before dedupe); fills *coords_out (2 per), *bits_out
 * (NW per), allocated with malloc; caller frees via free_buf. */
long maximal_cells(int n, const double *cx, const double *cy, const double *cc, const double *ss,
                   double m, double tol, int nthreads, double **coords_out, uint64_t **bits_out,
                   long *nvert_out) {
    n_sq = n; CX = cx; CY = cy; CC = cc; SS = ss; M = m; TOL = tol;
    NW = (n + 63) / 64;
    hcell = 0.75; ncell = (int)ceil(m / hcell);
    cell_start = calloc(ncell * ncell + 1, sizeof(int)); cell_list = malloc(sizeof(int) * n);
    int *cnt = calloc(ncell * ncell, sizeof(int));
    for (int i = 0; i < n; i++) cnt[cell_of(cx[i]) * ncell + cell_of(cy[i])]++;
    for (int c = 0; c < ncell * ncell; c++) cell_start[c + 1] = cell_start[c] + cnt[c];
    memset(cnt, 0, sizeof(int) * ncell * ncell);
    for (int i = 0; i < n; i++) { int c = cell_of(cx[i]) * ncell + cell_of(cy[i]); cell_list[cell_start[c] + cnt[c]++] = i; }
    free(cnt);
    /* corners */
    pt *corn = malloc(sizeof(pt) * 4 * n);
    for (int i = 0; i < n; i++) {
        double c = cc[i], s = ss[i];
        double sx[4] = {-0.5, 0.5, 0.5, -0.5}, sy[4] = {-0.5, -0.5, 0.5, 0.5};
        for (int k = 0; k < 4; k++) { corn[4 * i + k].x = cx[i] + sx[k] * c - sy[k] * s; corn[4 * i + k].y = cy[i] + sx[k] * s + sy[k] * c; }
    }
#ifdef _OPENMP
    if (nthreads > 0) omp_set_num_threads(nthreads);
    int T = omp_get_max_threads();
#else
    int T = 1;
#endif
    outbuf *obs = calloc(T, sizeof(outbuf));
    long nvert = 0;
#pragma omp parallel reduction(+:nvert)
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        uint64_t *S = malloc(sizeof(uint64_t) * NW);
        int *members = malloc(sizeof(int) * n), *others = malloc(sizeof(int) * n);
        pt *poly = malloc(sizeof(pt) * 128), *tmp = malloc(sizeof(pt) * 128);
#pragma omp for schedule(dynamic, 8)
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < 4; k++) { nvert++; process(corn[4 * i + k].x, corn[4 * i + k].y, &obs[tid], S, members, others, poly, tmp); }
            int cxi = cell_of(cx[i]), cyi = cell_of(cy[i]);
            for (int a = cxi - 2; a <= cxi + 2; a++) for (int b = cyi - 2; b <= cyi + 2; b++) {
                if (a < 0 || b < 0 || a >= ncell || b >= ncell) continue;
                int cid = a * ncell + b;
                for (int kk = cell_start[cid]; kk < cell_start[cid + 1]; kk++) {
                    int j = cell_list[kk];
                    if (j <= i) continue;
                    double dx = cx[i] - cx[j], dy = cy[i] - cy[j];
                    if (dx * dx + dy * dy > 2.0 + 1e-9) continue;
                    for (int ea = 0; ea < 4; ea++) for (int eb = 0; eb < 4; eb++) {
                        pt X;
                        if (seg_isect(corn[4 * i + ea], corn[4 * i + (ea + 1) % 4], corn[4 * j + eb], corn[4 * j + (eb + 1) % 4], 1e-9, &X)) {
                            nvert++; process(X.x, X.y, &obs[tid], S, members, others, poly, tmp);
                        }
                    }
                }
            }
        }
        free(S); free(members); free(others); free(poly); free(tmp);
    }
    long total = 0;
    for (int t = 0; t < T; t++) total += obs[t].cnt;
    double *coords = malloc(sizeof(double) * 2 * (total ? total : 1));
    uint64_t *bits = malloc(sizeof(uint64_t) * NW * (total ? total : 1));
    long pos = 0;
    for (int t = 0; t < T; t++) {
        if (obs[t].cnt) {
            memcpy(coords + 2 * pos, obs[t].coords, sizeof(double) * 2 * obs[t].cnt);
            memcpy(bits + (size_t)NW * pos, obs[t].bits, sizeof(uint64_t) * NW * obs[t].cnt);
        }
        pos += obs[t].cnt; free(obs[t].coords); free(obs[t].bits);
    }
    free(obs); free(corn); free(cell_start); free(cell_list);
    *coords_out = coords; *bits_out = bits; *nvert_out = nvert;
    return total;
}

void free_buf(void *p) { free(p); }
