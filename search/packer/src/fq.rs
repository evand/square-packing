//! fq: fast smooth quench for packings of n unit squares in a square (screening, f64).
//!
//! Formulation (min s subject to non-overlap), smooth by construction:
//! * every nearby pair (i, j) carries its own separating line: unit normal u = (cos phi, sin phi) and offset d, measured
//!   from the midpoint M of the two centres; constraints  u.(p - M) - d <= 0  for the 4 corners p of i and
//!   d - u.(q - M) <= 0  for the 4 corners q of j;
//! * every corner of every square inside [0, s]^2 (4 constraints per corner).
//! Each constraint is smooth in (x, y, theta, phi, d, s), so the PHR augmented Lagrangian
//!     L = s + sum_k (mu/2) max(0, g_k + lam_k/mu)^2 - lam_k^2/(2 mu)
//! is C^1: no max-over-axes kink (the reason L-BFGS on SAT penetration stalls at face-to-face contacts).
//! Inner: L-BFGS over (x, y, theta)_i, s, (phi, d)_pairs.  Outer: multiplier update, mu increase, Verlet pair list.
//! End: exact SAT check and repair by uniform scaling, so the written packing is feasible in f64.
//!
//! File format as packer: "n s", then "x y theta_degrees" per square.
//!
//!   fq quench --in A --out B [--loosen 1.02] [--kick SIGMA --seed K] [--mu0 10] [--v]
//!   fq check --in A           (min pair gap, min wall clearance)
//!   --dump-lp DIR: append every polish LP to DIR/lp-<pid>.jsonl (rows/columns with stable keys, Clarabel time; lp_bench.py)

use std::fs;
use std::time::Instant;
use std::sync::atomic::{AtomicU8, AtomicUsize, Ordering};
use std::sync::OnceLock;

// --dump-lp DIR: append every polish LP (rows, columns with stable keys, Clarabel time) to DIR/lp-<pid>.jsonl
static DUMP_LP: OnceLock<String> = OnceLock::new();
static LP_TAG: AtomicU8 = AtomicU8::new(0); // 0 step, 1 second-order correction, 2 flip search, 3 final force network
static LP_IT: AtomicUsize = AtomicUsize::new(0);
static IDENT_FPS: std::sync::Mutex<Vec<u64>> = std::sync::Mutex::new(Vec::new()); // --ident: load-network fingerprint per accepted step

const H: f64 = 0.5;
const SQRT2: f64 = std::f64::consts::SQRT_2;
const CS: [(f64, f64); 4] = [(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)];

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        let mut r = Rng(seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(0x632BE59BD9B4E019) | 1);
        for _ in 0..8 { r.next(); }
        r
    }
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn f(&mut self) -> f64 { (self.next() >> 11) as f64 / (1u64 << 53) as f64 }
    fn gauss(&mut self) -> f64 {
        let u1 = self.f().max(1e-300);
        let u2 = self.f();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------------------------------------------------------- geometry helpers
#[inline]
fn corner_off(c: f64, s: f64, k: usize) -> (f64, f64) {
    let (sx, sy) = CS[k];
    (H * (c * sx - s * sy), H * (s * sx + c * sy))
}

/// SAT gap of two unit squares (largest separating gap over the 4 face normals; negative = overlap depth).
fn sat_gap(xi: f64, yi: f64, ti: f64, xj: f64, yj: f64, tj: f64) -> f64 {
    let mut best = f64::NEG_INFINITY;
    let (dx, dy) = (xj - xi, yj - yi);
    for &t in &[ti, ti + std::f64::consts::FRAC_PI_2, tj, tj + std::f64::consts::FRAC_PI_2] {
        let (ax, ay) = (t.cos(), t.sin());
        let dc = (dx * ax + dy * ay).abs();
        let ri = H * ((ti - t).cos().abs() + (ti - t).sin().abs());
        let rj = H * ((tj - t).cos().abs() + (tj - t).sin().abs());
        best = best.max(dc - ri - rj);
    }
    best
}

/// Face separators of a pair (faces 0, 1 of i, then 0, 1 of j; normals oriented i -> j): (gap, phi, d).
fn face_seps(xi: f64, yi: f64, ti: f64, xj: f64, yj: f64, tj: f64) -> Vec<(f64, f64, f64)> {
    let (mx, my) = ((xi + xj) / 2.0, (yi + yj) / 2.0);
    let mut out = Vec::with_capacity(4);
    for &t0 in &[ti, ti + std::f64::consts::FRAC_PI_2, tj, tj + std::f64::consts::FRAC_PI_2] {
        let mut t = t0;
        if (xj - xi) * t.cos() + (yj - yi) * t.sin() < 0.0 { t += std::f64::consts::PI; }
        let (ux, uy) = (t.cos(), t.sin());
        let (mut ihi, mut jlo) = (f64::NEG_INFINITY, f64::INFINITY);
        for k in 0..4 {
            let (ox, oy) = corner_off(ti.cos(), ti.sin(), k);
            ihi = ihi.max(ux * (xi + ox - mx) + uy * (yi + oy - my));
            let (ox, oy) = corner_off(tj.cos(), tj.sin(), k);
            jlo = jlo.min(ux * (xj + ox - mx) + uy * (yj + oy - my));
        }
        out.push((jlo - ihi, t, (jlo + ihi) / 2.0));
    }
    out
}

/// Best initial separator for a pair: face normal (or centre line) with the largest gap, oriented i -> j.
fn init_sep(xi: f64, yi: f64, ti: f64, xj: f64, yj: f64, tj: f64) -> (f64, f64) {
    let (mx, my) = ((xi + xj) / 2.0, (yi + yj) / 2.0);
    let mut cands = vec![ti, ti + std::f64::consts::FRAC_PI_2, tj, tj + std::f64::consts::FRAC_PI_2];
    cands.push((yj - yi).atan2(xj - xi));
    let (mut bphi, mut bd, mut bg) = (0.0, 0.0, f64::NEG_INFINITY);
    for &t0 in &cands {
        for &t in &[t0, t0 + std::f64::consts::PI] {
            let (ux, uy) = (t.cos(), t.sin());
            if (xj - xi) * ux + (yj - yi) * uy < 0.0 { continue; }
            let (mut ihi, mut jlo) = (f64::NEG_INFINITY, f64::INFINITY);
            for k in 0..4 {
                let (ox, oy) = corner_off(ti.cos(), ti.sin(), k);
                ihi = ihi.max(ux * (xi + ox - mx) + uy * (yi + oy - my));
                let (ox, oy) = corner_off(tj.cos(), tj.sin(), k);
                jlo = jlo.min(ux * (xj + ox - mx) + uy * (yj + oy - my));
            }
            if jlo - ihi > bg { bg = jlo - ihi; bphi = t; bd = (jlo + ihi) / 2.0; }
        }
    }
    (bphi, bd)
}

// ---------------------------------------------------------------- problem
#[derive(Clone)]
struct Pair { i: usize, j: usize, phi: f64, d: f64, lam: [f64; 8], face: usize }

struct Prob {
    n: usize,
    pairs: Vec<Pair>,
    wlam: Vec<f64>, // 16 per square: corner k, wall w (x>=0, x<=s, y>=0, y<=s)
    mu: f64,
    fix_s: bool, // fixed side: pure overlap penalty, no s objective, d/ds = 0
}

#[inline]
fn term(g: f64, lam: f64, mu: f64) -> (f64, f64) {
    let t = g + lam / mu;
    if t > 0.0 { (0.5 * mu * t * t - lam * lam / (2.0 * mu), mu * t) } else { (-lam * lam / (2.0 * mu), 0.0) }
}

impl Prob {
    // variable layout: 3n (x, y, theta), then s, then 2 per pair (phi, d)
    fn pack(&self, x: &[f64], s: f64) -> Vec<f64> {
        let mut v = x.to_vec();
        v.push(s);
        for p in &self.pairs { v.push(p.phi); v.push(p.d); }
        v
    }
    fn unpack(&mut self, v: &[f64]) {
        let o = 3 * self.n + 1;
        for (k, p) in self.pairs.iter_mut().enumerate() { p.phi = v[o + 2 * k]; p.d = v[o + 2 * k + 1]; }
    }

    /// Augmented Lagrangian value and gradient.
    fn eval(&self, v: &[f64], g: &mut [f64]) -> f64 {
        let n = self.n;
        let mu = self.mu;
        for q in g.iter_mut() { *q = 0.0; }
        let s = v[3 * n];
        let mut f = if self.fix_s { 0.0 } else { s };
        g[3 * n] = if self.fix_s { 0.0 } else { 1.0 };
        let mut cs = vec![(0.0, 0.0); n];
        for i in 0..n { let t = v[3 * i + 2]; cs[i] = (t.cos(), t.sin()); }
        // walls
        for i in 0..n {
            let (c, sn) = cs[i];
            let (x, y) = (v[3 * i], v[3 * i + 1]);
            for k in 0..4 {
                let (ox, oy) = corner_off(c, sn, k);
                let (px, py) = (x + ox, y + oy);
                let gs = [-px, px - s, -py, py - s];
                for w in 0..4 {
                    let (fv, dg) = term(gs[w], self.wlam[16 * i + 4 * k + w], mu);
                    f += fv;
                    if dg != 0.0 {
                        // dpx/dtheta = -oy, dpy/dtheta = ox
                        match w {
                            0 => { g[3 * i] -= dg; g[3 * i + 2] += dg * oy; }
                            1 => { g[3 * i] += dg; g[3 * i + 2] -= dg * oy; g[3 * n] -= dg; }
                            2 => { g[3 * i + 1] -= dg; g[3 * i + 2] -= dg * ox; }
                            _ => { g[3 * i + 1] += dg; g[3 * i + 2] += dg * ox; g[3 * n] -= dg; }
                        }
                    }
                }
            }
        }
        // pairs
        let o = 3 * n + 1;
        for (pk, p) in self.pairs.iter().enumerate() {
            let (i, j) = (p.i, p.j);
            let phi = v[o + 2 * pk];
            let d = v[o + 2 * pk + 1];
            let (ux, uy) = (phi.cos(), phi.sin());
            let (hx, hy) = ((v[3 * i] - v[3 * j]) / 2.0, (v[3 * i + 1] - v[3 * j + 1]) / 2.0);
            let (mut gx, mut gti, mut gtj, mut gphi, mut gd) = (0.0, 0.0, 0.0, 0.0, 0.0);
            let (ci, si) = cs[i];
            let (cj, sj) = cs[j];
            for k in 0..4 {
                // corner of i: r = h + o_i;  g = u.r - d
                let (ox, oy) = corner_off(ci, si, k);
                let (rx, ry) = (hx + ox, hy + oy);
                let (fv, dg) = term(ux * rx + uy * ry - d, p.lam[k], mu);
                f += fv;
                if dg != 0.0 {
                    gx += dg;
                    gti += dg * (-ux * oy + uy * ox);
                    gphi += dg * (-uy * rx + ux * ry);
                    gd -= dg;
                }
                // corner of j: r' = -h + o_j;  g = d - u.r'
                let (ox, oy) = corner_off(cj, sj, k);
                let (rx, ry) = (-hx + ox, -hy + oy);
                let (fv, dg2) = term(d - (ux * rx + uy * ry), p.lam[4 + k], mu);
                f += fv;
                if dg2 != 0.0 {
                    gx += dg2;
                    gtj -= dg2 * (-ux * oy + uy * ox);
                    gphi -= dg2 * (-uy * rx + ux * ry);
                    gd += dg2;
                }
            }
            // both corner kinds: dg/dXi = u/2, dg/dXj = -u/2 (gx holds the summed dg)
            if gx != 0.0 {
                g[3 * i] += 0.5 * gx * ux; g[3 * i + 1] += 0.5 * gx * uy;
                g[3 * j] -= 0.5 * gx * ux; g[3 * j + 1] -= 0.5 * gx * uy;
            }
            g[3 * i + 2] += gti;
            g[3 * j + 2] += gtj;
            g[o + 2 * pk] += gphi;
            g[o + 2 * pk + 1] += gd;
        }
        if self.fix_s { g[3 * n] = 0.0; }
        f
    }

    /// Max constraint violation, and multiplier update lam <- max(0, lam + mu g).
    fn viol_update(&mut self, v: &[f64], update: bool) -> f64 {
        let n = self.n;
        let mu = self.mu;
        let s = v[3 * n];
        let mut mv: f64 = 0.0;
        for i in 0..n {
            let t = v[3 * i + 2];
            let (c, sn) = (t.cos(), t.sin());
            for k in 0..4 {
                let (ox, oy) = corner_off(c, sn, k);
                let (px, py) = (v[3 * i] + ox, v[3 * i + 1] + oy);
                let gs = [-px, px - s, -py, py - s];
                for w in 0..4 {
                    mv = mv.max(gs[w]);
                    if update { let l = &mut self.wlam[16 * i + 4 * k + w]; *l = (*l + mu * gs[w]).max(0.0); }
                }
            }
        }
        let o = 3 * n + 1;
        for (pk, p) in self.pairs.iter_mut().enumerate() {
            let (i, j) = (p.i, p.j);
            let (phi, d) = (v[o + 2 * pk], v[o + 2 * pk + 1]);
            let (ux, uy) = (phi.cos(), phi.sin());
            let (hx, hy) = ((v[3 * i] - v[3 * j]) / 2.0, (v[3 * i + 1] - v[3 * j + 1]) / 2.0);
            let (ti, tj) = (v[3 * i + 2], v[3 * j + 2]);
            for k in 0..4 {
                let (ox, oy) = corner_off(ti.cos(), ti.sin(), k);
                let gi = ux * (hx + ox) + uy * (hy + oy) - d;
                let (ox, oy) = corner_off(tj.cos(), tj.sin(), k);
                let gj = d - (ux * (-hx + ox) + uy * (-hy + oy));
                mv = mv.max(gi).max(gj);
                if update {
                    p.lam[k] = (p.lam[k] + mu * gi).max(0.0);
                    p.lam[4 + k] = (p.lam[4 + k] + mu * gj).max(0.0);
                }
            }
        }
        mv
    }

    /// Rebuild the pair list at cutoff; keep separators and multipliers of surviving pairs.
    fn rebuild(&mut self, x: &[f64], cut: f64) {
        let n = self.n;
        let mut old = std::collections::HashMap::new();
        for p in self.pairs.drain(..) { old.insert((p.i, p.j), p); }
        let c2 = cut * cut;
        for i in 0..n {
            for j in i + 1..n {
                let (dx, dy) = (x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1]);
                if dx * dx + dy * dy >= c2 { continue; }
                if let Some(p) = old.remove(&(i, j)) { self.pairs.push(p); continue; }
                let (phi, d) = init_sep(x[3 * i], x[3 * i + 1], x[3 * i + 2], x[3 * j], x[3 * j + 1], x[3 * j + 2]);
                self.pairs.push(Pair { i, j, phi, d, lam: [0.0; 8], face: 99 });
            }
        }
    }
}


// ---------------------------------------------------------------- SLP polish (trust region, lifted formulation)
#[derive(Clone, Copy)]
enum Spec { Wall(usize, usize, usize), PI(usize, usize), PJ(usize, usize) }
struct Con { g: f64, nz: Vec<(usize, f64)>, sp: Spec }

impl Prob {
    /// Constraints with g > -delta and their sparse gradients (variable layout as pack()).
    fn near(&self, v: &[f64], delta: f64) -> Vec<Con> {
        let n = self.n;
        let si = 3 * n;
        let s = v[si];
        let mut out = Vec::new();
        for i in 0..n {
            let t = v[3 * i + 2];
            let (c, sn) = (t.cos(), t.sin());
            for k in 0..4 {
                let (ox, oy) = corner_off(c, sn, k);
                let (px, py) = (v[3 * i] + ox, v[3 * i + 1] + oy);
                if -px > -delta { out.push(Con { g: -px, nz: vec![(3 * i, -1.0), (3 * i + 2, oy)], sp: Spec::Wall(i, k, 0) }); }
                if px - s > -delta { out.push(Con { g: px - s, nz: vec![(3 * i, 1.0), (3 * i + 2, -oy), (si, -1.0)], sp: Spec::Wall(i, k, 1) }); }
                if -py > -delta { out.push(Con { g: -py, nz: vec![(3 * i + 1, -1.0), (3 * i + 2, -ox)], sp: Spec::Wall(i, k, 2) }); }
                if py - s > -delta { out.push(Con { g: py - s, nz: vec![(3 * i + 1, 1.0), (3 * i + 2, ox), (si, -1.0)], sp: Spec::Wall(i, k, 3) }); }
            }
        }
        let o = 3 * n + 1;
        for (pk, p) in self.pairs.iter().enumerate() {
            let (i, j) = (p.i, p.j);
            let (phi, d) = (v[o + 2 * pk], v[o + 2 * pk + 1]);
            let (ux, uy) = (phi.cos(), phi.sin());
            let (hx, hy) = ((v[3 * i] - v[3 * j]) / 2.0, (v[3 * i + 1] - v[3 * j + 1]) / 2.0);
            let (ti, tj) = (v[3 * i + 2], v[3 * j + 2]);
            let xs = [(3 * i, 0.5 * ux), (3 * i + 1, 0.5 * uy), (3 * j, -0.5 * ux), (3 * j + 1, -0.5 * uy)];
            // separator attached to a face: phi moves with the owner's theta
            let phicol = match p.face { 0 | 1 => 3 * i + 2, 2 | 3 => 3 * j + 2, _ => o + 2 * pk };
            for k in 0..4 {
                let (ox, oy) = corner_off(ti.cos(), ti.sin(), k);
                let (rx, ry) = (hx + ox, hy + oy);
                let g = ux * rx + uy * ry - d;
                if g > -delta {
                    let mut nz = xs.to_vec();
                    nz.push((3 * i + 2, -ux * oy + uy * ox));
                    nz.push((phicol, -uy * rx + ux * ry));
                    nz.push((o + 2 * pk + 1, -1.0));
                    out.push(Con { g, nz, sp: Spec::PI(pk, k) });
                }
                let (ox, oy) = corner_off(tj.cos(), tj.sin(), k);
                let (rx, ry) = (-hx + ox, -hy + oy);
                let g = d - (ux * rx + uy * ry);
                if g > -delta {
                    let mut nz = xs.to_vec();
                    nz.push((3 * j + 2, ux * oy - uy * ox));
                    nz.push((phicol, uy * rx - ux * ry));
                    nz.push((o + 2 * pk + 1, 1.0));
                    out.push(Con { g, nz, sp: Spec::PJ(pk, k) });
                }
            }
        }
        out
    }
}

impl Prob {
    /// Apply the face tie to a trial point: phi of a face-attached separator moves with its owner's theta.
    fn tie(&self, v0: &[f64], v: &mut [f64]) {
        let o = 3 * self.n + 1;
        for (pk, p) in self.pairs.iter().enumerate() {
            let own = match p.face { 0 | 1 => p.i, 2 | 3 => p.j, _ => continue };
            v[o + 2 * pk] = v0[o + 2 * pk] + (v[3 * own + 2] - v0[3 * own + 2]);
        }
    }
    fn g_at(&self, v: &[f64], sp: Spec) -> f64 {
        let n = self.n;
        match sp {
            Spec::Wall(i, k, w) => {
                let t = v[3 * i + 2];
                let (ox, oy) = corner_off(t.cos(), t.sin(), k);
                let (px, py, s) = (v[3 * i] + ox, v[3 * i + 1] + oy, v[3 * n]);
                [-px, px - s, -py, py - s][w]
            }
            Spec::PI(pk, k) | Spec::PJ(pk, k) => {
                let p = &self.pairs[pk];
                let (i, j) = (p.i, p.j);
                let o = 3 * n + 1;
                let (phi, d) = (v[o + 2 * pk], v[o + 2 * pk + 1]);
                let (ux, uy) = (phi.cos(), phi.sin());
                let (hx, hy) = ((v[3 * i] - v[3 * j]) / 2.0, (v[3 * i + 1] - v[3 * j + 1]) / 2.0);
                if let Spec::PI(..) = sp {
                    let t = v[3 * i + 2];
                    let (ox, oy) = corner_off(t.cos(), t.sin(), k);
                    ux * (hx + ox) + uy * (hy + oy) - d
                } else {
                    let t = v[3 * j + 2];
                    let (ox, oy) = corner_off(t.cos(), t.sin(), k);
                    d - (ux * (-hx + ox) + uy * (-hy + oy))
                }
            }
        }
    }
}

/// One LP: min w_s + kappa tau  s.t.  grad g . w - tau <= -g/R (near constraints), |w| <= 1, tau >= 0.
/// Only columns that occur in some row (and s) are LP variables; the rest of w is 0.  Returns (w, tau, row duals).
fn slp_lp(nv: usize, si: usize, cons: &[Con], r: f64, kappa: f64, hq: &[(usize, usize, f64)], kp: &Prob) -> Option<(Vec<f64>, f64, Vec<f64>)> {
    use clarabel::algebra::*;
    use clarabel::solver::*;
    let mut col = vec![usize::MAX; nv];
    let mut used = Vec::new();
    let mut take = |k: usize, col: &mut Vec<usize>| { if col[k] == usize::MAX { col[k] = used.len(); used.push(k); } };
    take(si, &mut col);
    for c in cons { for &(k, _) in &c.nz { take(k, &mut col); } }
    for &(a, b, _) in hq { take(a, &mut col); take(b, &mut col); }
    let nu = used.len();
    let nx = nu + 1; // + tau
    let (mut ii, mut jj, mut vv) = (Vec::new(), Vec::new(), Vec::new());
    let mut b = Vec::new();
    let mut row = 0;
    for c in cons {
        for &(k, a) in &c.nz { ii.push(row); jj.push(col[k]); vv.push(a); }
        ii.push(row); jj.push(nu); vv.push(-1.0);
        b.push(-c.g / r);
        row += 1;
    }
    for k in 0..nu {
        ii.push(row); jj.push(k); vv.push(1.0); b.push(1.0); row += 1;
        ii.push(row); jj.push(k); vv.push(-1.0); b.push(1.0); row += 1;
    }
    ii.push(row); jj.push(nu); vv.push(-1.0); b.push(0.0); row += 1;
    let a = CscMatrix::new_from_triplets(row, nx, ii, jj, vv);
    let p = if hq.is_empty() { CscMatrix::<f64>::zeros((nx, nx)) } else {
        // SQP: objective + (r/2) w^T H w in the scaled variables (H convexified by the caller), upper triangle
        let (mut pi, mut pj, mut pv) = (Vec::new(), Vec::new(), Vec::new());
        for &(a, b, val) in hq {
            let (ca, cb) = (col[a].min(col[b]), col[a].max(col[b]));
            pi.push(ca); pj.push(cb); pv.push(r * val);
        }
        CscMatrix::new_from_triplets(nx, nx, pi, pj, pv)
    };
    let mut q = vec![0.0; nx];
    q[col[si]] = 1.0;
    q[nu] = kappa;
    let cones = [NonnegativeConeT(row)];
    let settings = DefaultSettingsBuilder::default().verbose(false).max_iter(200).build().unwrap();
    let t0 = Instant::now();
    let mut solver = DefaultSolver::new(&p, &q, &a, &b, &cones, settings).ok()?;
    solver.solve();
    if std::env::var("FQ_LPT").is_ok() { eprintln!("lp rows {row} cols {nx} iters {} sec {:.3}", solver.info.iterations, t0.elapsed().as_secs_f64()); }
    if let Some(dir) = DUMP_LP.get() { dump_lp(dir, kp, &used, cons, &b, kappa, r, row, nx, &solver, t0.elapsed().as_secs_f64()); }
    match solver.solution.status {
        SolverStatus::Solved | SolverStatus::AlmostSolved => {
            let x = &solver.solution.x;
            let mut w = vec![0.0; nv];
            for (c, &k) in used.iter().enumerate() { w[k] = x[c]; }
            Some((w, x[nu].max(0.0), solver.solution.z[..cons.len()].to_vec()))
        }
        _ => None,
    }
}

fn row_key(p: &Prob, v: &[f64], c: &Con) -> (u32, u32, u8, u8) {
    match c.sp {
        Spec::Wall(i, k, w) => (i as u32, 1_000_000 + w as u32, 99, k as u8),
        Spec::PI(pk, k) => { let q = &p.pairs[pk]; (q.i as u32, q.j as u32, canon_face(v, q), k as u8) }
        Spec::PJ(pk, k) => { let q = &p.pairs[pk]; (q.i as u32, q.j as u32, canon_face(v, q), 4 + k as u8) }
    }
}

/// Load-bearing contacts: rows with LP dual > 1e-6 max dual (interior point = maximal-support dual = force network).
fn load_keys(p: &Prob, v: &[f64], cons: &[Con], z: &[f64]) -> Vec<(u32, u32, u8, u8)> {
    let zm = z.iter().cloned().fold(0.0, f64::max);
    if std::env::var("FQ_ZHIST").is_ok() {
        let mut h = [0usize; 16];
        for &zz in z { let e = if zz <= 0.0 { 15 } else { ((-(zz / zm).log10()).floor() as usize).min(15) }; h[e] += 1; }
        eprintln!("zhist (decades below max) {:?}", h);
    }
    // contact level: drop the corner index (which corners of a face contact carry load is not basin-invariant)
    let mut out: Vec<_> = cons.iter().zip(z).filter(|(_, &zz)| zz > 1e-6 * zm).map(|(c, _)| { let k = row_key(p, v, c); (k.0, k.1, k.2, 0u8) }).collect();
    out.sort();
    out.dedup();
    out
}

/// Active contacts at v (g > -tol): keys (i, j or 1_000_000 + wall, face (99 = wall), corner (0-3 of i, 4-7 of j)).
fn contact_keys(p: &Prob, v: &[f64], tol: f64) -> Vec<(u32, u32, u8, u8)> {
    let mut out: Vec<(u32, u32, u8, u8)> = p.near(v, tol).iter().map(|c| match c.sp {
        Spec::Wall(i, k, w) => (i as u32, 1_000_000 + w as u32, 99, k as u8),
        Spec::PI(pk, k) => { let q = &p.pairs[pk]; (q.i as u32, q.j as u32, canon_face(v, q), k as u8) }
        Spec::PJ(pk, k) => { let q = &p.pairs[pk]; (q.i as u32, q.j as u32, canon_face(v, q), 4 + k as u8) }
    }).collect();
    out.sort();
    out
}

/// Face label independent of tie-breaks: the lowest face index parallel to the chosen face; 98 when two non-parallel faces
/// are both within 1e-7 of contact (corner-corner: either branch).
fn canon_face(v: &[f64], q: &Pair) -> u8 {
    let (i, j) = (q.i, q.j);
    let f = face_seps(v[3 * i], v[3 * i + 1], v[3 * i + 2], v[3 * j], v[3 * j + 1], v[3 * j + 2]);
    let par = |a: usize, b: usize| {
        let dd = (f[a].1 - f[b].1).rem_euclid(2.0 * std::f64::consts::PI);
        dd.min(2.0 * std::f64::consts::PI - dd) < 1e-6
    };
    let c = q.face.min(3);
    let best = (0..4).map(|k| f[k].0).fold(f64::NEG_INFINITY, f64::max);
    if (0..4).any(|k| !par(k, c) && f[k].0 > best - 1e-7 && f[k].0 > -1e-7) && f[c].0 > -1e-7 { return 98; }
    (0..4).find(|&k| par(k, c)).unwrap() as u8
}

fn fingerprint(keys: &[(u32, u32, u8, u8)]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &(a, b, c, d) in keys {
        for x in [a as u64, b as u64, c as u64, d as u64] { h ^= x; h = h.wrapping_mul(0x100000001b3); }
    }
    h
}

/// Identification statistics of an accepted step from v0 to v1: (max corner displacement over squares, smallest slack of
/// near-but-inactive constraints (g in (-0.05, -tol]), number of active, symmetric difference with `prev` active set).
fn ident(p: &Prob, v0: &[f64], v1: &[f64], prev: &mut Vec<(u32, u32, u8, u8)>, tol: f64) -> (f64, f64, usize, usize) {
    let n = p.n;
    let mut dm: f64 = 0.0;
    for i in 0..n {
        let d = (v1[3 * i] - v0[3 * i]).hypot(v1[3 * i + 1] - v0[3 * i + 1]) + H * SQRT2 * (v1[3 * i + 2] - v0[3 * i + 2]).abs();
        dm = dm.max(d);
    }
    let slack = p.near(v1, 0.05).iter().filter(|c| c.g <= -tol).map(|c| -c.g).fold(f64::INFINITY, f64::min);
    let act = contact_keys(p, v1, tol);
    let a: std::collections::HashSet<_> = act.iter().collect();
    let b: std::collections::HashSet<_> = prev.iter().collect();
    let diff = a.symmetric_difference(&b).count();
    let na = act.len();
    *prev = act;
    (dm, slack, na, diff)
}

/// Hessian of sum_k z_k g_k at v (rows from `cons`, multipliers `z`), with face-tied separators (phi -> owner theta),
/// convexified by diagonal dominance (H_aa >= sum_b |H_ab|).  Returns the upper triangle (a <= b) in original indices.
fn hess_lagr(p: &Prob, v: &[f64], cons: &[Con], z: &[f64]) -> Vec<(usize, usize, f64)> {
    use std::collections::HashMap;
    let n = p.n;
    let o = 3 * n + 1;
    let zm = z.iter().cloned().fold(0.0, f64::max);
    let mut h: HashMap<(usize, usize), f64> = HashMap::new();
    let mut add = |a: usize, b: usize, val: f64, h: &mut HashMap<(usize, usize), f64>| {
        if a == b { *h.entry((a, a)).or_insert(0.0) += val; }
        else { *h.entry((a, b)).or_insert(0.0) += val; *h.entry((b, a)).or_insert(0.0) += val; }
    };
    for (c, &zk) in cons.iter().zip(z) {
        if !(zk > 1e-9 * zm) { continue; }
        match c.sp {
            Spec::Wall(i, k, w) => {
                let th = v[3 * i + 2];
                let (ox, oy) = corner_off(th.cos(), th.sin(), k);
                let d2 = [ox, -ox, oy, -oy][w];
                add(3 * i + 2, 3 * i + 2, zk * d2, &mut h);
            }
            Spec::PI(pk, k) | Spec::PJ(pk, k) => {
                let pr = &p.pairs[pk];
                let (i, j) = (pr.i, pr.j);
                let own = match pr.face { 0 | 1 => i, _ => j };
                let phc = 3 * own + 2;
                let map = |a: usize| if a == usize::MAX { phc } else { a };
                const PHI: usize = usize::MAX;
                let phi = v[o + 2 * pk];
                let (ux, uy) = (phi.cos(), phi.sin());
                let (dx, dy) = (-uy, ux); // u'
                let (hx, hy) = ((v[3 * i] - v[3 * j]) / 2.0, (v[3 * i + 1] - v[3 * j + 1]) / 2.0);
                let mut e: Vec<(usize, usize, f64)> = Vec::new();
                let isi = matches!(c.sp, Spec::PI(..));
                let (sq_i, sgn) = if isi { (i, 1.0) } else { (j, -1.0) };
                let th = v[3 * sq_i + 2];
                let (ox, oy) = corner_off(th.cos(), th.sin(), k);
                let (rx, ry) = if isi { (hx + ox, hy + oy) } else { (-hx + ox, -hy + oy) };
                let ts = 3 * sq_i + 2;
                e.push((ts, ts, -sgn * (ux * ox + uy * oy)));            // theta theta
                e.push((PHI, PHI, -sgn * (ux * rx + uy * ry)));          // phi phi
                e.push((PHI, ts, sgn * (dx * (-oy) + dy * ox)));         // phi theta
                for (col, val) in [(3 * i, dx / 2.0), (3 * i + 1, dy / 2.0), (3 * j, -dx / 2.0), (3 * j + 1, -dy / 2.0)] {
                    e.push((PHI, col, val));
                }
                for (a, b, val) in e {
                    let (ma, mb) = (map(a), map(b));
                    if a != b && ma == mb { add(ma, ma, 2.0 * zk * val, &mut h); } else { add(ma, mb, zk * val, &mut h); }
                }
            }
        }
    }
    // convexify: diagonal dominance
    let mut off: HashMap<usize, f64> = HashMap::new();
    for (&(a, b), &val) in &h { if a != b { *off.entry(a).or_insert(0.0) += val.abs(); } }
    let mut out = Vec::new();
    let mut diag: HashMap<usize, f64> = HashMap::new();
    for (&(a, b), &val) in &h { if a == b { diag.insert(a, val); } }
    for (&a, &s) in &off { let d = diag.entry(a).or_insert(0.0); if *d < s { *d = s; } }
    for (&a, &d) in &diag { if d > 0.0 { out.push((a, a, d)); } }
    for (&(a, b), &val) in &h { if a < b && val != 0.0 { out.push((a, b, val)); } }
    out
}

struct POpt { sqp: bool, rmax: f64, r0: f64, rmin: f64, maxit: usize, cc_tol: f64, flip_top: usize, ident: bool, verbose: bool, stag_w: usize, stag_tol: f64 }

/// Snap each pair's separator to a face branch: `choice[pk]` = Some(face) forces that face while it stays ambiguous,
/// None = the face with the largest gap.  Returns the ambiguous pairs (>= 2 faces within cc_tol of contact and of the
/// best) with their candidate faces; forced choices that are no longer ambiguous are cleared.
fn dump_lp(dir: &str, p: &Prob, used: &[usize], cons: &[Con], b: &[f64], kappa: f64, r: f64, nrow: usize, nx: usize,
           solver: &clarabel::solver::DefaultSolver<f64>, sec: f64) {
    use std::io::Write;
    let n = p.n;
    let o = 3 * n + 1;
    let ck = |k: usize| -> String {
        if k < 3 * n { format!("{}{}", ["x", "y", "t"][k % 3], k / 3) } else if k == 3 * n { "s".into() } else {
            let q = &p.pairs[(k - o) / 2];
            format!("{}{}_{}", if (k - o) % 2 == 0 { "f" } else { "d" }, q.i, q.j)
        }
    };
    let rk = |c: &Con| -> String {
        match c.sp {
            Spec::Wall(i, k, w) => format!("w{i}_{k}_{w}"),
            Spec::PI(pk, k) => format!("a{}_{}_{k}", p.pairs[pk].i, p.pairs[pk].j),
            Spec::PJ(pk, k) => format!("b{}_{}_{k}", p.pairs[pk].i, p.pairs[pk].j),
        }
    };
    let mut col = std::collections::HashMap::new();
    for (c, &k) in used.iter().enumerate() { col.insert(k, c); }
    let mut s = String::with_capacity(64 * cons.len());
    s.push_str(&format!("{{\"n\":{n},\"tag\":{},\"it\":{},\"r\":{:e},\"kappa\":{:e},\"nrow_total\":{nrow},\"nx\":{nx},\"cols\":[",
        LP_TAG.load(Ordering::Relaxed), LP_IT.load(Ordering::Relaxed), r, kappa));
    for (c, &k) in used.iter().enumerate() { if c > 0 { s.push(','); } s.push_str(&format!("\"{}\"", ck(k))); }
    s.push_str("],\"rows\":[");
    for (ri, c) in cons.iter().enumerate() {
        if ri > 0 { s.push(','); }
        s.push_str(&format!("[\"{}\",{:e},[", rk(c), b[ri]));
        for (m, &(k, a)) in c.nz.iter().enumerate() { if m > 0 { s.push(','); } s.push_str(&format!("[{},{:e}]", col[&k], a)); }
        s.push_str("]]");
    }
    let x = &solver.solution.x;
    let obj = x[col[&(3 * n)]] + kappa * x[nx - 1];
    s.push_str(&format!("],\"clarabel\":{{\"status\":\"{:?}\",\"iters\":{},\"sec\":{:e},\"solve_time\":{:e},\"obj\":{:e}}}}}\n",
        solver.solution.status, solver.info.iterations, sec, solver.info.solve_time, obj));
    let path = format!("{dir}/lp-{}.jsonl", std::process::id());
    let mut f = fs::OpenOptions::new().create(true).append(true).open(path).expect("dump-lp");
    f.write_all(s.as_bytes()).unwrap();
}

fn snap(p: &mut Prob, x: &[f64], choice: &mut [Option<usize>], cc_tol: f64) -> Vec<(usize, Vec<usize>)> {
    let mut amb = Vec::new();
    for (pk, pr) in p.pairs.iter_mut().enumerate() {
        let (i, j) = (pr.i, pr.j);
        let f = face_seps(x[3 * i], x[3 * i + 1], x[3 * i + 2], x[3 * j], x[3 * j + 1], x[3 * j + 2]);
        let best = (0..4).max_by(|&a, &b| f[a].0.partial_cmp(&f[b].0).unwrap()).unwrap();
        // candidate branches: faces near contact and near the best, distinct in direction (parallel faces = one branch)
        let mut cand: Vec<usize> = vec![best];
        for q in 0..4 {
            if q == best || !(f[q].0 > f[best].0 - cc_tol && f[q].0 > -cc_tol) { continue; }
            let par = cand.iter().any(|&c| {
                let dd = (f[q].1 - f[c].1).rem_euclid(2.0 * std::f64::consts::PI);
                dd.min(2.0 * std::f64::consts::PI - dd) < 1e-6
            });
            if !par { cand.push(q); }
        }
        let mut c = best;
        if let Some(fc) = choice[pk] {
            if cand.contains(&fc) { c = fc; } else { choice[pk] = None; }
        }
        pr.phi = f[c].1;
        pr.d = f[c].2;
        pr.face = c;
        if cand.len() >= 2 {
            if std::env::var("FQ_DEBUG").is_ok() && amb.len() < 6 {
                eprintln!("amb pair {i} {j}: th {:.6} {:.6} faces {:?}", x[3*i+2].to_degrees(), x[3*j+2].to_degrees(), f.iter().map(|q| (q.0, q.1.to_degrees())).collect::<Vec<_>>());
            }
            amb.push((pk, cand));
        }
    }
    amb
}

/// Trust-region SLP on the lifted problem with face-branch separators (snapped each iteration), and a greedy branch
/// search over ambiguous (corner-corner-like) pairs when the step stalls.
/// Returns (s, x, iterations, final R, flips accepted).
fn polish(p: &mut Prob, s0: f64, x0: &[f64], o: &POpt, trace: &mut Vec<(usize, f64, f64, usize, f64, f64, usize, usize)>, p_load: &mut Vec<(u32, u32, u8, u8)>) -> (f64, Vec<f64>, usize, f64, usize) {
    let tp = Instant::now();
    let mut prev_act: Vec<(u32, u32, u8, u8)> = Vec::new();
    let n = p.n;
    let mut x = x0.to_vec();
    let mut s = s0;
    let mut r = o.r0;
    let kappa = 3.0 * s0;
    let mut it = 0;
    let mut flips = 0;
    p.rebuild(&x, SQRT2 + 0.1);
    let mut choice: Vec<Option<usize>> = vec![None; p.pairs.len()];
    let mut keys: Vec<(usize, usize)> = p.pairs.iter().map(|q| (q.i, q.j)).collect();
    let mut xref = x.clone();
    let merit = |p: &mut Prob, v: &[f64]| -> f64 { v[3 * p.n] + kappa * p.viol_update(v, false).max(0.0) };
    // LP at the current point with the current branch choice: (pred, w, tau, ncons)
    let try_lp = |p: &mut Prob, x: &[f64], s: f64, choice: &mut [Option<usize>], r: f64, hq: &[(usize, usize, f64)]| -> Option<(f64, Vec<f64>, Vec<Con>, f64, Vec<f64>)> {
        snap(p, x, choice, o.cc_tol);
        let v0 = p.pack(x, s);
        let m0 = merit(p, &v0);
        let cons = p.near(&v0, 5.0 * r + 1e-9);
        let (w, tau, z) = slp_lp(v0.len(), 3 * p.n, &cons, r, kappa, hq, p)?;
        Some((m0 - (s + r * w[3 * p.n] + kappa * r * tau), w, cons, m0, z))
    };
    let mut stalls = 0;
    let mut hq: Vec<(usize, usize, f64)> = Vec::new();
    let mut hist: Vec<f64> = Vec::new();
    while it < o.maxit {
        it += 1;
        hist.push(s);
        if o.stag_w > 0 && hist.len() > o.stag_w && hist[hist.len() - 1 - o.stag_w] - s < o.stag_tol { break; }
        let mut dmax: f64 = 0.0;
        for i in 0..n { dmax = dmax.max((x[3 * i] - xref[3 * i]).hypot(x[3 * i + 1] - xref[3 * i + 1])); }
        if dmax > 0.04 {
            // rebuild, carrying branch choices by key
            let old: std::collections::HashMap<(usize, usize), Option<usize>> = keys.iter().cloned().zip(choice.iter().cloned()).collect();
            p.rebuild(&x, SQRT2 + 0.1);
            keys = p.pairs.iter().map(|q| (q.i, q.j)).collect();
            choice = keys.iter().map(|k| *old.get(k).unwrap_or(&None)).collect();
            hq.clear();
            xref = x.clone();
        }
        LP_IT.store(it, Ordering::Relaxed); LP_TAG.store(0, Ordering::Relaxed);
        let Some((pred, w, mut cons, m0, zd)) = try_lp(p, &x, s, &mut choice, r, &hq) else { r /= 4.0; if r < o.rmin { break; } continue; };
        let small = pred < 1e-14 * s.max(1.0) || r < o.rmin;
        if small {
            // branch search: flip one ambiguous pair at a time, keep the best predicted decrease
            let amb = snap(p, &x, &mut choice, o.cc_tol);
            let rt = r.max(1e-6);
            // only pairs whose rows carry multiplier can gain from a flip: rank by total dual, try the top few
            let mut wt = vec![0.0; p.pairs.len()];
            for (c, &zz) in cons.iter().zip(zd.iter()) {
                if let Spec::PI(pk, _) | Spec::PJ(pk, _) = c.sp { wt[pk] += zz; }
            }
            let zmax = wt.iter().cloned().fold(0.0, f64::max);
            let mut amb: Vec<(usize, Vec<usize>)> = amb.into_iter().filter(|(pk, _)| wt[*pk] > 1e-6 * zmax).collect();
            amb.sort_by(|a, b| wt[b.0].partial_cmp(&wt[a.0]).unwrap());
            amb.truncate(o.flip_top);
            let mut best: Option<(f64, usize, usize)> = None;
            for (pk, cand) in &amb {
                if best.is_some() { break; }
                let (pk, cur) = (*pk, choice[*pk]);
                for &c in cand {
                    if Some(c) == cur { continue; }
                    let mut ch = choice.clone();
                    ch[pk] = Some(c);
                    LP_TAG.store(2, Ordering::Relaxed);
                    if let Some((pr2, ..)) = try_lp(p, &x, s, &mut ch, rt, &hq) {
                        if pr2 > 1e-12 * s && best.map_or(true, |b| pr2 > b.0) { best = Some((pr2, pk, c)); }
                    }
                    if best.is_some() { break; }
                }
            }
            if o.verbose { eprintln!("  stall at it {it}: R {r:.1e} pred {pred:.1e}, {} ambiguous pairs, best flip {:?}", amb.len(), best.map(|b| b.0)); }
            match best {
                Some((_, pk, c)) if stalls < 50 => { choice[pk] = Some(c); flips += 1; stalls += 1; r = rt; continue; }
                _ => break,
            }
        }
        let v0 = p.pack(&x, s);
        let nv = v0.len();
        let mut v1: Vec<f64> = (0..nv).map(|k| v0[k] + r * w[k]).collect();
        p.tie(&v0, &mut v1);
        let mut m1 = merit(p, &v1);
        let mut rho = (m0 - m1) / pred;
        let wmax = w.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
        if rho < 0.75 {
            // second-order correction: shift each row by its curvature error at the trial point
            for c in cons.iter_mut() {
                let lin: f64 = c.nz.iter().map(|&(k, a)| a * r * w[k]).sum();
                c.g = p.g_at(&v1, c.sp) - lin;
            }
            LP_TAG.store(1, Ordering::Relaxed);
            if let Some((w2, _, _)) = slp_lp(nv, 3 * n, &cons, r, kappa, &hq, p) {
                let mut v2: Vec<f64> = (0..nv).map(|k| v0[k] + r * w2[k]).collect();
                p.tie(&v0, &mut v2);
                let m2 = merit(p, &v2);
                if m2 < m1 { m1 = m2; v1 = v2; rho = (m0 - m1) / pred; }
            }
        }
        if o.verbose && (it < 20 || it % 10 == 0) { eprintln!("  slp {it}: R {r:.1e} s {:.13} pred {pred:.2e} rho {rho:.2}", s); }
        if rho >= 0.1 {
            x.copy_from_slice(&v1[..3 * n]);
            s = v1[3 * n];
            if o.sqp { hq = hess_lagr(p, &v1, &cons, &zd); }
            let (dm, sl, na, nd) = if o.ident {
                let (dm, sl, _, _) = ident(p, &v0, &v1, &mut Vec::new(), 1e-8);
                let lk = load_keys(p, &v0, &cons, &zd);
                let a: std::collections::HashSet<_> = lk.iter().collect();
                let b: std::collections::HashSet<_> = prev_act.iter().collect();
                let nd = a.symmetric_difference(&b).count();
                let na = lk.len();
                IDENT_FPS.lock().unwrap().push(fingerprint(&lk));
                prev_act = lk;
                (dm, sl, na, nd)
            } else { (0.0, 0.0, 0, 0) };
            trace.push((it, tp.elapsed().as_secs_f64(), m1, flips, dm, sl, na, nd));
            if rho > 0.5 && wmax > 0.5 { r = (2.0 * r).min(o.rmax); }
        } else {
            r /= 4.0;
        }
    }
    // leave separators consistent with x for the caller
    snap(p, &x, &mut choice, o.cc_tol);
    // final force network: one LP at small R at the final point
    LP_TAG.store(3, Ordering::Relaxed);
    if let Some((_, _, cons, _, z)) = try_lp(p, &x, s, &mut choice, r.max(1e-9).min(1e-6), &[]) {
        let v = p.pack(&x, s);
        *p_load = load_keys(p, &v, &cons, &z);
    }
    (s, x, it, r, flips)
}

// ---------------------------------------------------------------- L-BFGS
fn lbfgs(f: &mut dyn FnMut(&[f64], &mut [f64]) -> f64, v: &mut [f64], maxit: usize, gtol: f64, maxstep: f64) -> (f64, usize, f64) {
    let dim = v.len();
    let m = 12;
    let mut g = vec![0.0; dim];
    let mut e = f(v, &mut g);
    let mut nev = 1;
    let mut ss: Vec<Vec<f64>> = Vec::new();
    let mut ys: Vec<Vec<f64>> = Vec::new();
    let mut rho: Vec<f64> = Vec::new();
    let mut d = vec![0.0; dim];
    let mut xn = vec![0.0; dim];
    let mut gn = vec![0.0; dim];
    let mut alpha = vec![0.0; m];
    let mut stall = 0;
    let ginf = |g: &[f64]| g.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
    for _ in 0..maxit {
        if ginf(&g) < gtol { break; }
        for k in 0..dim { d[k] = -g[k]; }
        let h = ss.len();
        for k in (0..h).rev() {
            let a = rho[k] * dot(&ss[k], &d);
            alpha[k] = a;
            for q in 0..dim { d[q] -= a * ys[k][q]; }
        }
        if h > 0 {
            let gam = dot(&ss[h - 1], &ys[h - 1]) / dot(&ys[h - 1], &ys[h - 1]);
            for q in 0..dim { d[q] *= gam; }
        } else {
            let gm = ginf(&g);
            let sc = if gm > 0.0 { (1e-3 / gm).min(1.0) } else { 1.0 };
            for q in 0..dim { d[q] *= sc; }
        }
        for k in 0..h {
            let b = rho[k] * dot(&ys[k], &d);
            for q in 0..dim { d[q] += ss[k][q] * (alpha[k] - b); }
        }
        let mut gd = dot(&g, &d);
        if !(gd < 0.0) {
            ss.clear(); ys.clear(); rho.clear();
            let gm = ginf(&g);
            for q in 0..dim { d[q] = -g[q] * (1e-3 / gm).min(1.0); }
            gd = dot(&g, &d);
        }
        let dmax = ginf(&d);
        let mut step = if dmax > maxstep { maxstep / dmax } else { 1.0 };
        let mut ok = false;
        for _ in 0..50 {
            for q in 0..dim { xn[q] = v[q] + step * d[q]; }
            let en = f(&xn, &mut gn);
            nev += 1;
            if en <= e + 1e-4 * step * gd {
                let sv: Vec<f64> = (0..dim).map(|q| xn[q] - v[q]).collect();
                let yv: Vec<f64> = (0..dim).map(|q| gn[q] - g[q]).collect();
                let sy = dot(&sv, &yv);
                if sy > 1e-300 {
                    if ss.len() == m { ss.remove(0); ys.remove(0); rho.remove(0); }
                    rho.push(1.0 / sy); ss.push(sv); ys.push(yv);
                }
                if e - en <= 1e-15 * e.abs().max(1.0) { stall += 1; } else { stall = 0; }
                v.copy_from_slice(&xn);
                g.copy_from_slice(&gn);
                e = en;
                ok = true;
                break;
            }
            step *= 0.5;
        }
        if !ok {
            if ss.is_empty() { break; }
            ss.clear(); ys.clear(); rho.clear();
            stall += 1;
        }
        if stall > 10 { break; }
    }
    (e, nev, ginf(&g))
}

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x, y)| x * y).sum() }

// ---------------------------------------------------------------- quench
struct QOpt { mu0: f64, mu_max: f64, tol: f64, skin: f64, max_outer: usize, verbose: bool }

/// Returns (s, x, total evaluations, outer iterations, final violation before repair).
/// Minimise the pure overlap penalty at fixed side s (separators free, multipliers 0).  Returns (x, max violation, evals).
fn relax_fixed(p: &mut Prob, s: f64, x0: &[f64], skin: f64) -> (Vec<f64>, f64, usize) {
    let n = p.n;
    let cut = SQRT2 + skin;
    let mut x = x0.to_vec();
    p.fix_s = true;
    p.mu = 1.0;
    for l in p.wlam.iter_mut() { *l = 0.0; }
    p.rebuild(&x, cut);
    for pr in p.pairs.iter_mut() { pr.lam = [0.0; 8]; }
    let mut nev = 0;
    let mut viol = f64::INFINITY;
    for _round in 0..20 {
        let xref = x.clone();
        let mut v = p.pack(&x, s);
        let (_, ne, _) = { let pr = &*p; lbfgs(&mut |vv: &[f64], gg: &mut [f64]| pr.eval(vv, gg), &mut v, 4000, 1e-13, 0.05) };
        nev += ne;
        p.unpack(&v);
        x.copy_from_slice(&v[..3 * n]);
        viol = p.viol_update(&v, false);
        let mut dmax: f64 = 0.0;
        for i in 0..n { dmax = dmax.max((x[3 * i] - xref[3 * i]).hypot(x[3 * i + 1] - xref[3 * i + 1])); }
        if dmax > 0.5 * skin { p.rebuild(&x, cut); continue; }
        break;
    }
    p.fix_s = false;
    (x, viol, nev)
}

/// Shrink-hopping descent: grow until the fixed-side overlap problem is solved, then shrink by d (positions scaled about
/// the centre) while it stays solvable, halving d on failure (back to the last feasible state) down to dmin.
/// Returns (s, x) of the last feasible state, total evals, shrink steps.
fn shrink(s0: f64, x0: &[f64], d0: f64, dmin: f64, vtol: f64, verbose: bool) -> (f64, Vec<f64>, usize, usize) {
    let n = x0.len() / 3;
    let mut p = Prob { n, pairs: Vec::new(), wlam: vec![0.0; 16 * n], mu: 1.0, fix_s: true };
    let scale = |x: &[f64], s: f64, f: f64| -> Vec<f64> {
        let mut y = x.to_vec();
        for i in 0..n { y[3 * i] = s * f / 2.0 + (x[3 * i] - s / 2.0) * f; y[3 * i + 1] = s * f / 2.0 + (x[3 * i + 1] - s / 2.0) * f; }
        y
    };
    let mut s = s0;
    let mut x = x0.to_vec();
    let mut nev = 0;
    // grow to feasibility
    let mut feas: Option<(f64, Vec<f64>)> = None;
    for _ in 0..40 {
        let (xr, vi, ne) = relax_fixed(&mut p, s, &x, 0.3);
        nev += ne;
        if vi < vtol { feas = Some((s, xr)); break; }
        x = scale(&xr, s, 1.0 + d0);
        s *= 1.0 + d0;
    }
    let Some((mut sf, mut xf)) = feas else { return (f64::INFINITY, x, nev, 0); };
    let mut d = d0;
    let mut steps = 0;
    while d >= dmin && steps < 400 {
        steps += 1;
        let st = sf * (1.0 - d);
        let xt = scale(&xf, sf, 1.0 - d);
        let (xr, vi, ne) = relax_fixed(&mut p, st, &xt, 0.3);
        nev += ne;
        if vi < vtol { sf = st; xf = xr; if verbose { eprintln!("  shrink ok s {sf:.10} d {d:.1e}"); } }
        else { d /= 2.0; }
    }
    (sf, xf, nev, steps)
}

fn quench(s0: f64, x0: &[f64], o: &QOpt) -> (f64, Vec<f64>, usize, usize, f64, Prob) {
    let n = x0.len() / 3;
    let mut p = Prob { n, pairs: Vec::new(), wlam: vec![0.0; 16 * n], mu: o.mu0, fix_s: false };
    let cut = SQRT2 + o.skin;
    let mut x = x0.to_vec();
    let mut s = s0;
    p.rebuild(&x, cut);
    let mut xref = x.clone();
    let mut nev = 0;
    let mut prev_v = f64::INFINITY;
    let mut outer = 0;
    let mut viol = f64::INFINITY;
    while outer < o.max_outer {
        outer += 1;
        let mut v = p.pack(&x, s);
        let gtol = (1e-2 / p.mu.sqrt()).max(1e-11);
        let (_, ne, gi) = {
            let pr = &p;
            lbfgs(&mut |vv: &[f64], gg: &mut [f64]| pr.eval(vv, gg), &mut v, 20000, gtol, 0.05)
        };
        nev += ne;
        p.unpack(&v);
        x.copy_from_slice(&v[..3 * n]);
        s = v[3 * n];
        // pair list still valid?
        let mut dmax: f64 = 0.0;
        for i in 0..n {
            let (dx, dy) = (x[3 * i] - xref[3 * i], x[3 * i + 1] - xref[3 * i + 1]);
            dmax = dmax.max((dx * dx + dy * dy).sqrt());
        }
        if dmax > 0.5 * o.skin {
            p.rebuild(&x, cut);
            xref = x.clone();
            if o.verbose { eprintln!("  rebuild: {} pairs", p.pairs.len()); }
            continue; // re-solve with the new list before judging violation
        }
        viol = p.viol_update(&v, true);
        if o.verbose { eprintln!("outer {outer}: mu {:.1e} s {:.12} viol {:.2e} |g| {:.1e} evals {ne}", p.mu, s, viol, gi); }
        if viol < o.tol && p.mu >= 1e3 { break; }
        if viol > 0.25 * prev_v && p.mu < o.mu_max { p.mu *= 10.0; }
        prev_v = viol;
    }
    (s, x, nev, outer, viol, p)
}

/// Exact-ish (f64 SAT) feasibility: min pair gap and min wall clearance.
fn check(s: f64, x: &[f64]) -> (f64, f64) {
    let n = x.len() / 3;
    let mut mg = f64::INFINITY;
    for i in 0..n {
        for j in i + 1..n {
            let (dx, dy) = (x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1]);
            if dx * dx + dy * dy < 2.1 {
                mg = mg.min(sat_gap(x[3 * i], x[3 * i + 1], x[3 * i + 2], x[3 * j], x[3 * j + 1], x[3 * j + 2]));
            }
        }
    }
    let mut mw = f64::INFINITY;
    for i in 0..n {
        let t = x[3 * i + 2];
        let w = H * (t.cos().abs() + t.sin().abs());
        mw = mw.min(x[3 * i] - w).min(s - x[3 * i] - w).min(x[3 * i + 1] - w).min(s - x[3 * i + 1] - w);
    }
    (mg, mw)
}

/// Scale centres and side by f >= 1 until the f64 SAT check passes with margin.
fn repair(s: f64, x: &mut [f64], margin: f64) -> f64 {
    let n = x.len() / 3;
    // translate so the bounding box of the squares starts at the origin (walls x >= 0, y >= 0 tight or slack)
    let mut s = s;
    for _ in 0..60 {
        let (mg, mw) = check(s, x);
        let bad = (margin - mg).max(margin - mw);
        if bad <= 0.0 { return s; }
        let f = 1.0 + 2.0 * bad + 1e-15;
        for i in 0..n {
            x[3 * i] = (x[3 * i] - s / 2.0) * f + s * f / 2.0;
            x[3 * i + 1] = (x[3 * i + 1] - s / 2.0) * f + s * f / 2.0;
        }
        s *= f;
    }
    s
}

// ---------------------------------------------------------------- io
fn read_cfg(path: &str) -> (f64, Vec<f64>) {
    let txt = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let mut lines = txt.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let head: Vec<f64> = lines.next().unwrap().split_whitespace().map(|t| t.parse().unwrap()).collect();
    let (n, s) = (head[0] as usize, head[1]);
    let mut v = Vec::with_capacity(3 * n);
    for _ in 0..n {
        let t: Vec<f64> = lines.next().unwrap().split_whitespace().take(3).map(|t| t.parse().unwrap()).collect();
        v.push(t[0]); v.push(t[1]); v.push(t[2].to_radians());
    }
    (s, v)
}

fn write_cfg(path: &str, s: f64, v: &[f64]) {
    let n = v.len() / 3;
    let mut out = format!("{} {:.17}\n", n, s);
    for i in 0..n {
        out += &format!("{:.17} {:.17} {:.15}\n", v[3 * i], v[3 * i + 1], v[3 * i + 2].to_degrees().rem_euclid(90.0));
    }
    let tmp = format!("{path}.tmp");
    fs::write(&tmp, out).unwrap();
    fs::rename(&tmp, path).unwrap();
}

fn arg<T: std::str::FromStr>(a: &[String], k: &str, def: T) -> T {
    a.iter().rposition(|x| x == k).and_then(|i| a.get(i + 1)).and_then(|v| v.parse().ok()).unwrap_or(def)
}
fn sarg(a: &[String], k: &str) -> Option<String> { a.iter().rposition(|x| x == k).and_then(|i| a.get(i + 1)).cloned() }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let cmd = a.get(1).map(|s| s.as_str()).unwrap_or("");
    let inp = sarg(&a, "--in").expect("--in");
    if let Some(d) = sarg(&a, "--dump-lp") { fs::create_dir_all(&d).expect("--dump-lp dir"); DUMP_LP.set(d).ok(); }
    let (s, mut x) = read_cfg(&inp);
    match cmd {
        "check" => {
            let (mg, mw) = check(s, &x);
            println!("{} min_gap {:.3e} min_wall {:.3e}", s, mg, mw);
        }
        "quench" => {
            let t0 = Instant::now();
            let n = x.len() / 3;
            let sig: f64 = arg(&a, "--kick", 0.0);
            if sig > 0.0 {
                let mut r = Rng::new(arg(&a, "--seed", 1u64));
                for i in 0..n {
                    x[3 * i] += sig * r.gauss();
                    x[3 * i + 1] += sig * r.gauss();
                    x[3 * i + 2] += (20.0 * sig * r.gauss()).to_radians();
                }
            }
            let lo: f64 = arg(&a, "--loosen", 1.02);
            let s1 = s * lo;
            for i in 0..n {
                x[3 * i] = (x[3 * i] * lo).clamp(H, s1 - H);
                x[3 * i + 1] = (x[3 * i + 1] * lo).clamp(H, s1 - H);
            }
            let o = QOpt {
                mu0: arg(&a, "--mu0", 10.0), mu_max: arg(&a, "--mu-max", 1e2), tol: arg(&a, "--tol", 1e-7),
                skin: arg(&a, "--skin", 0.3), max_outer: arg(&a, "--outer", 60), verbose: a.iter().any(|t| t == "--v"),
            };
            let (mut sq, mut xq, nev, outer, viol, mut prob) = if a.iter().any(|t| t == "--no-alm") {
                // already (nearly) feasible: polish directly
                (s1, x.clone(), 0, 0, 0.0, Prob { n, pairs: Vec::new(), wlam: vec![0.0; 16 * n], mu: 1.0, fix_s: false })
            } else if a.iter().any(|t| t == "--shrink") {
                let (ss, xs, ne, st) = shrink(s1, &x, arg(&a, "--d0", 2e-3), arg(&a, "--dmin", 1e-6), arg(&a, "--vtol", 1e-10), o.verbose);
                if o.verbose { eprintln!("shrink: s {ss:.10} evals {ne} steps {st}"); }
                (ss, xs, ne, st, 0.0, Prob { n, pairs: Vec::new(), wlam: vec![0.0; 16 * n], mu: 1.0, fix_s: false })
            } else { quench(s1, &x, &o) };
            let t_alm = t0.elapsed().as_secs_f64();
            let mut xc = xq.clone();
            let s_coarse = repair(sq, &mut xc, 1e-13);
            let s_alm0 = sq;
            let mut pit = 0;
            let mut pr = 0.0;
            let mut flips = 0;
            let mut trace: Vec<(usize, f64, f64, usize, f64, f64, usize, usize)> = Vec::new();
            let mut load: Vec<(u32, u32, u8, u8)> = Vec::new();
            if !a.iter().any(|t| t == "--no-polish") {
                let po = POpt { sqp: a.iter().any(|t| t == "--sqp"), rmax: arg(&a, "--rmax", 4e-3), r0: arg(&a, "--r0", 1e-3), rmin: arg(&a, "--rmin", 1e-10), maxit: arg(&a, "--pit", 300), cc_tol: arg(&a, "--cc-tol", 1e-7), flip_top: arg(&a, "--flip-top", 8), ident: a.iter().any(|t| t == "--ident"), stag_w: arg(&a, "--stag-w", 15), stag_tol: arg(&a, "--stag-tol", 2e-9), verbose: o.verbose };
                let (s2, x2, it, rr, fl) = polish(&mut prob, sq, &xq, &po, &mut trace, &mut load);
                sq = s2; xq = x2; pit = it; pr = rr; flips = fl;
            }
            let (fp, ncont) = (fingerprint(&load), load.len());
            let sr = repair(sq, &mut xq, 1e-13);
            let (mg, mw) = check(sr, &xq);
            if let Some(out) = sarg(&a, "--out") { write_cfg(&out, sr, &xq); }
            println!("{{\"s\": {:.15}, \"s_alm\": {:.15}, \"viol\": {:.2e}, \"evals\": {}, \"outer\": {}, \"slp_it\": {}, \"slp_r\": {:.1e}, \"flips\": {}, \"min_gap\": {:.2e}, \"min_wall\": {:.2e}, \"sec\": {:.3}, \"s_coarse\": {:.15}, \"t_alm\": {:.3}, \"s_alm0\": {:.15}, \"trace\": [{}], \"fp\": \"{:016x}\", \"ncontacts\": {}{}}}",
                     sr, sq, viol, nev, outer, pit, pr, flips, mg, mw, t0.elapsed().as_secs_f64(), s_coarse, t_alm, s_alm0, trace.iter().map(|q| format!("[{},{:.4},{:.15},{},{:.3e},{:.3e},{},{}]", q.0, q.1, q.2, q.3, q.4, q.5, q.6, q.7)).collect::<Vec<_>>().join(","), fp, ncont, if a.iter().any(|t| t == "--ident") { format!(", \"fps\": [{}]", IDENT_FPS.lock().unwrap().iter().map(|h| format!("\"{:016x}\"", h)).collect::<Vec<_>>().join(",")) } else { String::new() });
        }
        "gradcheck" => {
            // random kick so many constraints are active, then compare analytic vs central differences
            let n = x.len() / 3;
            let mut r = Rng::new(7);
            for q in x.iter_mut() { *q += 0.05 * r.gauss(); }
            let mut p = Prob { n, pairs: Vec::new(), wlam: (0..16 * n).map(|_| r.f()).collect(), mu: 37.0, fix_s: false };
            p.rebuild(&x, SQRT2 + 0.4);
            for pr in p.pairs.iter_mut() { for l in pr.lam.iter_mut() { *l = r.f(); } pr.phi += 0.1 * r.gauss(); pr.d += 0.05 * r.gauss(); }
            let v = p.pack(&x, s * 0.97);
            let mut g = vec![0.0; v.len()];
            p.eval(&v, &mut g);
            let mut worst: f64 = 0.0;
            let mut gw = vec![0.0; v.len()];
            for k in 0..v.len() {
                let h = 1e-6;
                let mut a1 = v.clone(); a1[k] += h;
                let mut a2 = v.clone(); a2[k] -= h;
                let fd = (p.eval(&a1, &mut gw) - p.eval(&a2, &mut gw)) / (2.0 * h);
                worst = worst.max((fd - g[k]).abs() / (1.0 + g[k].abs()));
            }
            println!("dim {} pairs {} worst rel err {:.2e}", v.len(), p.pairs.len(), worst);
        }
        _ => { eprintln!("usage: fq quench|check --in FILE ..."); std::process::exit(2); }
    }
}
