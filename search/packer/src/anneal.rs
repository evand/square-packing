//! anneal: hard rounded squares in a square box, NPT compression MC with the shape morphing from disk to square
//! (10-08, Evan: start as a hot disk packing, anneal to squares as it cools).
//!
//! Shape: inner square of half-side h = 0.5 - r, Minkowski-rounded by radius r (r = 0.5: unit-diameter disk; r = 0: unit
//! square).  Overlap: Euclidean distance between the inner squares < 2r (r = 0: SAT).  Walls: inner corners >= r from
//! each wall of [0, s]^2.
//! MC: single-particle translate / rotate, isotropic box rescale (NPT, ln-area moves: accept exp(-bP dA + (N+1) ln(A'/A))),
//! shape step toward the scheduled r (accepted if overlap-free; if not, the smallest box growth <= 5 % that makes it
//! feasible: an annealing heuristic, not detailed balance).  Adaptive steps toward ~40 % acceptance.
//! Schedule over progress t in [0, 1] (sweeps): bP = bp0 (bp1/bp0)^t; r = 0.5 (1 - clamp((t - t0)/(t1 - t0)))^p.
//!
//!   anneal run --n 110 --seed 1 --sweeps 40000 --bp0 1 --bp1 3000 --t0 0.2 --t1 0.8 --p 1 --out X.txt
//!   anneal test        (distance routine vs brute force)
//! Output: packer format ("n s", then "x y theta_deg"), centres in [0, s]^2; a JSON summary on stdout.

use std::f64::consts::{FRAC_PI_2, PI, SQRT_2};

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        let mut r = Rng(seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(0x632BE59BD9B4E019) | 1);
        for _ in 0..8 { r.next(); }
        r
    }
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn f(&mut self) -> f64 { (self.next() >> 11) as f64 / (1u64 << 53) as f64 }
    fn u(&mut self) -> f64 { 2.0 * self.f() - 1.0 }
    fn idx(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
}

fn verts(x: f64, y: f64, t: f64, h: f64) -> [(f64, f64); 4] {
    let (c, s) = (t.cos(), t.sin());
    let mut v = [(0.0, 0.0); 4];
    for (k, &(a, b)) in [(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)].iter().enumerate() {
        v[k] = (x + h * (a * c - b * s), y + h * (a * s + b * c));
    }
    v
}

/// SAT gap of two squares of half-side h (max over the 4 face normals of the projected separation; < 0 = overlap).
fn sat_gap(xi: f64, yi: f64, ti: f64, xj: f64, yj: f64, tj: f64, h: f64) -> f64 {
    let mut best = f64::NEG_INFINITY;
    let (dx, dy) = (xj - xi, yj - yi);
    for &t in &[ti, ti + FRAC_PI_2, tj, tj + FRAC_PI_2] {
        let (ax, ay) = (t.cos(), t.sin());
        let dc = (dx * ax + dy * ay).abs();
        let ri = h * ((ti - t).cos().abs() + (ti - t).sin().abs());
        let rj = h * ((tj - t).cos().abs() + (tj - t).sin().abs());
        best = best.max(dc - ri - rj);
    }
    best
}

fn pt_seg(px: f64, py: f64, a: (f64, f64), b: (f64, f64)) -> f64 {
    let (ex, ey) = (b.0 - a.0, b.1 - a.1);
    let l2 = ex * ex + ey * ey;
    let t = if l2 > 0.0 { (((px - a.0) * ex + (py - a.1) * ey) / l2).clamp(0.0, 1.0) } else { 0.0 };
    (px - a.0 - t * ex).hypot(py - a.1 - t * ey)
}

/// Euclidean distance between two closed squares of half-side h (0 if they intersect).
fn sq_dist(xi: f64, yi: f64, ti: f64, xj: f64, yj: f64, tj: f64, h: f64) -> f64 {
    if sat_gap(xi, yi, ti, xj, yj, tj, h) <= 0.0 { return 0.0; }
    let (a, b) = (verts(xi, yi, ti, h), verts(xj, yj, tj, h));
    let mut d = f64::INFINITY;
    for k in 0..4 {
        for e in 0..4 {
            d = d.min(pt_seg(a[k].0, a[k].1, b[e], b[(e + 1) % 4]));
            d = d.min(pt_seg(b[k].0, b[k].1, a[e], a[(e + 1) % 4]));
        }
    }
    d
}

struct Sys { n: usize, x: Vec<f64>, s: f64, r: f64 }

impl Sys {
    /// does particle i at (x, y, t) overlap particle j (shape parameter r)?
    fn pair_ov(&self, x: f64, y: f64, t: f64, j: usize, r: f64) -> bool {
        let (xj, yj, tj) = (self.x[3 * j], self.x[3 * j + 1], self.x[3 * j + 2]);
        let h = 0.5 - r;
        let d2 = (x - xj).powi(2) + (y - yj).powi(2);
        if d2 < 1.0 - 1e-12 { return true; }                     // inscribed disks (radius 0.5) overlap
        let rc = h * SQRT_2 + r;
        if d2 >= 4.0 * rc * rc { return false; }                  // circumscribed disks apart
        if r <= 0.0 { return sat_gap(x, y, t, xj, yj, tj, h) < 0.0; }
        let g = sat_gap(x, y, t, xj, yj, tj, h);
        if g >= 2.0 * r { return false; }
        if g <= 0.0 { return true; }
        sq_dist(x, y, t, xj, yj, tj, h) < 2.0 * r
    }
    fn wall_ok(&self, x: f64, y: f64, t: f64, s: f64, r: f64) -> bool {
        let w = (0.5 - r) * (t.cos().abs() + t.sin().abs()) + r;
        x >= w && x <= s - w && y >= w && y <= s - w
    }
    fn one_ok(&self, i: usize, x: f64, y: f64, t: f64) -> bool {
        if !self.wall_ok(x, y, t, self.s, self.r) { return false; }
        (0..self.n).all(|j| j == i || !self.pair_ov(x, y, t, j, self.r))
    }
    /// full check of a configuration xs, side s, shape r
    fn all_ok(&self, xs: &[f64], s: f64, r: f64) -> bool {
        let tmp = Sys { n: self.n, x: xs.to_vec(), s, r };
        for i in 0..self.n {
            let (x, y, t) = (xs[3 * i], xs[3 * i + 1], xs[3 * i + 2]);
            if !tmp.wall_ok(x, y, t, s, r) { return false; }
            for j in i + 1..self.n { if tmp.pair_ov(x, y, t, j, r) { return false; } }
        }
        true
    }
}

fn area(r: f64) -> f64 { let h = 0.5 - r; 4.0 * h * h + 8.0 * h * r + PI * r * r }

fn arg<T: std::str::FromStr>(a: &[String], k: &str, def: T) -> T {
    a.iter().rposition(|x| x == k).and_then(|i| a.get(i + 1)).and_then(|v| v.parse().ok()).unwrap_or(def)
}

fn test() {
    let mut rng = Rng::new(3);
    let mut worst: f64 = 0.0;
    for trial in 0..3000 {
        let h = [0.5, 0.3, 0.1, 0.0][trial % 4];
        let (xi, yi, ti) = (0.0, 0.0, rng.f() * PI);
        let (xj, yj, tj) = (rng.u() * 2.0, rng.u() * 2.0, rng.f() * PI);
        let d = sq_dist(xi, yi, ti, xj, yj, tj, h);
        // brute force: dense boundary sampling of both squares (or the point) + inside test
        let samp = |x: f64, y: f64, t: f64| -> Vec<(f64, f64)> {
            let v = verts(x, y, t, h);
            let mut p = Vec::new();
            for e in 0..4 { for q in 0..400 { let s = q as f64 / 400.0; p.push((v[e].0 + s * (v[(e + 1) % 4].0 - v[e].0), v[e].1 + s * (v[(e + 1) % 4].1 - v[e].1))); } }
            p
        };
        let (pa, pb) = (samp(xi, yi, ti), samp(xj, yj, tj));
        let mut bd = f64::INFINITY;
        for a in &pa { for b in &pb { bd = bd.min((a.0 - b.0).hypot(a.1 - b.1)); } }
        let inside = sat_gap(xi, yi, ti, xj, yj, tj, h) <= 0.0;
        let err = if inside { d } else { (d - bd).abs() };
        // sampling resolution ~ 2h/400 * sqrt2
        if !inside && bd > 1e-9 && err > 2.0 * h / 400.0 + 1e-12 { println!("MISMATCH h {h} d {d} brute {bd}"); }
        if inside && d != 0.0 { println!("MISMATCH overlap h {h} d {d}"); }
        worst = worst.max(if inside { 0.0 } else { err });
    }
    println!("test done: worst |exact - brute| {worst:.2e} (brute resolution ~ {:.1e})", 2.0 * 0.5 / 400.0);
}

fn run(a: &[String]) {
    let n: usize = arg(a, "--n", 110);
    let seed: u64 = arg(a, "--seed", 1);
    let sweeps: usize = arg(a, "--sweeps", 40000);
    let (bp0, bp1): (f64, f64) = (arg(a, "--bp0", 1.0), arg(a, "--bp1", 3000.0));
    let (t0, t1, pw): (f64, f64, f64) = (arg(a, "--t0", 0.2), arg(a, "--t1", 0.8), arg(a, "--p", 1.0));
    let r_init: f64 = arg(a, "--r0", 0.5);
    let phi0: f64 = arg(a, "--phi0", 0.25);
    let out: String = arg(a, "--out", String::from("anneal.txt"));
    let mut rng = Rng::new(seed);
    let rsched = |t: f64| -> f64 {
        if t1 <= t0 { return if t >= t0 { 0.0 } else { r_init }; }
        r_init * (1.0 - ((t - t0) / (t1 - t0)).clamp(0.0, 1.0)).powf(pw)
    };
    let r_start = rsched(0.0);
    let s0 = (n as f64 * area(r_start) / phi0).sqrt();
    let mut sys = Sys { n, x: Vec::with_capacity(3 * n), s: s0, r: r_start };
    // random sequential placement
    let mut placed = 0;
    let mut tries = 0usize;
    while placed < n {
        tries += 1;
        assert!(tries < 10_000_000, "placement failed");
        let (x, y, t) = (rng.f() * s0, rng.f() * s0, rng.f() * FRAC_PI_2);
        sys.n = placed;
        if sys.one_ok(usize::MAX, x, y, t) { sys.x.extend_from_slice(&[x, y, t]); placed += 1; }
    }
    sys.n = n;
    let (mut dt, mut dr, mut dv) = (0.3f64, 0.3f64, 0.01f64);
    let (mut at, mut nt, mut av, mut nv) = (0usize, 0usize, 0usize, 0usize);
    let mut shape_fail = 0usize;
    let t_start = std::time::Instant::now();
    let report = sweeps / 10;
    for sw in 0..sweeps {
        let t = sw as f64 / sweeps as f64;
        let bp = bp0 * (bp1 / bp0).powf(t);
        // shape step toward the schedule
        let rt = rsched(t);
        if rt < sys.r - 1e-15 {
            let xs = sys.x.clone();
            if sys.all_ok(&xs, sys.s, rt) { sys.r = rt; }
            else {
                // smallest box growth (<= 5 %) that makes the new shape feasible (bisection on the scale factor)
                let feas = |f: f64, sys: &Sys| { let xf: Vec<f64> = xs.iter().enumerate().map(|(k, &v)| if k % 3 == 2 { v } else { v * f }).collect(); (sys.all_ok(&xf, sys.s * f, rt), xf) };
                let (ok, _) = feas(1.05, &sys);
                if ok {
                    let (mut lo, mut hi) = (1.0, 1.05);
                    for _ in 0..20 { let m = 0.5 * (lo + hi); if feas(m, &sys).0 { hi = m; } else { lo = m; } }
                    let (_, xf) = feas(hi, &sys);
                    sys.x = xf; sys.s *= hi; sys.r = rt;
                } else { shape_fail += 1; }
            }
        }
        for _ in 0..n {
            let i = rng.idx(n);
            let (x, y, th) = (sys.x[3 * i], sys.x[3 * i + 1], sys.x[3 * i + 2]);
            let (nx, ny, nth) = if rng.f() < 0.5 || sys.r >= 0.5 { (x + dt * rng.u(), y + dt * rng.u(), th) } else { (x, y, th + dr * rng.u()) };
            nt += 1;
            if sys.one_ok(i, nx, ny, nth) { sys.x[3 * i] = nx; sys.x[3 * i + 1] = ny; sys.x[3 * i + 2] = nth; at += 1; }
        }
        // box move (ln-area)
        {
            let dl = dv * rng.u();
            let f = (0.5 * dl).exp();
            let s1 = sys.s * f;
            let da = s1 * s1 - sys.s * sys.s;
            let lacc = -bp * da + (n as f64 + 1.0) * dl;
            nv += 1;
            if lacc >= 0.0 || rng.f() < lacc.exp() {
                let xf: Vec<f64> = sys.x.iter().enumerate().map(|(k, &v)| if k % 3 == 2 { v } else { v * f }).collect();
                if f >= 1.0 || sys.all_ok(&xf, s1, sys.r) { sys.x = xf; sys.s = s1; av += 1; }
            }
        }
        if sw % 50 == 49 {
            let ra = at as f64 / nt.max(1) as f64;
            let fa = if ra > 0.45 { 1.1 } else if ra < 0.35 { 0.9 } else { 1.0 };
            dt = (dt * fa).clamp(1e-7, 0.5); dr = (dr * fa).clamp(1e-7, 0.5);
            let rv = av as f64 / nv.max(1) as f64;
            dv = (dv * if rv > 0.45 { 1.1 } else if rv < 0.35 { 0.9 } else { 1.0 }).clamp(1e-9, 0.1);
            at = 0; nt = 0; av = 0; nv = 0;
        }
        if std::env::var("AN_V").is_ok() && report > 0 && sw % report == 0 {
            eprintln!("sw {sw} t {t:.2} bP {bp:.1} r {:.4} s {:.5} phi {:.4} dt {dt:.1e} dv {dv:.1e} shape_fail {shape_fail}", sys.r, sys.s, n as f64 * area(sys.r) / (sys.s * sys.s));
        }
    }
    // if the schedule did not reach r = 0, force it with box growth (report how much)
    let r_end = sys.r;
    let mut grow = 1.0;
    if sys.r > 0.0 {
        let xs = sys.x.clone();
        let mut f = 1.0;
        loop {
            let xf: Vec<f64> = xs.iter().enumerate().map(|(k, &v)| if k % 3 == 2 { v } else { v * f }).collect();
            if sys.all_ok(&xf, sys.s * f, 0.0) { sys.x = xf; sys.s *= f; sys.r = 0.0; grow = f; break; }
            f *= 1.01;
            if f > 2.0 { break; }
        }
    }
    let mut txt = format!("{} {:?}\n", n, sys.s);
    for i in 0..n { txt += &format!("{:?} {:?} {:?}\n", sys.x[3 * i], sys.x[3 * i + 1], (sys.x[3 * i + 2].to_degrees()).rem_euclid(90.0)); }
    std::fs::write(&out, txt).unwrap();
    println!("{{\"s\": {:.12}, \"r_end\": {:.4}, \"grow\": {:.4}, \"phi\": {:.5}, \"shape_fail\": {}, \"sec\": {:.2}}}",
             sys.s, r_end, grow, n as f64 / (sys.s * sys.s), shape_fail, t_start.elapsed().as_secs_f64());
}

// ---------------------------------------------------------------- regional melt (10-08)
fn verts_h(x: f64, y: f64, t: f64, h: f64) -> [(f64, f64); 4] { verts(x, y, t, h) }

/// SAT gap of two squares with half-sides hi, hj.
fn sat_gap2(xi: f64, yi: f64, ti: f64, hi: f64, xj: f64, yj: f64, tj: f64, hj: f64) -> f64 {
    let mut best = f64::NEG_INFINITY;
    let (dx, dy) = (xj - xi, yj - yi);
    for &t in &[ti, ti + FRAC_PI_2, tj, tj + FRAC_PI_2] {
        let (ax, ay) = (t.cos(), t.sin());
        let dc = (dx * ax + dy * ay).abs();
        let ri = hi * ((ti - t).cos().abs() + (ti - t).sin().abs());
        let rj = hj * ((tj - t).cos().abs() + (tj - t).sin().abs());
        best = best.max(dc - ri - rj);
    }
    best
}

fn sq_dist2(xi: f64, yi: f64, ti: f64, hi: f64, xj: f64, yj: f64, tj: f64, hj: f64) -> f64 {
    if sat_gap2(xi, yi, ti, hi, xj, yj, tj, hj) <= 0.0 { return 0.0; }
    let (a, b) = (verts_h(xi, yi, ti, hi), verts_h(xj, yj, tj, hj));
    let mut d = f64::INFINITY;
    for k in 0..4 { for e in 0..4 {
        d = d.min(pt_seg(a[k].0, a[k].1, b[e], b[(e + 1) % 4]));
        d = d.min(pt_seg(b[k].0, b[k].1, a[e], a[(e + 1) % 4]));
    } }
    d
}

/// Rounded squares with per-particle radius rr[i] (inner half-side 0.5 - rr[i]).
struct Sys2 { n: usize, x: Vec<f64>, s: f64, rr: Vec<f64> }
impl Sys2 {
    fn ov(&self, x: f64, y: f64, t: f64, ri: f64, j: usize) -> bool {
        let (xj, yj, tj, rj) = (self.x[3 * j], self.x[3 * j + 1], self.x[3 * j + 2], self.rr[j]);
        let d2 = (x - xj).powi(2) + (y - yj).powi(2);
        if d2 < 1.0 - 1e-12 { return true; }
        let (hi, hj) = (0.5 - ri, 0.5 - rj);
        let rc = hi * SQRT_2 + ri + hj * SQRT_2 + rj;
        if d2 >= rc * rc { return false; }
        let g = sat_gap2(x, y, t, hi, xj, yj, tj, hj);
        if ri + rj <= 0.0 { return g < 0.0; }
        if g >= ri + rj { return false; }
        if g <= 0.0 { return true; }
        sq_dist2(x, y, t, hi, xj, yj, tj, hj) < ri + rj
    }
    fn wall(&self, x: f64, y: f64, t: f64, ri: f64, s: f64) -> bool {
        let w = (0.5 - ri) * (t.cos().abs() + t.sin().abs()) + ri;
        x >= w && x <= s - w && y >= w && y <= s - w
    }
    fn ok(&self, i: usize, x: f64, y: f64, t: f64, ri: f64) -> bool {
        self.wall(x, y, t, ri, self.s) && (0..self.n).all(|j| j == i || !self.ov(x, y, t, ri, j))
    }
    fn all_ok(&self, xs: &[f64], s: f64) -> bool {
        let tmp = Sys2 { n: self.n, x: xs.to_vec(), s, rr: self.rr.clone() };
        (0..self.n).all(|i| { let (x, y, t) = (xs[3 * i], xs[3 * i + 1], xs[3 * i + 2]);
            tmp.wall(x, y, t, tmp.rr[i], s) && (i + 1..self.n).all(|j| !tmp.ov(x, y, t, tmp.rr[i], j)) })
    }
}

/// anneal melt --in F --out G --cx X --cy Y --rad R --rmax 0.1 --sweeps 4000 --bp 2000 --seed K
/// Squares within R of (cx, cy) get corner radius rmax (1 - d/R) (ramp up over t in [0, .25], hold to .55, down to 0 by
/// .85, then square); only squares within R + 1.5 move; NPT box moves at pressure bp.  Residual rounding at the end is
/// removed by the smallest uniform box growth.  Output packer format + JSON summary.
fn melt(a: &[String]) {
    let inp: String = arg(a, "--in", String::new());
    let out: String = arg(a, "--out", String::from("melt.txt"));
    let (cx, cy, rad): (f64, f64, f64) = (arg(a, "--cx", 0.0), arg(a, "--cy", 0.0), arg(a, "--rad", 2.5));
    let rmax: f64 = arg(a, "--rmax", 0.1);
    let sweeps: usize = arg(a, "--sweeps", 4000);
    let bp: f64 = arg(a, "--bp", 2000.0);
    let mut rng = Rng::new(arg(a, "--seed", 1u64));
    let txt = std::fs::read_to_string(&inp).expect("read --in");
    let mut it = txt.split_whitespace().map(|w| w.parse::<f64>().unwrap());
    let n = it.next().unwrap() as usize;
    let s0 = it.next().unwrap();
    let mut x = Vec::with_capacity(3 * n);
    for _ in 0..n { let (a1, b1, c1) = (it.next().unwrap(), it.next().unwrap(), it.next().unwrap()); x.extend_from_slice(&[a1, b1, c1.to_radians()]); }
    let w: Vec<f64> = (0..n).map(|i| (1.0 - (x[3 * i] - cx).hypot(x[3 * i + 1] - cy) / rad).max(0.0)).collect();
    let mobile: Vec<usize> = (0..n).filter(|&i| (x[3 * i] - cx).hypot(x[3 * i + 1] - cy) < rad + 1.5).collect();
    let mut sys = Sys2 { n, x, s: s0, rr: vec![0.0; n] };
    let prof = |t: f64| -> f64 { if t < 0.25 { t / 0.25 } else if t < 0.55 { 1.0 } else if t < 0.85 { 1.0 - (t - 0.55) / 0.3 } else { 0.0 } };
    let (mut dt, mut dr, mut dv) = (0.02f64, 0.02f64, 1e-3f64);
    let (mut at, mut nt, mut av, mut nv) = (0usize, 0usize, 0usize, 0usize);
    let t0 = std::time::Instant::now();
    for sw in 0..sweeps {
        let tt = sw as f64 / sweeps as f64;
        // shape: each particle toward its scheduled radius (growing the shape only where it fits)
        for &i in &mobile {
            let target = rmax * w[i] * prof(tt);
            if (target - sys.rr[i]).abs() < 1e-15 { continue; }
            let (xi, yi, ti) = (sys.x[3 * i], sys.x[3 * i + 1], sys.x[3 * i + 2]);
            if target > sys.rr[i] || sys.ok(i, xi, yi, ti, target) { sys.rr[i] = target; }
        }
        for _ in 0..mobile.len() {
            let i = mobile[rng.idx(mobile.len())];
            let (xi, yi, ti) = (sys.x[3 * i], sys.x[3 * i + 1], sys.x[3 * i + 2]);
            let (nx, ny, nth) = if rng.f() < 0.5 { (xi + dt * rng.u(), yi + dt * rng.u(), ti) } else { (xi, yi, ti + dr * rng.u()) };
            nt += 1;
            if sys.ok(i, nx, ny, nth, sys.rr[i]) { sys.x[3 * i] = nx; sys.x[3 * i + 1] = ny; sys.x[3 * i + 2] = nth; at += 1; }
        }
        if sw % 4 == 0 {
            let dl = dv * rng.u();
            let f = (0.5 * dl).exp();
            let s1 = sys.s * f;
            let lacc = -bp * (s1 * s1 - sys.s * sys.s) + (n as f64 + 1.0) * dl;
            nv += 1;
            if lacc >= 0.0 || rng.f() < lacc.exp() {
                let xf: Vec<f64> = sys.x.iter().enumerate().map(|(k, &v)| if k % 3 == 2 { v } else { v * f }).collect();
                if f >= 1.0 || sys.all_ok(&xf, s1) { sys.x = xf; sys.s = s1; av += 1; }
            }
        }
        if sw % 50 == 49 {
            let ra = at as f64 / nt.max(1) as f64;
            let fa = if ra > 0.45 { 1.1 } else if ra < 0.35 { 0.9 } else { 1.0 };
            dt = (dt * fa).clamp(1e-7, 0.3); dr = (dr * fa).clamp(1e-7, 0.3);
            let rv = av as f64 / nv.max(1) as f64;
            dv = (dv * if rv > 0.45 { 1.1 } else if rv < 0.35 { 0.9 } else { 1.0 }).clamp(1e-9, 0.05);
            at = 0; nt = 0; av = 0; nv = 0;
        }
    }
    let left = sys.rr.iter().filter(|&&r| r > 0.0).count();
    let mut grow = 1.0;
    if left > 0 {
        sys.rr = vec![0.0; n];
        let xs = sys.x.clone();
        let mut f = 1.0;
        loop {
            let xf: Vec<f64> = xs.iter().enumerate().map(|(k, &v)| if k % 3 == 2 { v } else { v * f }).collect();
            if sys.all_ok(&xf, sys.s * f) { sys.x = xf; sys.s *= f; grow = f; break; }
            f *= 1.001;
            if f > 1.5 { break; }
        }
    }
    let mut txt = format!("{} {:?}\n", n, sys.s);
    for i in 0..n { txt += &format!("{:?} {:?} {:?}\n", sys.x[3 * i], sys.x[3 * i + 1], (sys.x[3 * i + 2].to_degrees()).rem_euclid(90.0)); }
    std::fs::write(&out, txt).unwrap();
    println!("{{\"s\": {:.12}, \"s0\": {:.12}, \"mobile\": {}, \"rounded_left\": {}, \"grow\": {:.5}, \"sec\": {:.2}}}",
             sys.s, s0, mobile.len(), left, grow, t0.elapsed().as_secs_f64());
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(|s| s.as_str()) {
        Some("test") => test(),
        Some("run") => run(&a),
        Some("melt") => melt(&a),
        _ => eprintln!("usage: anneal test | anneal run --n N --seed K --sweeps S --bp0 --bp1 --t0 --t1 --p --out F"),
    }
}
