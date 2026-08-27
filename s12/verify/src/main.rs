//! EXACT verifier for weighted-unavoidable-set certificates in square packing.
//!
//! Certificate: atoms p_a = (X_a/D, Y_a/D) in C=[0,s]^2 with weights w_a = W_a/WD >= 0.
//! Claim to verify:  for EVERY unit square Q (any center, any angle) with Q subseteq C,
//!                   sum over atoms in the interior of Q of w_a  >=  1.
//!
//! Method (all arithmetic exact in i128):
//!   * angles theta_k = 2*arctan(k/N), k=0..K, covering [0,45deg]; cos/sin are RATIONAL:
//!         cos = (N^2-k^2)/(N^2+k^2),  sin = 2kN/(N^2+k^2).
//!   * for a gap delta_k = theta_{k+1}-theta_k, cos/sin of delta_k are rational too, and the
//!     unit square at ANY angle theta in [theta_k, theta_{k+1}] contains the concentric square
//!     of side sigma_k = 1/(cos delta_k + sin delta_k) at angle theta_k.
//!   * hence it suffices to check, for each k, EVERY placement of the sigma_k-square at angle
//!     theta_k whose center is admissible.  For fixed angle the covered weight is a step
//!     function of the center; we minimise it EXACTLY by sweeping the arrangement of
//!     breakpoints (q_a +- h) with an exact sliding-window sweep.
//!   * angles in [45deg,90deg) are covered by the D4 symmetry of the atom set (checked); for a
//!     non-symmetric set the same enumeration is run over all of [0,90deg] (k=0..N).
//!   * the admissible-centre box of a bin uses the bin's MINIMUM bounding-box width cos+sin
//!     (an endpoint, since cos+sin is unimodal on [0,90deg]).
//!   * input is validated: weights >= 0, points inside the container, s_den | s_num*D.
//!
//! Optional diagnostic mode (off unless the environment variable TIGHT_DUMP is set; when off,
//! behaviour and output are unchanged):
//!     TIGHT_DUMP=<path> TIGHT_THRESH=<integer>  verify CERT n N threads topk
//! writes every arrangement cell whose captured weight numerator (over W) is <= TIGHT_THRESH to
//! <path>, one line per cell, so that the NEAR-TIGHT pose set  {Q : weight(Q) <= thresh}  can be
//! analysed (search/tightset.py).  Each bin k gets a header line
//!     bin k N deg DEN cn sn g hh sg_n sg_d wm_n wm_d ntight area_tight area_box nwritten
//! (rotation cos = cn/g, sin = sn/g; deg = 2*atan(k/N) in degrees; every position numerator in
//! the cell lines is over the common denominator DEN = 2*sg_d*wm_d*D*g; the shrunk square has
//! side sg_n/sg_d and half-side hh/DEN; the admissible centre box is [L, s-L]^2 with
//! L = wm_n/(2 wm_d); ntight = number of cells with sum <= thresh in this bin, area_tight =
//! their total area clipped to the rotated centre box (f64 polygon clipping, container units^2),
//! area_box = (s - 2L)^2, nwritten = how many of the ntight cells follow), followed by cell lines
//!     c k a b c0 c1 sum cx cy area
//! meaning: at bin k, every centre with rotated-frame coordinates u0 in (a,b)/DEN, u1 in (c0,c1)/DEN
//! captures exactly weight sum/W with the sigma_k-shrunk square at angle theta_k; (cx,cy) is the
//! cell midpoint mapped back to container coordinates and area the cell's area clipped to the
//! centre box (both f64, for analysis only).  Cells are the raw arrangement cells: the integer
//! ranges may stick out of the admissible centre box (a superset -- clip against [L,s-L]^2
//! rotated by theta_k, as `area` does).  A strip with no atom captured at all (the EMPTY branch) is written as
//! one cell with sum 0 spanning the strip's whole y-range.  Since a true unit square at any angle
//! in the bin captures at least what the shrunk square captures, the dumped set is a SUPERSET of
//! the true near-tight set (the safe direction for a case analysis).
//! The arrangement has ~(2m)^2 cells per bin, so a loose threshold can select 10^9+ cells:
//! TIGHT_MAX=<n> (default 5_000_000) caps the cells written per run -- shared out equally as
//! TIGHT_MAX/(number of bins) per bin, the listed cells being a uniform reservoir sample of the
//! bin's tight cells when there are more -- and TIGHT_MAX=0 writes the bin statistics only.
//! Each bin line is followed by "h k i j area" lines: the clipped area of the bin's tight cells
//! whose midpoint falls in box (i, j) of a G x G grid over the container (G = TIGHT_GRID,
//! default 40), a measure-weighted spatial histogram that does not depend on the cap.
use std::env;
use std::fs;
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::Arc;
use std::thread;

#[derive(Clone)]
struct Cert { s_num: i128, s_den: i128, d: i128, wd: i128, atoms: Vec<(i128,i128,i128)> }

/// Malformed input is rejected with a message and exit status 2 -- never a panic, and never
/// anything that could be mistaken for a verdict (the words VERIFIED / NOT VERIFIED never appear).
fn die(msg: String) -> ! { eprintln!("ERROR: {}", msg); std::process::exit(2) }

fn tok(it: &mut std::str::SplitAsciiWhitespace, what: &str) -> i128 {
    match it.next() {
        None => die(format!("malformed certificate: missing {}", what)),
        Some(t) => t.parse::<i128>().unwrap_or_else(|_| die(format!("malformed certificate: {} is not an integer: {:?}", what, t))),
    }
}

/// Read and VALIDATE a certificate.  Beyond syntax, the checks below are part of the soundness
/// argument (see certificates/FORMAT.md and lean/Sqpack/Basic.lean, `packing_le_weight`):
///   * every weight must be >= 0  (hypothesis `hw` of the Lean theorem; a negative weight placed
///     where no square can see it would lower the total without affecting coverage),
///   * every point must lie in the closed container [0,s]^2,
///   * s*D must be an integer, i.e. s_den | s_num*D (the container side in atom units), which
///     the exact arithmetic below relies on,
///   * exactly m points, nothing after them.
fn read_cert(path: &str) -> Cert {
    let txt = match fs::read_to_string(path) { Ok(t) => t, Err(e) => die(format!("cannot read {}: {}", path, e)) };
    let mut it = txt.split_ascii_whitespace();
    let s_num = tok(&mut it, "s_num"); let s_den = tok(&mut it, "s_den");
    let d = tok(&mut it, "D"); let wd = tok(&mut it, "W");
    let m = tok(&mut it, "m");
    if s_num <= 0 || s_den <= 0 || d <= 0 || wd <= 0 { die("malformed certificate: s_num, s_den, D and W must be positive".to_string()); }
    if m < 0 { die("malformed certificate: m must be >= 0".to_string()); }
    if (s_num * d) % s_den != 0 { die(format!("s_den = {} does not divide s_num*D = {}: the container side must be an integer in atom units (coordinate denominator D must be a multiple of s_den)", s_den, s_num * d)); }
    let sd = s_num * d / s_den;                       // container side in atom units
    let mut atoms = Vec::new();
    for i in 1..=m {
        let x = tok(&mut it, &format!("X_{}", i)); let y = tok(&mut it, &format!("Y_{}", i)); let w = tok(&mut it, &format!("w_{}", i));
        if w < 0 { die(format!("point {} has negative weight {}: weights must be >= 0", i, w)); }
        if x < 0 || y < 0 || x > sd || y > sd { die(format!("point {} = ({}, {})/{} lies outside the container [0, {}/{}]^2", i, x, y, d, s_num, s_den)); }
        atoms.push((x, y, w));
    }
    if it.next().is_some() { die(format!("malformed certificate: trailing data after the {} points declared in the header", m)); }
    Cert{s_num,s_den,d,wd,atoms}
}

/// check the atom multiset is invariant under the dihedral group of the square [0,s]^2
fn check_symmetry(c: &Cert) -> bool {
    let sd = c.s_num * c.d / c.s_den;             // integer: validated in read_cert
    let mut v: Vec<(i128,i128,i128)> = c.atoms.clone(); v.sort();
    let f = |g: &dyn Fn(i128,i128)->(i128,i128)| -> bool {
        let mut u: Vec<(i128,i128,i128)> = c.atoms.iter().map(|&(x,y,w)| { let (a,b)=g(x,y); (a,b,w) }).collect();
        u.sort(); u == v
    };
    f(&|x,y| (sd-x, y)) && f(&|x,y| (x, sd-y)) && f(&|x,y| (y, x))
}

/// exact minimum covered weight over admissible centers, for one angle.
/// Rotation (cos,sin) = (cn,sn)/g.  Square side sigma = sg_n/sg_d.  Returns (min_num, min_den)
/// as weight numerator over c.wd, plus a witness.  All coordinates are handled over the common
/// denominator  M = d*g  for positions, and comparisons are cross-multiplied.
/// Near-tight-set dump state for one bin (TIGHT_DUMP mode; see the module doc).
struct Tight {
    thresh: i128,                                   // cells with sum <= thresh are "tight"
    cap: usize,                                     // at most this many cells are written per bin (reservoir sample)
    cells: Vec<(i128,i128,i128,i128,i128,f64,f64,f64)>, // (a, b, c0, c1, sum, cx, cy, area), positions over DEN
    count: i64,                                     // number of tight cells (written or not)
    area: f64,                                      // their area clipped to the admissible centre box, in container units^2
    grid: usize,                                    // G: the container is split into G x G boxes for the area histogram
    hist: Vec<f64>,                                 // G*G: clipped area of tight cells by midpoint box (row-major, ix*G+iy)
    s: f64,                                         // container side
    rng: u64,                                       // xorshift state for the reservoir (deterministic per bin)
}
/// Area (container units^2) of the u-space rectangle [a,b] x [c0,c1] (numerators over `den`)
/// clipped to the convex polygon `poly` (the admissible centre box rotated into u-space, ccw),
/// by Sutherland-Hodgman in f64.  Used only in TIGHT_DUMP mode.
fn clip_rect_area(a: i128, b: i128, c0: i128, c1: i128, poly: &[[i128;2]], den: f64) -> f64 {
    // shift by (a, c0) before converting to f64 so the cell keeps its full precision
    let sh = |p: [i128;2]| [ (p[0]-a) as f64 / den, (p[1]-c0) as f64 / den ];
    let mut cur: Vec<[f64;2]> = vec![[0.0,0.0], [(b-a) as f64/den, 0.0], [(b-a) as f64/den, (c1-c0) as f64/den], [0.0, (c1-c0) as f64/den]];
    let n = poly.len();
    for i in 0..n {
        let p = sh(poly[i]); let q = sh(poly[(i+1)%n]);
        let (ex, ey) = (q[0]-p[0], q[1]-p[1]);
        let side = |v: [f64;2]| ex*(v[1]-p[1]) - ey*(v[0]-p[0]);   // >= 0 : inside (left of the edge)
        let inp = std::mem::take(&mut cur);
        if inp.is_empty() { return 0.0; }
        let mut sp = inp[inp.len()-1]; let mut dp = side(sp);
        for &v in &inp {
            let dv = side(v);
            if dv >= 0.0 {
                if dp < 0.0 { let t = dp/(dp-dv); cur.push([sp[0]+t*(v[0]-sp[0]), sp[1]+t*(v[1]-sp[1])]); }
                cur.push(v);
            } else if dp >= 0.0 { let t = dp/(dp-dv); cur.push([sp[0]+t*(v[0]-sp[0]), sp[1]+t*(v[1]-sp[1])]); }
            sp = v; dp = dv;
        }
    }
    if cur.len() < 3 { return 0.0; }
    let mut s2 = 0.0;
    for i in 0..cur.len() { let p = cur[i]; let q = cur[(i+1)%cur.len()]; s2 += p[0]*q[1] - q[0]*p[1]; }
    s2.abs()/2.0
}

impl Tight {
    fn add(&mut self, cell: (i128,i128,i128,i128,i128,f64,f64,f64)) {
        let area = cell.7;
        self.count += 1; self.area += area;
        let g = self.grid;
        let ix = (((cell.5 / self.s) * g as f64).floor().max(0.0) as usize).min(g-1);
        let iy = (((cell.6 / self.s) * g as f64).floor().max(0.0) as usize).min(g-1);
        self.hist[ix*g + iy] += area;
        if self.cells.len() < self.cap { self.cells.push(cell); }
        else if self.cap > 0 {
            // reservoir sampling: keep each of the `count` cells with probability cap/count
            self.rng ^= self.rng << 13; self.rng ^= self.rng >> 7; self.rng ^= self.rng << 17;
            let j = (self.rng % (self.count as u64)) as usize;
            if j < self.cap { self.cells[j] = cell; }
        }
    }
}

/// `tight`: in dump mode every cell with sum <= thresh is counted (and appended, up to the cap)
/// as (a, b, c0, c1, sum, cx, cy) -- positions over DEN -- and an EMPTY strip is recorded (sum 0)
/// instead of returning early.  With `tight == None` the function is exactly the original.
fn min_cover_k(c: &Cert, cn: i128, sn: i128, g: i128, sg_n: i128, sg_d: i128, wm_n: i128, wm_d: i128, topk: usize,
               mut tight: Option<&mut Tight>) -> (i128, Vec<(i128,f64,f64)>) {
    // sg_n/sg_d = sigma (side of shrunken square), wm_n/wm_d = lower bound on min width
    let m = c.atoms.len();
    // Scaling.  Positions have denominator d.  Rotation has denominator g.  sigma = sg_n/sg_d.
    // Work with numerators over the common denominator  DEN = 2 * sg_d * wm_d * d * g.
    let k1: i128 = 2 * sg_d * wm_d;              // multiplier for rotated atom numerators
    let mut q: Vec<(i128,i128,i128)> = Vec::with_capacity(m);   // (q0, q1, weight)
    for &(x,y,w) in &c.atoms { q.push(((cn*x + sn*y)*k1, (-sn*x + cn*y)*k1, w)); }
    // h = sigma/2 = sg_n/(2 sg_d)  ->  over DEN:  h = sg_n * wm_d * d * g
    let hh: i128 = sg_n * wm_d * c.d * g;
    // admissible centre box [L,U]^2, L = wm/2 = wm_n/(2 wm_d), U = s - L
    // numerators over DEN/g  (i.e. over 2*sg_d*wm_d*d):
    let l_c: i128 = sg_d * wm_n * c.d;                                   // L * (2 sg_d wm_d d)
    let u_c: i128 = 2*sg_d*wm_d*c.d*c.s_num/c.s_den - sg_d*wm_n*c.d;     // U * (same)
    if u_c <= l_c { return (i128::MAX, Vec::new()); }
    // rotate the four corners into u-space: multiply by (cn,sn) -> numerators over DEN
    let corners: [[i128;2];4] = [[l_c,l_c],[u_c,l_c],[u_c,u_c],[l_c,u_c]];
    let rot: Vec<[i128;2]> = corners.iter().map(|p| [ cn*p[0] + sn*p[1], -sn*p[0] + cn*p[1] ]).collect();
    let px0 = rot.iter().map(|p| p[0]).min().unwrap();
    let px1 = rot.iter().map(|p| p[0]).max().unwrap();
    // x-breakpoints
    let mut bx: Vec<i128> = Vec::with_capacity(2*m+2);
    for t in &q { bx.push(t.0-hh); bx.push(t.0+hh); }
    bx.push(px0); bx.push(px1);
    bx.sort_unstable(); bx.dedup();
    let mut qs = q.clone(); qs.sort_unstable_by_key(|t| t.0);
    let qx: Vec<i128> = qs.iter().map(|t| t.0).collect();
    let mut best = i128::MAX;
    let mut heap: Vec<(i128,f64,f64)> = Vec::new();
    let mut den_f = 0.0f64;
    for w in 0..bx.len()-1 {
        let (a,b) = (bx[w], bx[w+1]);
        if b <= px0 || a >= px1 { continue; }
        // active: q0 in (b-h, a+h)   [strict: covered for every centre in the CLOSED cell]
        // CLOSED squares: atom counts iff |q0 - u0| <= h for every u0 in [a,b]
        let i0 = qx.partition_point(|&v| v <  b-hh);
        let i1 = qx.partition_point(|&v| v <= a+hh);
        // y-range of the admissible (rotated) square over the strip [a,b].
        // Computed in f64 and WIDENED: a superset is always sound (it can only add centres).
        let mut ylo_f = f64::INFINITY; let mut yhi_f = f64::NEG_INFINITY;
        let (af, bf) = (a as f64, b as f64);
        for i in 0..4 {
            let p = [rot[i][0] as f64, rot[i][1] as f64];
            let r = [rot[(i+1)%4][0] as f64, rot[(i+1)%4][1] as f64];
            if p[0] >= af && p[0] <= bf { if p[1]<ylo_f {ylo_f=p[1];} if p[1]>yhi_f {yhi_f=p[1];} }
            for &xc in &[af,bf] {
                let (x0,x1)=(p[0],r[0]);
                if (x0 <= xc && xc <= x1) || (x1 <= xc && xc <= x0) {
                    if x1 != x0 {
                        let y = p[1] + (r[1]-p[1])*(xc-x0)/(x1-x0);
                        if y<ylo_f {ylo_f=y;} if y>yhi_f {yhi_f=y;}
                    } else {
                        if p[1]<ylo_f {ylo_f=p[1];} if p[1]>yhi_f {yhi_f=p[1];}
                        if r[1]<ylo_f {ylo_f=r[1];} if r[1]>yhi_f {yhi_f=r[1];}
                    }
                }
            }
        }
        if !(ylo_f <= yhi_f) { continue; }
        let pad = 1e-9*(yhi_f-ylo_f).abs().max(1.0) + 1e-6*(hh as f64);
        let ylo = (ylo_f - pad).floor() as i128; let yhi = (yhi_f + pad).ceil() as i128;
        if ylo > yhi { continue; }
        if i1 <= i0 {
            if std::env::var("DBG").is_ok() {
                let dn = (2*sg_d*wm_d*c.d*g) as f64;
                eprintln!("  EMPTY strip u0 in [{:.6},{:.6}] (width {:.3e}), yrange [{:.6},{:.6}]",
                    (a as f64)/dn*(g as f64), (b as f64)/dn*(g as f64), ((b-a) as f64)/dn*(g as f64),
                    (ylo as f64)/dn*(g as f64), (yhi as f64)/dn*(g as f64));
            }
            let (u0,u1) = (((a+b)/2) as f64, ((ylo+yhi)/2) as f64);
            let dn = (2*sg_d*wm_d*c.d*g) as f64;
            let cx = ((cn as f64)*u0 - (sn as f64)*u1)/(g as f64)/dn;
            let cy = ((sn as f64)*u0 + (cn as f64)*u1)/(g as f64)/dn;
            if let Some(t) = tight.as_mut() {
                // dump mode: record the empty strip as one cell of sum 0 and keep sweeping
                if 0 <= t.thresh {
                    let area = clip_rect_area(a, b, ylo, yhi, &rot, dn);
                    t.add((a, b, ylo, yhi, 0, cx, cy, area));
                }
                if 0 < best { best = 0; }
                continue;
            }
            return (0, vec![(0,cx,cy)]);
        }                 // a whole strip of admissible centres covers nothing
        let mut ys: Vec<(i128,i128)> = qs[i0..i1].iter().map(|t| (t.1, t.2)).collect();
        ys.sort_unstable();
        let yv: Vec<i128> = ys.iter().map(|t| t.0).collect();
        let mut pre: Vec<i128> = Vec::with_capacity(ys.len()+1); pre.push(0);
        for t in &ys { pre.push(pre.last().unwrap() + t.1); }
        // y-breakpoints
        let mut by: Vec<i128> = Vec::with_capacity(2*ys.len()+2);
        for t in &ys { by.push(t.0-hh); by.push(t.0+hh); }
        by.push(ylo); by.push(yhi);
        by.sort_unstable(); by.dedup();
        for t in 0..by.len()-1 {
            let (c0,c1) = (by[t], by[t+1]);
            if c1 <= ylo || c0 >= yhi { continue; }
            let j0 = yv.partition_point(|&v| v <  c1-hh);
            let j1 = yv.partition_point(|&v| v <= c0+hh);
            let sum = if j1>j0 { pre[j1]-pre[j0] } else { 0 };
            if sum < best { best = sum;
                if std::env::var("DBG").is_ok() {
                    let dn = ((2*sg_d*wm_d*c.d) as f64)*(g as f64)/(g as f64);
                    let sc = (2*sg_d*wm_d*c.d*g) as f64;
                    let _ = dn;
                    eprintln!("  MIN cell: u0 in [{:.6},{:.6}], u1 in [{:.6},{:.6}]  yrange=[{:.6},{:.6}]  nx={} ny={} sum={}",
                        (a as f64)/sc, (b as f64)/sc, (c0 as f64)/sc, (c1 as f64)/sc,
                        (ylo as f64)/sc, (yhi as f64)/sc, i1-i0, j1-j0, sum);
                    eprintln!("    h = {:.6}", (hh as f64)/sc);
                    let mut lst: Vec<String> = ys.iter().map(|t| format!("{:.4}", (t.0 as f64)/sc)).collect();
                    lst.sort(); eprintln!("    active y: {}", lst.join(","));
                }
            }
            if let Some(t) = tight.as_mut() { if sum <= t.thresh {
                let dnf = (2*sg_d*wm_d*c.d*g) as f64;
                // area of the cell clipped to the rotated admissible centre box (f64 polygon clip)
                let area = clip_rect_area(a, b, c0, c1, &rot, dnf);
                let (u0,u1) = (((a+b)/2) as f64, ((c0+c1)/2) as f64);
                let cx = ((cn as f64)*u0 - (sn as f64)*u1)/(g as f64)/dnf;
                let cy = ((sn as f64)*u0 + (cn as f64)*u1)/(g as f64)/dnf;
                t.add((a, b, c0, c1, sum, cx, cy, area));
            } }
            if topk>0 && sum < c.wd {
                // witness centre: midpoint of the cell, rotated back to original coordinates
                let (u0,u1) = (((a+b)/2) as f64, ((c0+c1)/2) as f64);
                if den_f==0.0 { den_f = (2*sg_d*wm_d*c.d*g) as f64; }
                let cx = ((cn as f64)*u0 - (sn as f64)*u1)/(g as f64)/den_f;
                let cy = ((sn as f64)*u0 + (cn as f64)*u1)/(g as f64)/den_f;
                heap.push((sum,cx,cy));
                if heap.len()>4*topk { heap.sort_unstable_by_key(|t| t.0); heap.truncate(topk); }
            }
        }
    }
    heap.sort_unstable_by_key(|t| t.0); heap.truncate(topk);
    (best, heap)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 { die("usage: verify <certificate> <n> <N> [threads] [topk] [witness-file]".to_string()); }
    let cert = read_cert(&args[1]);
    let nn: i128 = args[2].parse().unwrap_or_else(|_| die(format!("n must be an integer, got {:?}", args[2])));   // claim: fewer than nn squares fit
    let bigN: i128 = args[3].parse().unwrap_or_else(|_| die(format!("N must be an integer, got {:?}", args[3])));  // angle parameter
    if nn < 1 || bigN < 1 { die("n and N must be >= 1".to_string()); }
    let threads: usize = if args.len()>4 { args[4].parse().unwrap_or_else(|_| die("threads must be an integer".to_string())) } else { 16 };
    let topk: usize = if args.len()>5 { args[5].parse().unwrap_or_else(|_| die("topk must be an integer".to_string())) } else { 0 };
    if threads < 1 { die("threads must be >= 1".to_string()); }
    // optional near-tight-set dump (see the module doc): TIGHT_DUMP=<path> TIGHT_THRESH=<int>
    let tight_thresh: Option<i128> = match env::var("TIGHT_DUMP") {
        Err(_) => None,
        Ok(_) => Some(env::var("TIGHT_THRESH").ok()
            .and_then(|s| s.trim().parse::<i128>().ok())
            .unwrap_or_else(|| die("TIGHT_DUMP is set but TIGHT_THRESH is missing or not an integer (weight numerator over W)".to_string()))),
    };
    // TIGHT_MAX: cap on the number of cells WRITTEN per run (default 5_000_000; 0 = per-bin
    // statistics only).  Cells beyond the cap are still counted and measured in the bin lines.
    let tight_cap: usize = match tight_thresh {
        None => 0,
        Some(_) => env::var("TIGHT_MAX").ok().map(|s| s.trim().parse::<usize>().unwrap_or_else(|_| die("TIGHT_MAX must be a non-negative integer".to_string()))).unwrap_or(5_000_000),
    };
    // TIGHT_GRID: resolution G of the per-bin G x G area histogram of tight cells (default 40).
    let tight_grid: usize = match tight_thresh {
        None => 1,
        Some(_) => env::var("TIGHT_GRID").ok().map(|s| s.trim().parse::<usize>().ok().filter(|&g| g >= 1).unwrap_or_else(|| die("TIGHT_GRID must be a positive integer".to_string()))).unwrap_or(40),
    };
    let tight_out: Option<Arc<std::sync::Mutex<std::io::BufWriter<fs::File>>>> = match (&tight_thresh, env::var("TIGHT_DUMP")) {
        (Some(th), Ok(path)) => {
            let f = fs::File::create(&path).unwrap_or_else(|e| die(format!("cannot create TIGHT_DUMP file {}: {}", path, e)));
            let mut w = std::io::BufWriter::new(f);
            use std::io::Write;
            let _ = writeln!(w, "# TIGHT_DUMP cert={} n={} N={} thresh={} W={} D={} s={}/{} max_cells={} grid={}", args[1], args[2], args[3], th, cert.wd, cert.d, cert.s_num, cert.s_den, tight_cap, tight_grid);
            let _ = writeln!(w, "# h k i j area   (area of this bin's tight cells whose midpoint lies in box [i,i+1)/G x [j,j+1)/G of the container, G=grid)");
            let _ = writeln!(w, "# bin k N deg DEN cn sn g hh sg_n sg_d wm_n wm_d ntight area_tight area_box nwritten   (cos=cn/g sin=sn/g; positions over DEN=2*sg_d*wm_d*D*g; half-side hh/DEN; centre box [L,s-L]^2, L=wm_n/(2 wm_d); ntight = cells with sum<=thresh, area_tight = their area clipped to the admissible centre box (container units^2, f64), area_box = (s-2L)^2, nwritten = cells listed below)");
            let _ = writeln!(w, "# c k a b c0 c1 sum cx cy area   (u0 in (a,b)/DEN, u1 in (c0,c1)/DEN capture sum/W; (cx,cy) = cell midpoint in container coords; area = clipped to the centre box)");
            Some(Arc::new(std::sync::Mutex::new(w)))
        }
        _ => None,
    };
    let tight_count = Arc::new(AtomicI64::new(0));
    let tight_written = Arc::new(AtomicI64::new(0));
    let out: Arc<std::sync::Mutex<Vec<(i128,f64,f64,f64)>>> = Arc::new(std::sync::Mutex::new(Vec::new()));
    // D4 symmetry lets us restrict angles to [0,45]; without it we cover [0,90).
    let sym = check_symmetry(&cert);
    println!("D4-symmetric atom set: {}  -> angles cover {}", sym, if sym {"[0,45] deg"} else {"[0,90) deg"});
    let total: i128 = cert.atoms.iter().map(|a| a.2).sum();
    println!("atoms={} total weight = {}/{} = {:.6}", cert.atoms.len(), total, cert.wd, total as f64/cert.wd as f64);
    println!("s = {}/{} = {:.6}", cert.s_num, cert.s_den, cert.s_num as f64/cert.s_den as f64);
    let weight_ok = total < nn * cert.wd;   // EXACT integer comparison
    if !weight_ok { println!("WEIGHT NOT < n={} (total*WD = {} vs {})", nn, total, nn*cert.wd); }
    // angle list k=0..K with t=k/N, need 2*arctan(K/N) >= 45deg  <=>  (K/N+1)^2 >= 2
    let mut kk: i128 = 0;
    if sym { while (kk + bigN)*(kk + bigN) < 2*bigN*bigN { kk += 1; } } else { kk = bigN; }
    println!("angles: k=0..{} (N={}), covering {}", kk, bigN, if sym {"[0,45deg]"} else {"[0,90deg]"});
    let bad = Arc::new(AtomicI64::new(0));
    let minw = Arc::new(std::sync::Mutex::new((i128::MAX, 0i128)));
    let cert = Arc::new(cert);
    let mut hs = Vec::new();
    let chunk = (kk as usize + threads) / threads;
    for t in 0..threads {
        let cert = cert.clone(); let bad = bad.clone(); let minw = minw.clone(); let out = out.clone();
        let tight_out = tight_out.clone(); let tight_count = tight_count.clone(); let tight_written = tight_written.clone();
        hs.push(thread::spawn(move || {
            for k in (t*chunk)..(((t+1)*chunk).min(kk as usize)) {
                let k = k as i128;
                let (c0,s0,g0) = (bigN*bigN - k*k, 2*k*bigN, bigN*bigN + k*k);
                let k2 = k+1;
                let (c1,s1,g1) = (bigN*bigN - k2*k2, 2*k2*bigN, bigN*bigN + k2*k2);
                // cos(delta), sin(delta) with delta = theta_{k+1}-theta_k, over g0*g1
                let cd = c0*c1 + s0*s1; let sd = c0*s1 - s0*c1;
                // sigma = 1/(cos d + sin d) = g0*g1/(cd+sd), ROUNDED DOWN to 1e-6 (sound: smaller square)
                const SCALE: i128 = 1_000_000;
                let sg_n = (g0*g1*SCALE)/(cd + sd); let sg_d = SCALE;
                // w_min = min over the bin of cos(theta)+sin(theta): unimodal on [0,90deg] with its
                // maximum at 45deg, so the minimum is at an endpoint theta_k or theta_{k+1}.
                // ROUNDED DOWN (sound: bigger centre box).  Using theta_k alone is wrong above 45deg.
                let wm_k0 = ((c0 + s0)*SCALE)/g0; let wm_k1 = ((c1 + s1)*SCALE)/g1;
                let wm_n = wm_k0.min(wm_k1); let wm_d = SCALE;
                let mut tstate = tight_thresh.map(|th| {
                    // per-bin share of the global cap (the cap is on the whole run; bins are of similar size)
                    let per_bin = if tight_cap == 0 { 0 } else { (tight_cap / (kk as usize + 1)).max(1) };
                    Tight { thresh: th, cap: per_bin, cells: Vec::new(), count: 0, area: 0.0,
                            grid: tight_grid, hist: vec![0.0; tight_grid*tight_grid],
                            s: (cert.s_num as f64)/(cert.s_den as f64), rng: 0x9E3779B97F4A7C15u64 ^ ((k as u64 + 1) * 0x2545F4914F6CDD1Du64) }
                });
                let (v, wit) = min_cover_k(&cert, c0, s0, g0, sg_n, sg_d, wm_n, wm_d, topk, tstate.as_mut());
                if let (Some(tw), Some(ts)) = (&tight_out, &tstate) {
                    use std::io::Write;
                    let den = 2*sg_d*wm_d*cert.d*g0;
                    let hh = sg_n*wm_d*cert.d*g0;
                    let deg = 2.0*((k as f64)/(bigN as f64)).atan().to_degrees();
                    let side_box = (cert.s_num as f64)/(cert.s_den as f64) - (wm_n as f64)/(wm_d as f64);
                    let area_box = if side_box > 0.0 { side_box*side_box } else { 0.0 };
                    let mut buf = format!("bin {} {} {:.9} {} {} {} {} {} {} {} {} {} {} {:.12e} {:.12e} {}\n",
                        k, bigN, deg, den, c0, s0, g0, hh, sg_n, sg_d, wm_n, wm_d, ts.count, ts.area, area_box, ts.cells.len());
                    for (idx, v) in ts.hist.iter().enumerate() {
                        if *v > 0.0 { buf.push_str(&format!("h {} {} {} {:.9e}\n", k, idx / ts.grid, idx % ts.grid, v)); }
                    }
                    for &(a,b,cc0,cc1,sum,cx,cy,area) in &ts.cells {
                        buf.push_str(&format!("c {} {} {} {} {} {} {:.9} {:.9} {:.9e}\n", k, a, b, cc0, cc1, sum, cx, cy, area));
                    }
                    tight_count.fetch_add(ts.count, Ordering::Relaxed);
                    tight_written.fetch_add(ts.cells.len() as i64, Ordering::Relaxed);
                    let mut w = tw.lock().unwrap();
                    w.write_all(buf.as_bytes()).unwrap_or_else(|e| die(format!("TIGHT_DUMP write failed: {}", e)));
                }
                if topk>0 { let th = 2.0*((k as f64)/(bigN as f64)).atan();
                    let mut o=out.lock().unwrap();
                    for (val,cx,cy) in wit { o.push((val, th, cx, cy)); } }
                let mut mg = minw.lock().unwrap();
                if v < mg.0 { *mg = (v, k); }
                drop(mg);
                if v < cert.wd { bad.fetch_add(1, Ordering::Relaxed); if topk==0 && bad.load(Ordering::Relaxed)<8 { println!("  FAIL at angle k={}: covered {}/{}", k, v, cert.wd); } }
            }
        }));
    }
    for h in hs { h.join().unwrap(); }
    if let Some(tw) = &tight_out {
        use std::io::Write;
        tw.lock().unwrap().flush().unwrap_or_else(|e| die(format!("TIGHT_DUMP flush failed: {}", e)));
        println!("TIGHT_DUMP: {} cells with captured weight <= {}/{}; wrote {} of them to {} (TIGHT_MAX={})",
            tight_count.load(Ordering::Relaxed), tight_thresh.unwrap(), cert.wd, tight_written.load(Ordering::Relaxed), env::var("TIGHT_DUMP").unwrap(), tight_cap);
    }
    if topk>0 {
        let mut o = out.lock().unwrap();
        o.sort_by(|a,b| a.0.cmp(&b.0));
        let mut f = String::new();
        for (v,th,cx,cy) in o.iter() { f.push_str(&format!("{} {} {} {}\n", *v as f64/cert.wd as f64, th, cx, cy)); }
        let outp = if args.len()>6 { args[6].clone() } else { "sep.txt".to_string() }; fs::write(&outp, f).unwrap();
        println!("wrote {} violated placements", o.len());
    }
    let mg = minw.lock().unwrap();
    println!("min covered weight over ALL placements = {}/{} = {:.6}  (at angle k={})", mg.0, cert.wd, mg.0 as f64/cert.wd as f64, mg.1);
    if bad.load(Ordering::Relaxed)==0 && weight_ok {
        println!("VERIFIED: every CLOSED unit square inside C covers weight >= 1, and total weight < {}.\n==> {} unit squares cannot be packed into any square of side < {}/{} = {:.9},\n    i.e.  s({}) >= {:.9}", nn, nn, cert.s_num, cert.s_den, cert.s_num as f64/cert.s_den as f64, nn, cert.s_num as f64/cert.s_den as f64);
    } else { println!("NOT VERIFIED"); }
}
