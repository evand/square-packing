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

/// Optional BRANCH data (see certificates/FORMAT.md, "Branch certificates").  The four corner
/// boxes are [0,r]^2 and its images under the symmetries of the container; a pose belongs to the
/// region iff its CENTRE lies in one of them.  The claim becomes: every closed unit square inside
/// C with centre in a corner box captures weight >= 1 + lam/W, every other one >= 1, and
/// total - lam*k < n*W.  Then no packing of n unit squares with EXACTLY k squares centred in the
/// corner boxes fits in any square of side < s (n <= sum_i w(S_i) - lam*k <= W - lam*k).
/// Since the centres of two interior-disjoint unit squares are >= 1 apart and the admissible
/// part [1/2, r]^2 of a corner box has diameter (r - 1/2)*sqrt(2) < 1 (checked exactly), each box
/// holds at most one square, so k ranges over 0..4 and the five branch certificates together
/// rule out every packing.  The region is D4-symmetric, so the [0,45deg] reduction is unaffected.
/// `lam[j]`, `k[j]` are per box (box 0 = [0,r]^2, 1 = [s-r,s]x[0,r], 2 = [0,r]x[s-r,s], 3 = [s-r,s]^2).
/// The one-number trailer form `lambda L / k K` means lam = [L;4] and K = total count (`kdot = L*K`);
/// the four-number form has k[j] in {0,1} and `kdot = sum_j lam[j]*k[j]`.  Unequal lambdas break the
/// D4 symmetry of the threshold, so the angle range is then the full [0,90deg).
#[derive(Clone)]
struct Region { r_num: i128, r_den: i128, lam: [i128;4], k: [i128;4], kdot: i128, single: bool, ktot: i128 }

#[derive(Clone)]
struct Cert { s_num: i128, s_den: i128, d: i128, wd: i128, atoms: Vec<(i128,i128,i128)>, region: Option<Region> }

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
    let region = match it.next() {
        None => None,
        Some("region") => {
            match it.next() { Some("corner") => (), other => die(format!("region trailer: unknown region kind {:?} (only 'corner' is supported)", other)) }
            let r_num = tok(&mut it, "r_num"); let r_den = tok(&mut it, "r_den");
            match it.next() { Some("lambda") => (), other => die(format!("region trailer: expected 'lambda', got {:?}", other)) }
            // one or four integers, up to the keyword `k`
            let mut lams: Vec<i128> = Vec::new();
            loop {
                match it.next() {
                    Some("k") => break,
                    Some(t) => lams.push(t.parse::<i128>().unwrap_or_else(|_| die(format!("region trailer: lambda value is not an integer: {:?}", t)))),
                    None => die("region trailer: expected 'k' after the lambda value(s)".to_string()),
                }
            }
            let mut ks: Vec<i128> = Vec::new();
            for t in it.by_ref() { ks.push(t.parse::<i128>().unwrap_or_else(|_| die(format!("malformed certificate: trailing data after the region trailer (got {:?})", t)))); }
            if r_num <= 0 || r_den <= 0 { die("region trailer: r must be positive".to_string()); }
            // one square per corner box needs (r - 1/2) sqrt2 < 1, i.e. (2r - 1)^2 < 2, exactly
            let t = 2 * r_num - r_den;
            if t * t >= 2 * r_den * r_den { die(format!("region trailer: r = {}/{} is too large: the admissible part of a corner box must have diameter < 1, i.e. (2r-1)^2 < 2", r_num, r_den)); }
            let (lam, k, kdot, single, ktot) = match (lams.len(), ks.len()) {
                (1, 1) => {
                    let (l, kk) = (lams[0], ks[0]);
                    if kk < 0 || kk > 4 { die(format!("region trailer: k = {} must be in 0..4 (one square per corner box)", kk)); }
                    ([l;4], [kk;4], l * kk, true, kk)
                }
                (4, 4) => {
                    let mut lam = [0i128;4]; let mut k = [0i128;4]; let mut kdot = 0; let mut ktot = 0;
                    for j in 0..4 {
                        if ks[j] < 0 || ks[j] > 1 { die(format!("region trailer: per-box k[{}] = {} must be 0 or 1", j, ks[j])); }
                        lam[j] = lams[j]; k[j] = ks[j]; kdot += lams[j] * ks[j]; ktot += ks[j];
                    }
                    (lam, k, kdot, false, ktot)
                }
                (a, b) => die(format!("region trailer: expected 1 or 4 lambda values and the same number of k values, got {} and {}", a, b)),
            };
            Some(Region { r_num, r_den, lam, k, kdot, single, ktot })
        }
        Some(t) => die(format!("malformed certificate: trailing data after the {} points declared in the header (got {:?})", m, t)),
    };
    Cert{s_num,s_den,d,wd,atoms,region}
}

/// u1-range of the convex polygon `poly` (in u-space) over the vertical strip u0 in [af, bf];
/// (INF, -INF) if they do not meet.  f64; callers widen the result.
fn yrange(poly: &[[f64;2]], af: f64, bf: f64) -> (f64, f64) {
    let mut ylo_f = f64::INFINITY; let mut yhi_f = f64::NEG_INFINITY;
    let n = poly.len();
    for i in 0..n {
        let p = poly[i]; let r = poly[(i+1)%n];
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
    (ylo_f, yhi_f)
}

/// Region membership of one arrangement cell, conservatively in both directions.
/// The cell is the u-space rectangle [u0a,u0b] x [u1a,u1b] (container units); its image in the
/// container is the convex hull of the four rotated corners, so it lies inside an axis-parallel box
/// iff its bounding box does, and it can meet the box only if its bounding box does.  Returns
/// (may_meet_a_corner_box, certainly_inside_a_corner_box); `pad` widens the first test and
/// shrinks the second, so f64 rounding can only make the verdict stricter.
/// Returns (bitmask of boxes the cell may meet, index of a box it certainly lies inside or 4).
fn region_flags(u0a: f64, u0b: f64, u1a: f64, u1b: f64, cnf: f64, snf: f64, gf: f64,
                boxes: &[[f64;2];4], rf: f64, pad: f64) -> (u8, usize) {
    let (mut x0, mut x1, mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for &(u0,u1) in &[(u0a,u1a),(u0b,u1a),(u0a,u1b),(u0b,u1b)] {
        let x = (cnf*u0 - snf*u1)/gf; let y = (snf*u0 + cnf*u1)/gf;
        if x<x0 {x0=x;} if x>x1 {x1=x;} if y<y0 {y0=y;} if y>y1 {y1=y;}
    }
    let mut may: u8 = 0; let mut inside: usize = 4;
    for (j, b) in boxes.iter().enumerate() {
        if x1 >= b[0]-pad && x0 <= b[0]+rf+pad && y1 >= b[1]-pad && y0 <= b[1]+rf+pad { may |= 1 << j; }
        if x0 >= b[0]+pad && x1 <= b[0]+rf-pad && y0 >= b[1]+pad && y1 <= b[1]+rf-pad { inside = j; }
    }
    (may, inside)
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
    let cur = clip_rect_poly(a, b, c0, c1, poly, den);
    if cur.len() < 3 { return 0.0; }
    let mut s2 = 0.0;
    for i in 0..cur.len() { let p = cur[i]; let q = cur[(i+1)%cur.len()]; s2 += p[0]*q[1] - q[0]*p[1]; }
    s2.abs()/2.0
}
/// The u-space rectangle [a,b] x [c0,c1] (numerators over `den`) clipped to the convex polygon
/// `poly` (ccw), by Sutherland-Hodgman in f64; vertices are returned RELATIVE to (a, c0), in
/// container units (the shift keeps the cell's full precision).  Empty if they do not meet.
fn clip_rect_poly(a: i128, b: i128, c0: i128, c1: i128, poly: &[[i128;2]], den: f64) -> Vec<[f64;2]> {
    let sh = |p: [i128;2]| [ (p[0]-a) as f64 / den, (p[1]-c0) as f64 / den ];
    let mut cur: Vec<[f64;2]> = vec![[0.0,0.0], [(b-a) as f64/den, 0.0], [(b-a) as f64/den, (c1-c0) as f64/den], [0.0, (c1-c0) as f64/den]];
    let n = poly.len();
    for i in 0..n {
        let p = sh(poly[i]); let q = sh(poly[(i+1)%n]);
        let (ex, ey) = (q[0]-p[0], q[1]-p[1]);
        let side = |v: [f64;2]| ex*(v[1]-p[1]) - ey*(v[0]-p[0]);   // >= 0 : inside (left of the edge)
        let inp = std::mem::take(&mut cur);
        if inp.is_empty() { return cur; }
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
    cur
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
               mut tight: Option<&mut Tight>) -> (i128, Vec<(i128,f64,f64,u8)>) {
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
    let rotf: Vec<[f64;2]> = rot.iter().map(|p| [p[0] as f64, p[1] as f64]).collect();
    let px0 = rot.iter().map(|p| p[0]).min().unwrap();
    let px1 = rot.iter().map(|p| p[0]).max().unwrap();
    // Region thresholds (branch certificates): f64 in container units, conservative both ways.
    let den_all = (2*sg_d*wm_d*c.d*g) as f64;        // DEN
    let sf = c.s_num as f64 / c.s_den as f64;
    let (cnf, snf, gf) = (cn as f64, sn as f64, g as f64);
    const RPAD: f64 = 1e-7;
    let reg = c.region.as_ref().map(|r| {
        let rf = r.r_num as f64 / r.r_den as f64;
        let boxes: [[f64;2];4] = [[0.0, 0.0], [sf - rf, 0.0], [0.0, sf - rf], [sf - rf, sf - rf]];   // lower-left corners
        let polys: Vec<[[f64;2];4]> = boxes.iter().map(|b| {
            let cs = [[b[0], b[1]], [b[0]+rf, b[1]], [b[0]+rf, b[1]+rf], [b[0], b[1]+rf]];
            let mut p = [[0.0;2];4];
            for i in 0..4 { p[i] = [(cnf*cs[i][0] + snf*cs[i][1])/gf, (-snf*cs[i][0] + cnf*cs[i][1])/gf]; }
            p
        }).collect();
        (rf, r.lam, boxes, polys)
    });
    // required extra weight (numerator over W) for a cell with flags (boxes it may meet, box it
    // lies inside): no contact -> 0; entirely inside box j -> lam[j]; straddling -> max over the
    // boxes met of max(lam[j], 0)
    let lams: [i128;4] = reg.as_ref().map(|r| r.1).unwrap_or([0;4]);
    let required = |may: u8, inside: usize| -> i128 {
        if may == 0 { return 0; }
        if inside < 4 { return lams[inside]; }
        let mut req = 0; for j in 0..4 { if may & (1 << j) != 0 && lams[j] > req { req = lams[j]; } } req
    };
    // x-breakpoints
    let mut bx: Vec<i128> = Vec::with_capacity(2*m+2);
    for t in &q { bx.push(t.0-hh); bx.push(t.0+hh); }
    bx.push(px0); bx.push(px1);
    bx.sort_unstable(); bx.dedup();
    let mut qs = q.clone(); qs.sort_unstable_by_key(|t| t.0);
    let qx: Vec<i128> = qs.iter().map(|t| t.0).collect();
    let mut best = i128::MAX;
    let mut heap: Vec<(i128,f64,f64,u8)> = Vec::new();
    let mut den_f = 0.0f64;
    let mut bands: Vec<(f64,f64)> = Vec::new();
    for w in 0..bx.len()-1 {
        let (a,b) = (bx[w], bx[w+1]);
        if b <= px0 || a >= px1 { continue; }
        // active: q0 in (b-h, a+h)   [strict: covered for every centre in the CLOSED cell]
        // CLOSED squares: atom counts iff |q0 - u0| <= h for every u0 in [a,b]
        let i0 = qx.partition_point(|&v| v <  b-hh);
        let i1 = qx.partition_point(|&v| v <= a+hh);
        // y-range of the admissible (rotated) square over the strip [a,b].
        // Computed in f64 and WIDENED: a superset is always sound (it can only add centres).
        let (af, bf) = (a as f64, b as f64);
        let (ylo_f, yhi_f) = yrange(&rotf, af, bf);
        if !(ylo_f <= yhi_f) { continue; }
        let pad = 1e-9*(yhi_f-ylo_f).abs().max(1.0) + 1e-6*(hh as f64);
        let ylo = (ylo_f - pad).floor() as i128; let yhi = (yhi_f + pad).ceil() as i128;
        if ylo > yhi { continue; }
        // u1-bands (container units, widened) where this strip meets a corner box; cells outside
        // every band cannot meet the region and need no per-cell test.
        bands.clear();
        if let Some((_, _, _, polys)) = &reg {
            for p in polys.iter() {
                let (lo, hi) = yrange(p, af/den_all, bf/den_all);
                if lo <= hi { bands.push((lo - RPAD, hi + RPAD)); }
            }
        }
        // region flags of a cell [a,b] x [c0,c1] (numerators over DEN)
        let flags = |c0: i128, c1: i128| -> (u8, usize) {
            match &reg {
                None => (0, 4),
                Some((rf, _, boxes, _)) => {
                    let (u1a, u1b) = (c0 as f64/den_all, c1 as f64/den_all);
                    if !bands.iter().any(|&(lo,hi)| u1b >= lo && u1a <= hi) { return (0, 4); }
                    region_flags(af/den_all, bf/den_all, u1a, u1b, cnf, snf, gf, boxes, *rf, RPAD)
                }
            }
        };
        // Witness choice for the LP: a pose strictly inside the cell clipped to the admissible
        // box, on the required side of the region boundary when the cell straddles it (flag 1:
        // inside a corner box, flag 0: outside all of them); the clipped midpoint otherwise.
        // Only the LP heuristics depend on this; the verdict does not.
        let pt_in_r = |u0: f64, u1: f64| -> Option<usize> {      // Some(j)=inside box j, Some(4)=outside all, None=ambiguous
            match &reg {
                None => Some(4),
                Some((rf, _, boxes, _)) => {
                    let x = (cnf*u0 - snf*u1)/gf; let y = (snf*u0 + cnf*u1)/gf;
                    let mut ins = 4usize; let mut near = false;
                    for (j, bb) in boxes.iter().enumerate() {
                        if x >= bb[0]+RPAD && x <= bb[0]+rf-RPAD && y >= bb[1]+RPAD && y <= bb[1]+rf-RPAD { ins = j; }
                        if x >= bb[0]-RPAD && x <= bb[0]+rf+RPAD && y >= bb[1]-RPAD && y <= bb[1]+rf+RPAD { near = true; }
                    }
                    if ins < 4 { Some(ins) } else if !near { Some(4) } else { None }
                }
            }
        };
        // want_in: the box index the witness should lie in, or 4 for "outside every box"
        let witness = |c0: i128, c1: i128, want_in: usize| -> (f64, f64) {
            // vertices of the cell clipped to the admissible polygon, absolute container units of u-space
            let (oa, oc) = (a as f64/den_all, c0 as f64/den_all);
            let verts: Vec<(f64,f64)> = clip_rect_poly(a, b, c0, c1, &rot, den_all).iter().map(|v| (v[0]+oa, v[1]+oc)).collect();
            let mut pick = if verts.len() >= 3 {
                let n = verts.len() as f64;
                (verts.iter().map(|v| v.0).sum::<f64>()/n, verts.iter().map(|v| v.1).sum::<f64>()/n)
            } else {
                let (cc0, cc1) = (c0.max(ylo), c1.min(yhi));
                (((a+b)/2) as f64/den_all, ((cc0+cc1)/2) as f64/den_all)
            };
            if reg.is_some() && verts.len() >= 3 && pt_in_r(pick.0, pick.1) != Some(want_in) {
                let (m0, m1) = pick;
                'outer: for &f in &[0.75f64, 0.5, 0.9, 0.99, 1.0] {
                    for &(x0, y0) in &verts {
                        let (u0, u1) = (f*x0 + (1.0-f)*m0, f*y0 + (1.0-f)*m1);
                        if pt_in_r(u0, u1) == Some(want_in) { pick = (u0, u1); break 'outer; }
                    }
                }
            }
            let (u0, u1) = pick;
            ((cnf*u0 - snf*u1)/gf, (snf*u0 + cnf*u1)/gf)
        };
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
            let (may, inside) = flags(ylo, yhi);
            let v = 0 - required(may, inside);
            let (cx, cy) = if tight.is_none() { witness(ylo, yhi, inside) } else { (cx, cy) };
            if let Some(t) = tight.as_mut() {
                // dump mode: record the empty strip as one cell of sum 0 and keep sweeping
                if 0 <= t.thresh {
                    let area = clip_rect_area(a, b, ylo, yhi, &rot, dn);
                    t.add((a, b, ylo, yhi, 0, cx, cy, area));
                }
                if v < best { best = v; }
                continue;
            }
            if v < best { best = v; }
            if v >= c.wd { continue; }          // only with lam <= -W and the strip inside the region: the corner is free
            return (v, vec![(v, cx, cy, if inside < 4 {(inside + 1) as u8} else {0})]);
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
            let (may, inside) = flags(c0, c1);
            let v = sum - required(may, inside);
            if v < best { best = v;
                if std::env::var("DBG").is_ok() {
                    let sc = (2*sg_d*wm_d*c.d*g) as f64;
                    eprintln!("  MIN cell: u0 in [{:.6},{:.6}], u1 in [{:.6},{:.6}]  yrange=[{:.6},{:.6}]  nx={} ny={} sum={} region=({},{})",
                        (a as f64)/sc, (b as f64)/sc, (c0 as f64)/sc, (c1 as f64)/sc,
                        (ylo as f64)/sc, (yhi as f64)/sc, i1-i0, j1-j0, sum, may, inside);
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
            if topk>0 && v < c.wd {
                // witness centre: midpoint of the cell, rotated back to original coordinates.
                // One witness per violated threshold: flag 1 = "needs 1 + lambda" (cell meets the
                // region), flag 0 = "needs 1" (cell is not entirely inside the region).
                // flag j+1 = "needs 1 + lambda_j" (cell meets box j), flag 0 = "needs 1" (cell not inside a box)
                for j in 0..4 { if may & (1 << j) != 0 && sum < c.wd + lams[j] { let (cx, cy) = witness(c0, c1, j); heap.push((sum - lams[j], cx, cy, (j + 1) as u8)); } }
                if inside == 4 && sum < c.wd { let (cx, cy) = witness(c0, c1, 4); heap.push((sum, cx, cy, 0)); }
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
    let out: Arc<std::sync::Mutex<Vec<(i128,f64,f64,f64,u8)>>> = Arc::new(std::sync::Mutex::new(Vec::new()));
    // D4 symmetry lets us restrict angles to [0,45]; without it we cover [0,90).  A branch
    // certificate with unequal per-box lambdas has a non-symmetric threshold: full range too.
    let sym_atoms = check_symmetry(&cert);
    let sym_lam = cert.region.as_ref().map(|r| r.lam.iter().all(|&l| l == r.lam[0])).unwrap_or(true);
    let sym = sym_atoms && sym_lam;
    println!("D4-symmetric atom set: {}{}  -> angles cover {}", sym_atoms, if sym_lam {""} else {" (but per-box lambdas differ)"}, if sym {"[0,45] deg"} else {"[0,90) deg"});
    let total: i128 = cert.atoms.iter().map(|a| a.2).sum();
    println!("atoms={} total weight = {}/{} = {:.6}", cert.atoms.len(), total, cert.wd, total as f64/cert.wd as f64);
    println!("s = {}/{} = {:.6}", cert.s_num, cert.s_den, cert.s_num as f64/cert.s_den as f64);
    let (lam, kreg) = match &cert.region {
        None => (0, 0),
        Some(r) => {
            if tight_thresh.is_some() { die("TIGHT_DUMP is not supported for branch certificates".to_string()); }
            if r.single {
                println!("BRANCH certificate: corner boxes [0, {}/{}]^2 = [0, {:.6}]^2 and images; lambda = {}/{} = {:+.7}; k = {} squares centred in the boxes",
                         r.r_num, r.r_den, r.r_num as f64/r.r_den as f64, r.lam[0], cert.wd, r.lam[0] as f64/cert.wd as f64, r.ktot);
                println!("  claim: squares centred in a corner box capture >= 1 + lambda, all others >= 1, and total - lambda*k < n");
            } else {
                println!("BRANCH certificate: corner boxes [0, {}/{}]^2 = [0, {:.6}]^2 and images, PER BOX: lambda = [{}] k = [{}]",
                         r.r_num, r.r_den, r.r_num as f64/r.r_den as f64,
                         r.lam.iter().map(|l| format!("{:+.7}", *l as f64/cert.wd as f64)).collect::<Vec<_>>().join(", "),
                         r.k.iter().map(|k| k.to_string()).collect::<Vec<_>>().join(", "));
                println!("  claim: squares centred in box j capture >= 1 + lambda_j, all others >= 1, and total - sum_j lambda_j*k_j < n");
            }
            (r.kdot, 1)
        }
    };
    let weight_ok = total - lam * kreg < nn * cert.wd;   // EXACT integer comparison (lam*kreg = sum_j lambda_j k_j)
    if !weight_ok { println!("WEIGHT NOT < n={} (total*WD - sum lambda*k = {} vs {})", nn, total - lam*kreg, nn*cert.wd); }
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
                    for (val,cx,cy,fl) in wit { o.push((val, th, cx, cy, fl)); } }
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
        for (v,th,cx,cy,fl) in o.iter() { f.push_str(&format!("{} {} {} {} {}\n", *v as f64/cert.wd as f64, th, cx, cy, fl)); }
        let outp = if args.len()>6 { args[6].clone() } else { "sep.txt".to_string() }; fs::write(&outp, f).unwrap();
        println!("wrote {} violated placements", o.len());
    }
    let mg = minw.lock().unwrap();
    if cert.region.is_some() {
        println!("min covered weight minus region threshold over ALL placements = {}/{} = {:.6}  (at angle k={})", mg.0, cert.wd, mg.0 as f64/cert.wd as f64, mg.1);
    } else {
        println!("min covered weight over ALL placements = {}/{} = {:.6}  (at angle k={})", mg.0, cert.wd, mg.0 as f64/cert.wd as f64, mg.1);
    }
    if bad.load(Ordering::Relaxed)==0 && weight_ok {
        match &cert.region {
            None => println!("VERIFIED: every CLOSED unit square inside C covers weight >= 1, and total weight < {}.\n==> {} unit squares cannot be packed into any square of side < {}/{} = {:.9},\n    i.e.  s({}) >= {:.9}", nn, nn, cert.s_num, cert.s_den, cert.s_num as f64/cert.s_den as f64, nn, cert.s_num as f64/cert.s_den as f64),
            Some(r) => if r.single {
                println!("VERIFIED: (branch k={}) every CLOSED unit square inside C centred in a corner box covers weight >= 1 + lambda, every other one >= 1, and total - lambda*k < {}.\n==> no packing of {} unit squares with exactly {} squares centred in the corner boxes [0, {}/{}]^2 (and images) fits in any square of side < {}/{} = {:.9}", r.ktot, nn, nn, r.ktot, r.r_num, r.r_den, cert.s_num, cert.s_den, cert.s_num as f64/cert.s_den as f64)
            } else {
                let kv: Vec<String> = r.k.iter().map(|k| k.to_string()).collect();
                println!("VERIFIED: (branch k={}) every CLOSED unit square inside C centred in corner box j covers weight >= 1 + lambda_j, every other one >= 1, and total - sum_j lambda_j*k_j < {}.\n==> no packing of {} unit squares with occupancy pattern ({}) of the corner boxes [0, {}/{}]^2 (and images; box j occupied iff a square is centred in it) fits in any square of side < {}/{} = {:.9}", kv.join(""), nn, nn, kv.join(","), r.r_num, r.r_den, cert.s_num, cert.s_den, cert.s_num as f64/cert.s_den as f64)
            },
        }
    } else { println!("NOT VERIFIED"); }
}
