//! `zmx2` --- an independent second exact checker for *mixed* covers (points + axis-parallel
//! segment densities) of `[0,s]^2`, format `tasks/line-cover/FORMAT.md` (mixed v1).
//!
//! Claim checked: for every closed unit square `Q subseteq [0,s]^2` (every centre, every angle),
//! `mu(Q) >= 1`, where `mu` = point masses + uniform-by-length masses on closed segments.
//!
//! Written from the statement only (no code or lemma write-up of `search/zm_mixed.py`).  Every
//! lemma used here is stated and proved in `search/ZMX2.md`; the section numbers below refer to it.
//! The existing `zmcheck` binary (`src/main.rs`) is untouched; this is a separate binary.
//!
//! Method (ZMX2.md sec 2-6): branch and bound over closed pose boxes
//! `B = [x0,x1] x [y0,y1] x [u0,u1]`, `u = tan(theta/2)`.  For each box a lower bound
//! `L(B) <= mu(Q)` for every *admissible* pose of `B` is computed as
//!   (points certainly inside `Q` for all poses of `B`)          -- exact integer test (Lemma P)
//! + (per axis-parallel grid line, or per pair of lines at distance 1, a lower bound on the
//!    captured segment mass)                                     -- Lemmas C, Z (pair), W (walls)
//! and the box is certified when `L(B) >= 1`, otherwise halved.
//!
//! Arithmetic: the point test and every mass comparison are exact integer computations.  The
//! chord-endpoint *geometry* (rational functions of u, with square roots at interior extrema)
//! is enclosed by IEEE-754 binary64 interval arithmetic with outward rounding (every result
//! widened by one ulp in the safe direction), then rounded outward to an integer grid of pitch
//! `1/(D 2^K)`; everything after that (the piecewise-linear masses and their exact minima) is
//! integer arithmetic again.  ZMX2.md sec 5 (Lemma R) explains why this is sound.
//!
//! Area densities (2026-09-30, search/ZMX2_AREA.md): polygon records that are axis-parallel
//! rectangles (Lemmas 1-3, K; covers without polygons take the old code path unchanged), the
//! theta = 0 mode `cert0` (Lemma Z0, exact), `--umin`, and the first-order bound at theta -> 0+
//! `--first-order` (Lemma U, off by default).

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, Write};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{mpsc, Arc};
use std::time::Instant;

type I = i128;

/// `--no-atoms`: keep points on grid lines as ordinary points (Lemma P only), for comparison.
static NO_ATOMS: AtomicBool = AtomicBool::new(false);
/// `--pair-points`: also make atom lines of every point abscissa/ordinate that has a partner
/// point at distance exactly 1 (sec 4.5); slower, needed for off-grid point pairs (s(32)).
static PAIR_POINTS: AtomicBool = AtomicBool::new(false);
/// `--sym-atoms`: also bound every box that the default bound does not certify with the mirrored
/// atom assignment (horizontal lines before vertical ones) and keep the larger bound (ZMX2.md
/// sec 4.9, Lemma A).  Off by default: the default behaviour (every shipped census) is unchanged.
static SYM_ATOMS: AtomicBool = AtomicBool::new(false);
/// `--mirror-only` (test mode): use the mirrored atom assignment instead of the default one.
static MIRROR_ONLY: AtomicBool = AtomicBool::new(false);
/// `--first-order`: boxes with u0 = 0 that the bound does not certify also get Lemma U (ZMX2_AREA.md sec 11).
static FIRST_ORDER: AtomicBool = AtomicBool::new(false);

fn gcd(a: I, b: I) -> I {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}
fn lcm(a: I, b: I) -> I {
    a / gcd(a, b) * b
}
fn die(msg: &str) -> ! {
    println!("ERROR: {}", msg);
    std::process::exit(2);
}

// =====================================================================================
// Outward-rounded binary64 intervals (ZMX2.md sec 5, Lemma R)
// =====================================================================================

#[derive(Clone, Copy, Debug)]
struct Iv {
    lo: f64,
    hi: f64,
}
#[inline]
fn dn(x: f64) -> f64 {
    debug_assert!(!x.is_nan());
    x.next_down()
}
#[inline]
fn up(x: f64) -> f64 {
    debug_assert!(!x.is_nan());
    x.next_up()
}
const INF: f64 = f64::INFINITY;

impl Iv {
    fn exact(x: f64) -> Iv {
        Iv { lo: x, hi: x }
    }
    /// enclosure of the rational n/d (d > 0); exact 0 stays exact.
    fn rat(n: I, d: I) -> Iv {
        assert!(d > 0);
        if n == 0 {
            return Iv::exact(0.0);
        }
        // |n|, |d| < 2^53: both convert to f64 exactly, so q is one correctly rounded quotient and
        // one ulp step each way would do (Lemma R); we widen by 4.
        assert!(n.abs() < (1 << 53) && d < (1 << 53), "rational too large for exact conversion");
        let q = (n as f64) / (d as f64);
        let (mut lo, mut hi) = (q, q);
        for _ in 0..4 {
            lo = dn(lo);
            hi = up(hi);
        }
        Iv { lo, hi }
    }
    fn add(self, o: Iv) -> Iv {
        let lo = self.lo + o.lo;
        let hi = self.hi + o.hi;
        assert!(!lo.is_nan() && !hi.is_nan(), "NaN in interval add");
        Iv { lo: dn(lo), hi: up(hi) }
    }
    fn neg(self) -> Iv {
        Iv { lo: -self.hi, hi: -self.lo }
    }
    fn sub(self, o: Iv) -> Iv {
        self.add(o.neg())
    }
    /// finite operands only
    fn mul(self, o: Iv) -> Iv {
        assert!(self.lo.is_finite() && self.hi.is_finite() && o.lo.is_finite() && o.hi.is_finite());
        if self.is_zero() || o.is_zero() {
            return Iv::exact(0.0); // an exact zero factor gives an exact zero
        }
        let p = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi];
        let mut lo = p[0];
        let mut hi = p[0];
        for &v in &p[1..] {
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
        Iv { lo: dn(lo), hi: up(hi) }
    }
    fn scale(self, k: f64) -> Iv {
        self.mul(Iv::exact(k))
    }
    fn recip_pos(self) -> Iv {
        assert!(self.lo > 0.0 && self.hi.is_finite());
        Iv { lo: dn(1.0 / self.hi), hi: up(1.0 / self.lo) }
    }
    fn div_pos(self, o: Iv) -> Iv {
        self.mul(o.recip_pos())
    }
    fn sqrt(self) -> Iv {
        assert!(self.hi >= 0.0);
        let lo = if self.lo <= 0.0 { 0.0 } else { dn(self.lo.sqrt()).max(0.0) };
        Iv { lo, hi: up(self.hi.sqrt()) }
    }
    fn hull(self, o: Iv) -> Iv {
        Iv { lo: self.lo.min(o.lo), hi: self.hi.max(o.hi) }
    }
    fn is_zero(self) -> bool {
        self.lo == 0.0 && self.hi == 0.0
    }
}

// ---------------------------------------------------------------- the chord-end functions
//
// Vertical line x = l, pose (c, u), d = c_x - l, C = (1-u^2)/(1+u^2), S = 2u/(1+u^2).
// With y measured from c_y the chord of the line in Q is I1 cap I2 (ZMX2.md sec 3, Lemma C):
//   I1 = [f1, f2],  f1 = (dC - 1/2)/S = (d-1/2)/(2u) - (d+1/2) u/2,
//                   f2 = (dC + 1/2)/S = (d+1/2)/(2u) + (1/2-d) u/2,
//   I2 = [g1, g2],  g1 = (-dS - 1/2)/C = -d T - R/2,  g2 = -d T + R/2,
//   T = 2u/(1-u^2) = tan(theta),  R = (1+u^2)/(1-u^2) = sec(theta).

/// h(u) = a/u + b u over the interval u (u.lo > 0), naive enclosure.
fn h_at(a: Iv, b: Iv, u: Iv) -> Iv {
    let t1 = if a.is_zero() { Iv::exact(0.0) } else { a.mul(u.recip_pos()) };
    t1.add(b.mul(u))
}

/// Enclosure of { a/u + b u : u in [u0,u1] } (a, b thin intervals, 0 <= u0 < u1), in the
/// extended reals; at u0 = 0 the value is the limit u -> 0+ (ZMX2.md Lemma H).
fn h_range(a: Iv, b: Iv, u0: Iv, u1: Iv, u0zero: bool) -> Iv {
    let full = Iv { lo: -INF, hi: INF };
    let v0 = if u0zero {
        if a.is_zero() {
            Iv::exact(0.0)
        } else if a.lo > 0.0 {
            Iv::exact(INF)
        } else if a.hi < 0.0 {
            Iv::exact(-INF)
        } else {
            return full;
        }
    } else {
        h_at(a, b, u0)
    };
    let v1 = h_at(a, b, u1);
    let mut r = v0.hull(v1);
    let pos = a.lo > 0.0 && b.lo > 0.0;
    let neg = a.hi < 0.0 && b.hi < 0.0;
    let mono = (a.lo >= 0.0 && b.hi <= 0.0) || (a.hi <= 0.0 && b.lo >= 0.0);
    if pos || neg {
        // one interior critical point u* = sqrt(a/b)
        let q = if pos { a.div_pos(b) } else { a.neg().div_pos(b.neg()) };
        let us = q.sqrt();
        let lo = us.lo.max(u0.lo);
        let hi = us.hi.min(u1.hi);
        if lo <= hi {
            if lo > 0.0 {
                r = r.hull(h_at(a, b, Iv { lo, hi }));
            } else {
                return full;
            }
        }
    } else if !mono {
        if u0zero || u0.lo <= 0.0 {
            return full;
        }
        r = r.hull(h_at(a, b, Iv { lo: u0.lo, hi: u1.hi }));
    }
    r
}

/// T(u) = 2u/(1-u^2) at a point (enclosure), 0 <= x < 1
fn t_pt(x: f64) -> Iv {
    let xi = Iv::exact(x);
    let den = Iv::exact(1.0).sub(xi.mul(xi));
    Iv::exact(2.0 * x).div_pos(den)
}
/// R(u) = (1+u^2)/(1-u^2) at a point (enclosure)
fn r_pt(x: f64) -> Iv {
    let xi = Iv::exact(x);
    let x2 = xi.mul(xi);
    Iv::exact(1.0).add(x2).div_pos(Iv::exact(1.0).sub(x2))
}
/// g(u) = -d T(u) + sg R(u)/2 over an interval u (T, R increasing on [0,1))
fn g_at(d: Iv, u: Iv, sg: f64) -> Iv {
    let t = Iv { lo: t_pt(u.lo).lo, hi: t_pt(u.hi).hi };
    let r = Iv { lo: r_pt(u.lo).lo, hi: r_pt(u.hi).hi };
    d.neg().mul(t).add(r.scale(0.5 * sg))
}
/// Enclosure of { g(u) : u in [u0,u1] }, 0 <= u0 < u1 < 1 (ZMX2.md Lemma G).
fn g_range(d: Iv, sg: f64, u0: Iv, u1: Iv) -> Iv {
    let v0 = g_at(d, u0, sg);
    let v1 = g_at(d, u1, sg);
    let mut r = v0.hull(v1);
    let crit = (sg > 0.0 && d.lo > 0.0) || (sg < 0.0 && d.hi < 0.0);
    if crit {
        let dd = if d.lo > 0.0 { d } else { d.neg() };
        let disc = Iv::exact(1.0).sub(dd.mul(dd).scale(4.0));
        if disc.hi >= 0.0 {
            let disc = Iv { lo: disc.lo.max(0.0), hi: disc.hi };
            let us = dd.scale(2.0).div_pos(Iv::exact(1.0).add(disc.sqrt()));
            let lo = us.lo.max(u0.lo);
            let hi = us.hi.min(u1.hi);
            if lo <= hi {
                r = r.hull(g_at(d, Iv { lo, hi }, sg));
            }
        }
    } else if d.lo < 0.0 && d.hi > 0.0 {
        r = r.hull(g_at(d, Iv { lo: u0.lo, hi: u1.hi }, sg));
    }
    r
}

// =====================================================================================
// The cover
// =====================================================================================

#[derive(Clone)]
struct Line {
    pos: i64,       // line position in 1/D units (x for vertical; -y for horizontal, see Lemma T)
    bps: Vec<i64>,  // breakpoints along the line, 1/D units, increasing
    cum: Vec<I>,    // F at bps, units 1/(W Lc)
    dens: Vec<I>,   // density on [bps[k], bps[k+1]], units 1/(W Lc) per 1/D of length
    atoms: Vec<(i64, I)>, // point masses on the line: (coordinate along, weight in units 1/(W Lc))
}

/// An area density: mass w/W spread uniformly over the closed axis-parallel rectangle
/// [x0,x1] x [y0,y1] (1/D units, x0 < x1, y0 < y1).  ZMX2_AREA.md.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Rect {
    x0: i64,
    x1: i64,
    y0: i64,
    y1: i64,
    w: I,
}

#[derive(Clone)]
struct Cover {
    s_num: I,
    s_den: I,
    d: I,        // coordinate denominator D
    w: I,        // weight denominator W
    sx: I,       // container side in 1/D units
    lc: I,       // lcm of segment lengths (1/D units)
    pts: Vec<(i64, i64, I)>, // aggregated points, weight numerators (units 1/W), w > 0 only
    vl: Vec<Line>, // vertical lines, sorted by pos
    hl: Vec<Line>, // horizontal lines in the rotated frame: pos = -y, coordinate along = x
    alt: Option<(Vec<Line>, Vec<Line>)>, // --sym-atoms: (vl, hl) under the mirrored assignment (sec 4.9)
    buckets: Vec<Vec<u32>>, // point index grid, cell 1/10 (in [0,s])
    nb: usize,
    raw_pts: Vec<(i64, i64, I)>,
    raw_segs: Vec<(i64, i64, i64, i64, I)>,
    rects: Vec<Rect>, // area densities on axis-parallel rectangles (ZMX2_AREA.md); empty = old behaviour
    total: (I, I), // total mass as a fraction (num, den)
    hash: u64,
}

fn tokens(text: &str) -> Vec<String> {
    let mut v = Vec::new();
    for line in text.lines() {
        let l = match line.find('#') {
            Some(p) => &line[..p],
            None => line,
        };
        for t in l.split_whitespace() {
            v.push(t.to_string());
        }
    }
    v
}

fn fnv(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn parse_cover(path: &str) -> Cover {
    let bytes = std::fs::read(path).unwrap_or_else(|e| die(&format!("cannot read {}: {}", path, e)));
    let text = String::from_utf8_lossy(&bytes).to_string();
    let tk = tokens(&text);
    let first = tk.first().cloned().unwrap_or_default();
    let mixed = first == "mixed";
    let mut it2 = tk.iter().skip(if mixed { 1 } else { 0 });
    let mut nx = |what: &str| -> I {
        let t = it2.next().unwrap_or_else(|| die(&format!("unexpected end of file reading {}", what)));
        t.parse::<I>().unwrap_or_else(|_| die(&format!("bad integer '{}' reading {}", t, what)))
    };
    if mixed {
        let ver = nx("version");
        if ver != 1 {
            die("only mixed format version 1 is supported");
        }
    }
    let s_num = nx("s_num");
    let s_den = nx("s_den");
    let d = nx("D");
    let w = nx("W");
    if s_num <= 0 || s_den <= 0 || d <= 0 || w <= 0 {
        die("header values must be positive");
    }
    // Input bounds (ZMX2_AUDIT.md F1): release builds do not check i128 overflow, so every input
    // integer is bounded here, before any arithmetic, such that all later sizes of ZMX2.md sec 2/5
    // hold for every accepted file (F * 2^30 * Lc < 2^110, grid integers and casts small).
    const HDR: I = 1 << 40;
    if s_num >= HDR || s_den >= HDR || w >= (1 << 50) {
        die("header value too large (s_num, s_den < 2^40, W < 2^50)");
    }
    if d >= (1 << 20) {
        die("D too large for this checker (limit 2^20)");
    }
    if (s_num * d) % s_den != 0 {
        die("s_den must divide s_num*D");
    }
    let sx = s_num * d / s_den;
    if sx >= (1 << 27) {
        die("s*D too large for this checker (limit 2^27)");
    }
    let inrange = |v: I| v >= 0 && v <= sx;
    let np = nx("np");
    if np < 0 {
        die("negative count");
    }
    let mut raw_pts = Vec::new();
    for _ in 0..np {
        let x = nx("point X");
        let y = nx("point Y");
        let wt = nx("point w");
        if !inrange(x) || !inrange(y) {
            die("point outside the container");
        }
        if wt < 0 {
            die("negative weight");
        }
        if wt >= (1 << 50) {
            die("weight too large (limit 2^50)");
        }
        raw_pts.push((x as i64, y as i64, wt));
    }
    let mut raw_segs = Vec::new();
    let mut raw_rects: Vec<Rect> = Vec::new();
    if mixed {
        let ns = nx("ns");
        if ns < 0 {
            die("negative count");
        }
        for _ in 0..ns {
            let x0 = nx("seg X0");
            let y0 = nx("seg Y0");
            let x1 = nx("seg X1");
            let y1 = nx("seg Y1");
            let wt = nx("seg w");
            if !inrange(x0) || !inrange(y0) || !inrange(x1) || !inrange(y1) {
                die("segment outside the container");
            }
            if wt < 0 {
                die("negative weight");
            }
            if wt >= (1 << 50) {
                die("weight too large (limit 2^50)");
            }
            if x0 == x1 && y0 == y1 {
                die("degenerate (zero-length) segment");
            }
            if x0 != x1 && y0 != y1 {
                die("unsupported: zmx2 handles axis-parallel segments only");
            }
            raw_segs.push((x0 as i64, y0 as i64, x1 as i64, y1 as i64, wt));
        }
        let npg = nx("npg");
        if npg < 0 {
            die("negative count");
        }
        // Polygons (ZMX2_AREA.md): only axis-parallel rectangles, given counter-clockwise.
        for _ in 0..npg {
            let k = nx("polygon k");
            let wt = nx("polygon w");
            if k != 4 {
                die("unsupported: zmx2 handles polygons that are axis-parallel rectangles only (k must be 4)");
            }
            let mut vs: Vec<(I, I)> = Vec::new();
            for _ in 0..4 {
                let x = nx("polygon X");
                let y = nx("polygon Y");
                if !inrange(x) || !inrange(y) {
                    die("polygon outside the container");
                }
                vs.push((x, y));
            }
            if wt < 0 {
                die("negative weight");
            }
            if wt >= (1 << 50) {
                die("weight too large (limit 2^50)");
            }
            let (xa, xb) = (vs.iter().map(|v| v.0).min().unwrap(), vs.iter().map(|v| v.0).max().unwrap());
            let (ya, yb) = (vs.iter().map(|v| v.1).min().unwrap(), vs.iter().map(|v| v.1).max().unwrap());
            if xa == xb || ya == yb {
                die("degenerate (zero-area) polygon");
            }
            // counter-clockwise corner cycle of [xa,xb] x [ya,yb]; the file's list must be a rotation of it
            let cyc = [(xa, ya), (xb, ya), (xb, yb), (xa, yb)];
            let ok = (0..4).any(|r| (0..4).all(|i| vs[i] == cyc[(i + r) % 4]));
            if !ok {
                die("unsupported: zmx2 handles polygons that are axis-parallel rectangles only (counter-clockwise)");
            }
            raw_rects.push(Rect { x0: xa as i64, x1: xb as i64, y0: ya as i64, y1: yb as i64, w: wt });
        }
    }
    if it2.next().is_some() {
        die("trailing tokens after the declared pieces");
    }
    // each weight < 2^50, and a file has far fewer than 2^70 pieces, so this sum cannot overflow
    let tsum: I =
        raw_pts.iter().map(|p| p.2).chain(raw_segs.iter().map(|s| s.4)).chain(raw_rects.iter().map(|r| r.w)).sum();
    if tsum >= (1 << 56) {
        die("total weight too large (limit 2^56)");
    }
    build_cover(s_num, s_den, d, w, sx, raw_pts, raw_segs, raw_rects, fnv(&bytes))
}

fn build_cover(
    s_num: I,
    s_den: I,
    d: I,
    w: I,
    sx: I,
    raw_pts: Vec<(i64, i64, I)>,
    raw_segs: Vec<(i64, i64, i64, i64, I)>,
    raw_rects: Vec<Rect>,
    hash: u64,
) -> Cover {
    // aggregate points
    let mut pm: BTreeMap<(i64, i64), I> = BTreeMap::new();
    for &(x, y, wt) in &raw_pts {
        *pm.entry((x, y)).or_insert(0) += wt;
    }
    // line positions: interior unit-grid lines and every segment's line (sec 4.5: atoms)
    let mut vpos: HashSet<i64> = HashSet::new();
    let mut hpos: HashSet<i64> = HashSet::new();
    let mut k = d;
    while k < sx {
        vpos.insert(k as i64);
        hpos.insert(k as i64);
        k += d;
    }
    for &(x0, y0, x1, _y1, wt) in &raw_segs {
        if wt == 0 {
            continue;
        }
        if x0 == x1 {
            vpos.insert(x0);
        } else {
            hpos.insert(y0);
        }
    }
    // ... and (--pair-points) every point abscissa X (ordinate Y) with a partner at X +- 1 (Y +- 1):
    // two parallel lines at distance 1 are what the pair lemma couples (Lemma Z).  A point is
    // assigned to the first that applies: vertical grid/segment line, horizontal grid/segment line,
    // vertical partner line, horizontal partner line; otherwise it stays an ordinary point.
    let mut vpart: HashSet<i64> = HashSet::new();
    let mut hpart: HashSet<i64> = HashSet::new();
    if PAIR_POINTS.load(Ordering::Relaxed) {
        let xs: HashSet<i64> = pm.keys().map(|k| k.0).chain(vpos.iter().copied()).collect();
        let ys: HashSet<i64> = pm.keys().map(|k| k.1).chain(hpos.iter().copied()).collect();
        let di = d as i64;
        for &x in &xs {
            if !vpos.contains(&x) && (xs.contains(&(x + di)) || xs.contains(&(x - di))) {
                vpart.insert(x);
            }
        }
        for &y in &ys {
            if !hpos.contains(&y) && (ys.contains(&(y + di)) || ys.contains(&(y - di))) {
                hpart.insert(y);
            }
        }
    }
    let no_atoms = NO_ATOMS.load(Ordering::Relaxed);
    // segments -> lines
    let mut lc: I = 1;
    for &(x0, y0, x1, y1, wt) in &raw_segs {
        if wt == 0 {
            continue;
        }
        let len = ((x1 - x0).abs() + (y1 - y0).abs()) as I;
        lc = lcm(lc, len);
        if lc > (1 << 24) {
            die("unsupported: lcm of segment lengths exceeds 2^24");
        }
    }
    let sets = AtomSets { vpos: &vpos, hpos: &hpos, vpart: &vpart, hpart: &hpart, no_atoms };
    // (--mirror-only, a test mode: the mirrored assignment alone, sec 4.9)
    let (vl, hl, pts) = build_lines(&pm, &raw_segs, &sets, lc, MIRROR_ONLY.load(Ordering::Relaxed));
    // --sym-atoms (ZMX2.md sec 4.9, Lemma A): the mirrored assignment rule (horizontal before
    // vertical at each rank).  It changes only which line a point is an atom of, never whether it
    // is an atom, so the ordinary points are the same (asserted).  Kept only if it differs.
    let alt = if SYM_ATOMS.load(Ordering::Relaxed) && !MIRROR_ONLY.load(Ordering::Relaxed) {
        let (vl2, hl2, pts2) = build_lines(&pm, &raw_segs, &sets, lc, true);
        assert!(pts2 == pts, "mirrored atom assignment changed the ordinary points");
        let same = |x: &Vec<Line>, y: &Vec<Line>| {
            x.len() == y.len() && x.iter().zip(y.iter()).all(|(a, b)| a.pos == b.pos && a.atoms == b.atoms)
        };
        if same(&vl2, &vl) && same(&hl2, &hl) {
            None
        } else {
            Some((vl2, hl2))
        }
    } else {
        None
    };
    // buckets of 1/10
    let nb = ((sx * 10 + d - 1) / d) as usize + 1;
    let mut buckets = vec![Vec::new(); nb * nb];
    for (i, &(x, y, _)) in pts.iter().enumerate() {
        let bx = ((x as I) * 10 / d) as usize;
        let by = ((y as I) * 10 / d) as usize;
        buckets[bx.min(nb - 1) * nb + by.min(nb - 1)].push(i as u32);
    }
    let tnum: I = raw_pts.iter().map(|p| p.2).sum::<I>()
        + raw_segs.iter().map(|s| s.4).sum::<I>()
        + raw_rects.iter().map(|r| r.w).sum::<I>();
    let mut rects: Vec<Rect> = raw_rects.iter().copied().filter(|r| r.w > 0).collect();
    rects.sort();
    let g = gcd(tnum, w).max(1);
    Cover {
        s_num,
        s_den,
        d,
        w,
        sx,
        lc,
        pts,
        vl,
        hl,
        alt,
        buckets,
        nb,
        raw_pts,
        raw_segs,
        rects,
        total: (tnum / g, w / g),
        hash,
    }
}

/// The line/atom position sets of a cover (sec 4.5).
struct AtomSets<'a> {
    vpos: &'a HashSet<i64>,
    hpos: &'a HashSet<i64>,
    vpart: &'a HashSet<i64>,
    hpart: &'a HashSet<i64>,
    no_atoms: bool,
}

/// Build the vertical and horizontal lines (densities + atoms) and the ordinary points.  A point
/// (x, y) is an atom of the first line in its rank list that exists, otherwise an ordinary point:
///   hfirst = false (default):  vertical grid/segment x, horizontal grid/segment y, vertical
///                              partner x, horizontal partner y;
///   hfirst = true (mirrored):  horizontal grid/segment y, vertical grid/segment x, horizontal
///                              partner y, vertical partner x.
/// Either rule makes every point an atom of at most one line, which is all Lemma DP needs, and
/// both make the same points atoms (sec 4.9).
fn build_lines(
    pm: &BTreeMap<(i64, i64), I>,
    raw_segs: &[(i64, i64, i64, i64, I)],
    st: &AtomSets,
    lc: I,
    hfirst: bool,
) -> (Vec<Line>, Vec<Line>, Vec<(i64, i64, I)>) {
    let mut pts: Vec<(i64, i64, I)> = Vec::new();
    // (orientation, pos) -> segments (a, b, w) and atoms (t, w) along the line; vertical: pos = x,
    // along = y; horizontal (rotated frame, Lemma T): pos = -y, along = x.
    let mut groups: BTreeMap<(u8, i64), (Vec<(i64, i64, I)>, Vec<(i64, I)>)> = BTreeMap::new();
    for (&(x, y), &wt) in pm {
        if wt <= 0 {
            continue;
        }
        let (vg, hg) = (st.vpos.contains(&x), st.hpos.contains(&y));
        let (vp, hp) = (st.vpart.contains(&x), st.hpart.contains(&y));
        // 0 = atom of the vertical line x, 1 = atom of the horizontal line y, 2 = ordinary point
        let o = if st.no_atoms {
            2
        } else if !hfirst {
            if vg {
                0
            } else if hg {
                1
            } else if vp {
                0
            } else if hp {
                1
            } else {
                2
            }
        } else if hg {
            1
        } else if vg {
            0
        } else if hp {
            1
        } else if vp {
            0
        } else {
            2
        };
        match o {
            0 => groups.entry((0, x)).or_default().1.push((y, wt * lc)),
            1 => groups.entry((1, -y)).or_default().1.push((x, wt * lc)),
            _ => pts.push((x, y, wt)),
        }
    }
    for &(x0, y0, x1, y1, wt) in raw_segs {
        if wt == 0 {
            continue;
        }
        if x0 == x1 {
            groups.entry((0, x0)).or_default().0.push((y0.min(y1), y0.max(y1), wt));
        } else {
            groups.entry((1, -y0)).or_default().0.push((x0.min(x1), x0.max(x1), wt));
        }
    }
    let mut vl = Vec::new();
    let mut hl = Vec::new();
    for ((o, pos), (segs, mut atoms)) in groups {
        let mut bs: Vec<i64> = segs.iter().flat_map(|&(a, b, _)| [a, b]).collect();
        bs.sort();
        bs.dedup();
        let nd = if bs.is_empty() { 0 } else { bs.len() - 1 };
        let mut dens = vec![0 as I; nd];
        for &(a, b, wt) in &segs {
            let len = (b - a) as I;
            let per = wt * (lc / len);
            let ia = bs.binary_search(&a).unwrap();
            let ib = bs.binary_search(&b).unwrap();
            for k in ia..ib {
                dens[k] += per;
            }
        }
        let mut cum = vec![0 as I; bs.len()];
        for k in 0..dens.len() {
            cum[k + 1] = cum[k] + dens[k] * ((bs[k + 1] - bs[k]) as I);
        }
        let tot: I = segs.iter().map(|&(_, _, wt)| wt * lc).sum();
        if !bs.is_empty() {
            assert_eq!(cum[cum.len() - 1], tot, "line cumulative mass mismatch");
        }
        atoms.sort();
        let line = Line { pos, bps: bs, cum, dens, atoms };
        if o == 0 {
            vl.push(line);
        } else {
            hl.push(line);
        }
    }
    vl.sort_by_key(|l| l.pos);
    hl.sort_by_key(|l| l.pos);
    (vl, hl, pts)
}

/// The cover reflected by y -> s - y (used by the unreduced sweep for theta in [45,90] deg).
fn reflect_y(c: &Cover) -> Cover {
    let sx = c.sx as i64;
    let p: Vec<_> = c.raw_pts.iter().map(|&(x, y, w)| (x, sx - y, w)).collect();
    let s: Vec<_> = c.raw_segs.iter().map(|&(x0, y0, x1, y1, w)| (x0, sx - y0, x1, sx - y1, w)).collect();
    let r: Vec<Rect> = c.rects.iter().map(|r| Rect { x0: r.x0, x1: r.x1, y0: sx - r.y1, y1: sx - r.y0, w: r.w }).collect();
    build_cover(c.s_num, c.s_den, c.d, c.w, c.sx, p, s, r, c.hash ^ 0x5bd1e995)
}

// ------------------------------------------------------------- D4 invariance of the measure

/// canonical form of a line density: maximal intervals of constant nonzero density
fn canon(l: &Line) -> Vec<(i64, i64, I)> {
    let mut v: Vec<(i64, i64, I)> = Vec::new();
    for k in 0..l.dens.len() {
        if l.dens[k] == 0 {
            continue;
        }
        if let Some(last) = v.last_mut() {
            if last.1 == l.bps[k] && last.2 == l.dens[k] {
                last.1 = l.bps[k + 1];
                continue;
            }
        }
        v.push((l.bps[k], l.bps[k + 1], l.dens[k]));
    }
    v
}

/// Exact check that mu is invariant under x -> s-x and x <-> y (ZMX2.md sec 7, Lemma S).
fn check_d4(c: &Cover) -> Result<(), String> {
    let sx = c.sx as i64;
    // all points (on-line atoms included), aggregated, zero weights dropped
    let mut pm: HashMap<(i64, i64), I> = HashMap::new();
    for &(x, y, w) in &c.raw_pts {
        *pm.entry((x, y)).or_insert(0) += w;
    }
    pm.retain(|_, w| *w != 0);
    for (&(x, y), &w) in &pm {
        for (gx, gy) in [(sx - x, y), (y, x)] {
            if pm.get(&(gx, gy)).copied().unwrap_or(0) != w {
                return Err(format!("point ({},{}) weight {} has no D4 image at ({},{})", x, y, w, gx, gy));
            }
        }
    }
    // vertical lines keyed by x, horizontal by y (hl.pos = -y); densities along the line
    let mut v: BTreeMap<i64, Vec<(i64, i64, I)>> = BTreeMap::new();
    let mut h: BTreeMap<i64, Vec<(i64, i64, I)>> = BTreeMap::new();
    for l in &c.vl {
        let cf = canon(l);
        if !cf.is_empty() {
            v.insert(l.pos, cf);
        }
    }
    for l in &c.hl {
        let cf = canon(l);
        if !cf.is_empty() {
            h.insert(-l.pos, cf);
        }
    }
    let refl = |f: &Vec<(i64, i64, I)>| -> Vec<(i64, i64, I)> {
        let mut r: Vec<_> = f.iter().map(|&(a, b, d)| (sx - b, sx - a, d)).collect();
        r.sort();
        r
    };
    // x -> s-x: vertical line x=a -> x = s-a, same density in y; horizontal line y=b -> itself,
    // density reflected in x.
    for (&a, f) in &v {
        if v.get(&(sx - a)) != Some(f) {
            return Err(format!("vertical line x={} density not matched by x={}", a, sx - a));
        }
    }
    for (&b, f) in &h {
        if &refl(f) != f {
            return Err(format!("horizontal line y={} density not symmetric under x->s-x", b));
        }
    }
    // x <-> y: vertical x=a <-> horizontal y=a, same density
    if v != h {
        return Err("vertical and horizontal line densities differ under x<->y".into());
    }
    // area densities (ZMX2_AREA.md sec 7): the multiset of weighted rectangles is invariant under both
    // generators (sufficient for the area measure to be invariant)
    let rs = c.rects.clone(); // sorted, w > 0
    let mut r1: Vec<Rect> = rs.iter().map(|r| Rect { x0: sx - r.x1, x1: sx - r.x0, ..*r }).collect();
    let mut r2: Vec<Rect> = rs.iter().map(|r| Rect { x0: r.y0, x1: r.y1, y0: r.x0, y1: r.x1, w: r.w }).collect();
    r1.sort();
    r2.sort();
    if r1 != rs || r2 != rs {
        return Err("area densities (rectangles) not invariant under x->s-x and x<->y".into());
    }
    Ok(())
}

// =====================================================================================
// Pose boxes
// =====================================================================================

#[derive(Clone, Copy, Debug)]
struct PBox {
    xn: [I; 2],
    yn: [I; 2],
    cl: u32, // centre denominator 10 * 2^cl
    un: [I; 2],
    ul: u32, // u denominator 8 * 2^ul
}
impl PBox {
    fn cd(&self) -> I {
        10 << self.cl
    }
    fn ud(&self) -> I {
        8 << self.ul
    }
    fn fstr(&self) -> String {
        let cd = self.cd() as f64;
        let ud = self.ud() as f64;
        format!(
            "x[{:.7},{:.7}] y[{:.7},{:.7}] u[{:.8},{:.8}] (th {:.4}-{:.4} deg)",
            self.xn[0] as f64 / cd,
            self.xn[1] as f64 / cd,
            self.yn[0] as f64 / cd,
            self.yn[1] as f64 / cd,
            self.un[0] as f64 / ud,
            self.un[1] as f64 / ud,
            2.0 * (self.un[0] as f64 / ud).atan().to_degrees(),
            2.0 * (self.un[1] as f64 / ud).atan().to_degrees()
        )
    }
    fn exact_str(&self) -> String {
        format!(
            "{}/{},{}/{},{}/{},{}/{},{}/{},{}/{}",
            self.xn[0],
            self.cd(),
            self.xn[1],
            self.cd(),
            self.yn[0],
            self.cd(),
            self.yn[1],
            self.cd(),
            self.un[0],
            self.ud(),
            self.un[1],
            self.ud()
        )
    }
}

// =====================================================================================
// The bound (ZMX2.md sec 2-6)
// =====================================================================================

const K: u32 = 30; // grid pitch 1/(D 2^K)

struct Ctx<'a> {
    cv: &'a Cover,
    g: I,      // D * 2^K
    gf: f64,   // same as f64 (exact)
    big: I,    // clamp value in grid units
    bigv: f64, // clamp value in length units
    target: I, // W * Lc * 2^K : the mass 1 in bound units
    ptw: I,    // Lc * 2^K : point weight numerator -> bound units
    // area densities (ZMX2_AREA.md sec 5): per rectangle r with weight w_r and area Pn_r (1/D^2 units)
    rect_tr: Vec<I>,  // floor(w_r Lc 2^K D^2 / Pn_r): the mass of area 1 of r, in bound units (rounded down)
    rect_wlc: Vec<I>, // w_r Lc
    rect_den: Vec<I>, // Pn_r 2^K: area phi*l (grid units phin, ln) has mass w_r Lc phin ln / (Pn_r 2^K)
}

impl<'a> Ctx<'a> {
    fn new(cv: &'a Cover) -> Ctx<'a> {
        let g = cv.d << K;
        let sf = (cv.sx as f64) / (cv.d as f64);
        let bigv = (4.0 * sf + 10.0).ceil();
        Ctx {
            cv,
            g,
            gf: g as f64,
            big: (bigv as I) * g,
            bigv,
            target: cv.w * cv.lc << K,
            ptw: cv.lc << K,
            rect_tr: cv
                .rects
                .iter()
                .map(|r| {
                    let pn = ((r.x1 - r.x0) as I) * ((r.y1 - r.y0) as I);
                    let num = r
                        .w
                        .checked_mul(cv.lc)
                        .and_then(|v| v.checked_mul(1 << K))
                        .and_then(|v| v.checked_mul(cv.d * cv.d))
                        .unwrap_or_else(|| die("area density: w Lc 2^K D^2 overflows i128"));
                    num / pn
                })
                .collect(),
            rect_wlc: cv
                .rects
                .iter()
                .map(|r| {
                    let v = r.w.checked_mul(cv.lc).unwrap_or_else(|| die("area density: w Lc overflows"));
                    // phin (grid units, < 2^(K+1) D) times this must stay far inside i128 (sec 5 of ZMX2_AREA.md)
                    if v >= (1 << 70) {
                        die("area density: w Lc too large (limit 2^70)");
                    }
                    v
                })
                .collect(),
            rect_den: cv.rects.iter().map(|r| (((r.x1 - r.x0) as I) * ((r.y1 - r.y0) as I)) << K).collect(),
        }
    }
    /// a grid integer <= v*G (up to the harmless clamp of Lemma R(iii))
    fn glo(&self, v: f64) -> I {
        assert!(!v.is_nan(), "NaN endpoint");
        if v <= -self.bigv {
            -self.big
        } else if v >= self.bigv {
            self.big
        } else {
            (dn(v * self.gf)).floor() as I
        }
    }
    /// a grid integer >= v*G (up to the clamp)
    fn ghi(&self, v: f64) -> I {
        assert!(!v.is_nan(), "NaN endpoint");
        if v <= -self.bigv {
            -self.big
        } else if v >= self.bigv {
            self.big
        } else {
            (up(v * self.gf)).ceil() as I
        }
    }
}

/// F(y) * 2^K for a grid integer y (units 1/(W Lc 2^K)); exact.
#[inline]
fn fval(l: &Line, y: I) -> I {
    let n = l.bps.len();
    if n == 0 {
        return 0;
    }
    let b0 = (l.bps[0] as I) << K;
    if y <= b0 {
        return 0;
    }
    let bl = (l.bps[n - 1] as I) << K;
    if y >= bl {
        return l.cum[n - 1] << K;
    }
    // largest k with bps[k]<<K <= y
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if ((l.bps[mid] as I) << K) <= y {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (l.cum[lo] << K) + l.dens[lo] * (y - ((l.bps[lo] as I) << K))
}

#[derive(Clone, Copy, Debug)]
struct Ends {
    zlo: I,  // zeta = L1 (lower end of I1), range
    zhi: I,  // also the upper bound L1_hi
    h1lo: I, // lower bound of H1 (upper end of I1)
    l2hi: I, // upper bound of L2
    h2lo: I, // lower bound of H2
}

/// Box data common to all lines, in the (possibly rotated) frame.
struct FrameBox {
    xn: [I; 2],
    xd: I,
    y0: Iv,
    y1: Iv,
    u0: Iv,
    u1: Iv,
    u0zero: bool,
    u0g: I,   // floor(u0 * G)
    q_lo: f64, // lower bound of q(u) = (1 - u^2 + 2u^3)/(2(1+u^2)) on [u0,u1] (Lemma W)
    wall_l: i64, // wall positions in this frame (1/D units)
    wall_r: i64,
    pb: PBox, // the box in this frame (exact data for the atom tests)
}

fn line_ends(ctx: &Ctx, fb: &FrameBox, pos: i64) -> Ends {
    let dd = ctx.cv.d;
    // d = x - l, as exact rationals over xd*D, at the two x endpoints
    let den = fb.xd * dd;
    let mut f1 = Iv { lo: INF, hi: -INF };
    let mut f2 = f1;
    let mut g1 = f1;
    let mut g2 = f1;
    for e in 0..2 {
        let dnum = fb.xn[e] * dd - (pos as I) * fb.xd; // d = dnum/den
        let dv = Iv::rat(dnum, den);
        let mm = Iv::rat(2 * dnum - den, 2 * den); // d - 1/2
        let pm = Iv::rat(2 * dnum + den, 2 * den); // d + 1/2
        let r1 = h_range(mm.scale(0.5), pm.scale(-0.5), fb.u0, fb.u1, fb.u0zero);
        let r2 = h_range(pm.scale(0.5), mm.scale(-0.5), fb.u0, fb.u1, fb.u0zero);
        f1 = f1.hull(r1);
        f2 = f2.hull(r2);
        g1 = g1.hull(g_range(dv, -1.0, fb.u0, fb.u1));
        g2 = g2.hull(g_range(dv, 1.0, fb.u0, fb.u1));
    }
    let zlo = ctx.glo(add_lo(fb.y0.lo, f1.lo));
    let mut zhi = ctx.ghi(add_hi(fb.y1.hi, f1.hi));
    let mut h1lo = ctx.glo(add_lo(fb.y0.lo, f2.lo));
    let l2hi = ctx.ghi(add_hi(fb.y1.hi, g1.hi));
    let h2lo = ctx.glo(add_lo(fb.y0.lo, g2.lo));
    // Lemma W (walls): line at distance exactly 1 from a wall, admissible poses only
    if pos - fb.wall_l == dd as i64 {
        h1lo = h1lo.max(ctx.glo(add_lo(fb.y0.lo, fb.q_lo)));
    }
    if fb.wall_r - pos == dd as i64 {
        zhi = zhi.min(ctx.ghi(add_hi(fb.y1.hi, -fb.q_lo)));
    }
    Ends { zlo, zhi, h1lo, l2hi, h2lo }
}
#[inline]
fn add_lo(a: f64, b: f64) -> f64 {
    let s = a + b;
    assert!(!s.is_nan());
    dn(s)
}
#[inline]
fn add_hi(a: f64, b: f64) -> f64 {
    let s = a + b;
    assert!(!s.is_nan());
    up(s)
}

/// An atom of a line in the current box: grid position along the line, weight (bound units),
/// and which of its four conditions hold for every (admissible) pose of the box (Lemma P per
/// condition, Lemma W at walls): [X <= 1/2 (H1), X >= -1/2 (L1), Y <= 1/2 (H2), Y >= -1/2 (L2)].
#[derive(Clone, Copy, Debug)]
struct Atom {
    y: I,
    w: I,
    c: [bool; 4],
}

fn single_bound(l: &Line, e: &Ends, at: &[Atom]) -> I {
    let lo = e.zhi.max(e.l2hi);
    let hi = e.h1lo.min(e.h2lo);
    let cont = if hi > lo { fval(l, hi) - fval(l, lo) } else { 0 };
    cont + at.iter().filter(|a| a.c.iter().all(|&c| c)).map(|a| a.w).sum::<I>()
}

/// Lemma Z: lines a (at l) and b (at l+1) share zeta = L1_a; H1_b = zeta + u.
/// Lemma Z: lines a (at l) and b (at l+1) share zeta = L1_a; H1_b = zeta + u >= zeta + u0.
/// Returns a lower bound of inf over zeta in [zlo, zhi] of the pair mass (sec 4.3, 4.5).
fn pair_bound(la: &Line, ea: &Ends, ata: &[Atom], lb: &Line, eb: &Ends, atb: &[Atom], u0g: I) -> I {
    let l_a = ea.l2hi;
    let h_a = ea.h1lo.min(ea.h2lo);
    let l_b = eb.zhi.max(eb.l2hi);
    let h2b = eb.h2lo;
    let h1b = eb.h1lo;
    let (zlo, zhi) = (ea.zlo, ea.zhi);
    if zlo > zhi {
        return single_bound(la, ea, ata) + single_bound(lb, eb, atb);
    }
    // atoms: unconditional (all four conditions certain), conditional-a (all but L1: counted iff
    // zeta <= y), conditional-b (all but H1: counted iff y <= zeta + u0)
    let mut uncond: I = 0;
    let mut ca: Vec<(I, I)> = Vec::new();
    let mut cb: Vec<(I, I)> = Vec::new();
    for a in ata {
        if a.c.iter().all(|&c| c) {
            uncond += a.w;
        } else if a.c[0] && a.c[2] && a.c[3] {
            ca.push((a.y, a.w));
        }
    }
    for b in atb {
        if b.c.iter().all(|&c| c) {
            uncond += b.w;
        } else if b.c[1] && b.c[2] && b.c[3] {
            cb.push((b.y, b.w));
        }
    }
    let fa_h = if h_a > l_a { fval(la, h_a) } else { 0 };
    let fb_l = fval(lb, l_b);
    let cont = |z: I| -> I {
        let alo = z.max(l_a);
        let a = if h_a > alo { fa_h - fval(la, alo) } else { 0 };
        let bhi = (z + u0g).max(h1b).min(h2b);
        let b = if bhi > l_b { fval(lb, bhi) - fb_l } else { 0 };
        a + b
    };
    // value at z, and the two one-sided limits (A right-limit / B left-limit)
    let eval = |z: I| -> I {
        let c = cont(z);
        let (mut a_at, mut a_plus, mut b_at, mut b_minus) = (0 as I, 0 as I, 0 as I, 0 as I);
        for &(y, w) in &ca {
            if y >= z {
                a_at += w;
                if y > z {
                    a_plus += w;
                }
            }
        }
        for &(y, w) in &cb {
            if y <= z + u0g {
                b_at += w;
                if y < z + u0g {
                    b_minus += w;
                }
            }
        }
        c + (a_at + b_minus).min(a_plus + b_at)
    };
    let mut best = eval(zlo).min(eval(zhi));
    let consider = |c: I, best: &mut I| {
        if c > zlo && c < zhi {
            let v = eval(c);
            if v < *best {
                *best = v;
            }
        }
    };
    for c in [l_a, h_a, h2b - u0g, h1b - u0g, l_b - u0g] {
        consider(c, &mut best);
    }
    for &bp in &la.bps {
        consider((bp as I) << K, &mut best);
    }
    for &bp in &lb.bps {
        consider(((bp as I) << K) - u0g, &mut best);
    }
    for &(y, _) in &ca {
        consider(y, &mut best);
    }
    for &(y, _) in &cb {
        consider(y - u0g, &mut best);
    }
    best + uncond
}

/// lower bound of the segment mass of one family of lines (DP over pairs, Lemma Z + sec 4.4)
/// which of the four conditions of point (px,py)/D hold for every pose of box b (exact, Lemma P)
fn pt_conds(d: I, px: I, py: I, b: &PBox) -> [bool; 4] {
    let cd = b.cd();
    let e = d * cd;
    let du = b.ud();
    let mut ok = [true; 4];
    for xi in 0..2 {
        for yi in 0..2 {
            let a = px * cd - b.xn[xi] * d;
            let bb = py * cd - b.yn[yi] * d;
            let qs = [
                (-2 * a - e, 4 * bb, 2 * a - e),
                (2 * a - e, -4 * bb, -2 * a - e),
                (-2 * bb - e, -4 * a, 2 * bb - e),
                (2 * bb - e, 4 * a, -2 * bb - e),
            ];
            for (k, &(p2, p1, p0)) in qs.iter().enumerate() {
                if ok[k] && !qmax_le0(p2, p1, p0, b.un[0], b.un[1], du) {
                    ok[k] = false;
                }
            }
        }
    }
    ok
}

/// the atoms of line l that can matter in this box, with their certain conditions
fn line_atoms(ctx: &Ctx, fb: &FrameBox, l: &Line) -> Vec<Atom> {
    if l.atoms.is_empty() {
        return Vec::new();
    }
    let dd = ctx.cv.d;
    let ylo = (fb.pb.yn[0] as f64) / (fb.pb.cd() as f64) - 0.75;
    let yhi = (fb.pb.yn[1] as f64) / (fb.pb.cd() as f64) + 0.75;
    let wall_l = l.pos - fb.wall_l == dd as i64;
    let wall_r = fb.wall_r - l.pos == dd as i64;
    let h1w = ctx.glo(add_lo(fb.y0.lo, fb.q_lo));
    let l1w = ctx.ghi(add_hi(fb.y1.hi, -fb.q_lo));
    let mut v = Vec::new();
    for &(t, w) in &l.atoms {
        let tf = t as f64 / dd as f64;
        if tf < ylo || tf > yhi {
            continue;
        }
        let mut c = pt_conds(dd, l.pos as I, t as I, &fb.pb);
        let yg = (t as I) << K;
        // Lemma W for atoms: H1 >= c_y + q (left wall), L1 <= c_y - q (right wall)
        if wall_l && yg <= h1w {
            c[0] = true;
        }
        if wall_r && yg >= l1w {
            c[1] = true;
        }
        v.push(Atom { y: yg, w: w << K, c });
    }
    v
}

/// lower bound of the segment + atom mass of one family of lines (DP over pairs, Lemma Z)
fn family_bound(ctx: &Ctx, fb: &FrameBox, lines: &[Line]) -> I {
    family_bound_ex(ctx, fb, lines, &[])
}

/// family_bound with the lines at the positions `skip` left out (they are bounded elsewhere, coupled
/// with an area cap: ZMX2_AREA.md Lemma K).  With `skip` empty this is exactly family_bound.
fn family_bound_ex(ctx: &Ctx, fb: &FrameBox, lines: &[Line], skip: &[i64]) -> I {
    let dd = ctx.cv.d;
    // active lines: a square meets a line only if |c_x - l| <= w/2 <= 0.7072
    let xlo = (fb.xn[0] as f64) / (fb.xd as f64) - 0.75;
    let xhi = (fb.xn[1] as f64) / (fb.xd as f64) + 0.75;
    let mut act: Vec<(usize, Ends, Vec<Atom>)> = Vec::new();
    for (i, l) in lines.iter().enumerate() {
        let p = l.pos as f64 / dd as f64;
        if p < xlo || p > xhi {
            continue;
        }
        if !skip.is_empty() && skip.contains(&l.pos) {
            continue;
        }
        act.push((i, line_ends(ctx, fb, l.pos), line_atoms(ctx, fb, l)));
    }
    let n = act.len();
    if n == 0 {
        return 0;
    }
    let debug = std::env::var("ZMX2_DEBUG").is_ok();
    // Lines pair only with a partner at distance exactly 1, so the active lines split into chains
    // of positions congruent mod 1 (sorted); the DP runs along each chain (Lemma DP).
    let mut chains: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
    for (k, (li, _, _)) in act.iter().enumerate() {
        chains.entry(lines[*li].pos.rem_euclid(dd as i64)).or_default().push(k);
    }
    let mut total: I = 0;
    for (_, ch) in chains {
        let m = ch.len();
        let mut best = vec![0 as I; m + 1]; // best[i] = best over the first i lines of the chain
        for i in 0..m {
            let (li, ei, ai) = &act[ch[i]];
            let sb = single_bound(&lines[*li], ei, ai);
            let mut b = best[i] + sb;
            if debug {
                let gf = ctx.g as f64;
                eprintln!(
                    "  line pos {:.4}: zeta [{:.6},{:.6}] h1lo {:.6} l2hi {:.6} h2lo {:.6} atoms {} single {:.6}",
                    lines[*li].pos as f64 / dd as f64,
                    ei.zlo as f64 / gf,
                    ei.zhi as f64 / gf,
                    ei.h1lo as f64 / gf,
                    ei.l2hi as f64 / gf,
                    ei.h2lo as f64 / gf,
                    ai.len(),
                    sb as f64 / ctx.target as f64
                );
            }
            if i >= 1 {
                let (lj, ej, aj) = &act[ch[i - 1]];
                if lines[*li].pos - lines[*lj].pos == dd as i64 {
                    let pb = pair_bound(&lines[*lj], ej, aj, &lines[*li], ei, ai, fb.u0g);
                    if debug {
                        eprintln!(
                            "  pair {:.3}-{:.3}: {:.6}",
                            lines[*lj].pos as f64 / dd as f64,
                            lines[*li].pos as f64 / dd as f64,
                            pb as f64 / ctx.target as f64
                        );
                    }
                    b = b.max(best[i - 1] + pb);
                }
            }
            best[i + 1] = b;
        }
        total += best[m];
    }
    total
}

struct BoxGeo {
    u0: Iv,
    u1: Iv,
    u0zero: bool,
    u0g: I,
    q_lo: f64,
    w_lo: f64, // lower bound of w(u) = cos + sin on [u0,u1]
}

fn box_geo(ctx: &Ctx, b: &PBox) -> BoxGeo {
    let ud = b.ud();
    let u0 = Iv::rat(b.un[0], ud);
    let u1 = Iv::rat(b.un[1], ud);
    let u0zero = b.un[0] == 0;
    let u0g = (b.un[0] * ctx.g).div_euclid(ud);
    // q lower bound: (1 - u1^2 + 2 u0^3) / (2 (1 + u1^2))
    let a = Iv::exact(u1.hi);
    let a2 = a.mul(a);
    let bb = Iv::exact(u0.lo.max(0.0));
    let b3 = bb.mul(bb).mul(bb);
    let num = Iv::exact(1.0).sub(a2).add(b3.scale(2.0));
    let den = Iv::exact(1.0).add(a2).scale(2.0);
    let q_lo = num.div_pos(den).lo;
    // w(u) = (1 + 2u - u^2)/(1 + u^2), unimodal on [0,1] with its max at sqrt2-1: min at ends
    let wpt = |x: Iv| -> f64 {
        let x2 = x.mul(x);
        Iv::exact(1.0).add(x.scale(2.0)).sub(x2).div_pos(Iv::exact(1.0).add(x2)).lo
    };
    let w_lo = wpt(u0).min(wpt(u1));
    BoxGeo { u0, u1, u0zero, u0g, q_lo, w_lo }
}

/// true if no pose of the box is admissible (Lemma E)
fn box_empty(ctx: &Ctx, b: &PBox, bg: &BoxGeo) -> bool {
    let cd = b.cd();
    let s = Iv::rat(ctx.cv.sx, ctx.cv.d);
    let half_w = bg.w_lo * 0.5; // exact halving (no underflow here)
    let rthr = s.sub(Iv::exact(half_w)).hi; // >= s - w_lo/2
    let x0 = Iv::rat(b.xn[0], cd);
    let x1 = Iv::rat(b.xn[1], cd);
    let y0 = Iv::rat(b.yn[0], cd);
    let y1 = Iv::rat(b.yn[1], cd);
    x1.hi < half_w || y1.hi < half_w || x0.lo > rthr || y0.lo > rthr
}

/// Lower bound of the line mass (segments + atoms), given `have` = the certain point weight
/// already counted: the default assignment's bound, and with `--sym-atoms` (sec 4.9, Lemma A),
/// if that does not reach the target, the larger of it and the mirrored assignment's bound.
fn segment_bound(ctx: &Ctx, b: &PBox, bg: &BoxGeo, have: I) -> I {
    let sb = segment_bound_lines(ctx, b, bg, &ctx.cv.vl, &ctx.cv.hl);
    if let Some((vl2, hl2)) = &ctx.cv.alt {
        if have + sb < ctx.target {
            return sb.max(segment_bound_lines(ctx, b, bg, vl2, hl2));
        }
    }
    sb
}

fn segment_bound_lines(ctx: &Ctx, b: &PBox, bg: &BoxGeo, vl: &[Line], hl: &[Line]) -> I {
    let (fv, fh) = frames(ctx, b, bg);
    family_bound(ctx, &fv, vl) + family_bound(ctx, &fh, hl)
}

/// The two frames of a box: vertical lines (identity) and horizontal lines (rotated, Lemma T).
fn frames(ctx: &Ctx, b: &PBox, bg: &BoxGeo) -> (FrameBox, FrameBox) {
    let cd = b.cd();
    let sxd = ctx.cv.sx as i64;
    // vertical lines: frame = identity; walls at 0 and s
    let fv = FrameBox {
        xn: b.xn,
        xd: cd,
        y0: Iv::rat(b.yn[0], cd),
        y1: Iv::rat(b.yn[1], cd),
        u0: bg.u0,
        u1: bg.u1,
        u0zero: bg.u0zero,
        u0g: bg.u0g,
        q_lo: bg.q_lo,
        wall_l: 0,
        wall_r: sxd,
        pb: *b,
    };
    // horizontal lines: rotated frame (x', y') = (-y, x) (Lemma T); walls at -s and 0
    let fh = FrameBox {
        xn: [-b.yn[1], -b.yn[0]],
        xd: cd,
        y0: Iv::rat(b.xn[0], cd),
        y1: Iv::rat(b.xn[1], cd),
        u0: bg.u0,
        u1: bg.u1,
        u0zero: bg.u0zero,
        u0g: bg.u0g,
        q_lo: bg.q_lo,
        wall_l: -sxd,
        wall_r: 0,
        pb: PBox { xn: [-b.yn[1], -b.yn[0]], yn: b.xn, cl: b.cl, un: b.un, ul: b.ul },
    };
    (fv, fh)
}

// =====================================================================================
// Area densities on axis-parallel rectangles (search/ZMX2_AREA.md)
// =====================================================================================

/// Lower bound of all the mass of `Q` except the ordinary points (which `have` already counts):
/// the lines (sec 4) and, if the cover has area densities, the rectangles (ZMX2_AREA.md sec 5).
/// Without rectangles this is exactly `segment_bound` (same code path).
fn mass_bound(ctx: &Ctx, b: &PBox, bg: &BoxGeo, have: I) -> I {
    let v = if ctx.cv.rects.is_empty() { segment_bound(ctx, b, bg, have) } else { area_bound(ctx, b, bg, have) };
    if have + v < ctx.target && FIRST_ORDER.load(Ordering::Relaxed) && 16 * b.un[1] <= b.ud() {
        // Lemma U on the box extended down to u = 0 (its poses include those of b; `have` = the points certain
        // on b, valid for the poses of b, which are all we bound here)
        let fo = if b.un[0] == 0 {
            first_order_bound(ctx, b, bg, have)
        } else {
            let mut b0 = *b;
            b0.un[0] = 0;
            let bg0 = box_geo(ctx, &b0);
            first_order_bound(ctx, &b0, &bg0, have)
        };
        if let Some(fo) = fo {
            // first_order_bound includes the point weight `have`
            return v.max(fo - have);
        }
    }
    v
}

/// enclosures of C(u) = (1-u^2)/(1+u^2) and S(u) = 2u/(1+u^2) at a thin u in [0,1)
fn cs_at(u: Iv) -> (Iv, Iv) {
    let u2 = u.mul(u);
    let n = Iv::exact(1.0).add(u2);
    (Iv::exact(1.0).sub(u2).div_pos(n), u.scale(2.0).div_pos(n))
}

/// enclosure of e(u) = C(u) - S(u) = (1 - 2u - u^2)/(1+u^2) at a thin u
fn e_at(u: Iv) -> Iv {
    let u2 = u.mul(u);
    Iv::exact(1.0).sub(u.scale(2.0)).sub(u2).div_pos(Iv::exact(1.0).add(u2))
}

/// upper bound of w(u) = C + S on [u0,u1] (w increases on [0, sqrt2-1], decreases after; max sqrt 2)
fn w_hi(bg: &BoxGeo) -> f64 {
    let wpt = |x: Iv| -> f64 {
        let x2 = x.mul(x);
        Iv::exact(1.0).add(x.scale(2.0)).sub(x2).div_pos(Iv::exact(1.0).add(x2)).hi
    };
    let mut w = wpt(bg.u0).max(wpt(bg.u1));
    // sqrt2 - 1 = 0.414213562373...: if [u0,u1] may contain it, the maximum sqrt 2 may be attained
    if bg.u0.lo <= 0.4142135624 && bg.u1.hi >= 0.4142135623 {
        w = w.max(up(std::f64::consts::SQRT_2));
    }
    w
}

/// Lemma 1 (exact containment of one side): h = hn/hd >= w(u)/2 for every u of the box, i.e.
/// (2h+1)u^2 - 2u + (2h-1) >= 0 on [u0,u1]; decided exactly (as Lemma P).
fn side_contained(hn: I, hd: I, b: &PBox) -> bool {
    if hn <= 0 {
        return false;
    }
    qmax_le0(-(2 * hn + hd), 2 * hd, -(2 * hn - hd), b.un[0], b.un[1], b.ud())
}

/// Lemma 2 (cap bound): an upper bound, valid for every pose of the box, of the area of the part of Q
/// beyond a line whose signed distance from the centre (positive: the centre is on the inner side)
/// is >= h.  Uses a(h, theta) decreasing in h and the exact piecewise formula of a(h, theta).
fn cap_ub(h: f64, bg: &BoxGeo, whi: f64) -> f64 {
    let m_hi_w = whi * 0.5; // >= M = w/2 on the box
    if h >= m_hi_w {
        return 0.0;
    }
    let m_lo_w = bg.w_lo * 0.5; // <= M
    let (c0, s0) = cs_at(bg.u0);
    let (c1, s1) = cs_at(bg.u1);
    // C decreasing and S increasing on [0,1)
    let (c_lo, c_hi, s_lo, s_hi) = (c1.lo, c0.hi, s0.lo.max(0.0), s1.hi);
    let p_lo = c_lo.min(s_lo).max(0.0); // p = min(C,S)
    let p_hi = c_hi.min(s_hi);
    let pp_lo = c_lo.max(s_lo); // P = max(C,S) >= 1/sqrt2
    let pp_hi = c_hi.max(s_hi);
    // m = |C - S|/2; C - S is decreasing in u
    let e0 = e_at(bg.u0);
    let e1 = e_at(bg.u1);
    let (m_lo, m_hi) = if e1.lo >= 0.0 {
        (e1.lo * 0.5, e0.hi * 0.5)
    } else if e0.hi <= 0.0 {
        (-e0.hi * 0.5, -e1.lo * 0.5)
    } else {
        (0.0, e0.hi.max(-e1.lo) * 0.5)
    };
    let mut a: f64 = 0.0;
    // (ii) m <= h <= M: a = s^2/(2pP), s = M - h in [0, p]; a <= s/(2P), and <= s^2/(2 p P)
    if h >= m_lo {
        let sp = Iv::exact(m_hi_w).sub(Iv::exact(h)).hi.min(p_hi).max(0.0);
        let two_pp = Iv::exact(pp_lo).scale(2.0);
        let mut v = Iv::exact(sp).div_pos(two_pp).hi;
        if p_lo > 0.0 {
            let v2 = Iv::exact(sp).mul(Iv::exact(sp)).div_pos(two_pp.mul(Iv::exact(p_lo))).hi;
            v = v.min(v2);
        }
        a = a.max(v);
    }
    // (iii) |h| <= m: a = 1/2 - h/P
    if h <= m_hi && h >= -m_hi {
        let v = if h >= 0.0 {
            Iv::exact(0.5).sub(Iv::exact(h).div_pos(Iv::exact(pp_hi))).hi
        } else {
            Iv::exact(0.5).add(Iv::exact(-h).div_pos(Iv::exact(pp_lo))).hi
        };
        a = a.max(v);
    }
    // (iv) -M <= h <= -m: a = 1 - (M+h)^2/(2pP), M + h in [0, p];  (v) h <= -M: a = 1
    if h <= -m_lo {
        let q = Iv::exact(m_lo_w).add(Iv::exact(h)).lo.max(0.0);
        let v = if q > 0.0 {
            let den = Iv::exact(p_hi).mul(Iv::exact(pp_hi)).scale(2.0);
            Iv::exact(1.0).sub(Iv::exact(q).mul(Iv::exact(q)).div_pos(den)).hi
        } else {
            1.0
        };
        a = a.max(v);
    }
    a.min(1.0)
}

/// Lemma K (a line coupled with the area cap beyond it).  `fb` is the frame, the line is at `pos`
/// (frame units), the rectangle `ri` lies on the side of the line that contains the centre of Q for
/// every pose of the box (`left`: the rectangle is to the right of the line, i.e. this is its left
/// side).  Returns a lower bound, valid at every pose of the box, of
///     nu_line(Q) - rho_ri * area(Q beyond the line)        (in bound units).
fn coupled_side(ctx: &Ctx, fb: &FrameBox, line: Option<&Line>, pos: i64, left: bool, ri: usize, whi: f64, bg: &BoxGeo) -> I {
    let dd = ctx.cv.d;
    let den = fb.xd * dd;
    // Lemma C/M enclosures of f1, f2, g1, g2 over the box (as line_ends)
    let mut f1 = Iv { lo: INF, hi: -INF };
    let mut f2 = f1;
    let mut g1 = f1;
    let mut g2 = f1;
    for e in 0..2 {
        let dnum = fb.xn[e] * dd - (pos as I) * fb.xd;
        let dv = Iv::rat(dnum, den);
        let mm = Iv::rat(2 * dnum - den, 2 * den);
        let pm = Iv::rat(2 * dnum + den, 2 * den);
        f1 = f1.hull(h_range(mm.scale(0.5), pm.scale(-0.5), fb.u0, fb.u1, fb.u0zero));
        f2 = f2.hull(h_range(pm.scale(0.5), mm.scale(-0.5), fb.u0, fb.u1, fb.u0zero));
        g1 = g1.hull(g_range(dv, -1.0, fb.u0, fb.u1));
        g2 = g2.hull(g_range(dv, 1.0, fb.u0, fb.u1));
    }
    // distance from the centre to the line, on the rectangle's side: >= h0 (exact numerators over den)
    let (hn_min, hn_max) = if left {
        (fb.xn[0] * dd - (pos as I) * fb.xd, fb.xn[1] * dd - (pos as I) * fb.xd)
    } else {
        ((pos as I) * fb.xd - fb.xn[1] * dd, (pos as I) * fb.xd - fb.xn[0] * dd)
    };
    assert!(hn_min >= 0, "coupled side with the centre beyond the line");
    let h0 = Iv::rat(hn_min, den).lo;
    // depth t = w/2 - h <= whi/2 - h0; cap <= l * phi, phi = t - min(t, p)/2, p = min(C,S) (Lemma K (a))
    let t = Iv::exact(whi * 0.5).sub(Iv::exact(h0)).hi.max(0.0);
    let (c1, _) = cs_at(fb.u1);
    let (_, s0) = cs_at(fb.u0);
    let p_min = c1.lo.min(s0.lo).max(0.0);
    let phi = Iv::exact(t).sub(Iv::exact(t.min(p_min) * 0.5)).hi.max(0.0);
    let phin = ctx.ghi(phi).max(0);
    let num = ctx.rect_wlc[ri].checked_mul(phin).expect("Lemma K: kappa overflow");
    let dk = ctx.rect_den[ri];
    let kappa = (num + dk - 1).div_euclid(dk); // ceil: bound units per grid unit of chord length
    // chord-end ranges (u > 0 poses, and the u -> 0+ limits), clipped to the extent of Q
    let mh = whi * 0.5;
    let lo_min_f = add_lo(fb.y0.lo, f1.lo.max(g1.lo));
    let lo_max_f = add_hi(fb.y1.hi, f1.hi.max(g1.hi));
    let hi_min_f = add_lo(fb.y0.lo, f2.lo.min(g2.lo));
    let hi_max_f = add_hi(fb.y1.hi, f2.hi.min(g2.hi));
    let ext_lo = add_lo(fb.y0.lo, -mh);
    let ext_hi = add_hi(fb.y1.hi, mh);
    let mut ilo = (ctx.glo(lo_min_f.max(ext_lo)), ctx.ghi(lo_max_f.min(ext_hi)));
    let mut ihi = (ctx.glo(hi_min_f.max(ext_lo)), ctx.ghi(hi_max_f.min(ext_hi)));
    // Lemma K (d), exact clip of the chord ends: g1 = -dT - R/2 = -dT - 1/2 - u^2/(1-u^2) >= -1/2 when
    // -d >= u1/(2(1-u1^2)) (as -dT >= -2du and u^2/(1-u^2) <= u u1/(1-u1^2)), so lo = c + max(f1, g1) >= c - 1/2 >=
    // y0 - 1/2 (a grid point; at u = 0 the limit chord starts at c - 1/2); symmetrically hi <= y1 + 1/2 when
    // d >= u1/(2(1-u1^2))
    {
        let u1v = Iv::exact(fb.u1.hi);
        let thr = u1v.div_pos(Iv::exact(1.0).sub(u1v.mul(u1v)).scale(2.0)).hi; // >= u1/(2(1-u1^2))
        let (dmin, dmax) = (Iv::rat(fb.xn[0] * dd - (pos as I) * fb.xd, den).lo, Iv::rat(fb.xn[1] * dd - (pos as I) * fb.xd, den).hi);
        let g = ctx.g;
        let pb = &fb.pb;
        if dmax <= -thr {
            let y0g = fdiv(pb.yn[0] * g, pb.cd()) - g / 2; // <= y0 - 1/2
            ilo.0 = ilo.0.max(y0g);
            ihi.0 = ihi.0.max(y0g);
        }
        if dmin >= thr {
            let y1g = cdiv(pb.yn[1] * g, pb.cd()) + g / 2; // >= y1 + 1/2
            ilo.1 = ilo.1.min(y1g);
            ihi.1 = ihi.1.min(y1g);
        }
    }
    // the chord is non-empty at every pose of the box (Lemma K (c))
    let far0 = fb.u0zero && 2 * hn_max > den; // a theta = 0 pose with |d| > 1/2 (empty chord)
    let nonempty = lo_max_f < hi_min_f && !far0;
    let gval = |y: I| -> I {
        let fv = match line {
            Some(l) => fval(l, y),
            None => 0,
        };
        fv - kappa.checked_mul(y).expect("Lemma K: kappa*y overflow")
    };
    // candidates: the range ends and the breakpoints inside (G is piecewise linear between them)
    let (lo_all, hi_all) = (ilo.0.min(ihi.0), ilo.1.max(ihi.1));
    let mut cands: Vec<I> = vec![ilo.0, ilo.1, ihi.0, ihi.1];
    if let Some(l) = line {
        for &bp in &l.bps {
            let y = (bp as I) << K;
            if y > lo_all && y < hi_all {
                cands.push(y);
            }
        }
    }
    cands.sort();
    cands.dedup();
    // inf over lo in ilo, hi in ihi, lo <= hi of G(hi) - G(lo): scan hi upwards with the running max of
    // G over the lo-candidates <= hi (Lemma K (b))
    let mut runmax: Option<I> = None;
    let mut best: Option<I> = None;
    for &c in &cands {
        let g = gval(c);
        if c >= ilo.0 && c <= ilo.1 {
            runmax = Some(runmax.map_or(g, |r| r.max(g)));
        }
        if c >= ihi.0 && c <= ihi.1 {
            if let Some(r) = runmax {
                let v = g - r;
                best = Some(best.map_or(v, |b| b.min(v)));
            }
        }
    }
    let mut val = match best {
        None => 0, // no non-empty chord is possible
        Some(v) => {
            if nonempty {
                v
            } else {
                v.min(0)
            }
        }
    };
    // atoms of the line certainly in Q (all four conditions)
    if let Some(l) = line {
        for a in line_atoms(ctx, fb, l) {
            if a.c.iter().all(|&c| c) {
                val += a.w;
            }
        }
    }
    let _ = bg;
    val
}

/// Lemma 3 (b): a lower bound of area(Q(c, u) cap {x <= ax (left) or >= ax, y <= ay (bottom) or >= ay}) at
/// one pose (c given by interval enclosures of a rational point, u thin), by clipping Q with interval
/// arithmetic.  None if a vertex cannot be classified (then no bound is used).
fn quad_area_lb(cx: Iv, cy: Iv, u: Iv, ax: Iv, ay: Iv, left: bool, bottom: bool) -> Option<f64> {
    let (c, s) = cs_at(u);
    let h = 0.5;
    let mut poly: Vec<(Iv, Iv)> = [(-h, -h), (h, -h), (h, h), (-h, h)]
        .iter()
        .map(|&(a, b)| (cx.add(c.scale(a)).sub(s.scale(b)), cy.add(s.scale(a)).add(c.scale(b))))
        .collect();
    // keep f(p) <= 0, f = sg (coord - lim)
    let clip = |poly: Vec<(Iv, Iv)>, axis: usize, sg: f64, lim: Iv| -> Option<Vec<(Iv, Iv)>> {
        let f = |p: &(Iv, Iv)| -> Iv {
            let v = if axis == 0 { p.0 } else { p.1 };
            v.sub(lim).scale(sg)
        };
        let n = poly.len();
        let fs: Vec<Iv> = poly.iter().map(f).collect();
        let mut out = Vec::new();
        for k in 0..n {
            let (fp, fq) = (fs[k], fs[(k + 1) % n]);
            let inp = fp.hi <= 0.0;
            let outp = fp.lo > 0.0;
            if !inp && !outp {
                return None;
            }
            if inp {
                out.push(poly[k]);
            }
            let inq = fq.hi <= 0.0;
            let outq = fq.lo > 0.0;
            if !inq && !outq {
                return None;
            }
            if (inp && outq && fp.hi < 0.0) || (outp && inq && fq.hi < 0.0) {
                // crossing point p + t (q - p), t = fp / (fp - fq) in (0, 1)
                let (p, q) = (poly[k], poly[(k + 1) % n]);
                let den = fp.sub(fq);
                let t = if den.lo > 0.0 {
                    fp.div_pos(den)
                } else if den.hi < 0.0 {
                    fp.neg().div_pos(den.neg())
                } else {
                    return None;
                };
                let t = Iv { lo: t.lo.max(0.0), hi: t.hi.min(1.0) };
                out.push((p.0.add(t.mul(q.0.sub(p.0))), p.1.add(t.mul(q.1.sub(p.1)))));
            }
        }
        Some(out)
    };
    poly = clip(poly, 0, if left { 1.0 } else { -1.0 }, ax)?;
    if poly.len() < 3 {
        return Some(0.0);
    }
    poly = clip(poly, 1, if bottom { 1.0 } else { -1.0 }, ay)?;
    if poly.len() < 3 {
        return Some(0.0);
    }
    let n = poly.len();
    let mut a = Iv::exact(0.0);
    for k in 0..n {
        let (p, q) = (poly[k], poly[(k + 1) % n]);
        a = a.add(p.0.mul(q.1).sub(q.0.mul(p.1)));
    }
    Some((a.lo * 0.5).max(0.0))
}

/// The bound with area densities (ZMX2_AREA.md sec 5): the larger of
///   config 0: sum_r rho_r * max(0, 1 - sum_sides cap_ub)            + the line bound (sec 4), and
///   config 1: sum_r rho_r * (1 - sum_{uncoupled sides} cap_ub)
///             + sum_{coupled sides} Lemma K + the line bound without the coupled lines.
fn area_bound(ctx: &Ctx, b: &PBox, bg: &BoxGeo, have: I) -> I {
    let cv = ctx.cv;
    let cd = b.cd();
    let dd = cv.d;
    let hd = cd * dd;
    let whi = w_hi(bg);
    let nr = cv.rects.len();
    // sides: 0 left (x = x0), 1 right (x = x1), 2 bottom (y = y0), 3 top (y = y1); hn/hd = distance from
    // the centre to the side line, inner side positive, minimised over the box
    let mut hn = vec![[0 as I; 4]; nr];
    let mut cont = vec![[false; 4]; nr];
    let mut cap = vec![[0.0f64; 4]; nr];
    for (ri, r) in cv.rects.iter().enumerate() {
        hn[ri] = [
            b.xn[0] * dd - (r.x0 as I) * cd,
            (r.x1 as I) * cd - b.xn[1] * dd,
            b.yn[0] * dd - (r.y0 as I) * cd,
            (r.y1 as I) * cd - b.yn[1] * dd,
        ];
        // a side on (or beyond) the container wall holds for every admissible pose (Lemma 1 (b))
        let at_wall = [r.x0 <= 0, r.x1 as I >= cv.sx, r.y0 <= 0, r.y1 as I >= cv.sx];
        for k in 0..4 {
            cont[ri][k] = at_wall[k] || side_contained(hn[ri][k], hd, b);
            if !cont[ri][k] {
                cap[ri][k] = cap_ub(Iv::rat(hn[ri][k], hd).lo, bg, whi);
            }
        }
    }
    // rho_r * (1 - capsum), rounded down, in bound units (may be negative); exact when capsum = 0 (Lemma 1)
    // Lemma 3 (corners): area(Q cap quadrant) >= that of the inscribed axis-parallel square c + [-r, r]^2,
    // r = 1/(2w) >= 1/(2 whi), for each of the four corner quadrants of each rectangle
    let rin = Iv::exact(1.0).div_pos(Iv::exact(whi).scale(2.0)).lo;
    let bx = [Iv::rat(b.xn[0], cd), Iv::rat(b.xn[1], cd), Iv::rat(b.yn[0], cd), Iv::rat(b.yn[1], cd)];
    let mut corner = vec![0.0f64; nr];
    for (ri, r) in cv.rects.iter().enumerate() {
        let rx0 = Iv::rat(r.x0 as I, dd);
        let rx1 = Iv::rat(r.x1 as I, dd);
        let ry0 = Iv::rat(r.y0 as I, dd);
        let ry1 = Iv::rat(r.y1 as I, dd);
        let cl = |v: f64| v.max(0.0).min(2.0 * rin);
        // overlap lengths of [c - r, c + r] with {x <= x0_R} (left), {x >= x1_R} (right), etc., minimised
        let lx = [cl(rx0.sub(bx[1]).add(Iv::exact(rin)).lo), cl(bx[0].add(Iv::exact(rin)).sub(rx1).lo)];
        let ly = [cl(ry0.sub(bx[3]).add(Iv::exact(rin)).lo), cl(bx[2].add(Iv::exact(rin)).sub(ry1).lo)];
        let mut acc = 0.0f64;
        for (i, &a) in lx.iter().enumerate() {
            for (j, &bb) in ly.iter().enumerate() {
                // the quadrant beyond corner (i: 0 left / 1 right, j: 0 bottom / 1 top)
                let sq = if a > 0.0 && bb > 0.0 { dn(a * bb) } else { 0.0 };
                // Lemma 3 (b): exact clip at the pessimal centre and the mid angle, minus the theta-Lipschitz term
                let ax = if i == 0 { rx0 } else { rx1 };
                let ay = if j == 0 { ry0 } else { ry1 };
                let cxp = if i == 0 { bx[1] } else { bx[0] };
                let cyp = if j == 0 { bx[3] } else { bx[2] };
                let mut v = sq;
                let um = Iv::rat(b.un[0] + b.un[1], 2 * b.ud());
                if let Some(ar) = quad_area_lb(cxp, cyp, um, ax, ay, i == 0, j == 0) {
                    // |theta - theta_m| <= u1 - u0 and |d area / d theta| <= 2 sqrt 2
                    let lip = Iv::rat(b.un[1] - b.un[0], b.ud()).mul(Iv::exact(std::f64::consts::SQRT_2)).scale(2.0);
                    v = v.max(Iv::exact(ar).sub(lip).lo);
                }
                if v > 0.0 {
                    acc = dn(acc + v).max(0.0);
                }
            }
        }
        corner[ri] = acc;
    }
    // rho_r * (1 - capsum + corner), rounded down, in bound units (may be negative); exact when capsum = 0
    let area_units = |ri: usize, capsum: f64| -> I {
        let tr = ctx.rect_tr[ri];
        if capsum == 0.0 {
            return tr;
        }
        // tr <= rho (exact) < tr + 1: the loss is scaled by an upper bound, the gain by a lower bound, so the
        // result is <= rho (1 - capsum + corner) also when that is negative (config 1)
        let gain = if corner[ri] > 0.0 { dn(corner[ri] * dn(tr as f64)).floor().max(0.0) as I } else { 0 };
        tr - (up(capsum * up((tr + 1) as f64)).ceil() as I) + gain
    };
    // sum of cap bounds, rounded up (an exact 0 stays 0)
    let capsum = |ri: usize, use_k: &dyn Fn(usize) -> bool| -> f64 {
        let mut acc = 0.0f64;
        for k in 0..4 {
            let c = cap[ri][k];
            if use_k(k) && c > 0.0 {
                acc = if acc == 0.0 { c } else { up(acc + c) };
            }
        }
        acc
    };
    // config 0
    let mut area0: I = 0;
    for ri in 0..nr {
        area0 += area_units(ri, capsum(ri, &|_| true)).max(0);
    }
    let lines0 = segment_bound(ctx, b, bg, have + area0);
    let tot0 = area0 + lines0;
    let debug = std::env::var("ZMX2_DEBUG").is_ok();
    if debug {
        eprintln!("  area: contained {:?} caps {:?} area0 {} lines0 {} (unit {})", cont, cap, area0, lines0, ctx.target);
    }
    if have + tot0 >= ctx.target {
        return tot0;
    }
    // config 1: couple every side whose line keeps the centre on the rectangle's side for the whole box and
    // that is not exactly contained; one coupling per line
    let (fv, fh) = frames(ctx, b, bg);
    let mut skip_v: Vec<i64> = Vec::new();
    let mut skip_h: Vec<i64> = Vec::new();
    let mut coupled: Vec<(usize, usize, i64)> = Vec::new(); // (rect, side, frame pos)
    for (ri, r) in cv.rects.iter().enumerate() {
        for k in 0..4 {
            if cont[ri][k] || hn[ri][k] < 0 {
                continue;
            }
            // frame position: vertical family x; horizontal family (rotated frame, Lemma T) -y
            let (pos, list) = match k {
                0 => (r.x0, &mut skip_v),
                1 => (r.x1, &mut skip_v),
                2 => (-r.y0, &mut skip_h),
                _ => (-r.y1, &mut skip_h),
            };
            if list.contains(&pos) {
                continue;
            }
            list.push(pos);
            coupled.push((ri, k, pos));
        }
    }
    if coupled.is_empty() {
        return tot0;
    }
    let mut area1: I = 0;
    for ri in 0..nr {
        let is_coupled = |k: usize| coupled.iter().any(|c| c.0 == ri && c.1 == k);
        area1 += area_units(ri, capsum(ri, &|k| !is_coupled(k)));
    }
    let line_sets: Vec<(&[Line], &[Line])> = match &cv.alt {
        None => vec![(&cv.vl[..], &cv.hl[..])],
        Some((vl2, hl2)) => vec![(&cv.vl[..], &cv.hl[..]), (&vl2[..], &hl2[..])],
    };
    let mut best1: Option<I> = None;
    for (vl, hl) in line_sets {
        let mut t = area1;
        for &(ri, k, pos) in &coupled {
            // left side (k = 0) and top side (k = 3, rotated frame): rectangle to the right of the line
            let (fb, lines, left) = match k {
                0 => (&fv, vl, true),
                1 => (&fv, vl, false),
                2 => (&fh, hl, false),
                _ => (&fh, hl, true),
            };
            let line = lines.iter().find(|l| l.pos == pos);
            let cs_v = coupled_side(ctx, fb, line, pos, left, ri, whi, bg);
            if debug {
                eprintln!("  coupled side rect {} side {} pos {}: {}", ri, k, pos, cs_v);
            }
            t += cs_v;
        }
        t += family_bound_ex(ctx, &fv, vl, &skip_v) + family_bound_ex(ctx, &fh, hl, &skip_h);
        best1 = Some(best1.map_or(t, |x| x.max(t)));
    }
    tot0.max(best1.unwrap())
}


// ------------------------------------------------------------------ points (Lemma P)

/// max over u in [U0/Du, U1/Du] of p2 u^2 + p1 u + p0 is <= 0  (exact)
#[inline]
fn qmax_le0(p2: I, p1: I, p0: I, u0: I, u1: I, du: I) -> bool {
    let val = |u: I| p2 * u * u + p1 * u * du + p0 * du * du;
    if val(u0) > 0 || val(u1) > 0 {
        return false;
    }
    if p2 < 0 {
        let m = -2 * p2;
        let v = p1 * du;
        if u0 * m < v && v < u1 * m && 4 * p2 * p0 - p1 * p1 < 0 {
            return false;
        }
    }
    true
}

/// point (px,py)/D lies in Q for every pose of the box (exact; Lemma P)
fn pt_certain_in(d: I, px: I, py: I, b: &PBox) -> bool {
    let cd = b.cd();
    let e = d * cd;
    let du = b.ud();
    for xi in 0..2 {
        for yi in 0..2 {
            let a = px * cd - b.xn[xi] * d;
            let bb = py * cd - b.yn[yi] * d;
            let qs = [
                (-2 * a - e, 4 * bb, 2 * a - e),
                (2 * a - e, -4 * bb, -2 * a - e),
                (-2 * bb - e, -4 * a, 2 * bb - e),
                (2 * bb - e, 4 * a, -2 * bb - e),
            ];
            for (p2, p1, p0) in qs {
                if !qmax_le0(p2, p1, p0, b.un[0], b.un[1], du) {
                    return false;
                }
            }
        }
    }
    true
}

/// float classification of a point against a box: 0 = never in Q (drop), 1 = maybe certain-in
/// (run the exact test), 2 = undecided (keep).  Only used to *skip* work; never certifies.
fn pt_float_class(px: f64, py: f64, x: [f64; 2], y: [f64; 2], u: [f64; 2]) -> u8 {
    let mut all_in = true;
    let mut maxes = [f64::NEG_INFINITY; 4];
    let mut mins = [f64::INFINITY; 4];
    for xi in 0..2 {
        for yi in 0..2 {
            let a = px - x[xi];
            let b = py - y[yi];
            let qs = [
                (-2.0 * a - 1.0, 4.0 * b, 2.0 * a - 1.0),
                (2.0 * a - 1.0, -4.0 * b, -2.0 * a - 1.0),
                (-2.0 * b - 1.0, -4.0 * a, 2.0 * b - 1.0),
                (2.0 * b - 1.0, 4.0 * a, -2.0 * b - 1.0),
            ];
            for (k, &(p2, p1, p0)) in qs.iter().enumerate() {
                let f = |t: f64| p2 * t * t + p1 * t + p0;
                let mut mx = f(u[0]).max(f(u[1]));
                let mut mn = f(u[0]).min(f(u[1]));
                if p2 != 0.0 {
                    let v = -p1 / (2.0 * p2);
                    if v > u[0] && v < u[1] {
                        mx = mx.max(f(v));
                        mn = mn.min(f(v));
                    }
                }
                if mx > 1e-9 {
                    all_in = false;
                }
                maxes[k] = maxes[k].max(mx);
                mins[k] = mins[k].min(mn);
            }
        }
    }
    for k in 0..4 {
        if mins[k] > 1e-9 {
            return 0;
        }
    }
    if all_in {
        1
    } else {
        2
    }
}

// =====================================================================================
// Float evaluation (diagnostics only: never used for a verdict)
// =====================================================================================

fn fline_mass(l: &Line, d: f64, lo: f64, hi: f64, wlc: f64) -> f64 {
    // F in floats
    let f = |y: f64| -> f64 {
        let yy = y * d;
        let n = l.bps.len();
        if yy <= l.bps[0] as f64 {
            return 0.0;
        }
        if yy >= l.bps[n - 1] as f64 {
            return l.cum[n - 1] as f64;
        }
        let mut k = 0;
        while k + 1 < n && (l.bps[k + 1] as f64) <= yy {
            k += 1;
        }
        l.cum[k] as f64 + l.dens[k] as f64 * (yy - l.bps[k] as f64)
    };
    let mut m = 0.0;
    for &(t, w) in &l.atoms {
        let tf = t as f64 / d;
        if tf >= lo - 1e-12 && tf <= hi + 1e-12 {
            m += w as f64 / wlc;
        }
    }
    if hi <= lo || l.bps.is_empty() {
        m
    } else {
        m + (f(hi) - f(lo)) / wlc
    }
}

/// float mu(Q(x,y,u)), u in [0,1)
fn float_mass(cv: &Cover, x: f64, y: f64, u: f64) -> f64 {
    let n = 1.0 + u * u;
    let c = (1.0 - u * u) / n;
    let s = 2.0 * u / n;
    let d = cv.d as f64;
    let mut m = 0.0;
    let wf = cv.w as f64;
    for &(px, py, w) in &cv.pts {
        let a = px as f64 / d - x;
        let b = py as f64 / d - y;
        let xx = a * c + b * s;
        let yy = -a * s + b * c;
        if xx.abs() <= 0.5 + 1e-12 && yy.abs() <= 0.5 + 1e-12 {
            m += w as f64 / wf;
        }
    }
    let wlc = (cv.w * cv.lc) as f64;
    let chord = |dd: f64, cy: f64| -> (f64, f64) {
        // vertical line at signed offset d = c_x - l
        let (l1, h1) = if s > 1e-15 {
            ((dd * c - 0.5) / s, (dd * c + 0.5) / s)
        } else if dd.abs() <= 0.5 + 1e-12 {
            (-1e9, 1e9)
        } else {
            (1e9, -1e9)
        };
        let l2 = (-dd * s - 0.5) / c;
        let h2 = (-dd * s + 0.5) / c;
        (cy + l1.max(l2), cy + h1.min(h2))
    };
    for l in &cv.vl {
        let (lo, hi) = chord(x - l.pos as f64 / d, y);
        m += fline_mass(l, d, lo, hi, wlc);
    }
    // horizontal: rotated frame c' = (-y, x)
    for l in &cv.hl {
        let (lo, hi) = chord(-y - l.pos as f64 / d, x);
        m += fline_mass(l, d, lo, hi, wlc);
    }
    for r in &cv.rects {
        let area = (r.x1 - r.x0) as f64 * (r.y1 - r.y0) as f64 / (d * d);
        m += float_rect_area(x, y, c, s, r, d) / area * (r.w as f64 / wf);
    }
    m
}

/// float area of Q(x, y; C, S) cap rectangle r (diagnostics only): clip the square by four half-planes
fn float_rect_area(x: f64, y: f64, c: f64, s: f64, r: &Rect, d: f64) -> f64 {
    let mut poly: Vec<(f64, f64)> = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
        .iter()
        .map(|&(a, b)| (x + c * a - s * b, y + s * a + c * b))
        .collect();
    // keep sg * coord(axis) <= lim
    let clip = |poly: Vec<(f64, f64)>, axis: usize, sg: f64, lim: f64| -> Vec<(f64, f64)> {
        let f = |p: &(f64, f64)| sg * (if axis == 0 { p.0 } else { p.1 }) - lim;
        let mut out = Vec::new();
        for k in 0..poly.len() {
            let (p, q) = (poly[k], poly[(k + 1) % poly.len()]);
            let (fp, fq) = (f(&p), f(&q));
            if fp <= 0.0 {
                out.push(p);
            }
            if (fp < 0.0 && fq > 0.0) || (fp > 0.0 && fq < 0.0) {
                let t = fp / (fp - fq);
                out.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
            }
        }
        out
    };
    poly = clip(poly, 0, -1.0, -(r.x0 as f64) / d);
    poly = clip(poly, 0, 1.0, r.x1 as f64 / d);
    poly = clip(poly, 1, -1.0, -(r.y0 as f64) / d);
    poly = clip(poly, 1, 1.0, r.y1 as f64 / d);
    if poly.len() < 3 {
        return 0.0;
    }
    let mut a = 0.0;
    for k in 0..poly.len() {
        let (p, q) = (poly[k], poly[(k + 1) % poly.len()]);
        a += p.0 * q.1 - q.0 * p.1;
    }
    0.5 * a.abs()
}

fn admissible_f(cv: &Cover, x: f64, y: f64, u: f64) -> bool {
    let n = 1.0 + u * u;
    let w = (1.0 - u * u + 2.0 * u) / n;
    let s = cv.sx as f64 / cv.d as f64;
    x >= w / 2.0 && x <= s - w / 2.0 && y >= w / 2.0 && y <= s - w / 2.0
}

/// float minimum of mu over admissible sample poses of a box (grid + pseudo-random)
fn float_min_box(cv: &Cover, b: &PBox, nsamp: usize) -> (f64, f64, f64, f64) {
    let cd = b.cd() as f64;
    let ud = b.ud() as f64;
    let (x0, x1) = (b.xn[0] as f64 / cd, b.xn[1] as f64 / cd);
    let (y0, y1) = (b.yn[0] as f64 / cd, b.yn[1] as f64 / cd);
    let (u0, u1) = (b.un[0] as f64 / ud, b.un[1] as f64 / ud);
    let mut best = (f64::INFINITY, 0.0, 0.0, 0.0);
    let mut seed: u64 = 0x9e3779b97f4a7c15 ^ (b.xn[0] as u64).wrapping_mul(31) ^ (b.un[0] as u64);
    let mut rnd = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let try_pose = |x: f64, y: f64, u: f64, best: &mut (f64, f64, f64, f64)| {
        if admissible_f(cv, x, y, u) {
            let m = float_mass(cv, x, y, u);
            if m < best.0 {
                *best = (m, x, y, u);
            }
        }
    };
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let x = x0 + (x1 - x0) * i as f64 / 2.0;
                let y = y0 + (y1 - y0) * j as f64 / 2.0;
                let u = u0 + (u1 - u0) * k as f64 / 2.0;
                try_pose(x, y, u, &mut best);
            }
        }
    }
    for _ in 0..nsamp {
        let x = x0 + (x1 - x0) * rnd();
        let y = y0 + (y1 - y0) * rnd();
        let u = u0 + (u1 - u0) * rnd();
        try_pose(x, y, u, &mut best);
    }
    best
}

// =====================================================================================
// Search
// =====================================================================================

#[derive(Clone)]
struct Settings {
    depth: u32,
    node_cap: usize,
    kappa: f64,
    uncert_cap: usize,
    tight: f64, // dump certified leaves with bound < (1 + tight) (diagnostics for the harness)
}

struct RootResult {
    id: usize,
    boxes: usize,
    cert: usize,
    empty: usize,
    uncert: Vec<PBox>,
    maxdepth: u32,
    capped: bool,
    ms: u128,
    tight: Vec<(PBox, I)>,
}

struct Node {
    b: PBox,
    depth: u32,
    inw: I,        // certain-in point weight (bound units)
    cand: Vec<u32>, // undecided points
}

fn initial_candidates(cv: &Cover, b: &PBox) -> Vec<u32> {
    let cd = b.cd() as f64;
    let d = cv.d as f64;
    let reach = 0.7072;
    let x0 = b.xn[0] as f64 / cd - reach;
    let x1 = b.xn[1] as f64 / cd + reach;
    let y0 = b.yn[0] as f64 / cd - reach;
    let y1 = b.yn[1] as f64 / cd + reach;
    let nb = cv.nb as i64;
    let bx0 = ((x0 * 10.0).floor() as i64 - 1).max(0);
    let bx1 = ((x1 * 10.0).floor() as i64 + 1).min(nb - 1);
    let by0 = ((y0 * 10.0).floor() as i64 - 1).max(0);
    let by1 = ((y1 * 10.0).floor() as i64 + 1).min(nb - 1);
    let mut v = Vec::new();
    for bx in bx0..=bx1 {
        for by in by0..=by1 {
            for &i in &cv.buckets[(bx * nb + by) as usize] {
                let (px, py, _) = cv.pts[i as usize];
                let (px, py) = (px as f64 / d, py as f64 / d);
                if px >= x0 && px <= x1 && py >= y0 && py <= y1 {
                    v.push(i);
                }
            }
        }
    }
    v
}

/// Evaluate a node: update point classes; return (bound, inw, remaining candidates).
fn eval_node(ctx: &Ctx, b: &PBox, inw: I, cand: &[u32]) -> (I, I, Vec<u32>, bool) {
    let cv = ctx.cv;
    let bg = box_geo(ctx, b);
    if box_empty(ctx, b, &bg) {
        return (0, inw, Vec::new(), true);
    }
    let cd = b.cd() as f64;
    let ud = b.ud() as f64;
    let d = cv.d as f64;
    let xf = [b.xn[0] as f64 / cd, b.xn[1] as f64 / cd];
    let yf = [b.yn[0] as f64 / cd, b.yn[1] as f64 / cd];
    let uf = [b.un[0] as f64 / ud, b.un[1] as f64 / ud];
    let mut inw = inw;
    let mut rest = Vec::with_capacity(cand.len());
    for &i in cand {
        let (px, py, w) = cv.pts[i as usize];
        match pt_float_class(px as f64 / d, py as f64 / d, xf, yf, uf) {
            0 => {}
            1 => {
                if pt_certain_in(cv.d, px as I, py as I, b) {
                    inw += w * ctx.ptw;
                } else {
                    rest.push(i);
                }
            }
            _ => rest.push(i),
        }
    }
    if inw >= ctx.target {
        return (inw, inw, rest, false);
    }
    let sb = mass_bound(ctx, b, &bg, inw);
    (inw + sb, inw, rest, false)
}

fn split(b: &PBox, kappa: f64) -> (PBox, PBox) {
    let cd = b.cd() as f64;
    let ud = b.ud() as f64;
    let wx = (b.xn[1] - b.xn[0]) as f64 / cd;
    let wy = (b.yn[1] - b.yn[0]) as f64 / cd;
    let wu = (b.un[1] - b.un[0]) as f64 / ud * kappa;
    let mut b = *b;
    if wu >= wx && wu >= wy {
        if (b.un[0] + b.un[1]) % 2 != 0 {
            b.un[0] *= 2;
            b.un[1] *= 2;
            b.ul += 1;
        }
        let m = (b.un[0] + b.un[1]) / 2;
        let mut l = b;
        let mut r = b;
        l.un[1] = m;
        r.un[0] = m;
        (l, r)
    } else {
        let dimx = wx >= wy;
        let need = if dimx { (b.xn[0] + b.xn[1]) % 2 != 0 } else { (b.yn[0] + b.yn[1]) % 2 != 0 };
        if need {
            for k in 0..2 {
                b.xn[k] *= 2;
                b.yn[k] *= 2;
            }
            b.cl += 1;
        }
        let mut l = b;
        let mut r = b;
        if dimx {
            let m = (b.xn[0] + b.xn[1]) / 2;
            l.xn[1] = m;
            r.xn[0] = m;
        } else {
            let m = (b.yn[0] + b.yn[1]) / 2;
            l.yn[1] = m;
            r.yn[0] = m;
        }
        (l, r)
    }
}

fn run_root(ctx: &Ctx, id: usize, root: PBox, st: &Settings) -> RootResult {
    let t0 = Instant::now();
    let mut res = RootResult {
        id,
        boxes: 0,
        cert: 0,
        empty: 0,
        uncert: Vec::new(),
        maxdepth: 0,
        capped: false,
        ms: 0,
        tight: Vec::new(),
    };
    let tight_thr = ctx.target + ((ctx.target as f64) * st.tight) as I;
    let cand0 = initial_candidates(ctx.cv, &root);
    let mut stack = vec![Node { b: root, depth: 0, inw: 0, cand: cand0 }];
    while let Some(nd) = stack.pop() {
        res.boxes += 1;
        res.maxdepth = res.maxdepth.max(nd.depth);
        if res.boxes > st.node_cap {
            res.capped = true;
            res.uncert.push(nd.b);
            for n in stack.drain(..) {
                res.uncert.push(n.b);
            }
            break;
        }
        let (bound, inw, rest, empty) = eval_node(ctx, &nd.b, nd.inw, &nd.cand);
        if empty {
            res.empty += 1;
            continue;
        }
        if bound >= ctx.target {
            res.cert += 1;
            if bound < tight_thr && res.tight.len() < 2000 {
                res.tight.push((nd.b, bound));
            }
            continue;
        }
        if nd.depth >= st.depth || nd.b.cl >= 28 || nd.b.ul >= 28 {
            res.uncert.push(nd.b);
            if res.uncert.len() >= st.uncert_cap {
                res.capped = true;
                for n in stack.drain(..) {
                    res.uncert.push(n.b);
                }
                break;
            }
            continue;
        }
        let (l, r) = split(&nd.b, st.kappa);
        stack.push(Node { b: r, depth: nd.depth + 1, inw, cand: rest.clone() });
        stack.push(Node { b: l, depth: nd.depth + 1, inw, cand: rest });
    }
    res.ms = t0.elapsed().as_millis();
    res
}

// =====================================================================================
// Lemma U: first order at theta = 0 (search/ZMX2_AREA.md sec 11)
// =====================================================================================
//
// For a box with u0 = 0 every pose P = (c, u) is written as P = (b + sigma (w(u) - 1)/2, u) with a *base* b
// (sigma = +1 / -1 per axis when the box reaches the low / high wall, 0 otherwise; then b is admissible
// at theta = 0).  Every piece of mu is bounded below at P by (its theta -> 0+ value at the base) + u * (a
// slope bound), and the infimum over the base of the sum of the base values is taken jointly (as Lemma
// Z0).  Line pairs at germs (Lemma Z) enter through inf over zeta, with their losses paid from the slope.
// Since 2026-10-01 (sec 11.5-11.7): everything is evaluated per cell of the base grid (rates, breakpoint
// corrections and the pair forms uniform on the cell, so the cell bound is concave and its minimum is at the
// cell's corners), and the poses near a germ pair can be split by its cone |eta| <= K u (Lemma V): inside, the
// perpendicular coordinate follows c = g + kappa u from the single base g; outside, the pair is uncut and one
// of its lines is an ordinary window line.

/// smallest and largest density of line l on the open grid interval (lo, hi) (0 off the support); F on
/// [e, e'] only depends on the density on (e, e')
fn dens_range(l: &Line, lo: I, hi: I) -> (I, I) {
    let n = l.bps.len();
    if n == 0 || hi <= lo {
        return (0, 0);
    }
    let mut mn: Option<I> = None;
    let mut mx: Option<I> = None;
    let mut upd = |v: I| {
        mn = Some(mn.map_or(v, |m: I| m.min(v)));
        mx = Some(mx.map_or(v, |m: I| m.max(v)));
    };
    let b0 = (l.bps[0] as I) << K;
    let bl = (l.bps[n - 1] as I) << K;
    if lo < b0 || hi > bl {
        upd(0);
    }
    for k in 0..n - 1 {
        let a = (l.bps[k] as I) << K;
        let bb = (l.bps[k + 1] as I) << K;
        if hi > a && lo < bb {
            upd(l.dens[k]);
        }
    }
    (mn.unwrap_or(0), mx.unwrap_or(0))
}

/// One line in Lemma U: its classification and the affine end bounds  lo <= b - 1/2 + al u,
/// hi >= b + 1/2 + be u  (b = base along-coordinate), when they hold.
struct ULine<'a> {
    l: &'a Line,
    lo_ok: bool,
    hi_ok: bool,
    al: f64,   // upper bound of the lo slope (per unit u, length units), path shift included
    be: f64,   // lower bound of the hi slope, path shift included
    al_g: f64, // the g-end parts alone (no path shift): g1 <= -1/2 + al_g u, g2 >= 1/2 + be_g u
    be_g: f64,
}

/// the path of one centre coordinate in Lemma U: `Path(sig)`: c = b + sig (w(u) - 1)/2 with a base b
/// (sec 11.2); `Cone`: c = g + kappa u with kappa in [k0, k1] and the single base g (the cone of a germ
/// pair, Lemma V, sec 11.5)
#[derive(Clone, Copy, Debug)]
enum AxisMode {
    Path(i32),
    Cone { g: I, k0: f64, k1: f64 },
}
/// one piece of one centre axis: the path and the centre range (grid units) whose poses it covers
#[derive(Clone, Copy, Debug)]
struct AxisSpec {
    mode: AxisMode,
    cr: (I, I),
}
/// bounds (upper, lower) of (c - b)/u along the path of an axis, u in (0, u1]; wlo = (1-u1)/(1+u1^2) rounded down
fn shift_rates(m: AxisMode, wlo: f64) -> (f64, f64) {
    match m {
        AxisMode::Path(1) => (1.0, wlo),
        AxisMode::Path(-1) => (-wlo, -1.0),
        AxisMode::Path(_) => (0.0, 0.0),
        AxisMode::Cone { k0, k1, .. } => (k1, k0),
    }
}

/// the end bounds of a line at frame position `pos` over the box (frame `fb`), u in [0, u1] (Lemma U (a));
/// `sh` = (upper, lower) bound of the shift rate of the along coordinate (shift_rates)
fn u_line<'a>(ctx: &Ctx, fb: &FrameBox, l: &'a Line, sh: (f64, f64), u1: f64) -> ULine<'a> {
    let dd = ctx.cv.d;
    let den = fb.xd * dd;
    let pos = l.pos;
    let mut f1 = Iv { lo: INF, hi: -INF };
    let mut f2 = f1;
    let mut dmin = INF;
    let mut dmax = -INF;
    for e in 0..2 {
        let dnum = fb.xn[e] * dd - (pos as I) * fb.xd;
        let dv = Iv::rat(dnum, den);
        dmin = dmin.min(dv.lo);
        dmax = dmax.max(dv.hi);
        let mm = Iv::rat(2 * dnum - den, 2 * den);
        let pm = Iv::rat(2 * dnum + den, 2 * den);
        f1 = f1.hull(h_range(mm.scale(0.5), pm.scale(-0.5), fb.u0, fb.u1, fb.u0zero));
        f2 = f2.hull(h_range(pm.scale(0.5), mm.scale(-0.5), fb.u0, fb.u1, fb.u0zero));
    }
    let u1i = Iv::exact(u1);
    let one_m = Iv::exact(1.0).sub(u1i.mul(u1i)); // 1 - u1^2 > 0
    // g1 = -dT - R/2 <= -1/2 + al_g u ;  g2 = -dT + R/2 >= 1/2 + be_g u   (T in [2u, 2u/(1-u1^2)], R >= 1)
    let al_g = if dmin >= 0.0 { Iv::exact(dmin).scale(-2.0).hi } else { Iv::exact(-dmin).scale(2.0).div_pos(one_m).hi };
    let be_g = if dmax <= 0.0 { Iv::exact(-dmax).scale(2.0).lo } else { Iv::exact(dmax).scale(-2.0).div_pos(one_m).lo };
    // domination of the f-ends (Lemma C): lo = c + max(f1, g1) <= c + g1-bound if f1 <= -1/2 + min(0, al_g u1)
    let lo_dom = f1.hi <= Iv::exact(-0.5).add(Iv::exact(al_g.min(0.0)).mul(u1i)).lo;
    let hi_dom = f2.lo >= Iv::exact(0.5).add(Iv::exact(be_g.max(0.0)).mul(u1i)).hi;
    // Lemma W (admissible poses): line at distance 1 from the low wall: H1 >= c + q, q >= 1/2 - u1 u;
    // from the high wall: L1 <= c - q
    let w_lo_wall = pos - fb.wall_l == dd as i64;
    let w_hi_wall = fb.wall_r - pos == dd as i64;
    // (an end that is not ok keeps its g-bound: in a pair that end is max(zeta, g-end) resp. min(zeta + u, g-end))
    let (lo_ok, al) = if lo_dom {
        (true, al_g)
    } else if w_hi_wall {
        (true, al_g.max(u1))
    } else {
        (false, al_g)
    };
    let (hi_ok, be) = if hi_dom {
        (true, be_g)
    } else if w_lo_wall {
        (true, be_g.min(-u1))
    } else {
        (false, be_g)
    };
    ULine { l, lo_ok, hi_ok, al: Iv::exact(al).add(Iv::exact(sh.0)).hi, be: Iv::exact(be).add(Iv::exact(sh.1)).lo, al_g, be_g }
}

/// densities left and right of the breakpoints of l strictly inside (lo, hi): Some((p, rho_left, rho_right))
/// if there is exactly one, None if none; Err if more
fn one_bp(l: &Line, lo: I, hi: I) -> Result<Option<(I, I, I)>, ()> {
    let n = l.bps.len();
    let mut found: Option<(I, I, I)> = None;
    for k in 0..n {
        let p = (l.bps[k] as I) << K;
        if p > lo && p < hi {
            if found.is_some() {
                return Err(());
            }
            let rl = if k == 0 { 0 } else { l.dens[k - 1] };
            let rr = if k + 1 == n { 0 } else { l.dens[k] };
            found = Some((p, rl, rr));
        }
    }
    Ok(found)
}

/// the sweep of an end moving at rate k for u <= u1, in grid units, rounded up: >= |k| u1 G
fn sweep(ctx: &Ctx, k: f64, u1: f64) -> I {
    ctx.ghi((k.abs() * u1 * 1.0000001).max(0.0)) + 1
}

/// Lemma U (b) for one end (sec 11.2, 11.6): the end e = b + off (b in [c0, c1], grid) moves to e + k u, u in
/// [0, u1].  hi = true (the hi end of a chord): F(e + k u) - F(e) >= k rho u G - corr(b);  hi = false (a lo end):
/// F(e + k u) - F(e) <= k rho u G + corr(b), with corr(b) >= 0 returned as a `Corr` term (its negative).  rho is the
/// smallest (gain) or largest (loss) density on the swept range; with `use_corr`, an end that moves the adverse
/// way across exactly one breakpoint p uses the density beyond p and corr(b) = (excess)+ min(dist(e, p)+, dmax).
fn end_rate<'a>(ctx: &Ctx, l: &Line, c0: I, c1: I, off: I, k: f64, u1: f64, hi: bool, use_corr: bool) -> (I, Option<UTerm<'a>>) {
    let sw = sweep(ctx, k, u1);
    let (e0, e1) = (c0 + off, c1 + off);
    let gain = if hi { k >= 0.0 } else { k <= 0.0 };
    if gain {
        let (lo, hi_) = if k >= 0.0 { (e0, e1 + sw) } else { (e0 - sw, e1) };
        return (dens_range(l, lo, hi_).0, None);
    }
    let (lo, hi_) = if k >= 0.0 { (e0, e1 + sw) } else { (e0 - sw, e1) };
    if use_corr {
        if let Ok(Some((p, rl, rr))) = one_bp(l, lo, hi_) {
            return if k < 0.0 {
                // moving left from the right piece (rr) into the left one (rl):
                // int_{e - d}^{e} rho <= rl d + (rr - rl)+ min((e - p)+, d);  (e - p)+ = (b - (p - off))+
                (rl, if rr > rl { Some(UTerm::Corr { coef: rr - rl, pb: p - off, sg: -1, dmax: sw }) } else { None })
            } else {
                // moving right from rl into rr: int_e^{e + d} rho <= rr d + (rl - rr)+ min((p - e)+, d)
                (rr, if rl > rr { Some(UTerm::Corr { coef: rl - rr, pb: p - off, sg: 1, dmax: sw }) } else { None })
            };
        }
    }
    (dens_range(l, lo, hi_).1, None)
}

/// k rho G * u1 rounded down (bound units), k rho G = the rate of an end per unit u
fn rate_u1(k: f64, rho: I, gf: f64, u1: f64) -> f64 {
    Iv::exact(k).mul(Iv::exact(rho as f64)).mul(Iv::exact(gf)).mul(Iv::exact(u1)).lo
}

/// Lemma U (b) on a base cell [c0, c1]: the window F(b + 1/2 + be u) - F(b - 1/2 + al u) >= W(b) + u s - corr(b):
/// returns s * u1 (bound units, rounded down) and the correction terms
fn u_slope<'a>(ctx: &Ctx, l: &Line, c0: I, c1: I, al: f64, be: f64, u1: f64, use_corr: bool) -> (f64, Vec<UTerm<'a>>) {
    let half = ctx.g / 2;
    let (rh, ch) = end_rate(ctx, l, c0, c1, half, be, u1, true, use_corr);
    let (rl, cl) = end_rate(ctx, l, c0, c1, -half, al, u1, false, use_corr);
    // (the lo-end term enters with a minus sign: rate_u1(-al, rl) <= -(al rl G u1), as rate_u1 rounds down)
    let s = Iv::exact(rate_u1(be, rh, ctx.gf, u1)).add(Iv::exact(rate_u1(-al, rl, ctx.gf, u1))).lo;
    (s, ch.into_iter().chain(cl).collect())
}

/// the sweep (grid units) of the zeta-dependent ends of a germ pair: a's lo end (rate al_a+) and b's hi end (rate
/// min(1, be_b))
fn pair_sweep(ctx: &Ctx, sl: (f64, f64, f64, f64), u1: f64) -> I {
    sweep(ctx, sl.0.max(0.0), u1).max(sweep(ctx, sl.3.min(1.0), u1))
}

/// Lemma U (c), per base cell (sec 11.2 (c), 11.6): a germ pair, a at l with its lo end cut by zeta, b at
/// l + 1 with its hi end cut by zeta + u; base along-coordinate b in the cell [c0, c1] (grid), H = b + 1/2,
/// L = b - 1/2.  For every zeta: a-term(u) >= max(0, X_a + u r_a), b-term(u) >= max(0, X_b + u r_b) with
/// X_a = F_a(H) - F_a(max(zeta, L)), X_b = F_b(min(zeta, H)) - F_b(L).  The zeta-axis is cut at the candidates
/// (b +- 1/2 and the breakpoints of both lines), whose order is the same at every b inside the cell; on each
/// zeta-cell the rates r and the choice affine/0 of each term are fixed for the whole base cell, so the bound
/// at each zeta-cell end is jointly affine in (b, u), and their minimum Pi(b, u) is jointly concave.  Returns
/// Pi at (c0, 0), (c1, 0), (c0, u1), (c1, u1) (without the end corrections) and the correction terms of the
/// two base-anchored ends (a's hi end, b's lo end).
fn u_pair_cell<'a>(ctx: &Ctx, a: &Line, bl: &Line, c0: I, c1: I, sl: (f64, f64, f64, f64), u1: f64, use_corr: bool) -> ([I; 2], [I; 2], Vec<UTerm<'a>>) {
    let half = ctx.g / 2;
    let gf = ctx.gf;
    let (al_a, be_a, al_b, be_b) = sl;
    let mut corr: Vec<UTerm> = Vec::new();
    let (rha, ca) = end_rate(ctx, a, c0, c1, half, be_a, u1, true, use_corr);
    let (rlb, cb) = end_rate(ctx, bl, c0, c1, -half, al_b, u1, false, use_corr);
    corr.extend(ca);
    corr.extend(cb);
    let ra_hi = rate_u1(be_a, rha, gf, u1);
    let rb_lo = rate_u1(-al_b, rlb, gf, u1);
    // candidates: C(y) a constant, R(o) the moving point b + o
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Zc {
        C(I),
        R(I),
    }
    let at = |z: Zc, b: I| -> I {
        match z {
            Zc::C(y) => y,
            Zc::R(o) => b + o,
        }
    };
    // (the breakpoints and their translates by -+ the sweep of the zeta-dependent ends, so that a zeta-cell's
    // swept range crosses no breakpoint it does not need to; the base grid has every candidate -+ 1/2)
    let psw = pair_sweep(ctx, sl, u1);
    let mut zs: Vec<Zc> = vec![Zc::R(-half), Zc::R(half)];
    for l in [a, bl] {
        for &bp in &l.bps {
            for y in [(bp as I) << K, ((bp as I) << K) - psw, ((bp as I) << K) + psw] {
                if y >= c0 - 3 * half && y <= c1 + 3 * half && !zs.contains(&Zc::C(y)) {
                    zs.push(Zc::C(y));
                }
            }
        }
    }
    // the order inside the open cell, by the value at its midpoint (doubled); ties only if c0 = c1
    zs.sort_by_key(|z| match *z {
        Zc::C(y) => 2 * y,
        Zc::R(o) => c0 + c1 + 2 * o,
    });
    let xa = |z: I, b: I| -> I {
        let (h, l) = (b + half, (b - half).max(z));
        if h > l {
            fval(a, h) - fval(a, l)
        } else {
            0
        }
    };
    let xb = |z: I, b: I| -> I {
        let (h, l) = (z.min(b + half), b - half);
        if h > l {
            fval(bl, h) - fval(bl, l)
        } else {
            0
        }
    };
    let far = 8 * half; // stands for -inf / +inf in the ranges of the sentinel cells
    let n = zs.len();
    let cs = [c0, c1];
    // pass 1: per zeta-cell the values of X_a, X_b at its ends and both vertices, the activity, the rates
    struct ZCell {
        vals: Vec<(I, I, usize)>, // (X_a, X_b, vertex)
        act_a: bool,
        act_b: bool,
        ra: f64,
        rb: f64,
    }
    let mut zcells: Vec<ZCell> = Vec::new();
    let mut pstar = [I::MAX; 2]; // the exact limit P at the two vertices: min of X_a + X_b
    for k in 0..=n {
        let lft = if k == 0 { None } else { Some(zs[k - 1]) };
        let rgt = if k == n { None } else { Some(zs[k]) };
        let ends: Vec<Zc> = [lft, rgt].iter().flatten().copied().collect();
        let mut vals: Vec<(I, I, usize)> = Vec::new();
        for (ci, &c) in cs.iter().enumerate() {
            for &z in &ends {
                let zz = at(z, c);
                let v = (xa(zz, c), xb(zz, c), ci);
                pstar[ci] = pstar[ci].min(v.0 + v.1);
                vals.push(v);
            }
        }
        let act_a = vals.iter().any(|v| v.0 > 0);
        let act_b = vals.iter().any(|v| v.1 > 0);
        // the zeta-range of this cell over the base cell
        let zlo = lft.map_or(at(zs[0], c0) - far, |z| at(z, c0));
        let zhi = rgt.map_or(at(zs[n - 1], c1) + far, |z| at(z, c1));
        let mut ra = ra_hi;
        if act_a && al_a > 0.0 {
            // a's lo end max(zeta, L + al_a u) <= max(zeta, L) + al_a u: loss <= al_a rho_max u
            let (m0, m1) = (zlo.max(c0 - half), zhi.max(c1 - half));
            let rho = dens_range(a, m0, m1 + sweep(ctx, al_a, u1)).1;
            ra = Iv::exact(ra).add(Iv::exact(rate_u1(-al_a, rho, gf, u1))).lo;
        }
        let mut rb = rb_lo;
        if act_b {
            // b's hi end min(zeta + u, H + be_b u) >= min(zeta, H) + kb u, kb = min(1, be_b)
            let kb = be_b.min(1.0);
            let (h0, h1) = (zlo.min(c0 + half), zhi.min(c1 + half));
            let sw = sweep(ctx, kb, u1);
            let rho = if kb >= 0.0 { dens_range(bl, h0, h1 + sw).0 } else { dens_range(bl, h0 - sw, h1).1 };
            rb = Iv::exact(rb).add(Iv::exact(rate_u1(kb, rho, gf, u1))).lo;
        }
        zcells.push(ZCell { vals, act_a, act_b, ra, rb });
    }
    // pass 2: the form of each term on each zeta-cell, fixed for the whole base cell and all u: affine (only where
    // the term is active, i.e. X >= 0 is its unclamped value) or 0 (always valid); among the forms that keep the
    // u = 0 values >= P at both vertices (so the limit stays exact), the one with the best u = u1 values
    let val_u1 = |zc: &ZCell, fa: bool, fbb: bool, v: &(I, I, usize)| -> (I, f64) {
        let b0 = (if fa { v.0 } else { 0 }) + (if fbb { v.1 } else { 0 });
        let mut t = Iv::exact(b0 as f64);
        if fa {
            t = t.add(Iv::exact(zc.ra));
        }
        if fbb {
            t = t.add(Iv::exact(zc.rb));
        }
        (b0, t.lo)
    };
    let mut v0 = [I::MAX; 2];
    let mut v1 = [f64::INFINITY; 2];
    for zc in &zcells {
        let mut best: Option<(f64, bool, bool)> = None;
        for fa in [false, true] {
            for fbb in [false, true] {
                if (fa && !zc.act_a) || (fbb && !zc.act_b) {
                    continue;
                }
                let ok = zc.vals.iter().all(|v| val_u1(zc, fa, fbb, v).0 >= pstar[v.2]);
                if !ok && !(fa == zc.act_a && fbb == zc.act_b) {
                    continue;
                }
                let w = zc.vals.iter().map(|v| val_u1(zc, fa, fbb, v).1).fold(f64::INFINITY, f64::min);
                if best.map_or(true, |(bw, _, _)| w > bw) {
                    best = Some((w, fa, fbb));
                }
            }
        }
        let (_, fa, fbb) = best.unwrap();
        for v in &zc.vals {
            let (b0, t) = val_u1(zc, fa, fbb, v);
            v0[v.2] = v0[v.2].min(b0);
            v1[v.2] = v1[v.2].min(t);
        }
    }
    // values below 2^80: a relative 2^-50 margin covers the f64 roundings of base as f64, then the floor
    let fl = |x: f64| -> I { dn(x - x.abs() * 1e-15 - 1.0).floor() as I };
    (v0, [fl(v1[0]), fl(v1[1])], corr)
}

/// Lemma U (d): for c in the box and theta in [0, theta1]: a lower bound `st` of d area(Q cap R)/d theta (per
/// radian) = sum over the edges of Q of the integral of (-s) over the part inside R (Reynolds; s = position along
/// the edge from its midpoint, counter-clockwise), counting s > 0 where possibly inside and s < 0 where certainly
/// inside; and enclosures of N = integral over the inside part of the outward normal (d area/dc = N).  64 pieces
/// per edge.
fn u_area_slope(b: &PBox, r: &Rect, d: I, c1: f64, s1: f64) -> (f64, Iv, Iv) {
    let cd = b.cd();
    let cx = Iv::rat(b.xn[0], cd).hull(Iv::rat(b.xn[1], cd));
    let cy = Iv::rat(b.yn[0], cd).hull(Iv::rat(b.yn[1], cd));
    let cc = Iv { lo: c1, hi: 1.0 };
    let ss = Iv { lo: 0.0, hi: s1 };
    let (rx0, rx1, ry0, ry1) = (Iv::rat(r.x0 as I, d), Iv::rat(r.x1 as I, d), Iv::rat(r.y0 as I, d), Iv::rat(r.y1 as I, d));
    let n = 64;
    let mut tot = 0.0f64;
    let mut nx = Iv::exact(0.0);
    let mut ny = Iv::exact(0.0);
    let h = Iv::exact(0.5);
    for e in 0..4 {
        // rotated outward normal: bottom (S, -C), right (C, S), top (-S, C), left (-C, -S)
        let (nnx, nny) = match e {
            0 => (ss, cc.neg()),
            1 => (cc, ss),
            2 => (ss.neg(), cc),
            _ => (cc.neg(), ss.neg()),
        };
        for k in 0..n {
            let sa = -0.5 + k as f64 / n as f64;
            let sb = -0.5 + (k + 1) as f64 / n as f64; // exact binary fractions
            let s = Iv { lo: sa, hi: sb };
            // body-frame point n/2 + s tau: bottom (s, -1/2), right (1/2, s), top (-s, 1/2), left (-1/2, -s)
            let (vx, vy) = match e {
                0 => (s, h.neg()),
                1 => (h, s),
                2 => (s.neg(), h),
                _ => (h.neg(), s.neg()),
            };
            let px = cx.add(cc.mul(vx)).sub(ss.mul(vy));
            let py = cy.add(ss.mul(vx)).add(cc.mul(vy));
            let surely = px.lo >= rx0.hi && px.hi <= rx1.lo && py.lo >= ry0.hi && py.hi <= ry1.lo;
            let maybe = px.hi >= rx0.lo && px.lo <= rx1.hi && py.hi >= ry0.lo && py.lo <= ry1.hi;
            // integral of -s over [sa, sb] = (sa^2 - sb^2)/2 (exact in binary64 for these sa, sb)
            let integ = (sa * sa - sb * sb) * 0.5;
            if (sb <= 0.0 && surely) || (sa >= 0.0 && maybe) {
                tot = dn(tot + integ);
            }
            let ds = 1.0 / n as f64;
            if surely {
                nx = nx.add(nnx.scale(ds));
                ny = ny.add(nny.scale(ds));
            } else if maybe {
                nx = nx.add(nnx.scale(ds).hull(Iv::exact(0.0)));
                ny = ny.add(nny.scale(ds).hull(Iv::exact(0.0)));
            }
        }
    }
    (tot, nx, ny)
}

/// a u-independent term of the base function of one coordinate (Lemma U): a line window; an edge term
/// rho min(G, max(0, sg (b - e)) / s1) (sg = -1: left-type line, e = l + 1/2; sg = +1: right-type, e = l - 1/2);
/// a breakpoint correction
enum UTerm<'a> {
    Window(&'a Line),
    Edge { rho: I, e: I, sg: I },
    /// minus coef * min(max(0, sg (pb - b)), dmax)
    Corr { coef: I, pb: I, sg: I, dmax: I },
}

impl<'a> UTerm<'a> {
    fn eval(&self, b: I, half: I, g: I, s1: f64) -> I {
        match self {
            UTerm::Window(l) => fval(l, b + half) - fval(l, b - half),
            UTerm::Edge { rho, e, sg } => {
                let t = sg * (b - e);
                if t <= 0 || *rho == 0 {
                    return 0;
                }
                // chord >= min(1, tau/s1) (length), i.e. min(G, t/s1) grid units, rounded down
                let lenf = dn(t as f64 / up(s1));
                let len = if lenf >= g as f64 { g } else { lenf.floor().max(0.0) as I };
                rho * len
            }
            UTerm::Corr { coef, pb, sg, dmax } => -coef * (sg * (pb - b)).max(0).min(*dmax),
        }
    }
    /// kinks of the term (grid points where its slope changes)
    fn kinks(&self, half: I, g: I, s1: f64, out: &mut Vec<I>) {
        match self {
            UTerm::Window(l) => {
                for &bp in &l.bps {
                    let y = (bp as I) << K;
                    out.push(y - half);
                    out.push(y + half);
                }
            }
            UTerm::Edge { e, sg, .. } => {
                let w = up(g as f64 * up(s1)).ceil() as I;
                for k in [0, w - 1, w, w + 1] {
                    out.push(e + sg * k);
                }
            }
            UTerm::Corr { pb, sg, dmax, .. } => {
                out.push(*pb);
                out.push(*pb - sg * dmax);
            }
        }
    }
}

/// a germ pair in Lemma U: along axis, the two lines, their slopes (al_a, be_a, al_b, be_b), and the edge terms of
/// the two lines on the perpendicular axis (Lemma U (e)), used instead of the pair value when chosen
struct UPair<'a> {
    along: usize,
    a: &'a Line,
    bl: &'a Line,
    sl: (f64, f64, f64, f64),
    edges: Vec<UTerm<'a>>,
}

/// the data of one base cell of one axis for one choice of the corrections: values of the u-independent terms at
/// its two vertices (u = 0 and u = u1 differ only through the pairs), slope * u1, pair values (u = 0, u = u1)
struct UVar {
    base: [I; 2],
    su1: f64,
    pv0: Vec<[I; 2]>,
    pv1: Vec<[I; 2]>,
}

/// Lemma U for one choice of the axis pieces (paths and centre ranges) and the pair-region overrides
/// `ovr` = (family, frame position of line a, region: 1 = a is a window line, 2 = b is) of Lemma V.
/// See ZMX2_AREA.md sec 11.  The bound is the minimum over the 2-d base cells of the best over the
/// correction variants of each 1-d cell and the pair/edge choices (each a valid bound for the poses whose
/// base lies in that cell).
fn first_order_sig(ctx: &Ctx, b: &PBox, bg: &BoxGeo, inw: I, ax: [AxisSpec; 2], ovr: &[(usize, i64, u8)]) -> Option<I> {
    let cv = ctx.cv;
    let g = ctx.g;
    let half = g / 2;
    let cd = b.cd();
    let u1 = bg.u1.hi;
    let u1i = Iv::exact(u1);
    let s1 = cs_at(bg.u1).1.hi; // S(u1) >= p = S(u) on the box (theta <= 45 deg)
    let whi = w_hi(bg);
    let ext = ctx.ghi(Iv::exact(whi).sub(Iv::exact(1.0)).scale(0.5).hi).max(0); // >= (w - 1)/2 on the box
    let sg = cv.sx << K;
    let wlo = Iv::exact(1.0).sub(u1i).div_pos(Iv::exact(1.0).add(u1i.mul(u1i))).lo; // (w - 1)/2 >= wlo u
    // base ranges per axis (grid units): Path: b = c - sig (w(u) - 1)/2 (bases of admissible poses are >= 1/2
    // from the walls); Cone: the single base g
    let base = |s: &AxisSpec| -> (I, I) {
        match s.mode {
            AxisMode::Path(1) => ((s.cr.0 - ext).max(half), s.cr.1),
            AxisMode::Path(-1) => (s.cr.0, (s.cr.1 + ext).min(sg - half)),
            AxisMode::Path(_) => s.cr,
            AxisMode::Cone { g, .. } => (g, g),
        }
    };
    let br = [base(&ax[0]), base(&ax[1])];
    if br[0].0 > br[0].1 || br[1].0 > br[1].1 {
        return None;
    }
    let sh = [shift_rates(ax[0].mode, wlo), shift_rates(ax[1].mode, wlo)];
    let (fv, fh) = frames(ctx, b, bg);
    let debug = std::env::var("ZMX2_DEBUG").is_ok();
    let mut konst: I = inw; // points and certain atoms
    // per axis: window lines (line, al, be) and u-independent terms (edges)
    let mut wins: [Vec<(&Line, f64, f64)>; 2] = [Vec::new(), Vec::new()];
    let mut terms: [Vec<UTerm>; 2] = [Vec::new(), Vec::new()];
    let mut pairs: Vec<UPair> = Vec::new();
    for fam in 0..2 {
        // vertical lines: frame = identity, perp x (axis 0), along y (axis 1); horizontal lines: frame x' = -y,
        // perp axis 1, along x (axis 0)
        let (fb, lines, along, perp) = if fam == 0 { (&fv, &cv.vl, 1usize, 0usize) } else { (&fh, &cv.hl, 0usize, 1usize) };
        let xlo = (fb.xn[0] as f64) / (fb.xd as f64) - 0.75;
        let xhi = (fb.xn[1] as f64) / (fb.xd as f64) + 0.75;
        let mut us: Vec<ULine> = lines
            .iter()
            .filter(|l| {
                let p = l.pos as f64 / cv.d as f64;
                p >= xlo && p <= xhi
            })
            .map(|l| u_line(ctx, fb, l, sh[along], u1))
            .collect();
        // Lemma V regions: line a resp. b of a germ pair is a window line there
        for &(f, pa, reg) in ovr {
            if f != fam {
                continue;
            }
            for ul in us.iter_mut() {
                if reg == 1 && ul.l.pos == pa {
                    ul.lo_ok = true;
                }
                if reg == 2 && ul.l.pos == pa + cv.d as i64 {
                    ul.hi_ok = true;
                }
            }
        }
        // the perp axis's path slopes in the frame (cone slack of the edge terms, sec 11.5): frame perp = real
        // perp (vertical) or minus it (horizontal)
        let slack = |left: bool| -> I {
            match ax[perp].mode {
                AxisMode::Cone { k0, k1, .. } => {
                    let (q0, q1) = if fam == 0 { (k0, k1) } else { (-k1, -k0) };
                    let v = if left { Iv::exact(q1).sub(Iv::exact(wlo)).hi } else { Iv::exact(-q0).sub(Iv::exact(wlo)).hi };
                    if v > 0.0 {
                        ctx.ghi(Iv::exact(v).mul(u1i).hi) + 1
                    } else {
                        0
                    }
                }
                _ => 0,
            }
        };
        // edge term of a line whose lo (left-type) / hi (right-type) end is the cut by an edge of Q (Lemma U (e)):
        // its chord has length >= min(1, tau / S(u1)) with tau >= 1/2 + l - b_perp (left) or 1/2 + b_perp - l
        // (right) in the frame (minus the cone slack), at a density >= rho_min over where the chord can lie;
        // returned on the real perp axis
        let edge_term = |ul: &ULine, left: bool| -> Option<UTerm> {
            let dd = cv.d;
            let d0 = fb.xn[0] * dd - (ul.l.pos as I) * fb.xd;
            let d1 = fb.xn[1] * dd - (ul.l.pos as I) * fb.xd;
            if (left && d0 < 0) || (!left && d1 > 0) {
                return None; // the centre must stay on one side (|d| = d resp. -d)
            }
            let (ylo, yhi) = (ctx.glo(add_lo(fb.y0.lo, -whi * 0.5)), ctx.ghi(add_hi(fb.y1.hi, whi * 0.5)));
            let (rho, _) = dens_range(ul.l, ylo, yhi);
            if rho == 0 {
                return None;
            }
            let lg = (ul.l.pos as I) << K;
            let (e, sgn) = if left { (lg + half - slack(true), -1) } else { (lg - half + slack(false), 1) };
            // frame -> real perp coordinate: vertical identity; horizontal x' = -y
            Some(if fam == 0 { UTerm::Edge { rho, e, sg: sgn } } else { UTerm::Edge { rho, e: -e, sg: -sgn } })
        };
        let mut used = vec![false; us.len()];
        for i in 0..us.len() {
            let ul = &us[i];
            for a in line_atoms(ctx, fb, ul.l) {
                if a.c.iter().all(|&c| c) {
                    konst += a.w;
                }
            }
            if ul.lo_ok && ul.hi_ok {
                used[i] = true;
                wins[along].push((ul.l, ul.al, ul.be));
            }
        }
        // germ pairs (Lemma Z): a at pos with its lo end cut, b at pos + D with its hi end cut
        for i in 0..us.len() {
            if used[i] || !(us[i].hi_ok && !us[i].lo_ok) {
                continue;
            }
            if let Some(j) = (0..us.len()).find(|&j| !used[j] && us[j].l.pos == us[i].l.pos + cv.d as i64 && us[j].lo_ok && !us[j].hi_ok) {
                used[i] = true;
                used[j] = true;
                let (a, bl) = (&us[i], &us[j]);
                let edges: Vec<UTerm> = [edge_term(a, true), edge_term(bl, false)].into_iter().flatten().collect();
                pairs.push(UPair { along, a: a.l, bl: bl.l, sl: (a.al, a.be, bl.al, bl.be), edges });
            }
        }
        // unpaired edge lines
        for i in 0..us.len() {
            if used[i] {
                continue;
            }
            let ul = &us[i];
            let et = if ul.hi_ok && !ul.lo_ok {
                edge_term(ul, true)
            } else if ul.lo_ok && !ul.hi_ok {
                edge_term(ul, false)
            } else {
                None
            };
            if let Some(t) = et {
                terms[perp].push(t);
            }
        }
    }
    // rectangles: base area in the grid minimum; slope by Reynolds and the path motion (Lemma U (d))
    let dd = cv.d;
    let (c1i, s1i) = cs_at(bg.u1);
    // w'(t)/2 = (1 - 2t - t^2)/(1+t^2)^2 is decreasing on [0, u1]: c'(t) = sig w'(t)/2 in sig [w'(u1)/2, 1]
    let wp1 = Iv::exact(1.0).sub(u1i.scale(2.0)).sub(u1i.mul(u1i)).div_pos(Iv::exact(1.0).add(u1i.mul(u1i)).mul(Iv::exact(1.0).add(u1i.mul(u1i)))).lo;
    let cprime = |m: AxisMode| -> Iv {
        match m {
            AxisMode::Path(1) => Iv { lo: wp1, hi: 1.0 },
            AxisMode::Path(-1) => Iv { lo: -1.0, hi: -wp1 },
            AxisMode::Path(_) => Iv::exact(0.0),
            AxisMode::Cone { k0, k1, .. } => Iv { lo: k0, hi: k1 },
        }
    };
    let mut area_su1 = 0.0f64; // area slope * u1, bound units, rounded down
    for (ri, r) in cv.rects.iter().enumerate() {
        // enclosures over the paths: the centre box extended by ext towards the base (cone paths stay in the box)
        let mut eb = *b;
        let extl = Iv::exact(whi).sub(Iv::exact(1.0)).scale(0.5).hi;
        let kx = (extl * cd as f64).ceil() as I + 1;
        for (axi, s) in ax.iter().enumerate() {
            match s.mode {
                AxisMode::Path(1) => {
                    if axi == 0 {
                        eb.xn[0] -= kx
                    } else {
                        eb.yn[0] -= kx
                    }
                }
                AxisMode::Path(-1) => {
                    if axi == 0 {
                        eb.xn[1] += kx
                    } else {
                        eb.yn[1] += kx
                    }
                }
                _ => {}
            }
        }
        let (st, nx, ny) = u_area_slope(&eb, r, dd, c1i.lo, s1i.hi);
        // d area/dt = theta'(t) dA/dtheta + c'(t) . N,  theta' in [2/(1+u1^2), 2]
        let thp = Iv { lo: Iv::exact(2.0).div_pos(Iv::exact(1.0).add(u1i.mul(u1i))).lo, hi: 2.0 };
        let du = thp.mul(Iv::exact(st)).add(cprime(ax[0].mode).mul(nx)).add(cprime(ax[1].mode).mul(ny)).lo;
        if du == 0.0 {
            continue;
        }
        let tr = if du < 0.0 { up((ctx.rect_tr[ri] + 1) as f64) } else { dn(ctx.rect_tr[ri] as f64) };
        if debug {
            eprintln!("    U area dA/dtheta >= {:.6}, N_x {:?}, N_y {:?}, dA/du >= {:.6}", st, nx, ny, du);
        }
        area_su1 = Iv::exact(area_su1).add(Iv::exact(du).mul(Iv::exact(tr)).mul(u1i)).lo;
    }
    // the base grid per axis: kinks of every term, correction kinks, rectangle sides -+ 1/2, range ends
    let mut grid: [Vec<I>; 2] = [vec![br[0].0, br[0].1], vec![br[1].0, br[1].1]];
    for axi in 0..2 {
        let out = &mut grid[axi];
        for t in &terms[axi] {
            t.kinks(half, g, s1, out);
        }
        for p in pairs.iter() {
            for t in &p.edges {
                // edge terms of a pair live on its perp axis
                if p.along != axi {
                    t.kinks(half, g, s1, out);
                }
            }
        }
        for &(l, al, be) in &wins[axi] {
            for &bp in &l.bps {
                let y = (bp as I) << K;
                for v in [y - half, y + half, y + half - sweep(ctx, al, u1), y - half + sweep(ctx, be, u1)] {
                    out.push(v);
                }
            }
        }
        for p in pairs.iter().filter(|p| p.along == axi) {
            for (l, k, hi) in [(p.a, p.sl.1, true), (p.bl, p.sl.2, false)] {
                let sw = sweep(ctx, k, u1);
                let psw = pair_sweep(ctx, p.sl, u1);
                for ln in [p.a, p.bl] {
                    for &bp in &ln.bps {
                        let y = (bp as I) << K;
                        for yy in [y, y - psw, y + psw] {
                            out.push(yy - half);
                            out.push(yy + half);
                        }
                    }
                }
                for &bp in &l.bps {
                    let y = (bp as I) << K;
                    out.push(if hi { y - half + sw } else { y + half - sw });
                }
            }
        }
        for r in &cv.rects {
            let es = if axi == 0 { [r.x0, r.x1] } else { [r.y0, r.y1] };
            for e in es {
                let e = (e as I) << K;
                out.push(e - half);
                out.push(e + half);
            }
        }
        let (lo, hi) = br[axi];
        out.retain(|&v| v >= lo && v <= hi);
        out.sort();
        out.dedup();
    }
    // per axis and base cell: the correction variants (all off / all on / per line the better one)
    let npairs = pairs.len();
    let mut cells: [Vec<(usize, usize, Vec<UVar>)>; 2] = [Vec::new(), Vec::new()];
    for axi in 0..2 {
        let xs = &grid[axi];
        let idx: Vec<(usize, usize)> = if xs.len() == 1 { vec![(0, 0)] } else { (0..xs.len() - 1).map(|i| (i, i + 1)).collect() };
        for (i0, i1) in idx {
            let (c0, c1) = (xs[i0], xs[i1]);
            let cs = [c0, c1];
            let mut fixed = [0 as I; 2]; // windows and edges at the two vertices
            for t in &terms[axi] {
                for k in 0..2 {
                    fixed[k] += t.eval(cs[k], half, g, s1);
                }
            }
            // per window line: (off, on) = (su1, corr values at the vertices)
            let mut opts: Vec<[(f64, [I; 2]); 2]> = Vec::new();
            for &(l, al, be) in &wins[axi] {
                let w = UTerm::Window(l);
                for k in 0..2 {
                    fixed[k] += w.eval(cs[k], half, g, s1);
                }
                let mut o = [(0.0, [0 as I; 2]); 2];
                for (j, uc) in [false, true].iter().enumerate() {
                    let (s, corr) = u_slope(ctx, l, c0, c1, al, be, u1, *uc);
                    let mut cvv = [0 as I; 2];
                    for t in &corr {
                        for k in 0..2 {
                            cvv[k] += t.eval(cs[k], half, g, s1);
                        }
                    }
                    o[j] = (s, cvv);
                }
                opts.push(o);
            }
            // per pair along this axis: (off, on) = (pv0, pv1) with the corrections included
            let mut popts: Vec<[([I; 2], [I; 2]); 2]> = Vec::new();
            for p in pairs.iter() {
                if p.along != axi {
                    popts.push([([0; 2], [0; 2]); 2]);
                    continue;
                }
                let mut o = [([0 as I; 2], [0 as I; 2]); 2];
                for (j, uc) in [false, true].iter().enumerate() {
                    let (p0, p1, corr) = u_pair_cell(ctx, p.a, p.bl, c0, c1, p.sl, u1, *uc);
                    let mut cvv = [0 as I; 2];
                    for t in &corr {
                        for k in 0..2 {
                            cvv[k] += t.eval(cs[k], half, g, s1);
                        }
                    }
                    o[j] = ([p0[0] + cvv[0], p0[1] + cvv[1]], [p1[0] + cvv[0], p1[1] + cvv[1]]);
                }
                popts.push(o);
            }
            let worst = |s: f64, c: [I; 2]| -> f64 { (c[0].min(c[1]) as f64) + s.min(0.0) };
            let pworst = |v: &([I; 2], [I; 2])| -> I { v.0[0].min(v.0[1]).min(v.1[0]).min(v.1[1]) };
            let mut vars: Vec<UVar> = Vec::new();
            let mut seen: Vec<(Vec<bool>, Vec<bool>)> = Vec::new();
            for mode in 0..3 {
                let wch: Vec<bool> = opts
                    .iter()
                    .map(|o| match mode {
                        0 => false,
                        1 => true,
                        _ => worst(o[1].0, o[1].1) > worst(o[0].0, o[0].1),
                    })
                    .collect();
                let pch: Vec<bool> = popts
                    .iter()
                    .map(|o| match mode {
                        0 => false,
                        1 => true,
                        _ => pworst(&o[1]) > pworst(&o[0]),
                    })
                    .collect();
                // a choice that changes nothing is the same variant
                let wch: Vec<bool> = wch.iter().zip(opts.iter()).map(|(&c, o)| c && (o[1].1 != [0, 0] || o[1].0 != o[0].0)).collect();
                let pch: Vec<bool> = pch.iter().zip(popts.iter()).map(|(&c, o)| c && o[1] != o[0]).collect();
                if seen.contains(&(wch.clone(), pch.clone())) {
                    continue;
                }
                seen.push((wch.clone(), pch.clone()));
                let mut bs = fixed;
                let mut su1 = 0.0f64;
                for (o, &c) in opts.iter().zip(wch.iter()) {
                    let (s, cvv) = o[c as usize];
                    su1 = dn(su1 + s);
                    bs[0] += cvv[0];
                    bs[1] += cvv[1];
                }
                let pv0: Vec<[I; 2]> = popts.iter().zip(pch.iter()).map(|(o, &c)| o[c as usize].0).collect();
                let pv1: Vec<[I; 2]> = popts.iter().zip(pch.iter()).map(|(o, &c)| o[c as usize].1).collect();
                vars.push(UVar { base: bs, su1, pv0, pv1 });
            }
            cells[axi].push((i0, i1, vars));
        }
    }
    // edge alternatives of the pairs, on their perp axis, at every grid vertex
    let evals: Vec<Vec<I>> = pairs
        .iter()
        .map(|p| {
            let perp = 1 - p.along;
            grid[perp].iter().map(|&v| p.edges.iter().map(|t| t.eval(v, half, g, s1)).sum()).collect()
        })
        .collect();
    // pair choices: bit set = pair value, clear = edge terms
    let masks: Vec<u32> = if npairs <= 4 {
        (0..(1u32 << npairs)).collect()
    } else {
        let all = (1u32 << npairs.min(31)) - 1;
        let mut m = vec![all, 0];
        for k in 0..npairs.min(31) {
            m.push(all ^ (1 << k));
            m.push(1 << k);
        }
        m
    };
    let ov = |c: I, r0: I, r1: I| -> I { ((c + half).min(r1) - (c - half).max(r0)).max(0) };
    let area_at = |x: I, y: I| -> Option<I> {
        let mut v: I = 0;
        for (ri, r) in cv.rects.iter().enumerate() {
            let ox = ov(x, (r.x0 as I) << K, (r.x1 as I) << K);
            let oy = ov(y, (r.y0 as I) << K, (r.y1 as I) << K);
            if ox == 0 || oy == 0 {
                continue;
            }
            let pn = ((r.x1 - r.x0) as I) * ((r.y1 - r.y0) as I);
            v += ctx.rect_wlc[ri].checked_mul(ox)?.checked_mul(oy)? / (pn << K);
        }
        Some(v)
    };
    let mut gmin: Option<I> = None;
    let mut arg = (0usize, 0usize, 0u32);
    for (cx, &(ix0, ix1, ref vxs)) in cells[0].iter().enumerate() {
        for (cy, &(iy0, iy1, ref vys)) in cells[1].iter().enumerate() {
            let xsv = [ix0, ix1];
            let ysv = [iy0, iy1];
            let mut areas = [[0 as I; 2]; 2];
            for a in 0..2 {
                for bb in 0..2 {
                    areas[a][bb] = area_at(grid[0][xsv[a]], grid[1][ysv[bb]])?;
                }
            }
            let mut best: Option<(I, u32)> = None;
            for vx in vxs {
                for vy in vys {
                    let su1 = Iv::exact(vx.su1).add(Iv::exact(vy.su1)).add(Iv::exact(area_su1)).lo.floor() as I;
                    for &m in &masks {
                        let mut cell_min: Option<I> = None;
                        for a in 0..2 {
                            for bb in 0..2 {
                                let (mut s0, mut s1v) = (vx.base[a] + vy.base[bb], vx.base[a] + vy.base[bb] + su1);
                                for (k, p) in pairs.iter().enumerate() {
                                    if (m >> k) & 1 == 1 {
                                        let (vv, c) = if p.along == 0 { (vx, a) } else { (vy, bb) };
                                        s0 += vv.pv0[k][c];
                                        s1v += vv.pv1[k][c];
                                    } else {
                                        let e = if p.along == 0 { evals[k][ysv[bb]] } else { evals[k][xsv[a]] };
                                        s0 += e;
                                        s1v += e;
                                    }
                                }
                                let v = s0.min(s1v) + areas[a][bb];
                                cell_min = Some(cell_min.map_or(v, |x| x.min(v)));
                            }
                        }
                        let cmv = cell_min.unwrap();
                        if best.map_or(true, |(x, _)| cmv > x) {
                            best = Some((cmv, m));
                        }
                    }
                }
            }
            let (bv, bm) = best.unwrap();
            if gmin.map_or(true, |x| bv < x) {
                gmin = Some(bv);
                arg = (cx, cy, bm);
            }
        }
    }
    let tot = konst + gmin?;
    if debug && std::env::var("ZMX2_DEBUG_CELL").is_ok() {
        let (cx, cy, bm) = arg;
        let tu = ctx.target as f64;
        let (ix0, ix1, vxs) = &cells[0][cx];
        let (iy0, iy1, vys) = &cells[1][cy];
        for (nvx, vx) in vxs.iter().enumerate() {
            for (nvy, vy) in vys.iter().enumerate() {
                eprintln!("    variant x{} y{}: su1 x {:.3e} y {:.3e} area {:.3e}", nvx, nvy, vx.su1 / tu, vy.su1 / tu, area_su1 / tu);
                for (a, &ix) in [*ix0, *ix1].iter().enumerate() {
                    for (bb, &iy) in [*iy0, *iy1].iter().enumerate() {
                        let mut line = format!("      corner x {:.9} y {:.9}: base x {:.9} y {:.9} area {:.9}", grid[0][ix] as f64 / ctx.gf, grid[1][iy] as f64 / ctx.gf, vx.base[a] as f64 / tu, vy.base[bb] as f64 / tu, area_at(grid[0][ix], grid[1][iy]).unwrap() as f64 / tu);
                        for (k, p) in pairs.iter().enumerate() {
                            let (vv, c) = if p.along == 0 { (vx, a) } else { (vy, bb) };
                            let e = if p.along == 0 { evals[k][iy] } else { evals[k][ix] };
                            line += &format!(" | pair{} along {} P0 {:.9} P1 {:.9} E {:.9}", k, p.along, vv.pv0[k][c] as f64 / tu, vv.pv1[k][c] as f64 / tu, e as f64 / tu);
                        }
                        eprintln!("{}", line);
                    }
                }
            }
        }
        let _ = bm;
    }
    if debug {
        let (cx, cy, bm) = arg;
        let (ix0, ix1, _) = &cells[0][cx];
        let (iy0, iy1, _) = &cells[1][cy];
        eprintln!(
            "  lemma U {:?}/{:?} ovr {:?}: base x [{:.9},{:.9}] y [{:.9},{:.9}] pairs {} windows {}/{} area_su1 {:.3e}; worst cell x [{:.9},{:.9}] y [{:.9},{:.9}] mask {:b} -> {:.9}",
            ax[0].mode,
            ax[1].mode,
            ovr,
            br[0].0 as f64 / ctx.gf,
            br[0].1 as f64 / ctx.gf,
            br[1].0 as f64 / ctx.gf,
            br[1].1 as f64 / ctx.gf,
            npairs,
            wins[0].len(),
            wins[1].len(),
            area_su1 / ctx.target as f64,
            grid[0][*ix0] as f64 / ctx.gf,
            grid[0][*ix1] as f64 / ctx.gf,
            grid[1][*iy0] as f64 / ctx.gf,
            grid[1][*iy1] as f64 / ctx.gf,
            bm,
            tot as f64 / ctx.target as f64
        );
    }
    Some(tot)
}

/// a germ pair whose cone (Lemma V) can split the poses of the box: family, frame position of line a, the
/// germ coordinate g (real, grid units; the centre coordinate where both lines are at distance 1/2) and the
/// cone slopes (k0, k1) of that real axis
struct GermSplit {
    fam: usize,
    pos_a: i64,
    g: I,
    k0: f64,
    k1: f64,
}

/// Lemma V candidates: germ pairs of the box with the germ inside the box's centre range on their perp axis
fn germ_splits(ctx: &Ctx, b: &PBox, bg: &BoxGeo, cr: [(I, I); 2]) -> Vec<GermSplit> {
    let cv = ctx.cv;
    let u1 = bg.u1.hi;
    let u1i = Iv::exact(u1);
    let (fv, fh) = frames(ctx, b, bg);
    let mut out = Vec::new();
    let half = ctx.g / 2;
    for fam in 0..2 {
        let (fb, lines, perp) = if fam == 0 { (&fv, &cv.vl, 0usize) } else { (&fh, &cv.hl, 1usize) };
        let xlo = (fb.xn[0] as f64) / (fb.xd as f64) - 0.75;
        let xhi = (fb.xn[1] as f64) / (fb.xd as f64) + 0.75;
        let us: Vec<ULine> = lines
            .iter()
            .filter(|l| {
                let p = l.pos as f64 / cv.d as f64;
                p >= xlo && p <= xhi
            })
            .map(|l| u_line(ctx, fb, l, (0.0, 0.0), u1))
            .collect();
        for a in us.iter().filter(|a| a.hi_ok && !a.lo_ok) {
            if let Some(bl) = us.iter().find(|bl| bl.l.pos == a.l.pos + cv.d as i64 && bl.lo_ok && !bl.hi_ok) {
                // K_a = 1 + 2 max(0, -al_g(a)) u1, K_b = 1 + 2 max(0, be_g(b)) u1 (rounded up)
                let ka = Iv::exact(1.0).add(Iv::exact((-a.al_g).max(0.0)).scale(2.0).mul(u1i)).hi;
                let kb = Iv::exact(1.0).add(Iv::exact(bl.be_g.max(0.0)).scale(2.0).mul(u1i)).hi;
                let pa = (a.l.pos as I) << K;
                // frame germ g' = l_a + 1/2; real: vertical g = g', horizontal g_y = -g'
                let (g, k0, k1) = if fam == 0 { (pa + half, -ka, kb) } else { (-pa - half, -kb, ka) };
                let (c0, c1) = cr[perp];
                // eta = d_a - 1/2 must stay in (-1, 1) on the box: |c - g| < 1
                if g >= c0 && g <= c1 && c0 > g - ctx.g && c1 < g + ctx.g {
                    out.push(GermSplit { fam, pos_a: a.l.pos, g, k0, k1 });
                }
            }
        }
    }
    out
}

/// Lemma U (+ Lemma V): the first-order bound of a box with u0 = 0 and u1 <= 1/16 (None if not applicable).
/// `inw` = the certain point weight.  First the plain bound (best over the admissible path signs: forced at a
/// wall, free otherwise); if it does not reach the target, the poses are split by the cones of germ pairs (one
/// per axis, on axes whose path sign is free) and the bound is the minimum over the pieces of the best over the
/// path signs of each piece.
fn first_order_bound(ctx: &Ctx, b: &PBox, bg: &BoxGeo, inw: I) -> Option<I> {
    if b.un[0] != 0 || 16 * b.un[1] > b.ud() {
        return None;
    }
    let g = ctx.g;
    let half = g / 2;
    let cd = b.cd();
    let whi = w_hi(bg);
    let ext = ctx.ghi(Iv::exact(whi).sub(Iv::exact(1.0)).scale(0.5).hi).max(0);
    let sg = ctx.cv.sx << K;
    let cr = [(fdiv(b.xn[0] * g, cd), cdiv(b.xn[1] * g, cd)), (fdiv(b.yn[0] * g, cd), cdiv(b.yn[1] * g, cd))];
    let opts = |r: (I, I)| -> Vec<i32> {
        if r.0 < half + ext {
            vec![1] // near the low wall: the base must be admissible
        } else if r.1 > sg - half - ext {
            vec![-1]
        } else {
            vec![0, 1, -1]
        }
    };
    // the best over the path signs of the Path axes of one piece
    let piece = |px: (AxisMode, (I, I)), py: (AxisMode, (I, I)), ovr: &[(usize, i64, u8)]| -> Option<I> {
        let sx: Vec<AxisMode> = match px.0 {
            AxisMode::Path(_) => opts(px.1).into_iter().map(AxisMode::Path).collect(),
            m => vec![m],
        };
        let sy: Vec<AxisMode> = match py.0 {
            AxisMode::Path(_) => opts(py.1).into_iter().map(AxisMode::Path).collect(),
            m => vec![m],
        };
        let mut best: Option<I> = None;
        for &mx in &sx {
            for &my in &sy {
                let spec = [AxisSpec { mode: mx, cr: px.1 }, AxisSpec { mode: my, cr: py.1 }];
                if let Some(v) = first_order_sig(ctx, b, bg, inw, spec, ovr) {
                    best = Some(best.map_or(v, |x| x.max(v)));
                    if v >= ctx.target {
                        return best;
                    }
                }
            }
        }
        best
    };
    let plain = piece((AxisMode::Path(0), cr[0]), (AxisMode::Path(0), cr[1]), &[]);
    if plain.map_or(false, |v| v >= ctx.target) {
        return plain;
    }
    // Lemma V: split by germ cones (at most one per axis, on free axes); every subset of the splittable axes is a
    // partition of the poses, so the best over the subsets is valid
    let splits = germ_splits(ctx, b, bg, cr);
    let mut sp: [Option<&GermSplit>; 2] = [None, None];
    for axi in 0..2 {
        // vertical pairs cut the x axis (family 0), horizontal pairs the y axis (family 1)
        if opts(cr[axi]).len() == 3 {
            sp[axi] = splits.iter().find(|s| s.fam == axi);
        }
    }
    let mut worst_best: Option<I> = None;
    for subset in [1u32, 2, 3] {
        if (0..2).any(|axi| (subset >> axi) & 1 == 1 && sp[axi].is_none()) {
            continue;
        }
        let mut pieces: [Vec<((AxisMode, (I, I)), Option<(usize, i64, u8)>)>; 2] = [Vec::new(), Vec::new()];
        for axi in 0..2 {
            match sp[axi].filter(|_| (subset >> axi) & 1 == 1) {
                Some(s) => {
                    let (c0, c1) = cr[axi];
                    let fam = s.fam;
                    pieces[axi].push(((AxisMode::Cone { g: s.g, k0: s.k0, k1: s.k1 }, cr[axi]), None));
                    // vertical: c <= g is region 1 (a window), c >= g region 2; horizontal: c_y >= g region 1
                    let (r_lo, r_hi) = if fam == 0 { (1u8, 2u8) } else { (2u8, 1u8) };
                    pieces[axi].push(((AxisMode::Path(0), (c0, s.g)), Some((fam, s.pos_a, r_lo))));
                    pieces[axi].push(((AxisMode::Path(0), (s.g, c1)), Some((fam, s.pos_a, r_hi))));
                }
                None => pieces[axi].push(((AxisMode::Path(0), cr[axi]), None)),
            }
        }
        let mut worst: Option<I> = None;
        'outer: for &(px, ox) in &pieces[0] {
            for &(py, oy) in &pieces[1] {
                let ovr: Vec<(usize, i64, u8)> = [ox, oy].into_iter().flatten().collect();
                match piece(px, py, &ovr) {
                    Some(v) => worst = Some(worst.map_or(v, |x| x.min(v))),
                    None => {
                        worst = None;
                        break 'outer;
                    }
                }
                if worst.map_or(false, |w| w < ctx.target) {
                    break 'outer; // this subset cannot certify
                }
            }
        }
        if let Some(w) = worst {
            worst_best = Some(worst_best.map_or(w, |x| x.max(w)));
            if w >= ctx.target {
                break;
            }
        }
    }
    let worst = worst_best;
    match (plain, worst) {
        (Some(p), Some(w)) => Some(p.max(w)),
        (p, w) => p.or(w),
    }
}

// =====================================================================================
// theta = 0 exactly: the axis-parallel squares (ZMX2_AREA.md sec 6, `zmx2 cert0`)
// =====================================================================================
//
// At theta = 0, Q(c) = [c_x - 1/2, c_x + 1/2] x [c_y - 1/2, c_y + 1/2] (closed).  A centre box
// B = [x0,x1] x [y0,y1] (exact rationals, denominator cd = 10 2^cl) is bounded by Lemma Z0:
//   points certainly in Q for all c in B; every line with |c_perp - l| <= 1/2 on all of B, by the exact
//   infimum over c_along of its window F(c + 1/2) - F(c - 1/2); every rectangle by the exact infimum of
//   its overlap area, which is (inf over c_x of the x-overlap) * (inf over c_y of the y-overlap).

/// floor(n/d) and ceil(n/d), d > 0
fn fdiv(n: I, d: I) -> I {
    n.div_euclid(d)
}
fn cdiv(n: I, d: I) -> I {
    -((-n).div_euclid(d))
}

/// Lemma Z0: lower bound of mu(Q(c, 0)) for every c in the box (bound units); exact integers.
/// xn, yn: numerators over cd.  At theta = 0 the vertical lines in Q contribute V(c_y) = sum of windows
/// F(c_y + 1/2) - F(c_y - 1/2), the horizontal ones H(c_x), each rectangle rho ox(c_x) oy(c_y); all are
/// piecewise linear (the area piecewise bilinear) on the grid of candidate abscissae/ordinates, so the
/// infimum over the box is the minimum over that grid.
fn bound0(ctx: &Ctx, xn: [I; 2], yn: [I; 2], cd: I) -> I {
    let cv = ctx.cv;
    let d = cv.d;
    let g = ctx.g;
    let half = g / 2; // G = D 2^K is even
    let mut pts_w: I = 0;
    // a coordinate p/D is in [c - 1/2, c + 1/2] for every c in [n0/cd, n1/cd]  <=>  2 p cd >= 2 n1 D - D cd
    // and 2 p cd <= 2 n0 D + D cd
    let inall = |p: I, n: [I; 2]| -> bool { 2 * p * cd >= 2 * n[1] * d - d * cd && 2 * p * cd <= 2 * n[0] * d + d * cd };
    // points: ordinary ones (line atoms below)
    for &(px, py, w) in &cv.pts {
        if inall(px as I, xn) && inall(py as I, yn) {
            pts_w += w * ctx.ptw;
        }
    }
    // grid range of the box along each axis (outward if the box is not on the grid)
    let gx = [fdiv(xn[0] * g, cd), cdiv(xn[1] * g, cd)];
    let gy = [fdiv(yn[0] * g, cd), cdiv(yn[1] * g, cd)];
    let mut xs: Vec<I> = vec![gx[0], gx[1]];
    let mut ys: Vec<I> = vec![gy[0], gy[1]];
    let push = |v: &mut Vec<I>, c: I, r: [I; 2]| {
        if c > r[0] && c < r[1] {
            v.push(c);
        }
    };
    let mut vin: Vec<&Line> = Vec::new();
    let mut hin: Vec<&Line> = Vec::new();
    for l in &cv.vl {
        if inall(l.pos as I, xn) {
            for &(t, w) in &l.atoms {
                if inall(t as I, yn) {
                    pts_w += w << K;
                }
            }
            for &bp in &l.bps {
                let y = (bp as I) << K;
                push(&mut ys, y - half, gy);
                push(&mut ys, y + half, gy);
            }
            vin.push(l);
        }
    }
    for l in &cv.hl {
        // horizontal line y = -pos, along x
        if inall(-(l.pos as I), yn) {
            for &(t, w) in &l.atoms {
                if inall(t as I, xn) {
                    pts_w += w << K;
                }
            }
            for &bp in &l.bps {
                let x = (bp as I) << K;
                push(&mut xs, x - half, gx);
                push(&mut xs, x + half, gx);
            }
            hin.push(l);
        }
    }
    // rectangles: overlap length of [c - 1/2, c + 1/2] with [r0, r1] (grid units), kinks at r0 -+ 1/2, r1 -+ 1/2
    let ov = |c: I, r0: I, r1: I| -> I { ((c + half).min(r1) - (c - half).max(r0)).max(0) };
    for r in &cv.rects {
        for e in [r.x0, r.x1] {
            let e = (e as I) << K;
            push(&mut xs, e - half, gx);
            push(&mut xs, e + half, gx);
        }
        for e in [r.y0, r.y1] {
            let e = (e as I) << K;
            push(&mut ys, e - half, gy);
            push(&mut ys, e + half, gy);
        }
    }
    xs.sort();
    xs.dedup();
    ys.sort();
    ys.dedup();
    let vv: Vec<I> = ys.iter().map(|&y| vin.iter().map(|l| fval(l, y + half) - fval(l, y - half)).sum()).collect();
    let hv: Vec<I> = xs.iter().map(|&x| hin.iter().map(|l| fval(l, x + half) - fval(l, x - half)).sum()).collect();
    let mut best: Option<I> = None;
    for (i, &x) in xs.iter().enumerate() {
        for (j, &y) in ys.iter().enumerate() {
            let mut v = pts_w + hv[i] + vv[j];
            for (ri, r) in cv.rects.iter().enumerate() {
                let ox = ov(x, (r.x0 as I) << K, (r.x1 as I) << K);
                let oy = ov(y, (r.y0 as I) << K, (r.y1 as I) << K);
                if ox == 0 || oy == 0 {
                    continue;
                }
                // mass = w Lc ox oy / (Pn 2^K)  (ox, oy in grid units G = D 2^K), rounded down
                let pn = ((r.x1 - r.x0) as I) * ((r.y1 - r.y0) as I);
                let num = ctx.rect_wlc[ri].checked_mul(ox).and_then(|t| t.checked_mul(oy));
                v += match num {
                    Some(n) => n / (pn << K),
                    None => {
                        // (not reached for the k2m3 cover) a float lower bound, rounded down
                        let a = Iv::rat(ox, g).mul(Iv::rat(oy, g)).lo;
                        dn(a * dn(ctx.rect_tr[ri] as f64)).floor().max(0.0) as I
                    }
                };
            }
            best = Some(best.map_or(v, |b| b.min(v)));
        }
    }
    best.unwrap()
}

fn cmd_cert0(args: &[String], cv: Cover) {
    let d4 = has(args, "--d4");
    let full = has(args, "--full");
    if d4 == full {
        die("exactly one of --d4 / --full is required");
    }
    let depth: u32 = arg_val(args, "--depth").map(|s| s.parse().unwrap()).unwrap_or(40);
    let uncert_cap: usize = arg_val(args, "--uncert-cap").map(|s| s.parse().unwrap()).unwrap_or(200);
    let tenths = cv.sx * 10;
    if tenths % cv.d != 0 {
        die("container side must be a multiple of 1/10");
    }
    let ncell = (tenths / cv.d) as i64;
    if d4 && ncell % 2 != 0 {
        die("--d4 needs the container side to be a multiple of 1/5");
    }
    if d4 {
        match check_d4(&cv) {
            Ok(()) => println!("D4: measure invariant under x->s-x and x<->y (exact)"),
            Err(e) => die(&format!("--d4: cover is not D4-invariant: {}", e)),
        }
    }
    if ncell < 10 {
        die("container side must be >= 1");
    }
    // admissible centres at theta = 0: [1/2, s - 1/2]; cells of 1/10 from 5 to ncell - 6 (D4: to ncell/2 - 1)
    let cmax = if d4 { ncell / 2 } else { ncell - 5 };
    let ctx = Ctx::new(&cv);
    println!(
        "# zmx2 cert0 (theta = 0) hash={:016x} mode={} depth={}{}",
        cv.hash,
        if d4 { "d4" } else { "full" },
        depth,
        if cv.rects.is_empty() { String::new() } else { format!(" area={}", cv.rects.len()) }
    );
    let t0 = Instant::now();
    let (mut nbox, mut ncert, mut nunc, mut maxd) = (0usize, 0usize, 0usize, 0u32);
    let mut nroots = 0usize;
    let mut worst: Option<(I, [I; 2], [I; 2], I)> = None;
    for i in 5..cmax {
        for j in 5..cmax {
            nroots += 1;
            // stack of (xn, yn, cd, depth)
            let mut st: Vec<([I; 2], [I; 2], I, u32)> = vec![([i as I, i as I + 1], [j as I, j as I + 1], 10, 0)];
            while let Some((xn, yn, cd, dep)) = st.pop() {
                nbox += 1;
                maxd = maxd.max(dep);
                let bd = bound0(&ctx, xn, yn, cd);
                if bd >= ctx.target {
                    ncert += 1;
                    continue;
                }
                if dep >= depth {
                    nunc += 1;
                    if worst.as_ref().map_or(true, |w| bd < w.0) {
                        worst = Some((bd, xn, yn, cd));
                    }
                    if nunc <= 20 {
                        println!(
                            "UNCERT0 box {}/{},{}/{},{}/{},{}/{} bound {:.9}",
                            xn[0], cd, xn[1], cd, yn[0], cd, yn[1], cd,
                            bd as f64 / ctx.target as f64
                        );
                    }
                    if nunc >= uncert_cap {
                        break;
                    }
                    continue;
                }
                // halve the longer side (x on ties); refine the denominator when needed
                let (mut xn, mut yn, mut cd) = (xn, yn, cd);
                let splitx = xn[1] - xn[0] >= yn[1] - yn[0];
                let odd = if splitx { (xn[0] + xn[1]) % 2 != 0 } else { (yn[0] + yn[1]) % 2 != 0 };
                if odd {
                    xn = [2 * xn[0], 2 * xn[1]];
                    yn = [2 * yn[0], 2 * yn[1]];
                    cd *= 2;
                }
                if splitx {
                    let m = (xn[0] + xn[1]) / 2;
                    st.push(([m, xn[1]], yn, cd, dep + 1));
                    st.push(([xn[0], m], yn, cd, dep + 1));
                } else {
                    let m = (yn[0] + yn[1]) / 2;
                    st.push((xn, [m, yn[1]], cd, dep + 1));
                    st.push((xn, [yn[0], m], cd, dep + 1));
                }
            }
            if nunc >= uncert_cap {
                break;
            }
        }
        if nunc >= uncert_cap {
            break;
        }
    }
    println!(
        "done0 in {:.1}s: roots {}, boxes {}, certified {}, uncertified {}, max depth {}",
        t0.elapsed().as_secs_f64(),
        nroots,
        nbox,
        ncert,
        nunc,
        maxd
    );
    if let Some((bd, xn, yn, cd)) = worst {
        println!(
            "worst uncertified box {}/{},{}/{},{}/{},{}/{}: bound {:.12}",
            xn[0], cd, xn[1], cd, yn[0], cd, yn[1], cd,
            bd as f64 / ctx.target as f64
        );
    }
    if nunc > 0 {
        println!("NOT VERIFIED (theta = 0): {} uncertified boxes", nunc);
    } else if d4 {
        println!("VERIFIED-D4 (theta = 0): every closed axis-parallel unit square in [0,s]^2 has mu >= 1 (D4-reduced)");
    } else {
        println!("VERIFIED (theta = 0): every closed axis-parallel unit square in [0,s]^2 has mu >= 1");
    }
}

// =====================================================================================
// CLI
// =====================================================================================

fn parse_rat(s: &str) -> (I, I) {
    let s = s.trim();
    if let Some(p) = s.find('/') {
        let n: I = s[..p].trim().parse().unwrap_or_else(|_| die(&format!("bad rational {}", s)));
        let d: I = s[p + 1..].trim().parse().unwrap_or_else(|_| die(&format!("bad rational {}", s)));
        (n, d)
    } else {
        (s.parse().unwrap_or_else(|_| die(&format!("bad rational {}", s))), 1)
    }
}

/// box from rationals: denominators must be of the form 10*2^a (centre) and 8*2^b (u)
fn box_from_rats(v: &[(I, I)]) -> PBox {
    let to_level = |r: (I, I), base: I| -> (I, u32) {
        for l in 0..40u32 {
            let den = base << l;
            if (r.0 * den) % r.1 == 0 {
                return (r.0 * den / r.1, l);
            }
        }
        die("box coordinate denominators must divide 10*2^k (centre) / 8*2^k (u)");
    };
    let cs: Vec<(I, u32)> = v[..4].iter().map(|&r| to_level(r, 10)).collect();
    let cl = cs.iter().map(|c| c.1).max().unwrap();
    let cn: Vec<I> = cs.iter().map(|&(n, l)| n << (cl - l)).collect();
    let us: Vec<(I, u32)> = v[4..6].iter().map(|&r| to_level(r, 8)).collect();
    let ul = us.iter().map(|c| c.1).max().unwrap();
    let un: Vec<I> = us.iter().map(|&(n, l)| n << (ul - l)).collect();
    PBox { xn: [cn[0], cn[1]], yn: [cn[2], cn[3]], cl, un: [un[0], un[1]], ul }
}

fn arg_val<'a>(args: &'a [String], key: &str) -> Option<&'a str> {
    args.iter().position(|a| a == key).map(|i| args.get(i + 1).map(|s| s.as_str()).unwrap_or_else(|| die(&format!("{} needs a value", key))))
}
fn has(args: &[String], key: &str) -> bool {
    args.iter().any(|a| a == key)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage:\n  zmx2 cert FILE (--d4 | --full) [--threads T] [--depth N] [--node-cap N] [--kappa F]\n\
             \x20                [--log PATH] [--xlo i --xhi i --ylo j --yhi j --bins a-b] [--uncert-cap N]\n\
             \x20                [--pair-points] [--sym-atoms] [--no-atoms] [--mirror-only] [--first-order]\n\
             \x20                [--umin m]   (bin 0 starts at u = 2^-m/8; restricted region, no verdict)\n\
             \x20 zmx2 cert0 FILE (--d4 | --full) [--depth N]         (theta = 0 exactly: axis-parallel squares)\n\
             \x20 zmx2 box FILE --box x0,x1,y0,y1,u0,u1 [--refl]      (rationals; bound breakdown)\n\
             \x20 zmx2 fmass FILE --x X --y Y --u U                   (float mu, diagnostics)\n\
             \x20 zmx2 fscan FILE [--pitch P] [--ubins N]             (float landscape of the D4 region)\n\
             \x20 zmx2 d4 FILE                                        (exact D4 invariance check)\n\
             \x20 zmx2 info FILE"
        );
        std::process::exit(2);
    }
    let cmd = args[1].as_str();
    if has(&args, "--no-atoms") {
        NO_ATOMS.store(true, Ordering::Relaxed);
    }
    if has(&args, "--pair-points") {
        PAIR_POINTS.store(true, Ordering::Relaxed);
    }
    if has(&args, "--sym-atoms") {
        SYM_ATOMS.store(true, Ordering::Relaxed);
    }
    if has(&args, "--mirror-only") {
        MIRROR_ONLY.store(true, Ordering::Relaxed);
    }
    if has(&args, "--first-order") {
        FIRST_ORDER.store(true, Ordering::Relaxed);
    }
    let cv = parse_cover(&args[2]);
    match cmd {
        "info" => {
            println!(
                "s = {}/{}  D = {}  W = {}  points {} (aggregated {})  segments {}  vlines {}  hlines {}  Lc {}",
                cv.s_num,
                cv.s_den,
                cv.d,
                cv.w,
                cv.raw_pts.len(),
                cv.pts.len(),
                cv.raw_segs.len(),
                cv.vl.len(),
                cv.hl.len(),
                cv.lc
            );
            println!("total = {}/{} = {:.12}", cv.total.0, cv.total.1, cv.total.0 as f64 / cv.total.1 as f64);
            for r in &cv.rects {
                println!(
                    "area density on [{}/{}, {}/{}] x [{}/{}, {}/{}], mass {}/{}",
                    r.x0, cv.d, r.x1, cv.d, r.y0, cv.d, r.y1, cv.d, r.w, cv.w
                );
            }
            if SYM_ATOMS.load(Ordering::Relaxed) {
                let na = |v: &[Line]| v.iter().map(|l| l.atoms.len()).sum::<usize>();
                match &cv.alt {
                    None => println!("sym-atoms: the mirrored assignment is the same (no point changes line)"),
                    Some((vl2, hl2)) => println!(
                        "sym-atoms: atoms vertical/horizontal {}/{} (default), {}/{} (mirrored)",
                        na(&cv.vl),
                        na(&cv.hl),
                        na(vl2),
                        na(hl2)
                    ),
                }
            }
        }
        "d4" => match check_d4(&cv) {
            Ok(()) => println!("D4: measure invariant under x->s-x and x<->y (exact)"),
            Err(e) => die(&format!("--d4: cover is not D4-invariant: {}", e)),
        },
        "fmass" => {
            let x: f64 = arg_val(&args, "--x").unwrap().parse().unwrap();
            let y: f64 = arg_val(&args, "--y").unwrap().parse().unwrap();
            let u: f64 = arg_val(&args, "--u").unwrap().parse().unwrap();
            println!("{:.12}", float_mass(&cv, x, y, u));
        }
        "fscan" => {
            let pitch: f64 = arg_val(&args, "--pitch").map(|s| s.parse().unwrap()).unwrap_or(0.01);
            let nu: usize = arg_val(&args, "--ubins").map(|s| s.parse().unwrap()).unwrap_or(50);
            let skipwall = has(&args, "--skip-wall");
            let umax: f64 = arg_val(&args, "--umax").map(|s| s.parse().unwrap()).unwrap_or(0.5); // skip poses touching a wall
            let outside = has(&args, "--outside"); // (area covers) skip poses with Q inside a rectangle
            let s = cv.sx as f64 / cv.d as f64;
            let n = (s / 2.0 / pitch).round() as usize;
            let mut hist = [0usize; 12];
            let mut worst: Vec<(f64, f64, f64, f64)> = Vec::new();
            for k in 0..=nu {
                let u = umax * k as f64 / nu as f64;
                for i in 0..=n {
                    for j in 0..=n {
                        let (x, y) = (i as f64 * pitch, j as f64 * pitch);
                        if !admissible_f(&cv, x, y, u) {
                            continue;
                        }
                        if skipwall {
                            let hw = 0.5 * (1.0 - u * u + 2.0 * u) / (1.0 + u * u);
                            let sf = cv.sx as f64 / cv.d as f64;
                            if x - hw < 1e-9 || y - hw < 1e-9 || x + hw > sf - 1e-9 || y + hw > sf - 1e-9 {
                                continue;
                            }
                        }
                        if outside {
                            // skip poses with Q inside a rectangle (area 1 there): w/2 margin to every side
                            let hw = 0.5 * (1.0 - u * u + 2.0 * u) / (1.0 + u * u);
                            let dd = cv.d as f64;
                            if cv.rects.iter().any(|r| {
                                x - hw >= r.x0 as f64 / dd
                                    && x + hw <= r.x1 as f64 / dd
                                    && y - hw >= r.y0 as f64 / dd
                                    && y + hw <= r.y1 as f64 / dd
                            }) {
                                continue;
                            }
                        }
                        let m = float_mass(&cv, x, y, u);
                        let bin = (((m - 1.0) / 0.01).floor().max(0.0) as usize).min(11);
                        hist[bin] += 1;
                        worst.push((m, x, y, u));
                    }
                }
            }
            worst.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            println!("histogram of mu-1 in 1% bins (last = >=11%): {:?}", hist);
            for w in worst.iter().take(15) {
                println!("  mu {:.6} at ({:.4},{:.4}) u {:.5} th {:.3} deg", w.0, w.1, w.2, w.3, 2.0 * w.3.atan().to_degrees());
            }
        }
        "box" => {
            let bs = arg_val(&args, "--box").unwrap_or_else(|| die("--box needed"));
            let v: Vec<(I, I)> = bs.split(',').map(parse_rat).collect();
            if v.len() != 6 {
                die("--box needs 6 rationals");
            }
            let cvr;
            let cvu = if has(&args, "--refl") {
                cvr = reflect_y(&cv);
                &cvr
            } else {
                &cv
            };
            let b = box_from_rats(&v);
            let ctx = Ctx::new(cvu);
            let bg = box_geo(&ctx, &b);
            let empty = box_empty(&ctx, &b, &bg);
            let cand = initial_candidates(cvu, &b);
            let mut inw: I = 0;
            let mut nin = 0;
            for &i in &cand {
                let (px, py, w) = cvu.pts[i as usize];
                if pt_certain_in(cvu.d, px as I, py as I, &b) {
                    inw += w * ctx.ptw;
                    nin += 1;
                }
            }
            let sb = mass_bound(&ctx, &b, &bg, inw);
            let tot = inw + sb;
            println!("box {}  {}", b.exact_str(), b.fstr());
            println!("empty {}", empty);
            println!("unit {}", ctx.target);
            println!("points_in {} weight {}", nin, inw);
            println!("segments {}", sb);
            println!("bound {} = {:.9}", tot, tot as f64 / ctx.target as f64);
            println!("certified {}", !empty && tot >= ctx.target);
        }
        "boxes" => {
            // batch: one box per stdin line "x0,x1,y0,y1,u0,u1"; prints "bound unit empty"
            let cvr;
            let cvu = if has(&args, "--refl") {
                cvr = reflect_y(&cv);
                &cvr
            } else {
                &cv
            };
            let ctx = Ctx::new(cvu);
            let stdin = std::io::stdin();
            for line in stdin.lock().lines() {
                let line = line.unwrap();
                if line.trim().is_empty() {
                    continue;
                }
                let v: Vec<(I, I)> = line.split(',').map(parse_rat).collect();
                let b = box_from_rats(&v);
                let bg = box_geo(&ctx, &b);
                let empty = box_empty(&ctx, &b, &bg);
                let cand = initial_candidates(cvu, &b);
                let mut inw: I = 0;
                for &i in &cand {
                    let (px, py, w) = cvu.pts[i as usize];
                    if pt_certain_in(cvu.d, px as I, py as I, &b) {
                        inw += w * ctx.ptw;
                    }
                }
                let sb = mass_bound(&ctx, &b, &bg, inw);
                println!("{} {} {} {}", inw + sb, ctx.target, empty as u8, inw);
            }
        }
        "cert" => cmd_cert(&args, cv),
        "cert0" => cmd_cert0(&args, cv),
        "boxes0" => {
            // theta = 0 batch: one centre box per stdin line "x0,x1,y0,y1" (rationals with denominators
            // dividing 10 2^k); prints "bound unit" (Lemma Z0)
            let ctx = Ctx::new(&cv);
            let stdin = std::io::stdin();
            for line in stdin.lock().lines() {
                let line = line.unwrap();
                if line.trim().is_empty() {
                    continue;
                }
                let v: Vec<(I, I)> = line.split(',').map(parse_rat).collect();
                if v.len() != 4 {
                    die("boxes0 needs 4 rationals per line");
                }
                let mut v6 = v.clone();
                v6.push((0, 1));
                v6.push((1, 8));
                let b = box_from_rats(&v6);
                println!("{} {}", bound0(&ctx, b.xn, b.yn, b.cd()), ctx.target);
            }
        }
        _ => die("unknown command"),
    }
}

fn cmd_cert(args: &[String], cv: Cover) {
    let d4 = has(args, "--d4");
    let full = has(args, "--full");
    if d4 == full {
        die("exactly one of --d4 / --full is required");
    }
    let threads: usize = arg_val(args, "--threads").map(|s| s.parse().unwrap()).unwrap_or(1);
    let st = Settings {
        depth: arg_val(args, "--depth").map(|s| s.parse().unwrap()).unwrap_or(40),
        node_cap: arg_val(args, "--node-cap").map(|s| s.parse().unwrap()).unwrap_or(20_000_000),
        kappa: arg_val(args, "--kappa").map(|s| s.parse().unwrap()).unwrap_or(1.5),
        uncert_cap: arg_val(args, "--uncert-cap").map(|s| s.parse().unwrap()).unwrap_or(200),
        tight: arg_val(args, "--tight").map(|s| s.parse().unwrap()).unwrap_or(0.0),
    };
    let mut tightf = arg_val(args, "--dump-tight").map(|p| std::fs::File::create(p).unwrap());
    // container must be a multiple of 1/10 (1/5 under --d4)
    let tenths = cv.sx * 10;
    if tenths % cv.d != 0 {
        die("container side must be a multiple of 1/10");
    }
    let ncell = (tenths / cv.d) as i64;
    if d4 && ncell % 2 != 0 {
        die("--d4 needs the container side to be a multiple of 1/5");
    }
    if d4 {
        match check_d4(&cv) {
            Ok(()) => println!("D4: measure invariant under x->s-x and x<->y (exact)"),
            Err(e) => die(&format!("--d4: cover is not D4-invariant: {}", e)),
        }
    }
    let cmax = if d4 { ncell / 2 } else { ncell };
    let lim = |k: &str, dflt: i64| -> i64 { arg_val(args, k).map(|s| s.parse().unwrap()).unwrap_or(dflt) };
    let (xlo, xhi, ylo, yhi) = (lim("--xlo", 0), lim("--xhi", cmax - 1), lim("--ylo", 0), lim("--yhi", cmax - 1));
    let (blo, bhi) = match arg_val(args, "--bins") {
        Some(s) => {
            let p: Vec<i64> = s.split('-').map(|t| t.parse().unwrap()).collect();
            (p[0], p[1])
        }
        None => (0, 3),
    };
    let umin: Option<u32> = arg_val(args, "--umin").map(|s| s.parse().unwrap());
    if let Some(m) = umin {
        if m == 0 || m > 20 {
            die("--umin m needs 1 <= m <= 20");
        }
    }
    let restricted = xlo != 0 || xhi != cmax - 1 || ylo != 0 || yhi != cmax - 1 || blo != 0 || bhi != 3 || umin.is_some();
    let passes = if d4 { 1 } else { 2 };
    let covers: Vec<Cover> = if d4 { vec![cv.clone()] } else { vec![cv.clone(), reflect_y(&cv)] };
    let mut roots: Vec<(usize, usize, PBox)> = Vec::new(); // (id, pass, box)
    for pass in 0..passes {
        for i in 0..cmax {
            for j in 0..cmax {
                for k in 0..4 {
                    let id = ((pass as i64 * cmax + i) * cmax + j) as usize * 4 + k as usize;
                    if i < xlo || i > xhi || j < ylo || j > yhi || k < blo || k > bhi {
                        continue;
                    }
                    let mut b = PBox {
                        xn: [i as I, i as I + 1],
                        yn: [j as I, j as I + 1],
                        cl: 0,
                        un: [k as I, k as I + 1],
                        ul: 0,
                    };
                    if k == 0 {
                        if let Some(m) = umin {
                            // --umin m: bin 0 becomes u in [2^-m/8, 1/8] (theta < 2 atan(2^-m/8) left out)
                            b.un = [1, 1 << m];
                            b.ul = m;
                        }
                    }
                    roots.push((id, pass, b));
                }
            }
        }
    }
    let header = format!(
        "# zmx2 cert hash={:016x} mode={} atoms={} depth={} node_cap={} kappa={} K={} region=x{}-{},y{}-{},bins{}-{}{}",
        cv.hash,
        if d4 { "d4" } else { "full" },
        // (--sym-atoms appends "+sym"; without it the header is the same as before sec 4.9)
        format!(
            "{}{}",
            if NO_ATOMS.load(Ordering::Relaxed) {
                "none"
            } else if PAIR_POINTS.load(Ordering::Relaxed) {
                "pairpts"
            } else {
                "grid"
            },
            format!(
                "{}{}",
                if MIRROR_ONLY.load(Ordering::Relaxed) {
                    "+mirror"
                } else if SYM_ATOMS.load(Ordering::Relaxed) {
                    "+sym"
                } else {
                    ""
                },
                if FIRST_ORDER.load(Ordering::Relaxed) { "+fo" } else { "" }
            )
        ),
        st.depth,
        st.node_cap,
        st.kappa,
        K,
        xlo,
        xhi,
        ylo,
        yhi,
        blo,
        bhi,
        // (covers with area densities / --umin only: otherwise the header is the same as before)
        format!(
            "{}{}",
            if cv.rects.is_empty() { String::new() } else { format!(" area={}", cv.rects.len()) },
            match umin {
                None => String::new(),
                Some(m) => format!(" umin=2^-{}/8", m),
            }
        )
    );
    println!("{}", header);
    println!(
        "cover: s = {}/{}, {} points, {} segments, total {}/{} = {:.11}",
        cv.s_num,
        cv.s_den,
        cv.raw_pts.len(),
        cv.raw_segs.len(),
        cv.total.0,
        cv.total.1,
        cv.total.0 as f64 / cv.total.1 as f64
    );
    if !cv.rects.is_empty() {
        println!("area densities: {} rectangle(s) (ZMX2_AREA.md)", cv.rects.len());
    }
    // resume
    let mut done: HashMap<usize, String> = HashMap::new();
    let log_path = arg_val(args, "--log").map(|s| s.to_string());
    if let Some(p) = &log_path {
        if let Ok(f) = std::fs::File::open(p) {
            let rd = std::io::BufReader::new(f);
            let mut first = true;
            for line in rd.lines() {
                let line = line.unwrap();
                if first {
                    first = false;
                    if line != header {
                        die(&format!("log {} was written with different settings:\n{}\nvs\n{}", p, line, header));
                    }
                    continue;
                }
                if line.starts_with("ROOT ") {
                    let id: usize = line.split_whitespace().nth(1).unwrap().parse().unwrap();
                    done.insert(id, line.clone());
                }
            }
        }
    }
    let mut logf = log_path.as_ref().map(|p| {
        let exists = std::path::Path::new(p).exists() && std::fs::metadata(p).map(|m| m.len() > 0).unwrap_or(false);
        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(p).unwrap();
        if !exists {
            writeln!(f, "{}", header).unwrap();
        }
        f
    });
    let todo: Vec<(usize, usize, PBox)> = roots.iter().filter(|r| !done.contains_key(&r.0)).cloned().collect();
    println!("roots: {} in region, {} already done (log), {} to run, {} threads", roots.len(), done.len(), todo.len(), threads);
    let covers = Arc::new(covers);
    let todo = Arc::new(todo);
    let next = Arc::new(AtomicUsize::new(0));
    let stop = Arc::new(AtomicBool::new(false));
    let (tx, rx) = mpsc::channel::<RootResult>();
    let t0 = Instant::now();
    let mut handles = Vec::new();
    for _ in 0..threads {
        let covers = covers.clone();
        let todo = todo.clone();
        let next = next.clone();
        let tx = tx.clone();
        let st = st.clone();
        let stop = stop.clone();
        handles.push(std::thread::spawn(move || {
            let ctxs: Vec<Ctx> = covers.iter().map(Ctx::new).collect();
            loop {
                if stop.load(Ordering::Relaxed) {
                    break;
                }
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= todo.len() {
                    break;
                }
                let (id, pass, b) = todo[i];
                let r = run_root(&ctxs[pass], id, b, &st);
                tx.send(r).unwrap();
            }
        }));
    }
    drop(tx);
    let mut n_new = 0usize;
    let total_todo = todo.len();
    let mut new_lines: Vec<String> = Vec::new();
    for r in rx {
        n_new += 1;
        let (_, pass, root) = *todo.iter().find(|t| t.0 == r.id).unwrap();
        let line = format!(
            "ROOT {} pass {} root {} boxes {} cert {} empty {} uncert {} maxdepth {} capped {} ms {}",
            r.id,
            pass,
            root.exact_str(),
            r.boxes,
            r.cert,
            r.empty,
            r.uncert.len(),
            r.maxdepth,
            r.capped as u8,
            r.ms
        );
        let mut extra = String::new();
        for ub in r.uncert.iter().take(40) {
            let fm = float_min_box(&covers[pass], ub, 300);
            extra.push_str(&format!(
                "UNCERT {} pass {} box {} | {} | float-min {:.6} at ({:.6},{:.6},u {:.7})\n",
                r.id,
                pass,
                ub.exact_str(),
                ub.fstr(),
                fm.0,
                fm.1,
                fm.2,
                fm.3
            ));
        }
        if !extra.is_empty() {
            print!("{}", extra);
        }
        if let Some(f) = tightf.as_mut() {
            for (tb, bd) in &r.tight {
                writeln!(f, "{} {} {}", pass, tb.exact_str(), bd).unwrap();
            }
        }
        if let Some(f) = logf.as_mut() {
            f.write_all(extra.as_bytes()).unwrap();
            writeln!(f, "{}", line).unwrap();
            f.flush().unwrap();
        }
        new_lines.push(line.clone());
        if n_new % 25 == 0 || n_new == total_todo {
            println!(
                "progress {}/{} roots, {:.0}s",
                n_new,
                total_todo,
                t0.elapsed().as_secs_f64()
            );
        }
    }
    for h in handles {
        h.join().unwrap();
    }
    // census over all roots (old + new)
    let mut all: HashMap<usize, String> = done;
    for l in new_lines {
        let id: usize = l.split_whitespace().nth(1).unwrap().parse().unwrap();
        all.insert(id, l);
    }
    let mut boxes = 0usize;
    let mut unc = 0usize;
    let mut cert = 0usize;
    let mut empty = 0usize;
    let mut maxd = 0u32;
    let mut ms: u128 = 0;
    let mut missing = 0;
    for r in &roots {
        match all.get(&r.0) {
            None => missing += 1,
            Some(l) => {
                let f: Vec<&str> = l.split_whitespace().collect();
                let g = |key: &str| -> u128 {
                    let p = f.iter().position(|t| *t == key).unwrap();
                    f[p + 1].parse().unwrap()
                };
                boxes += g("boxes") as usize;
                cert += g("cert") as usize;
                empty += g("empty") as usize;
                unc += g("uncert") as usize;
                maxd = maxd.max(g("maxdepth") as u32);
                ms += g("ms");
            }
        }
    }
    println!(
        "done in {:.0}s wall: roots {} (missing {}), boxes {}, certified {}, empty {}, uncertified {}, max depth {}, cpu {:.0}s",
        t0.elapsed().as_secs_f64(),
        roots.len(),
        missing,
        boxes,
        cert,
        empty,
        unc,
        maxd,
        ms as f64 / 1000.0
    );
    if missing > 0 {
        println!("INCOMPLETE: {} roots not run", missing);
    } else if unc > 0 {
        println!("NOT VERIFIED: {} uncertified boxes", unc);
    } else if restricted {
        println!("REGION CLEAN (restricted region; no verdict for the whole cover)");
    } else if d4 {
        println!("VERIFIED-D4: every closed unit square in [0,s]^2 has mu >= 1 (D4-reduced sweep, all roots)");
    } else {
        println!("VERIFIED: every closed unit square in [0,s]^2 has mu >= 1 (unreduced sweep, all roots)");
    }
}
