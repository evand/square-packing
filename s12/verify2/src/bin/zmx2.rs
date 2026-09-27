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
    buckets: Vec<Vec<u32>>, // point index grid, cell 1/10 (in [0,s])
    nb: usize,
    raw_pts: Vec<(i64, i64, I)>,
    raw_segs: Vec<(i64, i64, i64, i64, I)>,
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
    if (s_num * d) % s_den != 0 {
        die("s_den must divide s_num*D");
    }
    let sx = s_num * d / s_den;
    if d >= (1 << 22) {
        die("D too large for this checker (limit 2^22)");
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
        raw_pts.push((x as i64, y as i64, wt));
    }
    let mut raw_segs = Vec::new();
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
            if x0 == x1 && y0 == y1 {
                die("degenerate (zero-length) segment");
            }
            if x0 != x1 && y0 != y1 {
                die("unsupported: zmx2 handles axis-parallel segments only");
            }
            raw_segs.push((x0 as i64, y0 as i64, x1 as i64, y1 as i64, wt));
        }
        let npg = nx("npg");
        if npg != 0 {
            die("unsupported: zmx2 does not handle polygons (npg must be 0)");
        }
    }
    if it2.next().is_some() {
        die("trailing tokens after the declared pieces");
    }
    build_cover(s_num, s_den, d, w, sx, raw_pts, raw_segs, fnv(&bytes))
}

fn build_cover(
    s_num: I,
    s_den: I,
    d: I,
    w: I,
    sx: I,
    raw_pts: Vec<(i64, i64, I)>,
    raw_segs: Vec<(i64, i64, i64, i64, I)>,
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
    let mut pts: Vec<(i64, i64, I)> = Vec::new();
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
    // (orientation, pos) -> segments (a, b, w) and atoms (t, w) along the line; vertical: pos = x,
    // along = y; horizontal (rotated frame, Lemma T): pos = -y, along = x.
    let mut groups: BTreeMap<(u8, i64), (Vec<(i64, i64, I)>, Vec<(i64, I)>)> = BTreeMap::new();
    for (&(x, y), &wt) in &pm {
        if wt <= 0 {
            continue;
        }
        if !no_atoms && vpos.contains(&x) {
            groups.entry((0, x)).or_default().1.push((y, wt * lc));
        } else if !no_atoms && hpos.contains(&y) {
            groups.entry((1, -y)).or_default().1.push((x, wt * lc));
        } else if !no_atoms && vpart.contains(&x) {
            groups.entry((0, x)).or_default().1.push((y, wt * lc));
        } else if !no_atoms && hpart.contains(&y) {
            groups.entry((1, -y)).or_default().1.push((x, wt * lc));
        } else {
            pts.push((x, y, wt));
        }
    }
    for &(x0, y0, x1, y1, wt) in &raw_segs {
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
    // buckets of 1/10
    let nb = ((sx * 10 + d - 1) / d) as usize + 1;
    let mut buckets = vec![Vec::new(); nb * nb];
    for (i, &(x, y, _)) in pts.iter().enumerate() {
        let bx = ((x as I) * 10 / d) as usize;
        let by = ((y as I) * 10 / d) as usize;
        buckets[bx.min(nb - 1) * nb + by.min(nb - 1)].push(i as u32);
    }
    let tnum: I = raw_pts.iter().map(|p| p.2).sum::<I>() + raw_segs.iter().map(|s| s.4).sum::<I>();
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
        buckets,
        nb,
        raw_pts,
        raw_segs,
        total: (tnum / g, w / g),
        hash,
    }
}

/// The cover reflected by y -> s - y (used by the unreduced sweep for theta in [45,90] deg).
fn reflect_y(c: &Cover) -> Cover {
    let sx = c.sx as i64;
    let p: Vec<_> = c.raw_pts.iter().map(|&(x, y, w)| (x, sx - y, w)).collect();
    let s: Vec<_> = c.raw_segs.iter().map(|&(x0, y0, x1, y1, w)| (x0, sx - y0, x1, sx - y1, w)).collect();
    build_cover(c.s_num, c.s_den, c.d, c.w, c.sx, p, s, c.hash ^ 0x5bd1e995)
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

fn segment_bound(ctx: &Ctx, b: &PBox, bg: &BoxGeo) -> I {
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
    family_bound(ctx, &fv, &ctx.cv.vl) + family_bound(ctx, &fh, &ctx.cv.hl)
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
    m
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
    let sb = segment_bound(ctx, b, &bg);
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
            let s = cv.sx as f64 / cv.d as f64;
            let n = (s / 2.0 / pitch).round() as usize;
            let mut hist = [0usize; 12];
            let mut worst: Vec<(f64, f64, f64, f64)> = Vec::new();
            for k in 0..=nu {
                let u = 0.5 * k as f64 / nu as f64;
                for i in 0..=n {
                    for j in 0..=n {
                        let (x, y) = (i as f64 * pitch, j as f64 * pitch);
                        if !admissible_f(&cv, x, y, u) {
                            continue;
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
            let sb = segment_bound(&ctx, &b, &bg);
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
                let sb = segment_bound(&ctx, &b, &bg);
                println!("{} {} {} {}", inw + sb, ctx.target, empty as u8, inw);
            }
        }
        "cert" => cmd_cert(&args, cv),
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
    let restricted = xlo != 0 || xhi != cmax - 1 || ylo != 0 || yhi != cmax - 1 || blo != 0 || bhi != 3;
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
                    let b = PBox {
                        xn: [i as I, i as I + 1],
                        yn: [j as I, j as I + 1],
                        cl: 0,
                        un: [k as I, k as I + 1],
                        ul: 0,
                    };
                    roots.push((id, pass, b));
                }
            }
        }
    }
    let header = format!(
        "# zmx2 cert hash={:016x} mode={} atoms={} depth={} node_cap={} kappa={} K={} region=x{}-{},y{}-{},bins{}-{}",
        cv.hash,
        if d4 { "d4" } else { "full" },
        if NO_ATOMS.load(Ordering::Relaxed) {
            "none"
        } else if PAIR_POINTS.load(Ordering::Relaxed) {
            "pairpts"
        } else {
            "grid"
        },
        st.depth,
        st.node_cap,
        st.kappa,
        K,
        xlo,
        xhi,
        ylo,
        yhi,
        blo,
        bhi
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
