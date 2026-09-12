//! `zmcheck` --- an independent exact checker for weighted *closed* covers of `[0,m]^2`.
//!
//! Claim checked (`search/ZEROMARGIN.md` sec 1, `certificates/FORMAT.md`):
//!
//! > for every closed unit square `Q` with `Q subseteq [0,m]^2`, at every centre and every
//! > angle, the total weight of the certificate points lying in `Q` (a point on `dQ` counts)
//! > is at least `1`.
//!
//! This is a second implementation of the statement that `search/zeromargin.py` verifies.  It was
//! written from the *lemmas* of `search/RUNG2.md` / `search/ZEROMARGIN.md`, not from that code,
//! and its primitive set is deliberately different: one polynomial class, one exact maximum
//! routine, one inference rule.  See `search/RUNG2_XCHECK.md` for the proofs.
//!
//! ## The one polynomial class
//!
//! Pose `(x, y, u)`, `u = tan(theta/2) in [0,1]` (so `theta in [0,90]` deg, which is all angles
//! since a square is 90-deg symmetric), `C = 1-u^2 >= 0`, `S = 2u >= 0`, `N = 1+u^2 > 0`,
//! `cos = C/N`, `sin = S/N`, `w(theta) = cos+sin = (C+S)/N`.  With `a = p_x-x`, `b = p_y-y`,
//!
//!     X = (aC + bS)/N,   Y = (-aS + bC)/N,   p in Q  <=>  |X| <= 1/2 and |Y| <= 1/2.
//!
//! Multiplying by `2N > 0` gives the four **violation polynomials** of `RUNG2.md` sec 6.2
//!
//!     G_{p,0} =  2aC + 2bS - N      G_{p,1} = -2aC - 2bS - N
//!     G_{p,2} = -2aS + 2bC - N      G_{p,3} =  2aS - 2bC - N
//!
//! and `p in Q(x,y,u)` iff all four are `<= 0`.  Admissibility (`RUNG2.md` sec 1) is the same
//! kind of object: `x >= w/2` iff `A_1 := C+S-2xN <= 0`, `x <= m-w/2` iff
//! `A_2 := 2xN-2mN+C+S <= 0`, and `A_3`, `A_4` likewise in `y`.  **Every** inequality the checker
//! ever reasons about therefore lies in the class
//!
//!     P4 = { F(x,y,u) : affine in (x,y), degree <= 4 in u, rational coefficients },
//!
//! which is closed under nonnegative combinations.  Two facts drive everything:
//!
//! * **Lemma D (exact maximum).**  For `F in P4` and a pose box
//!   `B = [x0,x1] x [y0,y1] x [u0,u1]`, `max_B F` is attained at one of the four corners of the
//!   centre rectangle (for each fixed `u`, `F` is affine in `(x,y)`), so
//!   `max_B F = max_{4 corners} max_{u in [u0,u1]} F`.  The inner maximum of a univariate
//!   polynomial of degree `<= 2` is exact (endpoints, plus the vertex when the leading
//!   coefficient is negative and the vertex is interior); for degree `3, 4` the Bernstein
//!   convex-hull bound on `[u0,u1]` (`RUNG2.md` Lemma C), refined by bisection, is a *sound
//!   over-estimate*, so a certification is never created by slack, only lost.
//!
//! * **Lemma I (inference).**  Let the node's hypotheses be `F_i <= 0` (`F_i in P4`).  If there
//!   are rationals `mu_i >= 0` with `max_B (F - sum_i mu_i F_i) <= 0`, then `F <= 0` at every
//!   pose of `B` satisfying the hypotheses.  (`F <= sum mu_i F_i <= 0`.)
//!
//! Lemma I with the admissibility hypotheses reproduces `CORE` (all `mu = 0`), `P1`
//! (`mu = 1`) and `ADM`/Lemma A (`mu = cos` resp. `sin`, i.e. `mu = C/N`, `S/N` --- which is why
//! degree 4 is needed at all), *and* `clip_bin` (`RUNG2.md` sec 4.6: the admissibility hypotheses
//! are part of every node, so a box whose bin is almost entirely inadmissible is not tested over
//! the inadmissible part).  Lemma I with a *branch* hypothesis `+-G_{q,k} <= 0` reproduces
//! `CHAIN`'s Lemmas F, G and H (sec 6.2): see `DISJ` below.
//!
//! ## The search
//!
//! Root boxes: `x, y` on a pitch of `1/10` over `[0,m]`, `u` in 8 bins of `[0,1]`; the full
//! admissible domain, no symmetry reduction.  A box is a leaf when
//!
//!   * `EMPTY`  --- some admissibility hypothesis `A_i` satisfies `max_B(-A_i) < 0`, i.e. `A_i > 0`
//!     on all of `B`: no pose of `B` is admissible, nothing to prove;
//!   * `ADM`    --- the points `p` whose four conditions are all certified by Lemma I from the
//!     admissibility hypotheses alone have total weight `>= 1` (a *monotone witness certificate*
//!     in the sense of `RUNG2.md` Theorem 1);
//!   * `DISJ`   --- a binary tree of sign splits on violation polynomials covers `B`, and every
//!     leaf of that tree has its own Lemma-I witness set of weight `>= 1`.  Theorem 1 says a
//!     cover of weight `< m^2 = 16` *must* be certified this way almost everywhere.
//!
//! otherwise the box is halved (longest scaled side, with an angle bias at the walls) until the
//! depth limit, below which it is reported `UNCERTIFIED` with its exact coordinates.
//!
//! `DISJ` in detail.  Branching on the sign of a violation polynomial `G_q` splits `B` into the
//! two **closed** halves `{G_q <= 0}` and `{G_q >= 0}`, which cover `B`; each becomes a
//! hypothesis of its child.  In `{G_q <= 0}` the point `q` gains its missing condition; in
//! `{G_q >= 0}` other points gain theirs by Lemma I with `mu > 0` (this is `RUNG2.md` Lemma G),
//! and a child whose hypotheses force `G_q = G_{q'} = 0` certifies *both* `q` and `q'` (Lemma H
//! of `RUNG2.md` becomes a special case of Lemma I, so no "empty region" test is needed).
//!
//! ## What is exact and what is float
//!
//! Every fact that certifies a box is an `i128` integer computation: the coefficients of all
//! polynomials are integers after clearing the per-box denominators (centre denominator `Dc`, a
//! power of two times the certificate's `D`; angle denominator `M`, a power of two), and the
//! univariate tests are integer comparisons.  Floats appear only as (i) the reach filter that
//! decides which points to test at all (a documented superset: a point is kept unless
//! `|p-c|_inf > 0.7072 >= w/2` for every centre of the box), and (ii) the heuristics that order
//! the candidate list, choose the split dimension and choose the branch polynomial.  None of
//! them can create a certification.  The integer magnitudes are bounded a priori by the box's
//! denominators; a box whose denominators could overflow `i128` is refused rather than trusted
//! (see `SAFE_BITS`), so no result depends on wrap-around.
//!
//! ## Usage
//!
//!     zmcheck cert  FILE [--depth D] [--threads T] [--nodisj] [--dump F] [--xlo A --xhi B]
//!     zmcheck pose  FILE --x NUM/DEN --y NUM/DEN --u NUM/DEN
//!     zmcheck box   FILE --box "x0,x1,y0,y1,u0,u1"     (rationals; one box, verbose)
//!
//! Verdicts, in the style of `tests/rejection_tests.sh`: `VERIFIED:` + exit 0 when 0 boxes are
//! uncertified; `NOT VERIFIED` + exit 0 when some are (a reasoned refusal); `ERROR:` + exit 2 for
//! a malformed certificate (no verdict word is printed).

use std::fmt::Write as _;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

type I = i128;

// ---------------------------------------------------------------- univariate integer polynomials
//
// A degree-<=4 polynomial in the *local* bin parameter t in [0,1], u = (v0 + h t)/M.  Using the
// local parameter (rather than u itself) is what keeps the coefficient magnitudes small: see
// SAFE_BITS.

const DEG: usize = 5;
type Pu = [I; DEG];

#[inline]
fn padd(a: &Pu, b: &Pu) -> Pu {
    let mut r = [0; DEG];
    for i in 0..DEG {
        r[i] = a[i] + b[i];
    }
    r
}
#[inline]
fn psub(a: &Pu, b: &Pu) -> Pu {
    let mut r = [0; DEG];
    for i in 0..DEG {
        r[i] = a[i] - b[i];
    }
    r
}
#[inline]
fn pscale(a: &Pu, k: I) -> Pu {
    let mut r = [0; DEG];
    for i in 0..DEG {
        r[i] = a[i] * k;
    }
    r
}
/// product of two polynomials; panics (in debug) if the product would exceed degree 4.
#[inline]
fn pmul(a: &Pu, b: &Pu) -> Pu {
    let mut r = [0; DEG];
    for i in 0..DEG {
        if a[i] == 0 {
            continue;
        }
        for j in 0..DEG {
            if b[j] == 0 {
                continue;
            }
            debug_assert!(i + j < DEG, "degree overflow in pmul");
            r[i + j] += a[i] * b[j];
        }
    }
    r
}
#[inline]
fn pdeg(a: &Pu) -> usize {
    let mut d = 0;
    for i in 0..DEG {
        if a[i] != 0 {
            d = i;
        }
    }
    d
}

/// `max_{t in [0,1]} q(t) <= 0` --- exact for `deg <= 2`, a sound over-estimate for `deg 3,4`.
fn max_le0(q: &Pu) -> bool {
    max_bound_le(q, false)
}
/// `max_{t in [0,1]} q(t) < 0` --- sound in the same direction (never says `< 0` when it is not).
fn max_lt0(q: &Pu) -> bool {
    max_bound_le(q, true)
}

fn max_bound_le(q: &Pu, strict: bool) -> bool {
    let d = pdeg(q);
    // endpoints first (exact in every degree)
    let at0 = q[0];
    let at1: I = q.iter().sum();
    if strict {
        if at0 >= 0 || at1 >= 0 {
            return false;
        }
    } else if at0 > 0 || at1 > 0 {
        return false;
    }
    if d <= 1 {
        return true; // affine: the maximum is at an endpoint
    }
    if d == 2 {
        let (c0, c1, c2) = (q[0], q[1], q[2]);
        if c2 < 0 && c1 > 0 && c1 + 2 * c2 < 0 {
            // interior vertex at t* = -c1/(2c2); max value = (4 c0 c2 - c1^2)/(4 c2), 4c2 < 0
            let lhs = 4 * c0 * c2;
            let rhs = c1 * c1;
            return if strict { lhs > rhs } else { lhs >= rhs };
        }
        return true;
    }
    // degree 3 or 4: Bernstein convex-hull bound of RUNG2.md Lemma C on [0,1], refined by
    // bisection.  b_i = sum_{j<=i} c_j C(i,j)/C(4,j); scaled by 12 = lcm(1,4,6) to stay integral.
    bern_le(q, strict, 0)
}

const BERN_REFINE: u32 = 2; // 4 subintervals; each level multiplies the coefficients by 16

fn bern_le(q: &Pu, strict: bool, level: u32) -> bool {
    // 12/C(4,j) for j = 0..4
    const F: [I; 5] = [12, 3, 2, 3, 12];
    const C: [[I; 5]; 5] = [
        [1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [1, 3, 3, 1, 0],
        [1, 4, 6, 4, 1],
    ];
    let mut ok = true;
    for i in 0..5 {
        let mut b: I = 0;
        for j in 0..=i {
            b += q[j] * C[i][j] * F[j];
        }
        if strict {
            if b >= 0 {
                ok = false;
                break;
            }
        } else if b > 0 {
            ok = false;
            break;
        }
    }
    if ok {
        return true;
    }
    if level >= BERN_REFINE {
        return false;
    }
    // split [0,1]: q(s/2) and q((1+s)/2), both scaled by 16 to stay integral
    let mut lo = [0 as I; DEG];
    for j in 0..DEG {
        lo[j] = q[j] << (4 - j as u32); // c_j * 2^(4-j)
    }
    let mut hi = [0 as I; DEG];
    for j in 0..DEG {
        if q[j] == 0 {
            continue;
        }
        let k = q[j] << (4 - j as u32);
        for i in 0..=j {
            hi[i] += k * C[j][i]; // (1+s)^j expanded
        }
    }
    // endpoint pre-checks are cheap and exact
    let ok_lo = {
        let a0 = lo[0];
        let a1: I = lo.iter().sum();
        if strict {
            a0 < 0 && a1 < 0
        } else {
            a0 <= 0 && a1 <= 0
        }
    };
    let ok_hi = {
        let a0 = hi[0];
        let a1: I = hi.iter().sum();
        if strict {
            a0 < 0 && a1 < 0
        } else {
            a0 <= 0 && a1 <= 0
        }
    };
    ok_lo && ok_hi && bern_le(&lo, strict, level + 1) && bern_le(&hi, strict, level + 1)
}

// ---------------------------------------------------------------------------- class P4 elements
//
// F(Ax, Ay, t) = ax(t) * Ax + ay(t) * Ay + c(t), where the centre is (Ax/Dc, Ay/Dc).

#[derive(Clone, Copy)]
struct F4 {
    ax: Pu,
    ay: Pu,
    c: Pu,
}

impl F4 {
    fn zero() -> F4 {
        F4 { ax: [0; DEG], ay: [0; DEG], c: [0; DEG] }
    }
    fn add(&self, o: &F4) -> F4 {
        F4 { ax: padd(&self.ax, &o.ax), ay: padd(&self.ay, &o.ay), c: padd(&self.c, &o.c) }
    }
    fn sub(&self, o: &F4) -> F4 {
        F4 { ax: psub(&self.ax, &o.ax), ay: psub(&self.ay, &o.ay), c: psub(&self.c, &o.c) }
    }
    fn scale(&self, k: I) -> F4 {
        F4 { ax: pscale(&self.ax, k), ay: pscale(&self.ay, k), c: pscale(&self.c, k) }
    }
    /// multiply by a polynomial in `t` alone (used with `Chat`, `Shat`, `Nhat`, all `>= 0`)
    fn pmulp(&self, p: &Pu) -> F4 {
        F4 { ax: pmul(&self.ax, p), ay: pmul(&self.ay, p), c: pmul(&self.c, p) }
    }
    fn at(&self, axv: I, ayv: I) -> Pu {
        let mut r = [0; DEG];
        for i in 0..DEG {
            r[i] = self.ax[i] * axv + self.ay[i] * ayv + self.c[i];
        }
        r
    }
}

// ------------------------------------------------------------------------------------- pose box

#[derive(Clone, Copy)]
struct Bx {
    dc: I,   // centre denominator: 1000 * 2^dxy
    ax0: I,
    ax1: I,
    ay0: I,
    ay1: I,
    m: I, // angle denominator: 8 * 2^du
    v0: I,
    h: I,
    dxy: u32,
    du: u32,
    depth: u32,
}

/// Pre-computed per-box polynomials `Chat = M^2 C`, `Shat = M^2 S`, `Nhat = M^2 N` in `t`,
/// the four admissibility hypotheses and their products with `Chat`, `Shat`, `Nhat` (the
/// multiplier candidates of Lemma I), plus the float data of the pre-screen.
struct BoxPoly {
    ch: Pu,
    sh: Pu,
    nh: Pu,
    a: [F4; 4],
    /// `a_mul[i][j] = A_i * (Chat, Shat, Nhat)[j]`
    a_mul: [[F4; 3]; 4],
    /// float samples of the bin: `(cos, sin, w)` at `u0`, the midpoint and `u1`
    samp: Vec<(f64, f64, f64)>,
    /// sampled *admissible* poses of the box: `(x, y, cos, sin)` at the four centre corners
    /// clipped to the admissible range, crossed with five angles of the bin.  Used only by the
    /// float pre-screens, which can lose a certification but never create one.
    poses: Vec<(f64, f64, f64, f64)>,
    x0f: f64,
    x1f: f64,
    y0f: f64,
    y1f: f64,
    wall: bool,
}

fn box_poly_full(b: &Bx, mm: I, wall: bool) -> BoxPoly {
    let mut bp = box_poly(b);
    let mut a = [F4::zero(); 4];
    let mut a_mul = [[F4::zero(); 3]; 4];
    for i in 0..4 {
        a[i] = apoly(i, mm, b, &bp);
        a_mul[i][0] = a[i].pmulp(&bp.ch);
        a_mul[i][1] = a[i].pmulp(&bp.sh);
        a_mul[i][2] = a[i].pmulp(&bp.nh);
    }
    bp.a = a;
    bp.a_mul = a_mul;
    let dcf = b.dc as f64;
    bp.x0f = b.ax0 as f64 / dcf;
    bp.x1f = b.ax1 as f64 / dcf;
    bp.y0f = b.ay0 as f64 / dcf;
    bp.y1f = b.ay1 as f64 / dcf;
    bp.wall = wall;
    let mf = b.m as f64;
    for s in 0..5 {
        let u = (b.v0 as f64 + b.h as f64 * (s as f64) / 4.0) / mf;
        let n = 1.0 + u * u;
        let c = (1.0 - u * u) / n;
        let si = 2.0 * u / n;
        bp.samp.push((c, si, c + si));
    }
    let mf = mm as f64;
    for &(c, si, w) in &bp.samp {
        let lo = w / 2.0;
        let hi = mf - w / 2.0;
        for &xc in &[bp.x0f, bp.x1f] {
            for &yc in &[bp.y0f, bp.y1f] {
                let x = xc.max(lo).min(hi);
                let y = yc.max(lo).min(hi);
                if x >= bp.x0f - 1e-12
                    && x <= bp.x1f + 1e-12
                    && y >= bp.y0f - 1e-12
                    && y <= bp.y1f + 1e-12
                    && lo <= hi
                {
                    bp.poses.push((x, y, c, si));
                }
            }
        }
    }
    bp
}

/// `G_{p,k}/(2N)` in floats at one pose --- the geometric quantity `X - 1/2` etc.
#[inline]
fn gval(k: u8, px: f64, py: f64, pose: (f64, f64, f64, f64)) -> f64 {
    let (x, y, c, s) = pose;
    let (a, b) = (px - x, py - y);
    match k {
        0 => a * c + b * s - 0.5,
        1 => -(a * c + b * s) - 0.5,
        2 => -a * s + b * c - 0.5,
        _ => a * s - b * c - 0.5,
    }
}

fn box_poly(b: &Bx) -> BoxPoly {
    let (m, v0, h) = (b.m, b.v0, b.h);
    let m2 = m * m;
    // (v0 + h t)^2 = v0^2 + 2 v0 h t + h^2 t^2
    let sq: Pu = [v0 * v0, 2 * v0 * h, h * h, 0, 0];
    let mut ch = [0; DEG];
    ch[0] = m2;
    let ch = psub(&ch, &sq);
    let mut nh = [0; DEG];
    nh[0] = m2;
    let nh = padd(&nh, &sq);
    let sh: Pu = [2 * m * v0, 2 * m * h, 0, 0, 0];
    BoxPoly {
        ch,
        sh,
        nh,
        a: [F4::zero(); 4],
        a_mul: [[F4::zero(); 3]; 4],
        samp: Vec::new(),
        poses: Vec::new(),
        x0f: 0.0,
        x1f: 0.0,
        y0f: 0.0,
        y1f: 0.0,
        wall: false,
    }
}

/// `Ghat_{p,k} = Dc * M^2 * G_{p,k}` as an element of `P4`.  `xs = Dc * p_x`, `ys = Dc * p_y`.
fn gpoly(k: usize, xs: I, ys: I, b: &Bx, bp: &BoxPoly) -> F4 {
    let (ch, sh, nh) = (&bp.ch, &bp.sh, &bp.nh);
    let n_dc = pscale(nh, b.dc);
    match k {
        0 => F4 {
            ax: pscale(ch, -2),
            ay: pscale(sh, -2),
            c: psub(&padd(&pscale(ch, 2 * xs), &pscale(sh, 2 * ys)), &n_dc),
        },
        1 => F4 {
            ax: pscale(ch, 2),
            ay: pscale(sh, 2),
            c: psub(&psub(&pscale(ch, -2 * xs), &pscale(sh, 2 * ys)), &n_dc),
        },
        2 => F4 {
            ax: pscale(sh, 2),
            ay: pscale(ch, -2),
            c: psub(&padd(&pscale(sh, -2 * xs), &pscale(ch, 2 * ys)), &n_dc),
        },
        _ => F4 {
            ax: pscale(sh, -2),
            ay: pscale(ch, 2),
            c: psub(&padd(&pscale(sh, 2 * xs), &pscale(ch, -2 * ys)), &n_dc),
        },
    }
}

/// The four admissibility hypotheses, scaled the same way (`Dc * M^2 * A_i`):
/// `A_1 = C+S-2xN`, `A_2 = 2xN-2mN+C+S`, `A_3`, `A_4` the same in `y`.  `mm = m` (container side).
fn apoly(i: usize, mm: I, b: &Bx, bp: &BoxPoly) -> F4 {
    let cs = padd(&bp.ch, &bp.sh);
    let cs_dc = pscale(&cs, b.dc);
    let n2 = pscale(&bp.nh, 2);
    let n2m = pscale(&bp.nh, 2 * mm * b.dc);
    match i {
        0 => F4 { ax: pscale(&n2, -1), ay: [0; DEG], c: cs_dc },
        1 => F4 { ax: n2, ay: [0; DEG], c: psub(&cs_dc, &n2m) },
        2 => F4 { ax: [0; DEG], ay: pscale(&n2, -1), c: cs_dc },
        _ => F4 { ax: [0; DEG], ay: n2, c: psub(&cs_dc, &n2m) },
    }
}

/// Lemma D: `max_B F <= 0` (resp. `< 0`), by the four corners of the centre rectangle.
fn fmax_le0(f: &F4, b: &Bx) -> bool {
    for &axv in &[b.ax0, b.ax1] {
        for &ayv in &[b.ay0, b.ay1] {
            if !max_le0(&f.at(axv, ayv)) {
                return false;
            }
        }
    }
    true
}
fn fmax_lt0(f: &F4, b: &Bx) -> bool {
    for &axv in &[b.ax0, b.ax1] {
        for &ayv in &[b.ay0, b.ay1] {
            if !max_lt0(&f.at(axv, ayv)) {
                return false;
            }
        }
    }
    true
}

// --------------------------------------------------------------------------------- certificate

struct Cert {
    m_num: I,
    m_den: I,
    d: I,
    w: I,
    xs: Vec<I>, // X numerators over D
    ys: Vec<I>,
    wt: Vec<I>, // weight numerators over W
    total: I,   // sum of wt
}

fn load(path: &str) -> Result<Cert, String> {
    let txt = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let mut it = txt.split_ascii_whitespace();
    let mut next = |what: &str| -> Result<I, String> {
        let t = it.next().ok_or_else(|| format!("truncated file: expected {what}"))?;
        t.parse::<I>().map_err(|_| format!("expected integer for {what}, got {t:?}"))
    };
    let m_num = next("s_num")?;
    let m_den = next("s_den")?;
    let d = next("D")?;
    let w = next("W")?;
    let n = next("point count")?;
    if m_num <= 0 || m_den <= 0 || d <= 0 || w <= 0 {
        return Err("header values must be positive integers".into());
    }
    if n < 0 {
        return Err("negative point count".into());
    }
    // the container side must be a multiple of the coordinate unit (FORMAT.md well-formedness)
    if (m_num * d) % m_den != 0 {
        return Err(format!("s_den={m_den} does not divide s_num*D = {m_num}*{d}"));
    }
    let side_units = m_num * d / m_den; // container side in units of 1/D
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let mut wt = Vec::new();
    let mut total: I = 0;
    for i in 0..n {
        let x = next(&format!("X of point {}", i + 1))?;
        let y = next(&format!("Y of point {}", i + 1))?;
        let ww = next(&format!("w of point {}", i + 1))?;
        if ww < 0 {
            return Err(format!("point {} has negative weight {ww}", i + 1));
        }
        if x < 0 || y < 0 || x > side_units || y > side_units {
            return Err(format!(
                "point {} at ({x},{y})/{d} is outside the container [0,{}]^2",
                i + 1,
                side_units
            ));
        }
        xs.push(x);
        ys.push(y);
        wt.push(ww);
        total += ww;
    }
    if it.next().is_some() {
        return Err("trailing data after the last point line".into());
    }
    Ok(Cert { m_num, m_den, d, w, xs, ys, wt, total })
}

// ------------------------------------------------------------------------------------- checker

struct Checker {
    cert: Cert,
    mm: I,       // container side (integer; the checker requires an integer side)
    depth_max: u32,
    disj: bool,
    /// grid index of the points: cell (i,j) of pitch 1/10 -> point indices
    grid: Vec<Vec<u32>>,
    gn: usize,
    theta_bias: f64,
    branch_cap: usize,
    node_cap: usize,
    sign_depth: u32,
    try_branches: usize,
    rank_exact: usize,
    fail_cap: usize,
    seed_cap: usize,
    heur_cap: usize,
    noscreen: bool,
}

const SAFE_BITS: u32 = 74; // guard on 4*du + dxy (see the module header)

#[derive(Default)]
struct Census {
    boxes: usize,
    adm: usize,
    disj: usize,
    empty: usize,
    uncert: usize,
    maxdepth: u32,
    unsafe_boxes: usize,
    uncert_list: Vec<String>,
    dump: String,
}

impl Census {
    fn merge(&mut self, o: Census) {
        self.boxes += o.boxes;
        self.adm += o.adm;
        self.disj += o.disj;
        self.empty += o.empty;
        self.uncert += o.uncert;
        self.unsafe_boxes += o.unsafe_boxes;
        self.maxdepth = self.maxdepth.max(o.maxdepth);
        self.uncert_list.extend(o.uncert_list);
        self.dump.push_str(&o.dump);
    }
}

/// Per-box working state.  `cmask[i]` has bit `k` set when condition `k` of point `cand[i]` is
/// certified over the whole box by Lemma I from the admissibility hypotheses alone; `cand` is
/// sorted, so a child can inherit its parent's bits (its poses are a subset, so the parent's
/// Lemma-I certificate is still a certificate).
struct Work {
    cand: Vec<u32>,
    cmask: Vec<u8>,
    tset: Vec<u32>, // points with cmask == 0b1111 (a monotone witness set)
    twt: I,
}

impl Checker {
    fn wtarget(&self) -> I {
        self.cert.w
    }

    fn build_grid(cert: &Cert, mm: I) -> (Vec<Vec<u32>>, usize) {
        let gn = (mm * 10) as usize;
        let mut grid = vec![Vec::new(); gn * gn];
        for (i, (&x, &y)) in cert.xs.iter().zip(cert.ys.iter()).enumerate() {
            let ci = ((x * 10 / cert.d) as usize).min(gn - 1);
            let cj = ((y * 10 / cert.d) as usize).min(gn - 1);
            grid[ci * gn + cj].push(i as u32);
        }
        (grid, gn)
    }

    /// Points that could possibly be captured somewhere in the box.  Float superset: a captured
    /// point has `|p-c|_inf <= w/2 <= sqrt(2)/2 = 0.70711 < 0.7072`.
    fn reach(&self, b: &Bx) -> Vec<u32> {
        const R: f64 = 0.7072;
        let dcf = b.dc as f64;
        let (x0, x1) = (b.ax0 as f64 / dcf - R, b.ax1 as f64 / dcf + R);
        let (y0, y1) = (b.ay0 as f64 / dcf - R, b.ay1 as f64 / dcf + R);
        let gi0 = (x0 * 10.0).floor().max(0.0) as usize;
        let gi1 = ((x1 * 10.0).ceil() as i64).max(0) as usize;
        let gj0 = (y0 * 10.0).floor().max(0.0) as usize;
        let gj1 = ((y1 * 10.0).ceil() as i64).max(0) as usize;
        let mut out = Vec::new();
        let df = self.cert.d as f64;
        for i in gi0..=gi1.min(self.gn - 1) {
            for j in gj0..=gj1.min(self.gn - 1) {
                for &pi in &self.grid[i * self.gn + j] {
                    let px = self.cert.xs[pi as usize] as f64 / df;
                    let py = self.cert.ys[pi as usize] as f64 / df;
                    if px >= x0 && px <= x1 && py >= y0 && py <= y1 && self.cert.wt[pi as usize] > 0
                    {
                        out.push(pi);
                    }
                }
            }
        }
        out
    }

    /// Does the wall possibly bind on this box?  (`w_hi = max w` over the bin, `<= sqrt 2`.)
    fn wall_touch(&self, b: &Bx) -> bool {
        let whi = 1.41422_f64;
        let dcf = b.dc as f64;
        (b.ax0 as f64 / dcf) < whi / 2.0
            || (b.ax1 as f64 / dcf) > self.mm as f64 - whi / 2.0
            || (b.ay0 as f64 / dcf) < whi / 2.0
            || (b.ay1 as f64 / dcf) > self.mm as f64 - whi / 2.0
    }

    /// The indices of the admissibility hypothesis and the eliminating multiplier for each
    /// condition (Lemma A of `RUNG2.md` sec 3.1 read as a Lemma-I multiplier choice):
    ///   k=0: x -> A_1 with C,  y -> A_3 with S       k=1: x -> A_2 with C,  y -> A_4 with S
    ///   k=2: x -> A_2 with S,  y -> A_3 with C       k=3: x -> A_1 with S,  y -> A_4 with C
    /// (`ch = 0`, `sh = 1`, `nh = 2` index `BoxPoly::a_mul`.)
    #[inline]
    fn wall_slots(k: usize) -> (usize, usize, usize, usize) {
        match k {
            0 => (0, 0, 2, 1),
            1 => (1, 0, 3, 1),
            2 => (1, 1, 2, 0),
            _ => (0, 1, 3, 0),
        }
    }

    /// Lemma I with the admissibility hypotheses only: is condition `k` of point `pi` certified
    /// over the whole box?  (`ADM` of `RUNG2.md` sec 3, which contains `CORE` and `P1`.)
    fn cond_adm(&self, pi: u32, k: usize, b: &Bx, bp: &BoxPoly) -> bool {
        let sc = b.dc / self.cert.d;
        let xs = self.cert.xs[pi as usize] * sc;
        let ys = self.cert.ys[pi as usize] * sc;
        let g = gpoly(k, xs, ys, b, bp);
        // mu = 0: the plain box bound (this is CORE)
        if fmax_le0(&g, b) {
            return true;
        }
        if !bp.wall {
            return false;
        }
        let (ix, mx, iy, my) = Self::wall_slots(k);
        let gn = g.pmulp(&bp.nh); // target Nhat * Ghat, so that mu = C/N, S/N, 1 are available
        for ox in 0..3 {
            // mu_x in {0, the eliminating multiplier, 1}
            let t1 = if ox == 0 { gn } else { gn.sub(&bp.a_mul[ix][if ox == 1 { mx } else { 2 }]) };
            for oy in 0..3 {
                let t2 =
                    if oy == 0 { t1 } else { t1.sub(&bp.a_mul[iy][if oy == 1 { my } else { 2 }]) };
                if fmax_le0(&t2, b) {
                    return true;
                }
            }
        }
        false
    }

    /// Float pre-screen (a documented superset of the exact test, so it can only *lose*
    /// certifications).  Bit `k` of the result is set when the *ideal* Lemma-A bound --- the
    /// exact maximum of condition `k` over the admissible poses of the box, which every Lemma-I
    /// multiplier choice over-estimates --- is positive at some sampled angle of the bin: then
    /// there is an admissible pose of the box that does not capture `p`, and no multiplier
    /// choice can certify condition `k`.  A sampled angle with an empty admissible centre range
    /// is skipped (that is `clip_bin`, `RUNG2.md` sec 4.6, in float form).
    fn screen(&self, px: f64, py: f64, bp: &BoxPoly) -> u8 {
        const TOL: f64 = 1e-9;
        if self.noscreen {
            return 0;
        }
        let mf = self.mm as f64;
        let mut fail = 0u8;
        for &(c, s, w) in &bp.samp {
            let axx = bp.x0f.max(w / 2.0);
            let bxx = bp.x1f.min(mf - w / 2.0);
            let ayy = bp.y0f.max(w / 2.0);
            let byy = bp.y1f.min(mf - w / 2.0);
            if axx > bxx || ayy > byy {
                continue; // no admissible pose at this angle
            }
            if (px - axx) * c + (py - ayy) * s - 0.5 > TOL {
                fail |= 1;
            }
            if -((px - bxx) * c + (py - byy) * s) - 0.5 > TOL {
                fail |= 2;
            }
            if -(px - bxx) * s + (py - ayy) * c - 0.5 > TOL {
                fail |= 4;
            }
            if (px - axx) * s - (py - byy) * c - 0.5 > TOL {
                fail |= 8;
            }
        }
        fail
    }

    fn prepare(&self, b: &Bx, bp: &BoxPoly, parent: Option<&Work>) -> Work {
        let (cand, inh): (Vec<u32>, Vec<u8>) = match parent {
            None => {
                let c = self.reach(b);
                let n = c.len();
                (c, vec![0u8; n])
            }
            Some(w) => {
                let mut c = Vec::with_capacity(w.cand.len());
                let mut m = Vec::with_capacity(w.cand.len());
                for (i, &pi) in w.cand.iter().enumerate() {
                    if self.in_reach(pi, b) {
                        c.push(pi);
                        m.push(w.cmask[i]);
                    }
                }
                (c, m)
            }
        };
        let df = self.cert.d as f64;
        let mut cmask = Vec::with_capacity(cand.len());
        let mut tset: Vec<u32> = Vec::new();
        let mut twt: I = 0;
        for (i, &pi) in cand.iter().enumerate() {
            let mut m = inh[i];
            if m != 0b1111 {
                let px = self.cert.xs[pi as usize] as f64 / df;
                let py = self.cert.ys[pi as usize] as f64 / df;
                let sf = self.screen(px, py, bp);
                for k in 0..4usize {
                    let bit = 1u8 << k;
                    if m & bit != 0 || sf & bit != 0 {
                        continue;
                    }
                    if self.cond_adm(pi, k, b, bp) {
                        m |= bit;
                    }
                }
            }
            cmask.push(m);
            if m == 0b1111 {
                tset.push(pi);
                twt += self.cert.wt[pi as usize];
            }
        }
        Work { cand, cmask, tset, twt }
    }

    fn in_reach(&self, pi: u32, b: &Bx) -> bool {
        const R: f64 = 0.7072;
        let dcf = b.dc as f64;
        let df = self.cert.d as f64;
        let px = self.cert.xs[pi as usize] as f64 / df;
        let py = self.cert.ys[pi as usize] as f64 / df;
        px >= b.ax0 as f64 / dcf - R
            && px <= b.ax1 as f64 / dcf + R
            && py >= b.ay0 as f64 / dcf - R
            && py <= b.ay1 as f64 / dcf + R
    }

    /// `EMPTY`: some admissibility hypothesis is violated by *every* pose of the box.
    fn empty(&self, b: &Bx, bp: &BoxPoly) -> bool {
        for i in 0..4 {
            // -A_i < 0 on all of B  <=>  A_i > 0 on all of B  =>  no admissible pose
            let neg = F4::zero().sub(&bp.a[i]);
            if fmax_lt0(&neg, b) {
                return true;
            }
        }
        false
    }
}

// --------------------------------------------------------------------------- disjunctive search

/// A branch hypothesis: the sign of `Ghat_{q,kq}`.  `le = true` means the hypothesis
/// `Ghat_{q,kq} <= 0`, `le = false` means `-Ghat_{q,kq} <= 0`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Hyp {
    q: u32,
    kq: u8,
    le: bool,
}

/// A node's certified-condition cover: the bitset over condition slots, plus a note of which
/// candidates' `le` masks the Lemma-K closure has already folded in (so the closure is not
/// re-applied from scratch at every node).
#[derive(Clone)]
struct Cov {
    m: Vec<u64>,
    added: Vec<bool>,
}

struct Disj<'a> {
    ck: &'a Checker,
    b: Bx,
    bp: &'a BoxPoly,
    /// points not in T, with the list of their failing conditions, and their float coordinates
    fails: Vec<(u32, Vec<u8>, f64, f64)>,
    /// branch candidates `(point, condition, x, y, slot)`; `slot = 4*fi+k` is the condition
    /// slot of the candidate's own (single) failing condition, i.e. the bit that says
    /// "`G_{q,k} <= 0` has been proved at this node"
    branches: Vec<(u32, u8, f64, f64, usize)>,
    /// float pre-screen, one bitset over condition slots `4*fi+k` per hypothesis (see
    /// `build_plaus`); `emask` is the same thing computed exactly, on demand.
    pmask: Vec<Vec<u64>>,
    emask: Vec<Option<Vec<u64>>>,
    /// bit `4*fi+k` set for every condition that `fails[fi]` still needs
    need: Vec<u64>,
    nw: usize,
    nb: usize,
    score: u8,
    noscreen: bool,
    dbg: bool,
    twt: I,
    target: I,
    nodes: usize,
    node_cap: usize,
    /// a forced prefix of branch decisions (a Lemma-K seed pair)
    forced: Vec<usize>,
    fail_reports: usize,
}

impl<'a> Disj<'a> {
    fn cond_under_raw(&self, p: u32, k: u8, h: Hyp) -> bool {
        let ck = self.ck;
        let b = &self.b;
        let bp = self.bp;
        let sc = b.dc / ck.cert.d;
        let gp = gpoly(k as usize, ck.cert.xs[p as usize] * sc, ck.cert.ys[p as usize] * sc, b, bp);
        let gq = gpoly(
            h.kq as usize,
            ck.cert.xs[h.q as usize] * sc,
            ck.cert.ys[h.q as usize] * sc,
            b,
            bp,
        );
        // hypothesis F_h <= 0 with F_h = +-Ghat_q; target Ghat_p - mu F_h <= 0, mu = n/d >= 0.
        //   h.le  : F_h =  Ghat_q  ->  test  d*Ghat_p - n*Ghat_q
        //   !h.le : F_h = -Ghat_q  ->  test  d*Ghat_p + n*Ghat_q
        for (li, &(dn, nn)) in [(1 as I, 1 as I), (2, 1), (1, 2)].iter().enumerate() {
            let comb = if h.le {
                gp.scale(dn).sub(&gq.scale(nn))
            } else {
                gp.scale(dn).add(&gq.scale(nn))
            };
            if fmax_le0(&comb, b) {
                return true;
            }
            // ... and with the wall multipliers of Lemma A on top (only for mu = 1: the extra
            // multipliers are what a wall needs, the extra mu values what a slanted cut needs)
            let _ = li;
            if !bp.wall {
                continue;
            }
            let (ix, mx, iy, my) = Checker::wall_slots(k as usize);
            let base = comb.pmulp(&bp.nh);
            for ox in 0..3 {
                let t1 =
                    if ox == 0 { base } else { base.sub(&bp.a_mul[ix][if ox == 1 { mx } else { 2 }]) };
                for oy in 0..3 {
                    if ox == 0 && oy == 0 {
                        continue;
                    }
                    let t2 =
                        if oy == 0 { t1 } else { t1.sub(&bp.a_mul[iy][if oy == 1 { my } else { 2 }]) };
                    if fmax_le0(&t2, b) {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Float pre-screen of `cond_under`: is there a *sampled admissible* pose of the box at which
    /// the hypothesis holds and condition `k` of `p` fails?  If so no Lemma-I certificate exists
    /// and the exact test is skipped.  (Sound in the certifying direction: a pose counts as a
    /// counterexample only when both signs are outside the float tolerance.)
    fn screen_under(&self, px: f64, py: f64, k: u8, qx: f64, qy: f64, kq: u8, le: bool) -> bool {
        const TOL: f64 = 1e-9;
        for &pose in &self.bp.poses {
            let gq = gval(kq, qx, qy, pose);
            let hyp_holds = if le { gq < -TOL } else { gq > TOL };
            if hyp_holds && gval(k, px, py, pose) > TOL {
                return false;
            }
        }
        true
    }

    /// Build the float pre-screen masks.  There is one bitset per hypothesis (branch candidate x
    /// sign) over the *condition slots* `4*fi + k`; bit set means "the screen cannot rule out
    /// that Lemma I certifies condition `k` of `fails[fi]` from this hypothesis".  The relation
    /// depends only on the box, not on the node, which is why it can be tabulated once.
    fn build_plaus(&mut self) {
        let nf = self.fails.len();
        let nb = self.branches.len();
        self.nb = nb;
        self.nw = (nf * 4 + 63) / 64;
        self.need = vec![0u64; self.nw];
        self.pmask = vec![vec![0u64; self.nw]; nb * 2];
        self.emask = vec![None; nb * 2];
        for fi in 0..nf {
            for &k in &self.fails[fi].1 {
                let b = fi * 4 + k as usize;
                self.need[b / 64] |= 1u64 << (b % 64);
            }
        }
        for fi in 0..nf {
            let (p, px, py) = (self.fails[fi].0, self.fails[fi].2, self.fails[fi].3);
            let ks = self.fails[fi].1.clone();
            for &k in &ks {
                let slot = fi * 4 + k as usize;
                for bi in 0..nb {
                    let (q, kq, qx, qy, _) = self.branches[bi];
                    if q == p && kq == k {
                        // the hypothesis G_{p,k} <= 0 certifies condition k of p outright
                        self.pmask[bi * 2][slot / 64] |= 1u64 << (slot % 64);
                        continue;
                    }
                    for (si, le) in [(0usize, true), (1, false)] {
                        if self.noscreen || self.screen_under(px, py, k, qx, qy, kq, le) {
                            self.pmask[bi * 2 + si][slot / 64] |= 1u64 << (slot % 64);
                        }
                    }
                }
            }
        }
    }

    /// The *exact* mask of one hypothesis, built on demand: bit `4*fi+k` is set iff Lemma I
    /// certifies condition `k` of `fails[fi]` from this hypothesis (and the admissibility
    /// hypotheses).  Only slots the float screen left open are tested.
    fn exact_mask(&mut self, hi: usize) -> &Vec<u64> {
        if self.emask[hi].is_none() {
            let nf = self.fails.len();
            let (bi, le) = (hi / 2, hi % 2 == 0);
            let (q, kq, _, _, _) = self.branches[bi];
            let mut m = vec![0u64; self.nw];
            for fi in 0..nf {
                let p = self.fails[fi].0;
                let ks = self.fails[fi].1.clone();
                for &k in &ks {
                    let slot = fi * 4 + k as usize;
                    if self.pmask[hi][slot / 64] & (1u64 << (slot % 64)) == 0 {
                        continue;
                    }
                    let ok = if q == p && kq == k {
                        le
                    } else {
                        self.cond_under_raw(p, k, Hyp { q, kq, le })
                    };
                    if ok {
                        m[slot / 64] |= 1u64 << (slot % 64);
                    }
                }
            }
            self.emask[hi] = Some(m);
        }
        self.emask[hi].as_ref().unwrap()
    }

    /// Total weight of `T` plus every point all of whose failing conditions are in `cover`.
    fn weight_of(&self, cover: &[u64]) -> I {
        let mut w = self.twt;
        for fi in 0..self.fails.len() {
            let mut all = true;
            for &k in &self.fails[fi].1 {
                let b = fi * 4 + k as usize;
                if cover[b / 64] & (1u64 << (b % 64)) == 0 {
                    all = false;
                    break;
                }
            }
            if all {
                w += self.ck.cert.wt[self.fails[fi].0 as usize];
            }
        }
        w
    }

    /// `w[0]` is the weight of the `G <= 0` child, `w[1]` that of the `G >= 0` child.
    #[inline]
    fn combine(&self, w: [I; 2]) -> I {
        match self.score {
            0 => w[0].min(w[1]),
            1 => w[1],
            _ => w[0] + w[1],
        }
    }

    fn or_into(dst: &mut [u64], src: &[u64]) {
        for i in 0..dst.len() {
            dst[i] |= src[i];
        }
    }

    #[inline]
    fn has(m: &[u64], slot: usize) -> bool {
        m[slot / 64] & (1u64 << (slot % 64)) != 0
    }

    /// **Lemma K closure.**  `active` lists the node's branch decisions as hypothesis indices
    /// `2*bi + si` (`si = 0` is `G_{q_bi} <= 0`, `si = 1` is `G_{q_bi} >= 0`).  Start from the
    /// union of their masks; then, whenever the node has the hypothesis `G_q >= 0` for a
    /// candidate `q` *and* the cover already proves `G_q <= 0` (its own slot is set --- which is
    /// exactly what `RUNG2.md` Lemma H establishes when two pivots cannot both be violated),
    /// `G_q = 0` holds at every pose of the node, so `G_q <= 0` may be *added as a hypothesis*
    /// and everything its `le` mask certifies becomes available too.  Iterated to a fixpoint.
    ///
    /// This is the step the earlier version was missing: it is what makes the witness-free
    /// quadrant of two independent sliding cuts (the interior tile poses) certifiable instead of
    /// merely "provably empty".
    fn cover(&mut self, active: &[usize], exact: bool) -> Vec<u64> {
        let mut cov = Cov { m: vec![0u64; self.nw], added: vec![false; self.branches.len()] };
        for i in 0..active.len() {
            cov = self.cover_child(&cov, &active[..=i], exact);
        }
        cov.m
    }

    /// One step of the closure: the parent's cover plus the mask of the hypothesis just pushed
    /// (`active.last()`), then the Lemma-K closure re-run over the node's `ge` hypotheses.
    fn cover_child(&mut self, parent: &Cov, active: &[usize], exact: bool) -> Cov {
        let mut cov = parent.clone();
        let hi = *active.last().unwrap();
        if exact {
            let m = self.exact_mask(hi).clone();
            Self::or_into(&mut cov.m, &m);
        } else {
            let m = self.pmask[hi].clone();
            Self::or_into(&mut cov.m, &m);
        }
        loop {
            let mut grew = false;
            for idx in 0..active.len() {
                let h = active[idx];
                if h % 2 == 0 {
                    continue; // already the `le` hypothesis
                }
                let bi = h / 2;
                if cov.added[bi] {
                    continue;
                }
                if Self::has(&cov.m, self.branches[bi].4) {
                    cov.added[bi] = true;
                    let le = bi * 2;
                    if exact {
                        let m = self.exact_mask(le).clone();
                        Self::or_into(&mut cov.m, &m);
                    } else {
                        let m = self.pmask[le].clone();
                        Self::or_into(&mut cov.m, &m);
                    }
                    grew = true;
                }
            }
            if !grew {
                break;
            }
        }
        cov
    }

    /// The chain order of `RUNG2.md` sec 6.2, in floats: within one kind `k`, sort the branch
    /// candidates by an estimate of `max_B G_{q,k}` ascending.  Lemma F wants `G_{q_1} <= ... <=
    /// G_{q_k}` on `B`, and this is the cheap proxy for it; the ordering is a heuristic, every
    /// certification is still the exact Lemma I.
    fn chain_order(&self) -> Vec<Vec<usize>> {
        let mut by_kind: Vec<Vec<(f64, usize)>> = vec![Vec::new(); 4];
        for bi in 0..self.branches.len() {
            let (_, kq, qx, qy, _) = self.branches[bi];
            let mut m = f64::NEG_INFINITY;
            for &pose in &self.bp.poses {
                m = m.max(gval(kq, qx, qy, pose));
            }
            by_kind[kq as usize].push((m, bi));
        }
        let mut out: Vec<Vec<usize>> = Vec::new();
        for g in by_kind.iter_mut() {
            if g.is_empty() {
                continue;
            }
            g.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            out.push(g.iter().map(|e| e.1).collect());
        }
        // the biggest families first: they are the sliding cuts
        out.sort_by_key(|v| -(v.len() as i64));
        out
    }

    /// Float-plausible Lemma-K pairs: `(bi, bj)` such that the screen cannot rule out that the
    /// hypothesis `G_{q_bj} >= 0` certifies `G_{q_bi} <= 0` (and symmetrically).  These are the
    /// seeds of the two-chain certificate; `search_seeded` forces them to the top of the tree.
    fn kpairs(&self) -> Vec<(usize, usize)> {
        let nb = self.branches.len();
        let mut out = Vec::new();
        for bi in 0..nb {
            for bj in 0..nb {
                if bi == bj {
                    continue;
                }
                if self.branches[bi].1 == self.branches[bj].1 {
                    continue; // same kind: a chain, not a product
                }
                if Self::has(&self.pmask[bj * 2 + 1], self.branches[bi].4) {
                    out.push((bi, bj));
                }
            }
        }
        out
    }

    fn search(&mut self, active: &mut Vec<usize>, pc: &Cov, ec: &Cov, depth: u32,
              regions: &mut usize) -> bool {
        self.nodes += 1;
        if self.nodes > self.node_cap {
            return false;
        }
        if self.weight_of(&pc.m) >= self.target && self.weight_of(&ec.m) >= self.target {
            *regions += 1;
            return true;
        }
        if depth == 0 {
            return false;
        }
        if active.len() < self.forced.len() {
            // a seeded prefix: no ranking needed, the pair to split on is fixed
            let bi = self.forced[active.len()];
            let mut ok = true;
            for si in 0..2 {
                let hi = bi * 2 + si;
                active.push(hi);
                let pcc = self.cover_child(pc, active, false);
                let ecc = self.cover_child(ec, active, true);
                let r = self.search(active, &pcc, &ecc, depth - 1, regions);
                active.pop();
                if !r {
                    ok = false;
                    break;
                }
            }
            if !ok && self.dbg && self.fail_reports < 6 {
                self.fail_reports += 1;
                let ecx = self.cover(active, true);
                let sg: String = active
                    .iter()
                    .map(|&h| {
                        format!("{}{}", if h % 2 == 0 { "-" } else { "+" }, self.branches[h / 2].0)
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                eprintln!(
                    "    FORCED FAIL depth={} w={:.6} signs=[{}]",
                    active.len(),
                    self.weight_of(&ecx) as f64 / self.ck.cert.w as f64,
                    sg
                );
            }
            return ok;
        }
        // rank the branch candidates: cheaply by the float masks, then exactly for the best few.
        // Heuristic only --- the acceptance above is the exact test.
        // Rank the branch candidates.  `score` selects the shape of disjunction being looked
        // for: `Min` (maximise the weaker child) finds balanced cuts, `Ge` (maximise the
        // "G >= 0" child) builds a *chain* -- the pivot whose far side is richest first, which
        // is exactly `RUNG2.md` sec 6.2's greedy chain -- and `Sum` is a compromise.  All three
        // are heuristics over the same exact primitive; `classify` tries them in turn.
        let mut scored: Vec<(I, usize)> = Vec::new();
        let mut tmp = vec![0u64; self.nw];
        for bi in 0..self.branches.len() {
            if active.contains(&(bi * 2)) || active.contains(&(bi * 2 + 1)) {
                continue;
            }
            let mut w = [0 as I; 2];
            for si in 0..2 {
                tmp.copy_from_slice(&pc.m);
                Self::or_into(&mut tmp, &self.pmask[bi * 2 + si]);
                w[si] = self.weight_of(&tmp);
            }
            scored.push((self.combine(w), bi));
        }
        scored.sort_by(|a, b| b.0.cmp(&a.0));
        scored.truncate(self.ck.rank_exact);
        for e in scored.iter_mut() {
            let bi = e.1;
            let mut w = [0 as I; 2];
            for si in 0..2 {
                let hm = self.exact_mask(bi * 2 + si).clone();
                tmp.copy_from_slice(&ec.m);
                Self::or_into(&mut tmp, &hm);
                w[si] = self.weight_of(&tmp);
            }
            e.0 = self.combine(w);
        }
        scored.sort_by(|a, b| b.0.cmp(&a.0));
        if self.dbg && active.is_empty() {
            for &(w, bi) in scored.iter().take(6) {
                eprintln!(
                    "    cand bi={} pt={} k={} exact-min={:.6}",
                    bi, self.branches[bi].0, self.branches[bi].1,
                    w as f64 / self.ck.cert.w as f64
                );
            }
        }
        let order: Vec<usize> = if active.len() < self.forced.len() {
            vec![self.forced[active.len()]]
        } else {
            scored.iter().take(self.ck.try_branches).map(|e| e.1).collect()
        };
        for bi in order {
            if active.contains(&(bi * 2)) || active.contains(&(bi * 2 + 1)) {
                continue;
            }
            let mut ok = true;
            for si in 0..2 {
                let hi = bi * 2 + si;
                active.push(hi);
                let pcc = self.cover_child(pc, active, false);
                let ecc = self.cover_child(ec, active, true);
                let r = self.search(active, &pcc, &ecc, depth - 1, regions);
                active.pop();
                if !r {
                    ok = false;
                    break;
                }
            }
            if ok {
                return true;
            }
        }
        false
    }
}

// ------------------------------------------------------------------------------- the box driver

enum Leaf {
    Empty,
    Adm(I),
    Disj(usize),
    Split,
    Unsafe,
}

impl Checker {
    fn classify(&self, b: &Bx, bp: &BoxPoly, wk: &Work) -> (Leaf, Option<String>) {
        if self.empty(b, bp) {
            return (Leaf::Empty, None);
        }
        if wk.twt >= self.wtarget() {
            return (Leaf::Adm(wk.twt), Some(format!("ADM w={} pts={}", wk.twt, wk.tset.len())));
        }
        if !self.disj {
            return (Leaf::Split, None);
        }
        // build the failing-condition lists and the branch candidates
        let df = self.cert.d as f64;
        let mut fails: Vec<(u32, Vec<u8>, f64, f64)> = Vec::new();
        for (i, &pi) in wk.cand.iter().enumerate() {
            let m = wk.cmask[i];
            if m == 0b1111 {
                continue;
            }
            let ks: Vec<u8> = (0..4u8).filter(|k| m & (1 << k) == 0).collect();
            if ks.len() <= 2 {
                fails.push((
                    pi,
                    ks,
                    self.cert.xs[pi as usize] as f64 / df,
                    self.cert.ys[pi as usize] as f64 / df,
                ));
            }
        }
        if fails.len() > self.fail_cap {
            // keep the heaviest swing points: the missing weight is at most what they carry, so
            // dropping the light tail can only lose a certification, never create one
            fails.sort_by_key(|f| -(self.cert.wt[f.0 as usize]));
            fails.truncate(self.fail_cap);
        }
        // an upper bound on what any disjunction could reach: T plus every swing point.  A
        // tighter float bound (every point reachable from *some* hypothesis) is applied after
        // the pre-screen is built, below.
        let mut ub = wk.twt;
        for f in &fails {
            ub += self.cert.wt[f.0 as usize];
        }
        if ub < self.wtarget() {
            return (Leaf::Split, None);
        }
        let mut branches: Vec<(u32, u8, f64, f64, usize)> = fails
            .iter()
            .enumerate()
            .filter(|(_, f)| f.1.len() == 1)
            .map(|(fi, f)| (f.0, f.1[0], f.2, f.3, fi * 4 + f.1[0] as usize))
            .collect();
        branches.sort_by_key(|b| -(self.cert.wt[b.0 as usize]));
        branches.truncate(self.branch_cap);
        if branches.is_empty() {
            return (Leaf::Split, None);
        }
        let mut d = Disj {
            ck: self,
            b: *b,
            bp,
            fails,
            branches,
            pmask: Vec::new(),
            emask: Vec::new(),
            need: Vec::new(),
            nw: 0,
            nb: 0,
            score: 0,
            noscreen: self.noscreen,
            dbg: std::env::var("ZM_DEBUG").is_ok(),
            twt: wk.twt,
            target: self.wtarget(),
            nodes: 0,
            node_cap: self.node_cap,
            forced: Vec::new(),
            fail_reports: 0,
        };
        d.build_plaus();
        {
            // float upper bound over *all* hypotheses at once: if even that cannot reach 1, no
            // sign-splitting tree can, so go straight to subdivision with no exact work
            let mut all = vec![0u64; d.nw];
            for hi in 0..d.pmask.len() {
                Disj::or_into(&mut all, &d.pmask[hi].clone());
            }
            let ubf = d.weight_of(&all);
            if std::env::var("ZM_UB").is_ok() {
                let mut alle = vec![0u64; d.nw];
                for hi in 0..d.pmask.len() {
                    let m = d.exact_mask(hi).clone();
                    Disj::or_into(&mut alle, &m);
                }
                eprintln!(
                    "  exact-ub over ALL hypotheses = {:.9}",
                    d.weight_of(&alle) as f64 / self.cert.w as f64
                );
            }
            if std::env::var("ZM_DEBUG").is_ok() {
                eprintln!(
                    "  DISJ: fails={} branches={} w(T)={:.9} float-ub={:.9}",
                    d.fails.len(),
                    d.branches.len(),
                    d.twt as f64 / self.cert.w as f64,
                    ubf as f64 / self.cert.w as f64
                );
            }
            if ubf < self.wtarget() {
                return (Leaf::Split, None);
            }
        }
        if std::env::var("ZM_DEBUG").is_ok() {
            eprintln!(
                "  DISJ: fails={} branches={} w(T)={}",
                d.fails.len(),
                d.branches.len(),
                d.twt
            );
        }
        let zero = Cov { m: vec![0u64; d.nw], added: vec![false; d.branches.len()] };
        let z2 = zero.clone();
        let mut total_nodes = 0usize;
        for score in 0..3u8 {
            let mut active: Vec<usize> = Vec::new();
            let mut regions = 0usize;
            d.score = score;
            d.nodes = 0;
            d.node_cap = self.node_cap.min(self.heur_cap);
            let ok = d.search(&mut active, &zero, &z2, self.sign_depth, &mut regions);
            total_nodes += d.nodes;
            if ok {
                let desc =
                    format!("DISJ regions={} nodes={} score={}", regions, total_nodes, score);
                return (Leaf::Disj(regions), Some(desc));
            }
        }
        // ... and then with a Lemma-K seed pair forced to the top of the tree.  Two sliding cuts
        // of different kinds (the interior tile poses) leave a quadrant with no witness at all,
        // and the only thing that discharges it is Lemma K on the pair of pivots that cannot
        // both be violated; a weight-greedy branch choice never looks for that pair, so the
        // pairs are enumerated (an O(k^2) float scan, confirmed exactly inside `cover`) and
        // each is tried as a forced prefix.
        let pairs = d.kpairs();
        if d.dbg {
            eprintln!("  DISJ: {} float-plausible Lemma-K pairs", pairs.len());
        }
        let mut seeds: Vec<(I, usize, usize)> = Vec::new();
        for &(bi, bj) in pairs.iter() {
            let act = vec![bi * 2 + 1, bj * 2 + 1];
            let w = {
                let c = d.cover(&act, false);
                d.weight_of(&c)
            };
            seeds.push((w, bi, bj));
        }
        seeds.sort_by(|a, b| b.0.cmp(&a.0));
        if d.dbg {
            for &(w, bi, bj) in seeds.iter().take(8) {
                let act = vec![bi * 2 + 1, bj * 2 + 1];
                let ce = d.cover(&act, true);
                eprintln!(
                    "    seed ({},{}) pts ({},{}) kinds ({},{}) float-cov={:.6} exact-cov={:.6}",
                    bi, bj, d.branches[bi].0, d.branches[bj].0,
                    d.branches[bi].1, d.branches[bj].1,
                    w as f64 / self.cert.w as f64,
                    d.weight_of(&ce) as f64 / self.cert.w as f64
                );
            }
        }
        seeds.truncate(self.ck_seed_cap());
        for &(_, bi, bj) in seeds.iter() {
            for fo in [vec![bi, bj], vec![bj, bi]] {
                let mut active: Vec<usize> = Vec::new();
                let mut regions = 0usize;
                d.score = 0;
                d.nodes = 0;
                d.node_cap = self.node_cap.min(self.heur_cap);
                d.forced = fo;
                let ok = d.search(&mut active, &zero, &z2, self.sign_depth, &mut regions);
                total_nodes += d.nodes;
                if ok {
                    let desc = format!(
                        "DISJ regions={} nodes={} seed=({},{})",
                        regions, total_nodes, bi, bj
                    );
                    return (Leaf::Disj(regions), Some(desc));
                }
            }
        }
        // ... and finally the full two-chain enumeration: interleave the chain orders of the two
        // largest kinds and force that whole sequence.  The recursion then visits exactly the
        // product regions `R_r x R'_s` of `RUNG2.md` sec 6.2, with the inconsistent and
        // witness-free ones discharged by the Lemma-K closure in `cover` instead of by a
        // separate emptiness test.
        let chains = d.chain_order();
        if !chains.is_empty() {
            let mut forced: Vec<usize> = Vec::new();
            if chains.len() == 1 {
                forced = chains[0].clone();
            } else {
                let (a, b) = (&chains[0], &chains[1]);
                let n = a.len().max(b.len());
                for i in 0..n {
                    if i < a.len() {
                        forced.push(a[i]);
                    }
                    if i < b.len() {
                        forced.push(b[i]);
                    }
                }
            }
            let dep = forced.len() as u32 + 4;
            let mut active: Vec<usize> = Vec::new();
            let mut regions = 0usize;
            d.score = 0;
            d.nodes = 0;
            d.node_cap = self.node_cap;
            d.forced = forced;
            let ok = d.search(&mut active, &zero, &z2, dep, &mut regions);
            total_nodes += d.nodes;
            if d.dbg {
                eprintln!("  DISJ: two-chain forced order, nodes={} ok={}", d.nodes, ok);
            }
            if ok {
                let desc = format!("DISJ regions={} nodes={} chains", regions, total_nodes);
                return (Leaf::Disj(regions), Some(desc));
            }
        }
        (Leaf::Split, None)
    }

    fn ck_seed_cap(&self) -> usize {
        self.seed_cap
    }

    fn split(&self, b: &Bx) -> (Bx, Bx) {
        let dcf = b.dc as f64;
        let sx = (b.ax1 - b.ax0) as f64 / dcf;
        let sy = (b.ay1 - b.ay0) as f64 / dcf;
        let bias = if self.wall_touch(b) { self.theta_bias } else { 1.0 };
        // d theta / d u = 2/(1+u^2) in (1,2]; weight the angle extent by 2 * bias
        let su = 2.0 * bias * (b.h as f64 / b.m as f64);
        if su >= sx && su >= sy {
            let mut l = *b;
            let mut r = *b;
            l.m = b.m * 2;
            r.m = b.m * 2;
            l.v0 = b.v0 * 2;
            l.h = b.h;
            r.v0 = b.v0 * 2 + b.h;
            r.h = b.h;
            l.du = b.du + 1;
            r.du = b.du + 1;
            l.depth = b.depth + 1;
            r.depth = b.depth + 1;
            (l, r)
        } else if sx >= sy {
            let mut l = *b;
            let mut r = *b;
            for c in [&mut l, &mut r] {
                c.dc = b.dc * 2;
                c.ay0 = b.ay0 * 2;
                c.ay1 = b.ay1 * 2;
                c.dxy = b.dxy + 1;
                c.depth = b.depth + 1;
            }
            let mid = b.ax0 + b.ax1;
            l.ax0 = b.ax0 * 2;
            l.ax1 = mid;
            r.ax0 = mid;
            r.ax1 = b.ax1 * 2;
            (l, r)
        } else {
            let mut l = *b;
            let mut r = *b;
            for c in [&mut l, &mut r] {
                c.dc = b.dc * 2;
                c.ax0 = b.ax0 * 2;
                c.ax1 = b.ax1 * 2;
                c.dxy = b.dxy + 1;
                c.depth = b.depth + 1;
            }
            let mid = b.ay0 + b.ay1;
            l.ay0 = b.ay0 * 2;
            l.ay1 = mid;
            r.ay0 = mid;
            r.ay1 = b.ay1 * 2;
            (l, r)
        }
    }

    fn rat(&self, n: I, d: I) -> String {
        format!("{n}/{d}")
    }

    fn describe(&self, b: &Bx) -> String {
        format!(
            "x[{},{}] y[{},{}] u[{},{}]",
            self.rat(b.ax0, b.dc),
            self.rat(b.ax1, b.dc),
            self.rat(b.ay0, b.dc),
            self.rat(b.ay1, b.dc),
            self.rat(b.v0, b.m),
            self.rat(b.v0 + b.h, b.m)
        )
    }

    fn run_box(&self, b: Bx, parent: Option<&Work>, cen: &mut Census, dump: bool) {
        cen.boxes += 1;
        cen.maxdepth = cen.maxdepth.max(b.depth);
        if 4 * b.du + b.dxy > SAFE_BITS {
            cen.unsafe_boxes += 1;
            cen.uncert += 1;
            if cen.uncert_list.len() < 80 {
                cen.uncert_list.push(format!("{} [REFUSED: i128 headroom]", self.describe(&b)));
            }
            return;
        }
        let bp = box_poly_full(&b, self.mm, self.wall_touch(&b));
        let wk = self.prepare(&b, &bp, parent);
        let (leaf, desc) = self.classify(&b, &bp, &wk);
        match leaf {
            Leaf::Empty => {
                cen.empty += 1;
                if dump {
                    let _ = writeln!(cen.dump, "EMPTY {}", self.describe(&b));
                }
            }
            Leaf::Adm(_) => {
                cen.adm += 1;
                if dump {
                    let _ = writeln!(
                        cen.dump,
                        "ADM {} {} witnesses={:?}",
                        self.describe(&b),
                        desc.unwrap_or_default(),
                        wk.tset
                    );
                }
            }
            Leaf::Disj(_) => {
                cen.disj += 1;
                if dump {
                    let _ =
                        writeln!(cen.dump, "DISJ {} {}", self.describe(&b), desc.unwrap_or_default());
                }
            }
            Leaf::Unsafe => unreachable!(),
            Leaf::Split => {
                if b.depth >= self.depth_max {
                    cen.uncert += 1;
                    if cen.uncert_list.len() < 80 {
                        cen.uncert_list.push(self.describe(&b));
                    }
                    if dump {
                        let _ = writeln!(cen.dump, "UNCERT {}", self.describe(&b));
                    }
                    return;
                }
                let (l, r) = self.split(&b);
                self.run_box(l, Some(&wk), cen, dump);
                self.run_box(r, Some(&wk), cen, dump);
            }
        }
    }
}

// ------------------------------------------------------------------------------ exact pose mode

/// Exact captured weight at one rational pose, with no boxes and no subdivision: a point is
/// captured iff all four `Ghat_{p,k} <= 0` at that pose.
fn pose_weight(ck: &Checker, xn: I, xd: I, yn: I, yd: I, un: I, ud: I) -> (I, usize, bool) {
    // put the pose in the box representation with a degenerate box: Dc = lcm-ish, M = ud, h = 0
    let d = ck.cert.d;
    let dc = xd * yd * d;
    let ax = xn * (dc / xd);
    let ay = yn * (dc / yd);
    let b = Bx {
        dc,
        ax0: ax,
        ax1: ax,
        ay0: ay,
        ay1: ay,
        m: ud,
        v0: un,
        h: 0,
        dxy: 0,
        du: 0,
        depth: 0,
    };
    let bp = box_poly_full(&b, ck.mm, true);
    let mut w: I = 0;
    let mut cnt = 0;
    let sc = dc / d;
    for i in 0..ck.cert.xs.len() {
        let xs = ck.cert.xs[i] * sc;
        let ys = ck.cert.ys[i] * sc;
        let mut inside = true;
        for k in 0..4 {
            let g = gpoly(k, xs, ys, &b, &bp);
            if g.at(ax, ay)[0] > 0 {
                inside = false;
                break;
            }
        }
        if inside {
            w += ck.cert.wt[i];
            cnt += 1;
        }
    }
    // admissibility
    let mut adm = true;
    for i in 0..4 {
        let a = bp.a[i];
        if a.at(ax, ay)[0] > 0 {
            adm = false;
        }
    }
    (w, cnt, adm)
}

// ------------------------------------------------------------------------------------------ CLI

fn parse_rat(s: &str) -> Result<(I, I), String> {
    if let Some((a, b)) = s.split_once('/') {
        let n: I = a.trim().parse().map_err(|_| format!("bad rational {s}"))?;
        let d: I = b.trim().parse().map_err(|_| format!("bad rational {s}"))?;
        if d <= 0 {
            return Err(format!("bad rational {s}"));
        }
        Ok((n, d))
    } else if s.contains('.') {
        let neg = s.starts_with('-');
        let t = s.trim_start_matches('-');
        let (ip, fp) = t.split_once('.').unwrap();
        let mut d: I = 1;
        for _ in 0..fp.len() {
            d *= 10;
        }
        let ipv: I = if ip.is_empty() { 0 } else { ip.parse().map_err(|_| "bad decimal")? };
        let fpv: I = if fp.is_empty() { 0 } else { fp.parse().map_err(|_| "bad decimal")? };
        let n = ipv * d + fpv;
        Ok((if neg { -n } else { n }, d))
    } else {
        let n: I = s.trim().parse().map_err(|_| format!("bad rational {s}"))?;
        Ok((n, 1))
    }
}

fn die(msg: &str) -> ! {
    eprintln!("ERROR: {msg}");
    std::process::exit(2);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: zmcheck cert FILE [--depth D] [--threads T] [--nodisj] [--dump F] \
             [--xlo A] [--xhi B] [--theta-bias K] [--sign-depth S] [--node-cap N]\n\
             \x20      zmcheck pose FILE --x R --y R --u R\n\
             \x20      zmcheck box  FILE --box x0,x1,y0,y1,u0,u1"
        );
        std::process::exit(2);
    }
    let mode = args[1].clone();
    let path = args[2].clone();
    let mut depth = 18u32;
    let mut threads = 8usize;
    let mut disj = true;
    let mut dumpf: Option<String> = None;
    let mut xlo: Option<f64> = None;
    let mut xhi: Option<f64> = None;
    let mut ylo: Option<f64> = None;
    let mut yhi: Option<f64> = None;
    let mut theta_bias = 4.0f64;
    let mut sign_depth = 14u32;
    let mut node_cap = 400000usize;
    let mut branch_cap = 160usize;
    let mut try_branches = 3usize;
    let mut rank_exact = 64usize;
    let mut fail_cap = 1024usize;
    let mut seed_cap = 0usize;
    let mut heur_cap = 3000usize;
    let mut px = None;
    let mut py = None;
    let mut pu = None;
    let mut boxspec: Option<String> = None;
    let mut i = 3;
    while i < args.len() {
        let a = args[i].as_str();
        let mut need = || -> String {
            i += 1;
            args.get(i).cloned().unwrap_or_else(|| die(&format!("missing value after {a}")))
        };
        match a {
            "--depth" => depth = need().parse().unwrap_or_else(|_| die("bad --depth")),
            "--threads" => threads = need().parse().unwrap_or_else(|_| die("bad --threads")),
            "--nodisj" => disj = false,
            "--dump" => dumpf = Some(need()),
            "--xlo" => xlo = Some(need().parse().unwrap_or_else(|_| die("bad --xlo"))),
            "--xhi" => xhi = Some(need().parse().unwrap_or_else(|_| die("bad --xhi"))),
            "--ylo" => ylo = Some(need().parse().unwrap_or_else(|_| die("bad --ylo"))),
            "--yhi" => yhi = Some(need().parse().unwrap_or_else(|_| die("bad --yhi"))),
            "--theta-bias" => {
                theta_bias = need().parse().unwrap_or_else(|_| die("bad --theta-bias"))
            }
            "--sign-depth" => {
                sign_depth = need().parse().unwrap_or_else(|_| die("bad --sign-depth"))
            }
            "--node-cap" => node_cap = need().parse().unwrap_or_else(|_| die("bad --node-cap")),
            "--heur-cap" => {
                heur_cap = need().parse().unwrap_or_else(|_| die("bad --heur-cap"))
            }
            "--seed-cap" => {
                seed_cap = need().parse().unwrap_or_else(|_| die("bad --seed-cap"))
            }
            "--fail-cap" => {
                fail_cap = need().parse().unwrap_or_else(|_| die("bad --fail-cap"))
            }
            "--rank-exact" => {
                rank_exact = need().parse().unwrap_or_else(|_| die("bad --rank-exact"))
            }
            "--try-branches" => {
                try_branches = need().parse().unwrap_or_else(|_| die("bad --try-branches"))
            }
            "--branch-cap" => {
                branch_cap = need().parse().unwrap_or_else(|_| die("bad --branch-cap"))
            }
            "--x" => px = Some(need()),
            "--y" => py = Some(need()),
            "--u" => pu = Some(need()),
            "--box" => boxspec = Some(need()),
            _ => die(&format!("unknown option {a}")),
        }
        i += 1;
    }

    let cert = match load(&path) {
        Ok(c) => c,
        Err(e) => die(&e),
    };
    if cert.m_num % cert.m_den != 0 {
        die("this checker requires an integer container side (s_num/s_den must be an integer)");
    }
    let mm = cert.m_num / cert.m_den;
    if mm < 1 || mm > 12 {
        die("container side out of the supported range 1..12");
    }
    if cert.xs.is_empty() {
        die("certificate has no points");
    }
    if cert.d % 10 != 0 {
        die("this checker requires D to be a multiple of 10 (the root grid pitch is 1/10)");
    }
    let (grid, gn) = Checker::build_grid(&cert, mm);
    let total = cert.total;
    let w = cert.w;
    let ck = Checker {
        cert,
        mm,
        depth_max: depth,
        disj,
        grid,
        gn,
        theta_bias,
        branch_cap,
        node_cap,
        sign_depth,
        try_branches,
        rank_exact,
        fail_cap,
        seed_cap,
        heur_cap,
        noscreen: std::env::var("ZM_NOSCREEN").is_ok(),
    };

    println!(
        "certificate {path}: container [0,{}]^2, {} points, D={} W={}",
        mm,
        ck.cert.xs.len(),
        ck.cert.d,
        w
    );
    println!(
        "total weight = {}/{} = {:.9}",
        total,
        w,
        total as f64 / w as f64
    );

    match mode.as_str() {
        "pose" => {
            let (xn, xd) = parse_rat(&px.unwrap_or_else(|| die("pose needs --x"))).unwrap();
            let (yn, yd) = parse_rat(&py.unwrap_or_else(|| die("pose needs --y"))).unwrap();
            let (un, ud) = parse_rat(&pu.unwrap_or_else(|| die("pose needs --u"))).unwrap();
            if un < 0 || un > ud {
                die("u must lie in [0,1]");
            }
            let (cw, cnt, adm) = pose_weight(&ck, xn, xd, yn, yd, un, ud);
            println!(
                "pose x={}/{} y={}/{} u={}/{} admissible={}",
                xn, xd, yn, yd, un, ud, adm
            );
            println!(
                "EXACT captured weight = {}/{} = {:.9}  ({} points){}",
                cw,
                w,
                cw as f64 / w as f64,
                cnt,
                if cw < w { "   *** VIOLATION: < 1 ***" } else { "   (OK: >= 1)" }
            );
            std::process::exit(if cw < w { 1 } else { 0 });
        }
        "box" => {
            let spec = boxspec.unwrap_or_else(|| die("box needs --box"));
            let p: Vec<&str> = spec.split(',').collect();
            if p.len() != 6 {
                die("--box wants x0,x1,y0,y1,u0,u1");
            }
            let r: Vec<(I, I)> = p.iter().map(|s| parse_rat(s).unwrap()).collect();
            // put everything over a common centre denominator and a power-of-two angle denominator
            let dc = {
                let mut d = ck.cert.d;
                for k in 0..4 {
                    d = lcm(d, r[k].1);
                }
                d
            };
            let m = lcm(r[4].1, r[5].1);
            let b = Bx {
                dc,
                ax0: r[0].0 * (dc / r[0].1),
                ax1: r[1].0 * (dc / r[1].1),
                ay0: r[2].0 * (dc / r[2].1),
                ay1: r[3].0 * (dc / r[3].1),
                m,
                v0: r[4].0 * (m / r[4].1),
                h: r[5].0 * (m / r[5].1) - r[4].0 * (m / r[4].1),
                dxy: 0,
                du: 0,
                depth: 0,
            };
            let bp = box_poly_full(&b, ck.mm, ck.wall_touch(&b));
            let wk = ck.prepare(&b, &bp, None);
            println!(
                "box {}: reach={} T={} w(T)={}/{} = {:.9}",
                ck.describe(&b),
                wk.cand.len(),
                wk.tset.len(),
                wk.twt,
                w,
                wk.twt as f64 / w as f64
            );
            let (leaf, desc) = ck.classify(&b, &bp, &wk);
            let name = match leaf {
                Leaf::Empty => "EMPTY",
                Leaf::Adm(_) => "ADM",
                Leaf::Disj(_) => "DISJ",
                Leaf::Split => "needs subdivision",
                Leaf::Unsafe => "REFUSED",
            };
            println!("verdict: {name}  {}", desc.unwrap_or_default());
        }
        "cert" => {
            // root boxes: pitch 1/10 in x and y over [0,m], 8 bins of [0,1] in u
            let n = (mm * 10) as usize;
            let mut roots: Vec<Bx> = Vec::new();
            for i in 0..n {
                let x0 = 100 * i as I;
                if let Some(v) = xlo {
                    if (x0 as f64) / 1000.0 + 1e-12 < v {
                        continue;
                    }
                }
                if let Some(v) = xhi {
                    if (x0 as f64) / 1000.0 > v + 1e-12 {
                        continue;
                    }
                }
                for j in 0..n {
                    let y0 = 100 * j as I;
                    if let Some(v) = ylo {
                        if (y0 as f64) / 1000.0 + 1e-12 < v {
                            continue;
                        }
                    }
                    if let Some(v) = yhi {
                        if (y0 as f64) / 1000.0 > v + 1e-12 {
                            continue;
                        }
                    }
                    for k in 0..8 {
                        roots.push(Bx {
                            dc: 1000,
                            ax0: x0,
                            ax1: x0 + 100,
                            ay0: 100 * j as I,
                            ay1: 100 * j as I + 100,
                            m: 8,
                            v0: k,
                            h: 1,
                            dxy: 0,
                            du: 0,
                            depth: 0,
                        });
                    }
                }
            }
            let partial = xlo.is_some() || xhi.is_some() || ylo.is_some() || yhi.is_some();
            if partial {
                println!("PARTIAL SWEEP: centre range restricted; the verdict can never be VERIFIED");
            }
            println!(
                "{} root boxes (pitch 1/10 in x,y; 8 bins of u in [0,1]); depth limit {}; disj {}; \
                 {} threads",
                roots.len(),
                depth,
                disj,
                threads
            );
            let t0 = std::time::Instant::now();
            let ckA = Arc::new(ck);
            let rootsA = Arc::new(roots);
            let idx = Arc::new(AtomicUsize::new(0));
            let acc = Arc::new(Mutex::new(Census::default()));
            let done = Arc::new(AtomicUsize::new(0));
            let gbox = Arc::new(AtomicUsize::new(0));
            let gunc = Arc::new(AtomicUsize::new(0));
            let dumping = dumpf.is_some();
            let mut hs = Vec::new();
            for _ in 0..threads.max(1) {
                let ck = ckA.clone();
                let roots = rootsA.clone();
                let idx = idx.clone();
                let acc = acc.clone();
                let done = done.clone();
                let gbox = gbox.clone();
                let gunc = gunc.clone();
                let nroot = roots.len();
                hs.push(std::thread::spawn(move || {
                    let mut local = Census::default();
                    let (mut pb, mut pu) = (0usize, 0usize);
                    loop {
                        let i = idx.fetch_add(1, Ordering::Relaxed);
                        if i >= nroot {
                            break;
                        }
                        ck.run_box(roots[i], None, &mut local, dumping);
                        gbox.fetch_add(local.boxes - pb, Ordering::Relaxed);
                        gunc.fetch_add(local.uncert - pu, Ordering::Relaxed);
                        pb = local.boxes;
                        pu = local.uncert;
                        let d = done.fetch_add(1, Ordering::Relaxed) + 1;
                        if d % 500 == 0 {
                            println!(
                                "  progress: roots {}/{}  boxes {}  uncertified {}  {:.0}s",
                                d,
                                nroot,
                                gbox.load(Ordering::Relaxed),
                                gunc.load(Ordering::Relaxed),
                                t0.elapsed().as_secs_f64()
                            );
                            use std::io::Write;
                            let _ = std::io::stdout().flush();
                        }
                    }
                    acc.lock().unwrap().merge(local);
                }));
            }
            for h in hs {
                h.join().unwrap();
            }
            let cen = acc.lock().unwrap();
            let secs = t0.elapsed().as_secs_f64();
            println!(
                "done in {:.0}s: boxes {}, max depth {}",
                secs, cen.boxes, cen.maxdepth
            );
            println!(
                "  leaves: ADM {}  DISJ {}  EMPTY {}  UNCERTIFIED {}",
                cen.adm, cen.disj, cen.empty, cen.uncert
            );
            if cen.unsafe_boxes > 0 {
                println!("  {} boxes REFUSED for i128 headroom", cen.unsafe_boxes);
            }
            if !cen.uncert_list.is_empty() {
                println!("  first uncertified boxes:");
                for s in cen.uncert_list.iter().take(40) {
                    println!("    {s}");
                }
            }
            if let Some(f) = dumpf {
                std::fs::write(&f, &cen.dump).unwrap_or_else(|e| die(&format!("dump: {e}")));
                println!("  leaf dump -> {f}");
            }
            if cen.uncert == 0 && !partial {
                println!(
                    "VERIFIED: every closed unit square in [0,{}]^2 captures weight >= 1; \
                     total weight {}/{} = {:.9}",
                    mm,
                    total,
                    w,
                    total as f64 / w as f64
                );
                std::process::exit(0);
            } else {
                println!(
                    "NOT VERIFIED: {} uncertified boxes{}",
                    cen.uncert,
                    if partial { " (partial sweep)" } else { "" }
                );
                std::process::exit(0);
            }
        }
        _ => die(&format!("unknown mode {mode}")),
    }
}

fn lcm(a: I, b: I) -> I {
    a / gcd(a, b) * b
}
fn gcd(a: I, b: I) -> I {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}
