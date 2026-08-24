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
//!   * angles in [45deg,90deg) are covered by the D4 symmetry of the atom set (checked).
use std::env;
use std::fs;
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::Arc;
use std::thread;

#[derive(Clone)]
struct Cert { s_num: i128, s_den: i128, d: i128, wd: i128, atoms: Vec<(i128,i128,i128)> }

fn read_cert(path: &str) -> Cert {
    let txt = fs::read_to_string(path).unwrap();
    let mut it = txt.split_ascii_whitespace().map(|x| x.parse::<i128>().unwrap());
    let s_num = it.next().unwrap(); let s_den = it.next().unwrap();
    let d = it.next().unwrap(); let wd = it.next().unwrap();
    let m = it.next().unwrap() as usize;
    let mut atoms = Vec::with_capacity(m);
    for _ in 0..m { let x=it.next().unwrap(); let y=it.next().unwrap(); let w=it.next().unwrap(); atoms.push((x,y,w)); }
    Cert{s_num,s_den,d,wd,atoms}
}

/// check the atom multiset is invariant under the dihedral group of the square [0,s]^2
fn check_symmetry(c: &Cert) -> bool {
    assert!(c.s_num * c.d % c.s_den == 0, "s*D must be an integer");
    let sd = c.s_num * c.d / c.s_den;
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
fn min_cover_k(c: &Cert, cn: i128, sn: i128, g: i128, sg_n: i128, sg_d: i128, wm_n: i128, wm_d: i128, topk: usize) -> (i128, Vec<(i128,f64,f64)>) {
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
    let cert = read_cert(&args[1]);
    let nn: i128 = args[2].parse().unwrap();          // claim: fewer than nn squares fit
    let bigN: i128 = args[3].parse().unwrap();        // angle parameter
    let threads: usize = if args.len()>4 { args[4].parse().unwrap() } else { 16 };
    let topk: usize = if args.len()>5 { args[5].parse().unwrap() } else { 0 };
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
    println!("angles: k=0..{} (N={}), covering [0,45deg]", kk, bigN);
    let bad = Arc::new(AtomicI64::new(0));
    let minw = Arc::new(std::sync::Mutex::new((i128::MAX, 0i128)));
    let cert = Arc::new(cert);
    let mut hs = Vec::new();
    let chunk = (kk as usize + threads) / threads;
    for t in 0..threads {
        let cert = cert.clone(); let bad = bad.clone(); let minw = minw.clone(); let out = out.clone();
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
                // w_min = cos(theta_k)+sin(theta_k) = (c0+s0)/g0, ROUNDED DOWN (sound: bigger centre box)
                let wm_n = ((c0 + s0)*SCALE)/g0; let wm_d = SCALE;
                let (v, wit) = min_cover_k(&cert, c0, s0, g0, sg_n, sg_d, wm_n, wm_d, topk);
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
