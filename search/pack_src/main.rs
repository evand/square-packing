// Unit-square packing: heuristic search for packings of n unit squares in an s x s square.
// Variables per square: (x, y, theta). Non-overlap via separating-axis penetration depth.
use std::env;
use std::sync::{Arc, Mutex};
use std::thread;

const PI2: f64 = std::f64::consts::FRAC_PI_2;

#[inline]
fn hw(a: f64) -> f64 { // half-width of unit square along direction at angle a relative to square axes
    (a.cos().abs() + a.sin().abs()) * 0.5
}
#[inline]
fn dhw(a: f64) -> f64 { // derivative of hw
    (-a.cos().signum() * a.sin() + a.sin().signum() * a.cos()) * 0.5
}

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self { Rng(seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407) | 1) }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        self.0 = x; x
    }
    fn f(&mut self) -> f64 { (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64 }
    fn gauss(&mut self) -> f64 {
        let u1 = self.f().max(1e-12); let u2 = self.f();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// energy and gradient; v layout [x0,y0,t0, x1,y1,t1, ...]
fn energy(v: &[f64], n: usize, s: f64, g: &mut [f64]) -> f64 {
    for x in g.iter_mut() { *x = 0.0; }
    let mut e = 0.0;
    // walls
    for i in 0..n {
        let (x, y, t) = (v[3*i], v[3*i+1], v[3*i+2]);
        let w = t.cos().abs() + t.sin().abs();
        let dw = -t.cos().signum() * t.sin() + t.sin().signum() * t.cos();
        // left
        let vl = w * 0.5 - x;
        if vl > 0.0 { e += vl*vl; g[3*i] += -2.0*vl; g[3*i+2] += 2.0*vl*0.5*dw; }
        let vr = x + w*0.5 - s;
        if vr > 0.0 { e += vr*vr; g[3*i] += 2.0*vr; g[3*i+2] += 2.0*vr*0.5*dw; }
        let vb = w*0.5 - y;
        if vb > 0.0 { e += vb*vb; g[3*i+1] += -2.0*vb; g[3*i+2] += 2.0*vb*0.5*dw; }
        let vt = y + w*0.5 - s;
        if vt > 0.0 { e += vt*vt; g[3*i+1] += 2.0*vt; g[3*i+2] += 2.0*vt*0.5*dw; }
    }
    // pairs
    for i in 0..n {
        for j in (i+1)..n {
            let (xi, yi, ti) = (v[3*i], v[3*i+1], v[3*i+2]);
            let (xj, yj, tj) = (v[3*j], v[3*j+1], v[3*j+2]);
            let dx = xj - xi; let dy = yj - yi;
            let phis = [ti, ti + PI2, tj, tj + PI2];
            let mut best = f64::INFINITY; let mut ba = 0usize;
            for (a, &ph) in phis.iter().enumerate() {
                let (c, sn) = (ph.cos(), ph.sin());
                let o = hw(ph - ti) + hw(ph - tj) - (dx*c + dy*sn).abs();
                if o < best { best = o; ba = a; }
            }
            if best > 0.0 {
                e += best*best;
                let ph = phis[ba];
                let (c, sn) = (ph.cos(), ph.sin());
                let dot = dx*c + dy*sn;
                let sd = if dot >= 0.0 { 1.0 } else { -1.0 };
                let k = 2.0*best;
                // position derivatives: o = hi+hj-|dot|
                g[3*i]   += k * ( sd*c);
                g[3*i+1] += k * ( sd*sn);
                g[3*j]   += k * (-sd*c);
                g[3*j+1] += k * (-sd*sn);
                // angle derivatives
                let ddot_dphi = -dx*sn + dy*c;
                if ba < 2 { // axis belongs to i: phi = ti + k*pi/2
                    g[3*i+2] += k * ( dhw(ph - tj) - sd*ddot_dphi );
                    g[3*j+2] += k * (-dhw(ph - tj));
                } else {    // axis belongs to j
                    g[3*j+2] += k * ( dhw(ph - ti) - sd*ddot_dphi );
                    g[3*i+2] += k * (-dhw(ph - ti));
                }
            }
        }
    }
    e
}

// L-BFGS with backtracking line search
fn minimize(v: &mut Vec<f64>, n: usize, s: f64, maxit: usize) -> f64 {
    let d = 3*n;
    let m = 10usize;
    let mut g = vec![0.0; d];
    let mut e = energy(v, n, s, &mut g);
    let mut sl: Vec<Vec<f64>> = Vec::new();
    let mut yl: Vec<Vec<f64>> = Vec::new();
    let mut rho: Vec<f64> = Vec::new();
    let mut q = vec![0.0; d];
    for _it in 0..maxit {
        let gn: f64 = g.iter().map(|a| a*a).sum::<f64>().sqrt();
        if e < 1e-24 || gn < 1e-14 { break; }
        // two-loop recursion
        q.copy_from_slice(&g);
        let k = sl.len();
        let mut alpha = vec![0.0; k];
        for idx in (0..k).rev() {
            let a = rho[idx] * sl[idx].iter().zip(q.iter()).map(|(x,y)| x*y).sum::<f64>();
            alpha[idx] = a;
            for t in 0..d { q[t] -= a * yl[idx][t]; }
        }
        let scale = if k > 0 {
            let sy: f64 = sl[k-1].iter().zip(yl[k-1].iter()).map(|(x,y)| x*y).sum();
            let yy: f64 = yl[k-1].iter().map(|x| x*x).sum();
            if yy > 0.0 { sy/yy } else { 1.0 }
        } else { 1.0/(1.0+gn) };
        for t in 0..d { q[t] *= scale; }
        for idx in 0..k {
            let b = rho[idx] * yl[idx].iter().zip(q.iter()).map(|(x,y)| x*y).sum::<f64>();
            for t in 0..d { q[t] += (alpha[idx]-b) * sl[idx][t]; }
        }
        // direction = -q
        let mut dir: Vec<f64> = q.iter().map(|x| -x).collect();
        let mut dg: f64 = dir.iter().zip(g.iter()).map(|(x,y)| x*y).sum();
        if dg > 0.0 { for t in 0..d { dir[t] = -g[t]; } dg = -gn*gn; }
        let mut step = 1.0;
        let mut vn = v.clone();
        let mut gn2 = vec![0.0; d];
        let mut e2 = e;
        let mut ok = false;
        for _ls in 0..40 {
            for t in 0..d { vn[t] = v[t] + step*dir[t]; }
            e2 = energy(&vn, n, s, &mut gn2);
            if e2 <= e + 1e-4*step*dg { ok = true; break; }
            step *= 0.5;
        }
        if !ok { break; }
        let mut sv = vec![0.0; d]; let mut yv = vec![0.0; d];
        for t in 0..d { sv[t] = vn[t]-v[t]; yv[t] = gn2[t]-g[t]; }
        let sy: f64 = sv.iter().zip(yv.iter()).map(|(x,y)| x*y).sum();
        if sy > 1e-20 {
            sl.push(sv); yl.push(yv); rho.push(1.0/sy);
            if sl.len() > m { sl.remove(0); yl.remove(0); rho.remove(0); }
        }
        v.copy_from_slice(&vn); g.copy_from_slice(&gn2); e = e2;
    }
    e
}

fn rand_config(rng: &mut Rng, n: usize, s: f64) -> Vec<f64> {
    let mut v = vec![0.0; 3*n];
    for i in 0..n {
        v[3*i] = 0.5 + rng.f()*(s-1.0);
        v[3*i+1] = 0.5 + rng.f()*(s-1.0);
        v[3*i+2] = if rng.f() < 0.35 { 0.0 } else { rng.f()*PI2 };
    }
    v
}

// try to achieve a feasible packing at side s; returns (best energy, config)
fn attempt(n: usize, s: f64, seed: u64, iters: usize) -> (f64, Vec<f64>) {
    let mut rng = Rng::new(seed);
    let mut best = rand_config(&mut rng, n, s);
    let mut be = minimize(&mut best, n, s, 600);
    let mut temp = 0.30;
    for it in 0..iters {
        let mut v = best.clone();
        // perturbation: shake a random subset
        let mode = rng.next_u64() % 4;
        match mode {
            0 => { for i in 0..n { v[3*i] += temp*rng.gauss(); v[3*i+1] += temp*rng.gauss(); v[3*i+2] += temp*rng.gauss(); } }
            1 => { let k = 1 + (rng.next_u64() as usize) % 3;
                   for _ in 0..k { let i = (rng.next_u64() as usize) % n;
                       v[3*i] = 0.5 + rng.f()*(s-1.0); v[3*i+1] = 0.5 + rng.f()*(s-1.0); v[3*i+2] = rng.f()*PI2; } }
            2 => { let i = (rng.next_u64() as usize) % n; v[3*i+2] += (rng.f()-0.5)*PI2; }
            _ => { for i in 0..n { v[3*i] += 0.5*temp*rng.gauss(); v[3*i+1] += 0.5*temp*rng.gauss(); } }
        }
        let e = minimize(&mut v, n, s, 500);
        if e < be { be = e; best = v; }
        if be < 1e-22 { break; }
        if it % 200 == 199 { temp = if temp > 0.05 { temp*0.7 } else { 0.30 }; }
    }
    (be, best)
}

fn main() {
    let a: Vec<String> = env::args().collect();
    let n: usize = a[1].parse().unwrap();
    let s: f64 = a[2].parse().unwrap();
    let trials: usize = a[3].parse().unwrap();
    let iters: usize = a[4].parse().unwrap();
    let threads: usize = if a.len() > 5 { a[5].parse().unwrap() } else { 8 };
    let seed0: u64 = if a.len() > 6 { a[6].parse().unwrap() } else { 12345 };
    let res = Arc::new(Mutex::new((f64::INFINITY, Vec::<f64>::new())));
    let mut hs = Vec::new();
    for t in 0..threads {
        let res = res.clone();
        hs.push(thread::spawn(move || {
            for k in 0..(trials/threads).max(1) {
                let seed = seed0 ^ ((t as u64) << 32) ^ (k as u64).wrapping_mul(0x9E3779B97F4A7C15);
                let (e, v) = attempt(n, s, seed, iters);
                let mut r = res.lock().unwrap();
                if e < r.0 { r.0 = e; r.1 = v; eprintln!("[t{} k{}] E={:.3e}", t, k, e); }
                if r.0 < 1e-22 { break; }
            }
        }));
    }
    for h in hs { h.join().unwrap(); }
    let r = res.lock().unwrap();
    println!("BEST_E {:.6e}", r.0);
    print!("CONFIG");
    for x in r.1.iter() { print!(" {:.15}", x); }
    println!();
}
