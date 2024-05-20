mod aos;

use std::arch::x86_64::_mm256_abs_epi8;
use crate::aos::{_SOFTENING, Double3, init_aos};
const THREAD_NUM: usize = 16;
const BODY_NUMBER: usize = 2048;
const SIZE: usize = BODY_NUMBER / THREAD_NUM;

fn sequential_ap_update(n: usize, m: &Vec<f64>, p: &mut Vec<Double3>, v: &mut Vec<Double3>, dt: f64){
    for i in 0 ..n {
        let mut f = Double3{x: 0.0f64, y:0.0f64, z:0.0f64};
        for j in 0..n {
            let s = p[j] - p[i];
            let d = s.dot() + _SOFTENING * _SOFTENING;
            let d_inv = 1.0f64 / d.sqrt();
            let d_inv3 = d_inv * d_inv * d_inv;

            f = f + s * d_inv3 * m[j];
        }
        v[i] = v[i] + f * dt;
    }

    for i in 0..n {
        p[i] = p[i] + v[i] * dt;
    }
}

fn sequential_ap_simulate(n: usize, dt: f64, t_end:f64){
    let (m, mut p, mut v) = init_aos(n);
    let mut t = 0.0f64;

    while t < t_end {
        sequential_ap_update(n, &m, &mut p, &mut v, dt);
        t += dt;
    }
    let ek = aos::ek(&m, &v);
    let ep = aos::ep(&m, &p);
    let etot = ek + ep;
    println!("Etot = {etot}");
}

fn main() {
    let n: usize = 2048;
    let dt = 0.01f64;
    let t_end = 1.0f64;

    use std::time::Instant;
    let now = Instant::now();

    sequential_ap_simulate(n, dt, t_end);

    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
}

