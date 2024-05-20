use std::f64::consts::PI;
use rand::prelude::ThreadRng;
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct Double3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl std::ops::Add for Double3 {
    type Output = Double3;
    fn add(self, other: Double3) -> Double3 {
        Double3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z
        }
    }
}

impl std::ops::Add<f64> for Double3 {
    type Output = Double3;
    fn add(self, other: f64) -> Double3 {
        Double3 {
            x: self.x + other,
            y: self.y + other,
            z: self.z + other
        }
    }
}

impl std::ops::Sub for Double3 {
    type Output = Double3;
    fn sub(self, other: Double3) -> Double3 {
        Double3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z
        }
    }
}

impl std::ops::Sub<f64> for Double3 {
    type Output = Double3;
    fn sub(self, other: f64) -> Double3 {
        Double3 {
            x: self.x - other,
            y: self.y - other,
            z: self.z - other
        }
    }
}

impl std::ops::Mul for Double3 {
    type Output = Double3;
    fn mul(self, other: Double3) -> Double3 {
        Double3 {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z
        }
    }
}

impl std::ops::Mul<f64> for Double3 {
    type Output = Double3;
    fn mul(self, other: f64) -> Double3 {
        Double3 {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other
        }
    }
}

impl std::ops::Div for Double3 {
    type Output = Double3;
    fn div(self, other: Double3) -> Double3 {
        Double3 {
            x: self.x / other.x,
            y: self.y / other.y,
            z: self.z / other.z
        }
    }
}

impl std::ops::Div<f64> for Double3 {
    type Output = Double3;
    fn div(self, other: f64) -> Double3 {
        Double3 {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other
        }
    }
}

impl Double3{
    pub fn dot(self) -> f64{
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn distance(self, other: Double3) -> f64{
        (
            (self.x - other.x) * (self.x - other.x) +
                (self.y - other.y) * (self.y - other.y) +
                (self.z - other.z) * (self.z - other.z)
        ).sqrt()
    }
}


pub const _SOFTENING:f64 = 0.025f64;
const _M:f64 = 1.0f64;

pub fn ep(m: &Vec<f64>, p: &Vec<Double3>) -> f64 {
    let mut epot = 0.0f64;
    for i in 0..m.len() {
        for j in i+1..m.len() {
            let d = p[i].distance(p[j]);
            epot += -1.0f64 * m[i] * m[j] / d;
        }
    }
    return epot;
}

pub fn ek(m: &Vec<f64>, v: &Vec<Double3>) -> f64 {
    let mut ekin = 0.0f64;
    for i in 0..m.len() {
        ekin += 0.5f64 * m[i] * v[i].dot();
    }
    return ekin;
}

pub fn scale_d3_array(m: &mut Vec<Double3>, scale: f64) {
    for i in 0..m.len(){
        m[i] = m[i] * scale;
    }
}

pub fn init_mass(n: usize) -> Vec<f64> {
    let m_i = 1.0f64 / (n as f64);
    let m: Vec<f64> = vec![m_i; n];
    return m;
}

pub fn init_pos(n: usize, rng: &mut ThreadRng) -> Vec<Double3> {
    let mut p: Vec<Double3> = vec![Double3 {x: 0.0f64, y: 0.0f64, z:0.0f64}; n];
    for i in 0..n {
        let r = rng.gen::<f64>();
        let x = (1.0f64 - 2.0f64 * rng.gen::<f64>()).acos();
        let y = rng.gen::<f64>() * 2.0f64 * PI;

        let pi_x = r * x.sin() * y.cos();
        let pi_y = r * x.sin() * y.sin();
        let pi_z = r * x.cos();
        p[i] = Double3 {x:pi_x, y: pi_y, z: pi_z};
    }
    return p;
}

pub fn init_vel(n: usize, rng: &mut ThreadRng) -> Vec<Double3> {
    let mut v: Vec<Double3> = vec![Double3 {x: 0.0f64, y: 0.0f64, z:0.0f64}; n];
    for i in 0..n {
        let vi_x = 1.0f64 - 2.0f64 * rng.gen::<f64>();
        let vi_y = 1.0f64 - 2.0f64 * rng.gen::<f64>();
        let vi_z = 1.0f64 - 2.0f64 * rng.gen::<f64>();
        v[i] =Double3 {x:vi_x, y: vi_y, z: vi_z};
    }
    return v;
}

pub fn move_to_center(m: &Vec<f64>, p: &mut Vec<Double3>, v: &mut Vec<Double3>) {

    let mut pt = Double3{x: 0.0f64, y:0.0f64, z:0.0f64};
    let mut vt = Double3{x: 0.0f64, y:0.0f64, z:0.0f64};

    for i in 0..m.len() {
        pt = pt + p[i] * m[i];
        vt = vt + v[i] * m[i];
    }

    pt = pt / _M;
    vt = vt / _M;

    for i in 0..m.len() {
        p[i] = p[i] - pt;
        v[i] = v[i] - vt;
    }
}

pub fn rescale_energy(m: &Vec<f64>, p: &mut Vec<Double3>, v: &mut Vec<Double3>) {
    let epot = ep(m, p);
    let ekin = ek(m, v);

    let virial_ratio = 0.5f64;
    let qv = (virial_ratio * epot.abs() / ekin).sqrt();

    scale_d3_array(v, qv);

    let beta = ((1.0f64 - virial_ratio) * epot / (epot + ekin)).abs();

    scale_d3_array(p, beta);
    scale_d3_array(v, 1.0f64 / beta.sqrt());

    let epot = ep(m, p);
    let beta = epot / -0.5f64;
    scale_d3_array(p, beta);
    scale_d3_array(v, 1.0f64 / beta.sqrt());
}

pub fn init_aos(n: usize) -> (Vec<f64>, Vec<Double3>, Vec<Double3>) {
    let mut rng = rand::thread_rng();

    let m = init_mass(n);

    let mut p = init_pos(n, &mut rng);

    let mut v = init_vel(n, &mut rng);

    move_to_center(&m, &mut p, &mut v);

    rescale_energy(&m, &mut p, &mut v);

    return (m, p, v);
}

