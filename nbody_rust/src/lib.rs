use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use ndarray::prelude::*;
use ndarray::Array;
use std::f64;



struct Simulation {
    x: Array1<f64>,
    y: Array1<f64>,
    z: Array1<f64>,
    mass: Array1<f64>,
    vx: Array1<f64>,
    vy: Array1<f64>,
    vz: Array1<f64>,
    g: f64,
    dt: f64,
    total_time: f64,
    time: f64,
    num_bodies: usize,
    x_history: Vec<Vec<f64>>,
    y_history: Vec<Vec<f64>>,
    z_history: Vec<Vec<f64>>,
    vx_history: Vec<Array1<f64>>,
    vy_history: Vec<Array1<f64>>,
    vz_history: Vec<Array1<f64>>,
    numba: bool,
    dark_matter: bool,
}

impl Simulation {
    fn new(x: Array1<f64>, y: Array1<f64>, z: Array1<f64>, mass: Array1<f64>, vx: Array1<f64>, vy: Array1<f64>, vz: Array1<f64>, dt: f64, total_time: f64, numba: bool) -> Simulation {
        Simulation {
            x: x,
            y: y,
            z: z,
            mass: mass,
            vx: vx,
            vy: vy,
            vz: vz,
            g: 6.67430E-11,
            dt: dt,
            total_time: total_time,
            time: 0.0,
            num_bodies: x.len(),
            x_history: Vec::new(),
            y_history: Vec::new(),
            z_history: Vec::new(),
            vx_history: Vec::new(),
            vy_history: Vec::new(),
            vz_history: Vec::new(),
            numba: numba,
            dark_matter: false,
        }
    }

    fn calculate_acceleration(&self, i: usize) -> (f64, f64, f64) {
        let mut ax = 0.0;
        let mut ay = 0.0;
        let mut az = 0.0;
        for j in 0..self.num_bodies {
            if j != i {
                let r = ((self.x[j] - self.x[i]).powi(2) + (self.y[j] - self.y[i]).powi(2) + (self.z[j] - self.z[i]).powi(2)).sqrt();
                ax += -self.g * self.mass[j] * (self.x[i] - self.x[j]) / r.powi(3);
                ay += -self.g * self.mass[j] * (self.y[i] - self.y[j]) / r.powi(3);
                az += -self.g * self.mass[j] * (self.z[i] - self.z[j]) / r.powi(3);
            }
        }
        (ax, ay, az)
    }

    fn integrate(&mut self, i: usize) {
        let (ax, ay, az) = self.calculate_acceleration(i);
        self.vx[i] += ax * self.dt / 2.0;
        self.vy[i] += ay * self.dt / 2.0;
        self.vz[i] += az * self.dt / 2.0;

        self.x[i] += self.vx[i] * self.dt;
        self.y[i] += self.vy[i] * self.dt;
        self.z[i] += self.vz[i] * self.dt;

        let (ax, ay, az) = self.calculate_acceleration(i);
        self.vx[i] += ax * self.dt / 2.0;
        self.vy[i] += ay * self.dt / 2.0;
        self.vz[i] += az * self.dt / 2.0;

        self.time += self.dt;
    }

    fn run_simulation(&mut self) {
        let t = (0..).map(|i| i as f64 * self.dt).take_while(|&x| x < self.total_time).collect::<Vec<_>>();
        for _ in 0..t.len() {
            for i in 0..self.num_bodies {
                self.integrate(i);
            }
            let x_com = self.x[0];
            let y_com = self.y[0];
            let z_com = self.z[0];
            self.x_history.push(self.x.to_vec());
            self.y_history.push(self.y.to_vec());
            self.z_history.push(self.z.to_vec());
            self.vx_history.push(self.vx.clone());
            self.vy_history.push(self.vy.clone());
            self.vz_history.push(self.vz.clone());
        }
    }
}

#[pymodule]
fn my_simulation(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfunction]
    fn py_simulation(
        _py: Python,
        x: Vec<f64>,
        y: Vec<f64>,
        z: Vec<f64>,
        mass: Vec<f64>,
        vx: Vec<f64>,
        vy: Vec<f64>,
        vz: Vec<f64>,
        dt: f64,
        total_time: f64,
        numba: bool,
    ) -> PyResult<Simulation> {
        Simulation::new(
            Array::from(x),
            Array::from(y),
            Array::from(z),
            Array::from(mass),
            Array::from(vx),
            Array::from(vy),
            Array::from(vz),
            dt,
            total_time,
            numba,
        )
        .map_err(PyErr)
    }

    Ok(())
}