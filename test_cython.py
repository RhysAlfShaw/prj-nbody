from inital_conditions import Galaxy
from matplotlib import pyplot as plt
import numpy as np
from nbody import Simulation

pc2m = 3.086e16

galaxy = Galaxy(
    N=10000,
    R=1000*pc2m,
    z=0.1*pc2m,
    sigma_R=0.1*pc2m,
    sigma_z=0.1*pc2m,
    bluge_height=0.1*pc2m,
    bludge_frac=0.2
)

X, Y, Z = galaxy.X, galaxy.Y, galaxy.Z
yrinsec = 3.154e7
dt = 1e2 * yrinsec
t_end = 1e6 * yrinsec
Vx = np.zeros_like(X).tolist()
Vy = np.zeros_like(Y).tolist()
Vz = np.zeros_like(Z).tolist()

sim = Simulation(X, Y, Z, galaxy.Mass, Vx,Vy,Vz, dt, t_end)
sim.run_simulation()