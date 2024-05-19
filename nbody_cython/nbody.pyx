import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h":
    double sqrt(double)

cdef double G = 6.67430E-11

cdef void calculate_acceleration(double[:] X, double[:] Y, double[:] Z, double[:] Mass, int num_bodies, int i, double[:] res) nogil:
    cdef:
        int j
        double r, ax = 0.0, ay = 0.0, az = 0.0
    for j in range(num_bodies):
        if j != i:
            with gil:
                r = sqrt((X[j] - X[i]) ** 2 + (Y[j] - Y[i]) ** 2 + (Z[j] - Z[i]) ** 2)
            ax += -G * Mass[j] * (X[i] - X[j]) / (r ** 3)
            ay += -G * Mass[j] * (Y[i] - Y[j]) / (r ** 3)
            az += -G * Mass[j] * (Z[i] - Z[j]) / (r ** 3)
    res[0] = ax
    res[1] = ay
    res[2] = az

cdef class Simulation:
    cdef:
        double[:] X, Y, Z, Mass, Vx, Vy, Vz
        double dt, total_time, time
        int num_bodies
        list X_history, Y_history, Z_history, Vx_history, Vy_history, Vz_history

    def __init__(self, X, Y, Z, Mass, Vx, Vy, Vz, dt, total_time):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Mass = Mass
        self.Vx = Vx
        self.Vy = Vy
        self.Vz = Vz
        self.dt = dt
        self.total_time = total_time
        self.time = 0
        self.num_bodies = len(X)
        self.X_history = []
        self.Y_history = []
        self.Z_history = []
        self.Vx_history = []
        self.Vy_history = []
        self.Vz_history = []

    cpdef void calculate_acceleration(self, int i, double[:] res):
        calculate_acceleration(self.X, self.Y, self.Z, self.Mass, self.num_bodies, i, res)

    cpdef void integrate(self, int i):
        cdef double[3] acc
        self.calculate_acceleration(i, acc)
        self.Vx[i] += acc[0] * self.dt / 2
        self.Vy[i] += acc[1] * self.dt / 2
        self.Vz[i] += acc[2] * self.dt / 2
        self.X[i] += self.Vx[i] * self.dt
        self.Y[i] += self.Vy[i] * self.dt
        self.Z[i] += self.Vz[i] * self.dt
        self.calculate_acceleration(i, acc)
        self.Vx[i] += acc[0] * self.dt / 2
        self.Vy[i] += acc[1] * self.dt / 2
        self.Vz[i] += acc[2] * self.dt / 2
        self.time += self.dt

    cpdef void run_simulation(self):
        cdef int i, j
        cdef double[:] t = np.arange(0, self.total_time, self.dt)
        for i in range(len(t)):
            self.time = t[i]
            for j in range(self.num_bodies):
                self.integrate(j)
            self.X_history.append(self.X[:].copy())
            self.Y_history.append(self.Y[:].copy())
            self.Z_history.append(self.Z[:].copy())
            self.Vx_history.append(self.Vx[:].copy())
            self.Vy_history.append(self.Vy[:].copy())
            self.Vz_history.append(self.Vz[:].copy())