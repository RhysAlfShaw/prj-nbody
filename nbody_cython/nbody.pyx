# distutils: language=c++
# cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cdef class Simulation:
    cdef public np.ndarray[DTYPE_t, ndim=1] X, Y, Z, Mass, Vx, Vy, Vz
    cdef public DTYPE_t[:] Vx_buf, Vy_buf, Vz_buf

    cdef public DTYPE_t G, dt, total_time, time, x_com, y_com, z_com
    cdef public int num_bodies
    cdef public list X_history, Y_history, Z_history, Vx_history, Vy_history, Vz_history
    cdef public bint dark_matter

    def __init__(self, X, Y, Z, Mass, Vx, Vy, Vz, dt, total_time):
        self.X = np.asarray(X, dtype=DTYPE)
        self.Y = np.asarray(Y, dtype=DTYPE)
        self.Z = np.asarray(Z, dtype=DTYPE)
        self.Mass = np.asarray(Mass, dtype=DTYPE)
        self.Vx = np.asarray(Vx, dtype=DTYPE)
        self.Vy = np.asarray(Vy, dtype=DTYPE)
        self.Vz = np.asarray(Vz, dtype=DTYPE)
        self.Vx_buf = self.Vx
        self.Vy_buf = self.Vy
        self.Vz_buf = self.Vz
        self.G = 6.67430E-11
        self.dt = dt
        self.total_time = total_time
        self.time = 0
        self.num_bodies = len(self.X)
        self.X_history = []
        self.Y_history = []
        self.Z_history = []
        self.Vx_history = []
        self.Vy_history = []
        self.Vz_history = []
        self.dark_matter = False

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef calculate_acceleration(self, int i):
        cdef DTYPE_t ax, ay, az, r
        cdef int j
        ax = 0.0
        ay = 0.0
        az = 0.0
        for j in range(self.num_bodies):
            if j != i:
                r = ((self.X[j] - self.X[i]) ** 2 + (self.Y[j] - self.Y[i]) ** 2 + (self.Z[j] - self.Z[i]) ** 2) ** 0.5
                ax += -self.G * self.Mass[j] * (self.X[i] - self.X[j]) / (r ** 3)
                ay += -self.G * self.Mass[j] * (self.Y[i] - self.Y[j]) / (r ** 3)
                az += -self.G * self.Mass[j] * (self.Z[i] - self.Z[j]) / (r ** 3)
        return np.array([ax, ay, az], dtype=DTYPE)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef integrate(self, int i):
        """
        Performs a leapfrog integration step for a single body.
        """
        cdef np.ndarray[DTYPE_t, ndim=1] acc = self.calculate_acceleration(i)
        # Update the velocities
        self.Vx_buf[i] += acc[0] * self.dt / 2
        self.Vy_buf[i] += acc[1] * self.dt / 2
        self.Vz_buf[i] += acc[2] * self.dt / 2

        # Update the positions
        self.X[i] += self.Vx_buf[i] * self.dt
        self.Y[i] += self.Vy_buf[i] * self.dt
        self.Z[i] += self.Vz_buf[i] * self.dt

        # Update the velocities again
        acc = self.calculate_acceleration(i)
        self.Vx_buf[i] += acc[0] * self.dt / 2
        self.Vy_buf[i] += acc[1] * self.dt / 2
        self.Vz_buf[i] += acc[2] * self.dt / 2

        self.time += self.dt

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def run_simulation(self):
        """
        Runs the simulation until the total time is reached.
        """
        cdef int i
        self.t = np.arange(0, self.total_time, self.dt)
        for i in range(len(self.t)):
            self.time = self.t[i]
            for i in range(self.num_bodies):
                self.integrate(i)
            # Compute the center of mass
            self.x_com = np.sum(self.X * self.Mass) / np.sum(self.Mass)
            self.y_com = np.sum(self.Y * self.Mass) / np.sum(self.Mass)
            self.z_com = np.sum(self.Z * self.Mass) / np.sum(self.Mass)
            # Store the current positions and velocities
            self.X_history.append(self.X[:].copy() - self.x_com)
            self.Y_history.append(self.Y[:].copy() - self.y_com)
            self.Z_history.append(self.Z[:].copy() - self.z_com)
            self.Vx_history.append(self.Vx[:].copy())
            self.Vy_history.append(self.Vy[:].copy())
            self.Vz_history.append(self.Vz[:].copy())