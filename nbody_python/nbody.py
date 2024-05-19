import numba
import numpy as np

@numba.jit(nopython=True)
def calculate_acceleration(X, Y, Z, Mass, G, num_bodies, i):
    ax = np.float64(0.0)
    ay = np.float64(0.0)
    az = np.float64(0.0)
    for j in range(num_bodies):
        if j != i:
            r = ((X[j] - X[i]) ** 2 + (Y[j] - Y[i]) ** 2 + (Z[j] - Z[i]) ** 2) ** 0.5
            ax += -G * Mass[j] * (X[i] - X[j]) / (r ** 3)
            ay += -G * Mass[j] * (Y[i] - Y[j]) / (r ** 3)
            az += -G * Mass[j] * (Z[i] - Z[j]) / (r ** 3)

    return np.float64(ax), np.float64(ay), np.float64(az)


class Simulation:

    def __init__(self, X, Y, Z, Mass, Vx, Vy, Vz, dt, total_time, numba=True):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Mass = Mass
        self.Vx = Vx
        self.Vy = Vy
        self.Vz = Vz
        self.G = 6.67430E-11
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
        self.numba = numba
        self.dark_matter = False

    def calculate_acceleration(self, i):
        if self.numba == True:
            ax, ay, az = calculate_acceleration(self.X, self.Y, self.Z, self.Mass, self.G, self.num_bodies, i)
        else:
            ax = 0.0
            ay = 0.0
            az = 0.0
            for j in range(self.num_bodies):
                if j != i:
                    r = ((self.X[j] - self.X[i]) ** 2 + (self.Y[j] - self.Y[i]) ** 2 + (self.Z[j] - self.Z[i]) ** 2) ** 0.5
                    ax += -self.G * self.Mass[j] * (self.X[i] - self.X[j]) / (r ** 3)
                    ay += -self.G * self.Mass[j] * (self.Y[i] - self.Y[j]) / (r ** 3)
                    az += -self.G * self.Mass[j] * (self.Z[i] - self.Z[j]) / (r ** 3)
        return ax, ay, az


    def integrate(self, i):
        """
        Performs a leapfrog integration step for a single body.
        """
        ax, ay, az = self.calculate_acceleration(i)

        self.Vx[i] += ax * self.dt / 2
        self.Vy[i] += ay * self.dt / 2
        self.Vz[i] += az * self.dt / 2

        self.X[i] += self.Vx[i] * self.dt
        self.Y[i] += self.Vy[i] * self.dt
        self.Z[i] += self.Vz[i] * self.dt

        ax, ay, az = self.calculate_acceleration(i)
        self.Vx[i] += ax * self.dt / 2
        self.Vy[i] += ay * self.dt / 2
        self.Vz[i] += az * self.dt / 2

        self.time += self.dt


    def run_simulation(self):
        """
        Runs the simulation until the total time is reached.
        """
        import numpy as np
        from tqdm import tqdm
        self.t = []
        self.t = np.arange(0, self.total_time, self.dt)
        for i in tqdm(range(len(self.t)),total=len(self.t), desc='Running Simulation'):
            self.time = self.t[i]
            for i in range(self.num_bodies):
                self.integrate(i)
            # Store the current positions and velocities
            X_com, Y_com, Z_com = self.X[0], self.Y[0], self.Z[0]
            self.X_history.append(self.X[:].copy()-X_com)
            self.Y_history.append(self.Y[:].copy()-Y_com)
            self.Z_history.append(self.Z[:].copy()-Z_com)
            self.Vx_history.append(self.Vx[:].copy())
            self.Vy_history.append(self.Vy[:].copy())
            self.Vz_history.append(self.Vz[:].copy())