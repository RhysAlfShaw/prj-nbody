import numpy as np

class Galaxy:
    """

    Class to create the initial conditions of a galaxy.

    """

    def __init__(self, N, R, z, sigma_R, sigma_z, M, Msol = 2E30,G = 6.67430e-11):
        """
        Parameters
        ----------
        N : int
            Number of particles in the disc.
        R : float
            Radius of the disc.
        z : float
            Height of the disc.
        sigma_R : float
            Standard deviation of the radial distribution.
        sigma_z : float
            Standard deviation of the vertical distribution.
        M : float
            Mass of the galaxy.
        G : float
            Gravitational constant.
        """
        self.N = N
        self.R = R
        self.z = z
        self.sigma_R = sigma_R
        self.sigma_z = sigma_z
        self.M = M*Msol
        self.G = G
        self.Mass_array = self.assign_mass()
        self.Pos_array = self.create_disc()
        self.Vel_array = self



    # create a 3d distribution of paticles that emulate a disc liek galaxy
    def create_disc(self):
        """
        Create a 3D distribution of particles that emulate a disc like galaxy.

        Parameters
        ----------
        N : int
            Number of particles in the disc.
        R : float
            Radius of the disc.
        z : float
            Height of the disc.
        sigma_R : float
            Standard deviation of the radial distribution.
        sigma_z : float
            Standard deviation of the vertical distribution.

        Returns
        -------
        ndarray
            3D distribution of particles.
        """
        # create the radial distribution
        R = np.random.normal(0, self.sigma_R, self.N)
        # create the vertical distribution
        z = np.random.normal(0, self.sigma_z, self.N)
        # create the azimuthal distribution
        phi = np.random.uniform(0, 2*np.pi, self.N)
        # create the 3D distribution
        x = R*np.cos(phi)
        y = R*np.sin(phi)
        self.Pos_array = np.array([x, y, z]).T


    # assign each particle a velocity based on its distance from the center following v = sqrt(GM/r)
    def assign_velocity(self):
        """
        Assign each particle a velocity based on its distance from the center following v = sqrt(GM/r).

        Parameters
        ----------
        disc : ndarray
            3D distribution of particles.
        M : float
            Mass of the galaxy.
        G : float
            Gravitational constant.

        Returns
        -------
        ndarray
            3D distribution of particles with assigned velocities.
        """
        # calculate the distance from the center
        r = np.sqrt(np.sum(self.Pos_array**2, axis=1))
        # calculate the velocity
        v = np.sqrt(self.G*self.Mass_array/r)
        # create the velocity vectors
        theta = np.arctan2(self.Pos_array[:,1], self.Pos_array[:,0])
        v_x = -v*np.sin(theta)
        v_y = v*np.cos(theta)
        v_z = np.zeros_like(v)
        self.Vel_array = np.array([v_x, v_y, v_z]).T



    # Assign masses based on the kroupa IMF
    def assign_mass(self):
        """
        Assign masses based on the Kroupa IMF.

        Parameters
        ----------
        N : int
            Number of particles.

        Returns
        -------
        ndarray
            Array of masses.
        """
        # create the masses
        m = np.zeros(self.N)
        for i in range(self.N):
            r = np.random.uniform(0, 1)
            if r < 0.08:
                m[i] = np.random.uniform(0.1, 0.5)
            elif r < 0.3:
                m[i] = np.random.uniform(0.5, 1)
            else:
                m[i] = np.random.uniform(1, 150)
        return m



    # create a plot of the distribution
    def plot_disc(self):
        """
        Create a plot of the 3D distribution of particles.

        Parameters
        ----------
        disc : ndarray
            3D distribution of particles.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.Pos_array[:,0], self.Pos_array[:,1], self.Pos_array[:,2])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()