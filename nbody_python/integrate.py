import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def gravitational_acceleration(r, m, G=6.67430e-11):
    """
    Calculates the gravitational acceleration on a particle due to a set of masses.

    Parameters:
    r (numpy.ndarray): Array of particle positions (n, 3).
    m (numpy.ndarray): Array of particle masses (n,).
    G (float): Gravitational constant.

    Returns:
    numpy.ndarray: Array of accelerations (n, 3).
    """
    n = len(r)
    a = np.zeros_like(r)
    for i in range(n):
        for j in range(n):
            if i != j:
                dr = r[j] - r[i]
                dist = np.linalg.norm(dr)
                a[i] += -G * m[j] * dr / (dist ** 3)
    return a

def leapfrog_galaxy(r0, v0, m, G, t0, t_end, dt):
    """
    Simulates the evolution of a galaxy using the Leapfrog integration method.

    Parameters:
    r0 (numpy.ndarray): Initial positions of particles (n, 3).
    v0 (numpy.ndarray): Initial velocities of particles (n, 3).
    m (numpy.ndarray): Masses of particles (n,).
    G (float): Gravitational constant.
    t0 (float): Initial time.
    t_end (float): Final time.
    dt (float): Time step.

    Returns:
    numpy.ndarray: Array of particle positions (n, 3, num_steps).
    numpy.ndarray: Array of particle velocities (n, 3, num_steps).
    numpy.ndarray: Array of times (num_steps,).

    """
    num_steps = int((t_end - t0) / dt)
    n = len(r0)
    r = np.zeros((n, 3, num_steps + 1))
    v = np.zeros((n, 3, num_steps + 1))
    t = np.linspace(t0, t_end, num_steps + 1)

    # Initial conditions
    r[:, :, 0] = r0
    v[:, :, 0] = v0

    # Perform Leapfrog integration
    for i in tqdm(range(1, num_steps + 1), desc="Simulating galaxy evolution"):
        a = gravitational_acceleration(r[:, :, i-1], m, G)
        v[:, :, i] = v[:, :, i-1] + dt * a
        r[:, :, i] = r[:, :, i-1] + dt * v[:, :, i]

    return r, v, t

# Example usage
# Initialize galaxy parameters
n_particles = 100
r0 = np.random.uniform(-1e12, 1e12, (n_particles, 3))
v0 = np.random.uniform(-1e5, 1e5, (n_particles, 3))
m = np.random.uniform(1e30, 1e32, n_particles)
G = 6.67430e-11
t0 = 0.0
t_end = 1e9  # 1 billion years
dt = 1e5  # 100,000 years

# Simulate galaxy evolution
r, v, t = leapfrog_galaxy(r0, v0, m, G, t0, t_end, dt)

# Plot the galaxy evolution
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(r[0, 0, 0], r[0, 1, 0], r[0, 2, 0], c='r', s=10)
for i in range(1, len(t)):
    ax.scatter(r[:, 0, i], r[:, 1, i], r[:, 2, i], c='b', s=1)
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('Galaxy Evolution')
plt.show()