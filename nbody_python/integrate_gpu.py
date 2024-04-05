import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from tqdm import tqdm

def calculate_acceleration_chunk(args):
    """
    Calculates the gravitational acceleration on a chunk of particles.
    """
    r, m, start, end, G = args
    a = np.zeros((end - start, 3))
    for i in range(start, end):
        for j in range(len(r)):
            if i != j:
                dr = r[j] - r[i]
                dist = np.linalg.norm(dr)
                a[i - start] += -G * m[j] * dr / (dist ** 3)
    return a

def gravitational_acceleration(r, m, G=6.67430e-11):
    """
    Calculates the gravitational acceleration on a particle due to a set of masses.
    """
    n = len(r)
    a = np.zeros_like(r)

    # Use multiprocessing to parallelize the acceleration calculation
    num_processes = mp.cpu_count()
    chunk_size = n // num_processes

    ctx = mp.get_context("spawn")
    with ctx.Pool() as pool:
        args = [(r, m, i * chunk_size, (i + 1) * chunk_size, G) for i in range(num_processes)]
        results = pool.map(calculate_acceleration_chunk, args)

    for i, res in enumerate(results):
        a[i * chunk_size:(i + 1) * chunk_size] = res

    return a

def leapfrog_galaxy(r0, v0, m, G, t0, t_end, dt):
    """
    Simulates the evolution of a galaxy using the Leapfrog integration method.
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
n_particles = 1000
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