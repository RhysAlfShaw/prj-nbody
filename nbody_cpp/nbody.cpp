#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

const double G = 6.67430e-11; // Gravitational constant
const double Msol = 1.989e30; // Solar mass
const double AU = 1.496e11;   // Astronomical unit

struct Body {
    double x, y, z;      // Position
    double vx, vy, vz;   // Velocity
    double mass;         // Mass
};

// Function to calculate the acceleration due to gravity
std::vector<double> acceleration(const std::vector<Body>& bodies, int i) {
    std::vector<double> acc(3, 0.0);
    for (int j = 0; j < bodies.size(); j++) {
        if (i != j) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dz = bodies[j].z - bodies[i].z;
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            double f = G * bodies[i].mass * bodies[j].mass / (r * r * r);
            acc[0] += f * dx / r;
            acc[1] += f * dy / r;
            acc[2] += f * dz / r;
        }
    }
    return acc;
}

int main() {
    std::vector<Body> bodies;
    std::string filename;
    std::cout << "Enter the input file name: ";
    std::cin >> filename;

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << "\n";
        return 1;
    }

    int n;
    file >> n;

    bodies.resize(n);

    double dt, Tend;
    for (int i = 0; i < n; i++) {
        file >> bodies[i].x >> bodies[i].y >> bodies[i].z
             >> bodies[i].vx >> bodies[i].vy >> bodies[i].vz
             >> bodies[i].mass >> dt >> Tend;
    }
    
    // Convert the initial positions and velocities to SI units
    for (int i = 0; i < n; i++) {
        bodies[i].x *= AU;
        bodies[i].y *= AU;
        bodies[i].z *= AU;
        bodies[i].mass *= Msol;
    }


    std::ofstream output_file("output.txt");

    double t = 0.0;
    while (t < Tend) {
        for (int i = 0; i < n; i++) {
            std::vector<double> acc = acceleration(bodies, i);
            bodies[i].vx += acc[0] * dt / 2.0;
            bodies[i].vy += acc[1] * dt / 2.0;
            bodies[i].vz += acc[2] * dt / 2.0;
        }

        for (int i = 0; i < n; i++) {
            bodies[i].x += bodies[i].vx * dt;
            bodies[i].y += bodies[i].vy * dt;
            bodies[i].z += bodies[i].vz * dt;
        }

        for (int i = 0; i < n; i++) {
            std::vector<double> acc = acceleration(bodies, i);
            bodies[i].vx += acc[0] * dt / 2.0;
            bodies[i].vy += acc[1] * dt / 2.0;
            bodies[i].vz += acc[2] * dt / 2.0;
        }

        t += dt;

        for (int i = 0; i < n; i++) {
            output_file << t << " " << bodies[i].x << " " << bodies[i].y << " " << bodies[i].z << " "
                        << bodies[i].vx << " " << bodies[i].vy << " " << bodies[i].vz << "\n";
        }
    }

    file.close();
    output_file.close();
    std::cout << "Simulation completed. Results saved in output.txt.\n";

    return 0;
}