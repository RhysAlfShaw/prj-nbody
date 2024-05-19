#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

namespace std;

const double G = 6.67430e-11; // Gravitational constant
const double Msol = 1.989e30; // Solar mass
const double AU = 1.496e11;   // Astronomical unit

struct BodyT {
    vector<double> x, y, z;      // Position
    vector<double> vx, vy, vz;   // Velocity
    vector<double> mass;         // Mass
    <vector<double> t;         // Acceleration
};

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

vector<double> integrate(int i, double dt, vector<Body>& bodies) {

    acc = acceleration(bodies, i);
    bodies[i].vx += acc[0] * dt / 2;
    bodies[i].vy += acc[1] * dt / 2;
    bodies[i].vz += acc[2] * dt / 2;

    bodies[i].x += bodies[i].vx * dt;
    bodies[i].y += bodies[i].vy * dt;
    bodies[i].z += bodies[i].vz * dt;

    acc = acceleration(bodies, i);
    bodies[i].vx += acc[0] * dt / 2;
    bodies[i].vy += acc[1] * dt / 2;
    bodies[i].vz += acc[2] * dt / 2;

    return bodies;
}

Body center_of_mass(const std::vector<Body>& bodies) {
    double x = 0, y = 0, z = 0, m = 0;
    for (int i = 0; i < bodies.size(); i++) {
        x += bodies[i].mass * bodies[i].x;
        y += bodies[i].mass * bodies[i].y;
        z += bodies[i].mass * bodies[i].z;
        m += bodies[i].mass;
    }
    return {x / m, y / m, z / m, 0, 0, 0, m};
}

struct Body run_simulation(const std::vector<Body>& bodies, double dt, double tmax) {
    double time = 0;
    for (double t = 0; t < tmax; t += dt) {
        for (int i = 0; i < bodies.size(); i++) {
            bodies = integrate(i, dt, bodies);
        }
        Body com = center_of_mass(bodies);

        // save the Body data into the BodyT struct
        BodyT bodyt;
        bodyt.x.push_back(com.x);
        bodyt.y.push_back(com.y);
        bodyt.z.push_back(com.z);
        bodyt.vx.push_back(com.vx);
        bodyt.vy.push_back(com.vy);
        bodyt.vz.push_back(com.vz);
        bodyt.mass.push_back(com.mass);
    }
    return bodies;
}


int main() {
    
    printf("Hello, World!\n");

    return 0;
}