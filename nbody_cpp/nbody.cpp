#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <list>

const double G = 6.67430e-11; // Gravitational constant
const double Msol = 1.989e30; // Solar mass
const double AU = 1.496e11;   // Astronomical unit

struct Body {
    double x, y, z;     // Position
    double vx, vy, vz;  // Velocity
    double mass;        // Mass
};

struct BodyT {
    std::vector<double> x, y, z;     // Position
    std::vector<double> vx, vy, vz;  // Velocity
    std::vector<double> mass;        // Mass
    std::vector<double> ax, ay, az;  // Acceleration
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

std::vector<Body> integrate(const std::vector<Body>& bodies, double dt) {
    std::vector<Body> new_bodies(bodies.size());
    for (int i = 0; i < bodies.size(); i++) {
        std::vector<double> acc = acceleration(bodies, i);
        new_bodies[i].vx = bodies[i].vx + acc[0] * dt;
        new_bodies[i].vy = bodies[i].vy + acc[1] * dt;
        new_bodies[i].vz = bodies[i].vz + acc[2] * dt;
        new_bodies[i].x = bodies[i].x + new_bodies[i].vx * dt;
        new_bodies[i].y = bodies[i].y + new_bodies[i].vy * dt;
        new_bodies[i].z = bodies[i].z + new_bodies[i].vz * dt;
        new_bodies[i].mass = bodies[i].mass;
        }
    return new_bodies;
}

Body center_of_mass(const std::vector<Body>& bodies) {
    double x = 0, y = 0, z = 0, m = 0;
    for (const auto& body : bodies) {
        x += body.mass * body.x;
        y += body.mass * body.y;
        z += body.mass * body.z;
        m += body.mass;
    }
    return {x / m, y / m, z / m, 0, 0, 0, m};
}

BodyT run_simulation(const std::vector<Body>& bodies, double dt, double tmax) {
    BodyT bodyt;
    std::vector<Body> current_bodies = bodies;
    double time = 0;
    while (time < tmax) {
        current_bodies = integrate(current_bodies, dt);
        Body com = center_of_mass(current_bodies);
        bodyt.x.push_back(com.x);
        bodyt.y.push_back(com.y);
        bodyt.z.push_back(com.z);
        bodyt.vx.push_back(com.vx);
        bodyt.vy.push_back(com.vy);
        bodyt.vz.push_back(com.vz);
        bodyt.mass.push_back(com.mass);
        std::vector<double> acc = acceleration(current_bodies, 0); // Assuming the first body is the center of mass
        bodyt.ax.push_back(acc[0]);
        bodyt.ay.push_back(acc[1]);
        bodyt.az.push_back(acc[2]);
        time += dt;
    }
    return bodyt;
}

extern "C" BodyT Simulate(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Z, const std::vector<double>& VX, const std::vector<double>& VY, const std::vector<double>& VZ, const std::vector<double>& MASS, double dt, double tmax) {
    // Set up initial conditions for the bodies
    std::vector<Body> bodies;
    for (int i = 0; i < X.size(); i++) {
        bodies.push_back({X[i], Y[i], Z[i], VX[i], VY[i], VZ[i], MASS[i]});
    }

    return run_simulation(bodies, dt, tmax);
}

int main () {
    // Set up initial conditions for the bodies
    std::vector<Body> bodies = {
        {0, 0, 0, 0, 1.61126479e+00, 0, 1*Msol}, // Sun
        {0.39*AU, 0, 0, 10237, 5.23688815e+04, 0, 1.6425e-07*Msol}, // Mercury
        {0.72*AU, 0, 0, 23240, -3.37179105e+04, 0, 2.4335e-06*Msol}, // Venus
        {1*AU, 0, 0, 11775, -3.41185838e+04, 0, 2.9860e-06*Msol}, // Earth
        {1.52*AU, 0, 0, -24169,1.87423599e+04 , 0, 3.1950e-07*Msol}, // Mars
        {5.2*AU, 0, 0, 14774,1.39683240e+04 , 0, 9.4900e-04*Msol}, // Jupiter
        {9.58*AU, 0, 0, -11509, 2.10079310e+02, 0, 2.8415e-04*Msol}, // Saturn
        {19.22*AU, 0, 0, -1558, 4.75195067e+03, 0, 4.3405e-05*Msol}, // Uranus
        {30.05*AU, 0, 0, -9918, 2.05761038e+04, 0, 5.1200e-05*Msol} // Neptune
    };

    // Run the simulation
    double yrinsec = 3.154e7;
    double dt = 0.01 * yrinsec; // Time step in years
    double tmax = 10 * yrinsec; // Total time in years

    BodyT bodyt = run_simulation(bodies, dt, tmax);

    // save the output struct to a file
    double time = 0;
    std::ofstream file("output.txt");
    for (int i = 0; i < bodyt.x.size()*9; i++) {
        // after 9 loops add dt to time
        if (i % 9 == 0) {
            time += dt;
        }
        file << time << " " << bodyt.x[i] << " " << bodyt.y[i] << " " << bodyt.z[i] << " " << bodyt.vx[i] << " " << bodyt.vy[i] << " " << bodyt.vz[i] << " " << bodyt.mass[i] << " " << bodyt.ax[i] << " " << bodyt.ay[i] << " " << bodyt.az[i] << std::endl;
    }
    file.close();
    
    return 0;
}

// to compile for terminal use: g++ -o nbody nbody.cpp