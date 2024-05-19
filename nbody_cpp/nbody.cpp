#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

const double G = 6.67430e-11; // Gravitational constant

const double yrinsec = 3.154e7;
const double dt = 0.01*yrinsec; // Time step

struct Body {
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double m;
};

std::vector<Body> bodies;

std::ofstream outFile;

void saveData() {
    for (const Body& body : bodies) {
        outFile << body.x << " " << body.y << " " << body.z << " "
                << body.vx << " " << body.vy << " " << body.vz << " "
                << body.m << " " << body.ax << " " << body.ay << " " << body.az <<"\n";
    }
}

void calculateAccelerations(int i) {
    for (int j = 0; j < bodies.size(); j++) {
        if (i != j) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dz = bodies[j].z - bodies[i].z;
            double r2 = (dx * dx) + (dy * dy) + (dz * dz);
            double r = std::sqrt(r2);
            double f = - G * bodies[j].m ;
            bodies[i].ax += f * dx / (r * r * r);
            bodies[i].ay += f * dy / (r * r * r);
            bodies[i].az += f * dz / (r * r * r);
        }
    }
}

void leapfrogIntegration() {
    for (int i = 0; i < bodies.size(); i++) {
        calculateAccelerations(i);
        bodies[i].vx += bodies[i].ax * dt * 0.5;
        bodies[i].vy += bodies[i].ay * dt * 0.5;
        bodies[i].vz += bodies[i].az * dt * 0.5;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    
        calculateAccelerations(i);
        bodies[i].vx += bodies[i].ax * dt * 0.5;
        bodies[i].vy += bodies[i].ay * dt * 0.5;
        bodies[i].vz += bodies[i].az * dt * 0.5;
    };

    saveData();
}



int main() {
    std::string outputFile = "output.txt";
    outFile.open(outputFile);

    double Msol = 1.989e30; // Solar mass
    double AU = 1.496e11; // Astronomical unit

    // bodies = {
    //     {0, 0, 0, 0, 1.61126479e+00, 0,0,0,0, 1*Msol}, // Sun
    //     // {0.39*AU, 0, 0, 10237, 5.23688815e+04, 0, 1.6425e-07*Msol}, // Mercury
    //     // {0.72*AU, 0, 0, 23240, -3.37179105e+04, 0, 2.4335e-06*Msol}, // Venus
    //      {1*AU, 0, 0, 11775, -3.41185838e+04, 0,0,0,0, 2.9860e-06*Msol}, // Earth
    //     // {1.52*AU, 0, 0, -24169,1.87423599e+04 , 0, 3.1950e-07*Msol}, // Mars
    //     // {5.2*AU, 0, 0, 14774,1.39683240e+04 , 0, 9.4900e-04*Msol}, // Jupiter
    //     // {9.58*AU, 0, 0, -11509, 2.10079310e+02, 0, 2.8415e-04*Msol}, // Saturn
    //     // {19.22*AU, 0, 0, -1558, 4.75195067e+03, 0, 4.3405e-05*Msol}, // Uranus
    //     // {30.05*AU, 0, 0, -9918, 2.05761038e+04, 0, 5.1200e-05*Msol} // Neptune
    // };

    // set initial conditions
    
    bodies.resize(2);

    bodies[0].x = 0;
    bodies[0].y = 0;
    bodies[0].z = 0;
    bodies[0].vx = 0;
    bodies[0].vy = 1.91126479e-01;
    bodies[0].vz = 0;
    bodies[0].m = 1*Msol;

    bodies[1].x = 1*AU;
    bodies[1].y = 0;
    bodies[1].z = 0;
    bodies[1].vx = 13942;
    bodies[1].vy = -3.34e+04;
    bodies[1].vz = 0;
    bodies[1].m = 2.9860e-06*Msol;





    int numSteps = 1000;
    for (int i = 0; i < numSteps; i++) {
        leapfrogIntegration();
        // ajust for center of mass
        double xcm = 0;
        double ycm = 0;
        double zcm = 0;
        double mtot = 0;
        for (const Body& body : bodies) {
            xcm += body.x * body.m;
            ycm += body.y * body.m;
            zcm += body.z * body.m;
            mtot += body.m;
        }
        xcm /= mtot;
        ycm /= mtot;
        zcm /= mtot;

        for (Body& body : bodies) {
            body.x -= xcm;
            body.y -= ycm;
            body.z -= zcm;
        }
    }

    outFile.close();

    return 0;
}
