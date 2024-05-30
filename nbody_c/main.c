#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>

#define G 6.67430e-11       // Gravitational constant
#define NUM_BODIES 9        // Number of celestial bodies (Sun + 8 planets)
#define DT 86400.0          // Time step in seconds (1 day)
#define SCALE 1e-10         // Scale factor for drawing
#define WINDOW_WIDTH 800    // Window width
#define WINDOW_HEIGHT 600   // Window height

typedef struct {
    double x, y, z;     // Position
    double vx, vy, vz;  // Velocity
    double mass;        // Mass
    Uint32 color;       // Color for rendering
} Body;

SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;


// Function to calculate the acceleration due to gravity
void calculate_acceleration(Body bodies[], int i, double ax[], double ay[], double az[]) {
    double dx, dy, dz, r, r3;
    ax[i] = ay[i] = az[i] = 0.0;

    for (int j = 0; j < NUM_BODIES; j++) {
        if (i != j) {
            dx = bodies[j].x - bodies[i].x;
            dy = bodies[j].y - bodies[i].y;
            dz = bodies[j].z - bodies[i].z;
            r = sqrt(dx * dx + dy * dy + dz * dz);
            r3 = r * r * r;
            ax[i] += G * bodies[j].mass * dx / r3;
            ay[i] += G * bodies[j].mass * dy / r3;
            az[i] += G * bodies[j].mass * dz / r3;
        }
    }
}

// Leapfrog integrator
void leapfrog_integrate(Body bodies[], double dt) {
    double ax[NUM_BODIES], ay[NUM_BODIES], az[NUM_BODIES];
    double vx_half[NUM_BODIES], vy_half[NUM_BODIES], vz_half[NUM_BODIES];

    // Calculate acceleration at current time
    for (int i = 0; i < NUM_BODIES; i++) {
        calculate_acceleration(bodies, i, ax, ay, az);
    }

    // Update velocity for half-step
    for (int i = 0; i < NUM_BODIES; i++) {
        vx_half[i] = bodies[i].vx + 0.5 * dt * ax[i];
        vy_half[i] = bodies[i].vy + 0.5 * dt * ay[i];
        vz_half[i] = bodies[i].vz + 0.5 * dt * az[i];
    }

    // Update position for full step
    for (int i = 0; i < NUM_BODIES; i++) {
        bodies[i].x += dt * vx_half[i];
        bodies[i].y += dt * vy_half[i];
        bodies[i].z += dt * vz_half[i];
    }

    // Calculate acceleration at new time
    for (int i = 0; i < NUM_BODIES; i++) {
        calculate_acceleration(bodies, i, ax, ay, az);
    }

    // Update velocity for full step
    for (int i = 0; i < NUM_BODIES; i++) {
        bodies[i].vx += dt * ax[i];
        bodies[i].vy += dt * ay[i];
        bodies[i].vz += dt * az[i];
    }
}

void draw_bodies(Body bodies[]) {
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
    SDL_RenderClear(renderer);

    for (int i = 0; i < NUM_BODIES; i++) {
        Uint8 r, g, b, a;
        SDL_GetRGBA(bodies[i].color, bodies[i].color & 0xFF000000 ? bodies[i].color : 0xFF000000, &r, &g, &b, &a);
        SDL_SetRenderDrawColor(renderer, r, g, b, a);
        SDL_RenderDrawPoint(renderer, WINDOW_WIDTH / 2 + bodies[i].x * SCALE, WINDOW_HEIGHT / 2 + bodies[i].y * SCALE);
    }

    SDL_RenderPresent(renderer);
}

int main(int argc, char* argv[]) {
    Body bodies[NUM_BODIES];

    // Initialize SDL
    SDL_Init(SDL_INIT_VIDEO);

    // Create window
    window = SDL_CreateWindow("Solar System Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_SHOWN);

    // Create renderer
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    // Initialize the bodies with their initial positions, velocities, masses, and colors
    // (You need to provide the initial conditions)
    // Example:
    bodies[0].x = 0.0;
    bodies[0].y = 0.0;
    bodies[0].z = 0.0;
    bodies[0].vx = 0.0;
    bodies[0].vy = 0.0;
    bodies[0].vz = 0.0;
    bodies[0].mass = 1.9885e30;
    bodies[0].color = 0xFFFF00FF; // Sun (yellow)

    // Simulation loop
    int num_steps = 365 * 100; // Simulate for 100 years
    for (int step = 0; step < num_steps; step++) {
        leapfrog_integrate(bodies, DT);
        draw_bodies(bodies);
        SDL_Delay(10); // Adjust the delay to control the animation speed
    }

    // Clean up
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}