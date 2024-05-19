# Nbody implementation Analysis.

Here we plan to create a simple Nbody implementation and analyse the preformance of the code in multiple languages.

We attempt, for each language, to implement a single threaded, multi-threaded and GPU version of the code. We will then compare the performance of the code in each language.

The languages we will be using are:
- c++
- Python
- Julia
- Fortran
- cython
- Rust?

## The Nbody problem

The Nbody problem is a simple problem in physics. It is the problem of predicting the motion of a group of celestial objects that interact with each other gravitationally. The problem is defined by the following equation:
```
F = G * m1 * m2 / r^2
```
Where F is the force between two objects, G is the gravitational constant, m1 and m2 are the masses of the two objects and r is the distance between the two objects.

The force between two objects is a vector, and the total force on an object is the sum of the forces between that object and all other objects.

The acceleration of an object is given by:
```
a = F / m
```
Where a is the acceleration of the object and m is the mass of the object.

The position of an object is given by:
```
r = r + v * dt + 0.5 * a * dt^2
```
Where r is the position of the object, v is the velocity of the object, a is the acceleration of the object, dt is the time step and t is the time.

The velocity of an object is given by:
```
v = v + 0.5 * (a + a_new) * dt
```
Where v is the velocity of the object, a is the acceleration of the object, a_new is the new acceleration of the object, dt is the time step and t is the time.

## The Nbody implementation

To preform the numerical interation of the Nbody problem we will use the leapfrog integration method. This method is a second order symplectic integrator that is commonly used in Nbody simulations.
the leapfrog integration method is defined by the following equations:
```
v = v + 0.5 * a * dt
r = r + v * dt
a = F / m
v = v + 0.5 * a * dt
```
Where v is the velocity of the object, r is the position of the object, a is the acceleration of the object, F is the force on the object and m is the mass of the object.