import numpy as np

class Planets:

    def __init__(self):
        self.X = []
        self.Y = []
        self.Z = []
        self.Mass = []
        self.Vx = []
        self.Vy = []
        self.Vz = []
        self.G = 6.67430E-11
        self.Mass = []

    def add_star(self,Mass):
        self.Mass.append(Mass)
        self.X.append(0)
        self.Y.append(0)
        self.Z.append(0)
        self.Vx.append(0)
        self.Vy.append(0)
        self.Vz.append(0)


    def add_planet(self,Mass,Oribital_radius,eccentricity,inclination):

        # turn the orbital semi-major axis into a x,y,z position, assume we are on the x-y plane
        self.X.append(Oribital_radius)
        self.Y.append(0)
        self.Z.append(0)

        # using star mass self.mass[0] and the orbital radius we can calculate the orbital velocity
        # v = sqrt(GM/r)
        self.Vx.append(0)
        self.Vy.append(np.sqrt(self.G*self.Mass[0]/Oribital_radius))
        self.Vz.append(0)
        self.Mass.append(Mass)