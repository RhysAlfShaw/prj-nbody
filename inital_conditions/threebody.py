import numpy as np

class ThreeBody:

    """

    Class to create a 3 body system with one planet, inspired by the 3-body problem.

    M1 is the 
    """

    def __init__(self, M1, M2, M3, Mp=3E-6, G= 6.67430e-11, Msol=2E30, Star_sep = 10, Planet_R=1, v_orbit=0.1):

        self.M1 = M1*Msol
        self.M2 = M2*Msol
        self.M3 = M3*Msol
        self.Mp = Mp*Msol
        self.Star_sep = Star_sep
        self.Planet_R = Planet_R
        self.G = G
        self.v_orbit = v_orbit

    def calculate_planet_orbit(self):
        """
        Calculate the orbit of the planet around star M1
        """

        F = self.G*self.M1*self.Mp/self.Star_sep**2
        a = F/self.Mp
        v = np.sqrt(a*self.Planet_R)
        return v

    def calculate_star_orbit(self):
        """
        Calculate the orbit of the star around the planet
        """
        F = self.G*self.M2*self.M3/self.Star_sep**2
        a = F/self.M2
        v = np.sqrt(a*self.Star_sep)
        return v
    
    def position_stars(self):

        # Star 1
        x1 = self.Star_sep   
        y1 = 0  
        z1 = 0  
        
        # Star 2
        x2 = self.Star_sep
        y2 = 0
        z2 = 0
        
        # Star 3
        x3 = 0
        y3 = 0
        z3 = 0

        # Planet
        x3p = self.Star_sep + self.Planet_R
        y3p = 0
        z3p = 0

        return x1, y1, x2, y2, x3, y3

    def velocity_stars(self):

        # Star 1
        v1x = 0
        v1y = self.v_orbit
        v1z = 0

        # Star 2
        v2x = 0
        v2y = self.v_orbit
        v2z = 0

        # Star 3
        v3x = 0
        v3y = self.v_orbit
        v3z = 0

        # Planet
        v3px = 0
        v3py = self.calculate_planet_orbit()
        v3pz = 0

        return v1x, v1y, v2x, v2y, v3x, v3y, v3px, v3py
        