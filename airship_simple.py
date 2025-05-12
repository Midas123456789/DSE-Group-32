import math
from scipy.optimize import fsolve, fmin
from shapely.measurement import length
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from Aerodynamic_Atmospheric.ISA_Calculator import ISA_Calculator


class Airship:

    def __init__(self, FR, startvolume, lobes, velocity, altitude, payload):
        """
        Initializes an Airship object.

        Parameters:
        - FR (float): Fineness ratio
        - AR (float): Aspect ratio
        - volume (float): Volume of the airship in cubic meters
        - length (float): Length of the airship in meters
        - n_eng (int): Number of engines
        - payload (float): Payload capacity in kilograms
        - altitude (float): Operating altitude in meters
        """
        # Variables for lobes
        NLs = {1: 2, 2: 2.25, 3: 2.4, 4: 2.5, 5: 2.54}
        self.FR = FR
        self.n_lobes = lobes
        self.pound = 0.45359237

        self.volume = startvolume
        self.p = 1.6075
        self.reference_volume = self.volume ** (2 / 3)
        self.velocity = velocity
        self.altitude = altitude

        self.density_sl = 0.002377  # slung/ft^2 # density at SL
        self.isa = ISA_Calculator(altitude=self.altitude * 0.3048, velocity=self.velocity)
        self.density = self.isa.results[self.altitude * 0.3048]['Density [kg/m³]'] / 515.35549

        self.temperature = self.isa.results[self.altitude * 0.3048]['Temperature [K]']
        self.pressure = self.isa.results[self.altitude * 0.3048]['Pressure [Pa]']
        self.molarmass = 4.0026*10e-3 /self.pound
        self.idealpresconstant = 8.31446261815324 / self.pound
        self.helium = (self.pressure*self.molarmass)/(self.molarmass*self.idealpresconstant)

        self.mu = 3.66 * 10 ** -7  # find out
        self.mu_cr = 0.0209 * self.isa.dynamic_viscosity(self.isa.results[self.altitude * 0.3048]["Temperature [K]"])

        self.NL = NLs[self.n_lobes]  # find out, it is a number according to number of lobes

        self.gas_density = 0.0646  # lb/ft3 for helium at sea level

        self.payload = payload


        # Parameters for hybrid n_lobes >1

    def geomertic_parameters(self):
        """
        Calculates the geometric parameters of the airship.

        """
        # Calculate basic geometric parameters
        volume = self.volume
        diameter_e = (6 * volume / (math.pi * self.FR)) ** (1 / 3)
        length = self.FR * diameter_e
        ratio = -0.0178 * self.n_lobes ** 2 + 0.361 * self.n_lobes + 0.575

        # Calculate the diameter of the airship
        diameter_c = diameter_e/ ratio
        width = (1 + self.n_lobes) * diameter_c / 2
        self.asp_ratio = 4 * self.w ** 2 / (math.pi * length * width)

        # Calculate wet surface area
        surface = math.pi * ((length ** self.p * width ** self.p + length ** self.p * diameter_c ** self.p
                                   + width ** self.p * diameter_c ** self.p) / 3) ** (1 / self.p)
        ratioperimenter = 1.122 - 0.1226 * (1 / self.n_lobes)
        wetsurface = ratioperimenter * self.surface


        return diameter_e, length, diameter_c, width, self.asp_ratio, surface, wetsurface


    def aerodynamic_properties(self, length, wetsurface, asp_ratio):
        """
        Calculates the aerodynamic properties of the airship.

        """

        self.q = 0.5 * self.density * (self.velocity ** 2)
        Re = self.density * self.velocity * length / self.mu_cr
        Cf = 0.455 / (math.log10(Re) ** 2.58)
        FFbody = 1 + 1.5 / self.FR ** 1.5 + 7 / self.FR ** 3
        self.CD0 = FFbody * Cf * wetsurface / self.reference_volume

        K = (-0.0146 * (1 / self.asp_ratio) ** 4 + 0.182 * (1 / self.asp_ratio) ** 3 - 0.514 * (1 / self.asp_ratio) ** 2 +
                  0.838 * (1 / self.asp_ratio) - 0.053)
        self.K = K / self.NL

        return self.CD0

    # task 13
    def buoyancy_lift(self):
        """
        Calculates the buoyant lift of the airship.

        """
        self.W_gas = self.gas_density * self.volume                     #density changes high up
        self.buoyancy = self.volume*(self.density - self.gas_density)   # self.W_gas * self.density_sigma
        self.BR = 0.7 #self.buoyancy/self.Wg
        return self.buoyancy

    def fuel_weight(self):
        self.fuel_total = 300
        return self.fuel_total

    def prelimanary_weight(self):

        self.wzf = self.buoyancy/self.BR - self.fuel_total
        self.woe = self.wzf - self.payload
        self.wg = self.wzf + self.fuel_total
        return self.wzf, self.woe, self.wg

    # Exercise/Line 19
    def calculate_lift(self):
        """
        Calculate the lift force of the airship based on the weight of the airship and the buoyancy force.
        Inputs: WHO (float): weight of airship in lbf
        Outputs: lift (float): lift force in lbf
        """
        self.L_aer0 = self.buoyancy / self.BR - self.buoyancy
        self.CL = self.L_aer0 / (self.q * self.reference_volume)
        self.CD = (self.CD0 + self.K * self.CL ** 2)

        return self.L_aer0


    def altitude_iteration(self):

        alt_step = 1000
        graph_df = [self.altitude, self.volume, self.wg, self.CL, self.CD]
        altitude = np.linspace(0, self.altitude, alt_step)

        for altitudes in altitude:
            airship = Airship(self.FR, self.volume, self.n_lobes, self.velocity, altitudes, self.payload)
            airship.complete()
            print(airship.q)
            row = [altitudes, airship.volume, airship.wg, airship.CL, airship.CD]

            graph_df = np.vstack((graph_df, row))




