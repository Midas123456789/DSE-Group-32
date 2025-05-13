import os
import sys
import math
from scipy.optimize import fsolve, fmin
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Setup path to import ISA_Calculator
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
isa_dir = os.path.join(parent_dir, 'Aerodynamic_Atmospheric', 'Data')
sys.path.insert(0, isa_dir)

from Aerodynamic_Atmospheric.ISA_Calculator import ISA_Calculator


class Airship:
    def __init__(self, FR, startvolume, lobes, velocity, altitude, payload):
        """
        Initializes an Airship object.
        """
        # Constants and parameters
        self.FR = FR
        self.n_lobes = lobes
        self.startvolume = startvolume
        self.volume = startvolume
        self.velocity = velocity
        self.altitude = altitude
        self.payload = payload

        self.p = 1.6075
        self.reference_volume = self.volume ** (2 / 3)

        # Conversion constants
        self.pound = 0.45359237
        self.feet = 3.2808399

        # Lobes scaling factor
        NLs = {1: 2, 2: 2.25, 3: 2.4, 4: 2.5, 5: 2.54}
        self.NL = NLs[self.n_lobes]

        # Atmospheric properties
        self.isa = ISA_Calculator(altitude=self.altitude * 0.3048, velocity=self.velocity)
        self.density = self.isa.results[self.altitude * 0.3048]['Density [kg/m³]'] / 515.35549
        self.temperature = self.isa.results[self.altitude * 0.3048]['Temperature [K]']
        self.pressure = self.isa.results[self.altitude * 0.3048]['Pressure [Pa]']

        self.mu = 3.66e-7
        self.mu_cr = 0.0209 * self.isa.dynamic_viscosity(self.temperature)

        # Helium properties
        self.molarmass = 4.0026e-3 / self.pound
        self.idealpresconstant = 8.31446261815324 / (self.pound * self.feet ** 2)
        self.gas_density = (self.pressure * self.molarmass) / (self.temperature * self.idealpresconstant)

        #self.gas_density = 0.0646  # lb/ft^3 (helium at sea level)

    def geomertic_parameters(self):
        """
        Calculates the geometric parameters of the airship.
        """
        diameter_e = (6 * self.volume / (math.pi * self.FR)) ** (1 / 3)
        self.length = self.FR * diameter_e
        ratio = -0.0178 * self.n_lobes ** 2 + 0.361 * self.n_lobes + 0.575

        diameter_c = diameter_e / ratio
        width = (1 + self.n_lobes) * diameter_c / 2
        self.asp_ratio = 4 * width ** 2 / (math.pi * self.length * width)

        surface = math.pi * ((self.length ** self.p * width ** self.p +
                              self.length ** self.p * diameter_c ** self.p +
                              width ** self.p * diameter_c ** self.p) / 3) ** (1 / self.p)

        ratioperimeter = 1.122 - 0.1226 * (1 / self.n_lobes)
        self.wetsurface = ratioperimeter * surface

        return diameter_e, self.length, diameter_c, width, self.asp_ratio, surface, self.wetsurface

    def aerodynamic_properties(self):
        """
        Calculates the aerodynamic properties of the airship.
        """
        self.q = 0.5 * self.density * (self.velocity ** 2)
        Re = self.density * self.velocity * self.length / self.mu_cr
        Cf = 0.455 / (math.log10(Re) ** 2.58)
        FFbody = 1 + 1.5 / self.FR ** 1.5 + 7 / self.FR ** 3
        self.CD0 = FFbody * Cf * self.wetsurface / self.reference_volume

        K = (-0.0146 * (1 / self.asp_ratio) ** 4 + 0.182 * (1 / self.asp_ratio) ** 3 -
             0.514 * (1 / self.asp_ratio) ** 2 + 0.838 * (1 / self.asp_ratio) - 0.053)
        self.K = K / self.NL

        return self.CD0

    def buoyancy_lift(self):
        """
        Calculates the buoyant lift of the airship.
        """
        self.W_gas = self.gas_density * self.volume
        self.buoyancy = self.volume * (self.density - self.gas_density)
        self.BR = 0.7
        return self.buoyancy

    def fuel_weight(self):
        self.fuel_total = 300  # constant fuel weight (can be made dynamic)
        return self.fuel_total

    def prelimanary_weight(self):
        self.wzf = self.buoyancy / self.BR - self.fuel_total
        self.woe = self.wzf - self.payload
        self.wg = self.wzf + self.fuel_total
        return self.wzf, self.woe, self.wg

    def calculate_lift(self):
        """
        Calculate the lift force of the airship.
        """
        self.L_aer0 = self.buoyancy / self.BR - self.buoyancy
        self.CL = self.L_aer0 / (self.q * self.reference_volume)
        self.CD = self.CD0 + self.K * self.CL ** 2
        return self.L_aer0

    def complete(self):
        self.geomertic_parameters()
        self.aerodynamic_properties()
        self.buoyancy_lift()
        self.fuel_weight()
        self.prelimanary_weight()
        self.calculate_lift()

    def altitude_iteration(self):
        """
        Performs iteration over different altitudes and plots results.
        """
        alt_step = 1000
        altitudes = np.linspace(0, self.altitude, alt_step)
        graph_data = []

        for alt in altitudes:
            # Update volume based on pressure ratio (ideal gas approximation)
            isa_local = ISA_Calculator(altitude=alt * 0.3048, velocity=self.velocity)
            pressure = isa_local.results[alt * 0.3048]['Pressure [Pa]']
            volume_new = self.startvolume * (self.pressure / pressure)

            airship = Airship(self.FR, volume_new, self.n_lobes, self.velocity, alt, self.payload)
            airship.complete()
            graph_data.append([alt, airship.volume, airship.wg, airship.CL, airship.CD])

        df = pd.DataFrame(graph_data, columns=['altitude', 'Volume', 'WeightG', 'CL', 'CD'])

        # Plotting
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))

        axs[0, 0].plot(df['Volume'], df['altitude'])
        axs[0, 0].set_title('Altitude vs Volume')
        axs[0, 0].set_xlabel('Volume [ft³]')
        axs[0, 0].set_ylabel('Altitude [ft]')

        axs[0, 1].plot(df['WeightG'], df['altitude'])
        axs[0, 1].set_title('Altitude vs Weight')
        axs[0, 1].set_xlabel('Weight [lb]')

        axs[1, 0].plot(df['CL'], df['altitude'])
        axs[1, 0].set_title('Altitude vs CL')
        axs[1, 0].set_xlabel('CL')

        axs[1, 1].plot(df['CD'], df['altitude'])
        axs[1, 1].set_title('Altitude vs CD')
        axs[1, 1].set_xlabel('CD')

        plt.tight_layout()
        plt.show()


# Run a sample airship calculation
if __name__ == "__main__":
    simplehybrid = Airship(FR=3, startvolume=2e8, lobes=3, velocity=34.5, altitude=40000, payload=1000)
    simplehybrid.complete()
    simplehybrid.altitude_iteration()
