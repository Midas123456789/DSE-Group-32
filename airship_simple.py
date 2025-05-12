import math
from scipy.optimize import fsolve, fmin
from Aerodynamic_Atmospheric.ISA_Calculator import ISA_Calculator


class Airship:

    def __init__(self, FR, volume, lobes, velocity, altitude, payload):
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

        self.volume = volume
        self.p = 1.6075
        self.reference_volume = self.volume ** (2 / 3)
        self.velocity = velocity
        self.altitude = altitude
        self.maxdensity = 0.001868
        self.density_sl = 0.002377  # slung/ft^2 # density at SL
        self.isa = ISA_Calculator(altitude=self.altitude * 0.3048, velocity=self.velocity)
        self.density = self.isa.results[self.altitude * 0.3048]['Density [kg/m³]'] / 515.35549
        self.density_sigma = self.density / self.density_sl
        self.mu = 3.66 * 10 ** -7  # find out
        self.mu_cr = 0.0209 * self.isa.dynamic_viscosity(self.isa.results[self.altitude * 0.3048]["Temperature [K]"])
        self.n_engines = 4
        self.NL = NLs[self.n_lobes]  # find out, it is a number according to number of lobes
        self.gas_density = 0.0646  # lb/ft3 for helium at sea level
        self.fuelres = 1251  # find out later
        self.efficienty_eng = 0.65
        self.range = range
        self.payload = payload


        # Parameters for hybrid n_lobes >1

    def geomertic_parameters(self):
        """
        Calculates the geometric parameters of the airship.

        """
        # Calculate basic geometric parameters
        self.de = (6 * self.volume / (math.pi * self.FR)) ** (1 / 3)
        self.length = self.FR * self.de
        self.ratio_dedc = -0.0178 * self.n_lobes ** 2 + 0.361 * self.n_lobes + 0.575

        # Calculate the diameter of the airship
        self.dc = self.de / self.ratio_dedc
        self.ht = self.dc  # diameter of the lobes
        self.w = (1 + self.n_lobes) * self.dc / 2
        self.AR = 4 * self.w ** 2 / (math.pi * self.length * self.w)

        # Calculate wet surface area
        self.surface = math.pi * ((self.length ** self.p * self.w ** self.p + self.length ** self.p * self.ht ** self.p
                                   + self.w ** self.p * self.ht ** self.p) / 3) ** (1 / self.p)
        self.ratioperimenter = 1.122 - 0.1226 * (1 / self.n_lobes)
        self.wetsurface = self.ratioperimenter * self.surface

        return self.de, self.length, self.dc, self.w, self.ht, self.AR, self.surface, self.wetsurface


    def aerodynamic_properties(self):
        """
        Calculates the aerodynamic properties of the airship.

        """

        self.q = 0.5 * self.density * (self.velocity ** 2)
        self.Re = self.density * self.velocity * self.length / self.mu_cr
        self.Cf = 0.455 / (math.log10(self.Re) ** 2.58)
        self.FFbody = 1 + 1.5 / self.FR ** 1.5 + 7 / self.FR ** 3
        self.CD0 = self.FFbody * self.Cf * self.wetsurface / self.reference_volume

        self.K = -0.0146 * (1 / self.AR) ** 4 + 0.182 * (1 / self.AR) ** 3 - 0.514 * (1 / self.AR) ** 2 + 0.838 * (
                    1 / self.AR) - 0.053
        self.K = self.K / self.NL

        return self.CD0, self.K

    # task 13
    def buoyant_lift(self):
        """
        Calculates the buoyant lift of the airship.

        """
        self.buoyancy = self.gas_density * self.volume * self.density_sigma
        self.BR = 0.7 #self.buoyancy/self.Wg
        return self.buoyancy

    def fuel_calculations(self):

        self.fuel_total = 300
        return self.fuel_total

    def prelimanary_weight(self):

        self.wzf = self.buoyancy / self.BR - self.fuelres
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
        self.density_sea = 0.002377
        self.CL = self.L_aer0 / (self.q * self.reference_volume)
        self.CD = (self.CD0 + self.K * self.CL ** 2)

        return self.L_aer0

    # Exercise/Line 22
    def engine(self):
        """
        Calculate the power required per engine at maximum velocity.
        Inputs: Vmax (float): maximum velocity of the airship in ft/s
                D_maxpower (float): drag force at maximum power in lbf
                n_eng (float): engine efficiency [dimensionless]
                NE (int): number of engines [#]
        Outputs: P_per_engine (float): power required per engine in hp

        """

        # 550 = 550 ft-lbf/s = 1 hp
        self.n = 10
        # self.D_maxpower = 13346
        self.P_hp_per_engine = ((self.CD*self.q*self.reference_volume) / (
                    self.n_engines * self.efficienty_eng)) / 550  # power required per engine in hp
        return self.P_hp_per_engine

    def complete(self):

        self.geomertic_parameters()
        self.aerodynamic_properties()
        self.buoyant_lift()
        self.fuel_calculations()
        self.prelimanary_weight()
        self.calculate_lift()
        self.engine()

        return

