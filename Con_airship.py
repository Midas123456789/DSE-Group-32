import math
from scipy.optimize import fsolve, fmin
from Aerodynamic_Atmospheric.ISA_Calculator import ISA_Calculator

class Con_Airship:
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
        self.FR = FR
        self.n_lobes = lobes
        self.volume = volume
        self.p = 1.6075     
        self.reference_volume = self.volume**(2/3)
        self.velocity = velocity
        self.altitude = altitude
        self.density = 0.00211              # problem for later density at
        self.maxdensity = 0.001868
        self.density_sl = 0.002377 #slung/ft^2 # density at SL
        self.isa = ISA_Calculator(altitude=self.altitude*0.3048,velocity=self.velocity)
        self.density = self.isa.results[self.altitude*0.3048]['Density [kg/m³]']/515.35549
        self.density_sigma = self.density/self.density_sl
        self.mu = 3.66*10**-7               # find out
        self.n_engines = 4
        self.NL = 2.4               # find out, it is a number according to number of lobes
        self.gas_density =  0.0646 #lb/ft3 for helium at sea level
        self.fuelres = 1251        #find out later
        self.efficienty_eng = 0.65

        self.payload = payload
        #self.length = length
        #self.n_eng = n_eng
        #self.payload = payload
        #self.altitude = altitude
        #self.volume2_3 = self.volume ** (2 / 3)

        #Parameters for hybrid n_lobes >1



    
    def geomertic_parameters(self):
        """
        Calculates the geometric parameters of the airship.

        """
        # Calculate basic geometric parameters
        self.de = (6*self.volume/(math.pi*self.FR))**(1/3)
        self.length = self.FR * self.de

        # Calculate the diameter of the airship
        self.AR = 4*self.de/(math.pi*self.length)

        # Calculate wet surface area
        self.surface = math.pi*((self.length**self.p*self.de**self.p + self.length**self.p*self.de**self.p + self.de**self.p*self.de**self.p)/3)**(1/self.p)

        return self.de, self.length, self.AR, self.surface
    
    def tailvolume(self):
        """
        Calculates the tail volume and surface area of the airship.

        """

        self.CHT = -0.0051*(1000000/self.volume)+0.0717
        self.CVT = -0.0049*(1000000/self.volume)+0.0641

        self.lengthtail = 0.38*self.length

        self.surface_ht = self.CHT*(self.reference_volume*self.length)/self.lengthtail
        self.surface_vt = self.CVT*(self.reference_volume*self.length)/self.lengthtail

        return self.CHT, self.CVT, self.lengthtail, self.surface_ht, self.surface_vt


    def aerodynamic_properties(self):
        """
        Calculates the aerodynamic properties of the airship.

        """

        self.q = 0.5 * self.density * (self.velocity**2)
        self.Re = self.density*self.velocity*self.length/self.mu
        self.Cf = 0.455 / (math.log10(self.Re)**2.58)
        self.FFbody = 1+ 1.5/self.FR**1.5+7/self.FR**3
        self.CD0 = self.FFbody * self.Cf * self.surface / self.reference_volume

        # Drag coefficient
        self.tctail = 0.15  # find the formula later
        self.ARtail = 1     # find aspect ratio later
        self.FFtail = 1 +1.2*self.tctail+100*self.tctail**4
        
        self.ctail = 0.5*((self.ARtail*self.surface_ht/2)**0.5 +(self.ARtail*self.surface_vt/2)**0.5)
        self.Retail = self.density*self.velocity*self.ctail/self.mu
        self.Cfetail = 0.455 / (math.log10(self.Retail)**2.58)

        self.surfacewettails = 2.2*(self.surface_ht + self.surface_vt)
        self.CD0tail = self.FFtail * self.Cfetail * (self.surfacewettails) / self.reference_volume

        # Calculate rest of CD0

        self.CD0_cab_gond = (0.108 *self.CD0 * self.reference_volume + 7.7) / self.reference_volume
        self.CD0_eng_nac = (self.n_engines * 4.25) / self.reference_volume
        self.CD0_eng_cooling = (self.n_engines * (2e-6 * self.volume + 4.1)) / self.reference_volume
        self.CD0_eng_mount = (0.044 * self.CD0 * self.reference_volume + 0.92) / self.reference_volume
        self.CD0_cables = (9.7e-6 * self.volume + 10.22) / self.reference_volume
        self.CD0_acls = 0.0002  # From the image
        self.CD0int = (4.78e-6 * self.volume)/self.reference_volume

        self.CD0 = self.CD0 + self.CD0tail + self.CD0_cab_gond + self.CD0_eng_nac + self.CD0_eng_cooling + self.CD0_eng_mount + self.CD0_cables + self.CD0_acls + self.CD0int

        self.K = -0.0146*(1/self.AR)**4+0.182*(1/self.AR)**3-0.514*(1/self.AR)**2+0.838*(1/self.AR)-0.053
        self.K = self.K/self.NL
        
        return self.CD0, self.CD0tail, self.K
    #task 13
    def buoyant_lift(self):
        """
        Calculates the buoyant lift of the airship.

        """
        self.buoyancy= self.gas_density*self.volume*self.density_sigma
        self.BR = 0.9   #self.buoyancy/self.Wg
        return self.buoyancy
    
    def prelimanary_weight(self):

        self.wzf = self.buoyancy/self.BR-self.fuelres
        self.woe = self.wzf-self.payload
        return self.wzf, self.woe

    def fuel_calculations(self):


        self.range = 1000
        self.BSFC = 0.48
        self.wland = self.wzf + self.fuelres
        self.wh1 = self.wland-self.buoyancy
        self.A = (326 * self.efficienty_eng) / (self.BSFC * (self.K * self.CD0) ** 0.5)
        self.B = self.q * self.reference_volume * (self.CD0/ self.K) ** 0.5
        self.wh0 = self.B*math.tan((self.range/self.A)+math.atan(self.wh1/self.B))
        self.fuel_total = self.wh0-self.wh1+self.fuelres
        if self.fuel_total < 0:
            self.fuel_total = 0
            print(f"Fuel is negative: {self.fuel_total} lb")
        return self.fuel_total

    def weight_calculations(self):
        self.wg = self.wzf+self.fuel_total
        self.br_takeoff = self.buoyancy/self.wg
        return self.wg, self.br_takeoff

    # Exercise/Line 19
    def calculate_lift(self):
        """
        Calculate the lift force of the airship based on the weight of the airship and the buoyancy force.
        Inputs: WHO (float): weight of airship in lbf
        Outputs: lift (float): lift force in lbf
        """
        self.L_a = self.wh0
        self.density_sea = 0.002377
        #self.K = 0.295
        self.vmax = 1.1*self.velocity
        self.qmax = 0.5 * self.density_sea * (self.vmax)**2  # dynamic pressure in lbf/ft²
        self.CL_maxpower = self.wh0 / (self.qmax * self.reference_volume)  # lift coefficient at maximum power
        self.D_maxpower = (self.CD0 + self.K * self.CL_maxpower ** 2) * self.qmax * self.reference_volume  # drag force at maximum power in lbf

        return self.L_a, self.qmax, self.CL_maxpower, self.D_maxpower

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
        self.n = 20
        #self.D_maxpower = 13346
        self.P_hp_per_engine = ((self.D_maxpower * self.vmax) / (self.n_engines * self.efficienty_eng)) / 550  # power required per engine in hp
        self.P = self.P_hp_per_engine * 550  # Convert power from hp to ft-lbf/s
        self.C_S = ((self.density_sea * self.velocity ** 5) / (self.P * self.n ** 2)) ** (1 / 5)  # speed coefficient [dimensionless]
        self.J = 0.156 * self.C_S ** 2 + 0.241 * self.C_S + 0.138  # propeller advance ratio [dimensionless]
        self.D_prop = self.velocity / (self.J * self.n)  # propeller diameter in ft
        self.eta_prop = 0.139 * self.C_S ** 3 - 0.749 * self.C_S ** 2 + 1.37 * self.C_S + 0.0115  # propeller efficiency [dimensionless]


        return self.P_hp_per_engine


    # Exercise/Line 27
    def further_weight(self):
        """
        Calculate the hull fabric load based on the safety factor, internal pressure, and radius of the hull.
        Inputs: FS (float): safety factor [dimensionless]
                P_int (float): internal pressure in psi
                height (float): height of the hull in ft
        Outputs: q_hull (float): hull fabric load in lb/in
                q for shear flow in lb/in

        """
        self.P_int = (1.2 * self.qmax + 0.0635 * self.de)  # internal pressure in lbf/ft²
        self.P_int = self.P_int / 144  # Convert from lbf/ft² to psi
        self.height_inches = self.de * 12  # Convert height from ft to inches
        # NOTE 1 ft = 12 inches
        self.FS = 4  # factor of safety
        self.q_hull = self.FS * self.P_int * (self.height_inches / 2)  # hull fabric load in psi

        material_strengths = {
            'Polyester (weave)': {'a': 0.0453, 'b': 1.962},
            'Vectran (weave)': {'a': 0.0141, 'b': 1.882},
            'Vectran (laminate)': {'a': 0.0085, 'b': 1.365},
            'Dyneema (laminate)': {'a': 0.0063, 'b': 0.889},
            'material': {'a': 0.0085, 'b': 1.365}
        }

        # Linear relationship calculation using coefficients a and b
        a = material_strengths['material']['a']
        b = material_strengths['material']['b']
        self.w_hull = a * self.q_hull + b

        # NOTE 16 oz = 1 lb, 1 yd² = 9 ft², 1.2 = manufacturing factor, 1.26 = attachements factor
        self.W_envelope = self.w_hull * self.surface * 1.2 * 1.26 / (16 * 9)  # conversion to lbf/ft²

        self.f_septum = a * (1.5 * self.q_hull) + b
        self.W_septum = 2*1.06*self.f_septum*0.75*math.pi*self.de*self.length/4/16/9
        self.W_body = self.W_envelope + self.W_septum

        return self.q_hull

    def ballonet(self):

        self.vball = self.volume*((1/self.density_sigma-1))
        self.nball = 6
        self.surfaceballonet = math.pi*(3*self.vball/math.pi/6)**(2/3)*self.nball
        self.W_ball = 0.035*self.surfaceballonet

        self.faf = 1.26
        self.W_ssf = (self.surface_ht+self.surface_vt)*self.faf*0.8
        self.W_cs = (self.surface_ht+self.surface_vt)*0.2
        self.W_tails = self.W_ssf + self.W_cs
        self.W_act = 1.15*(self.surface_ht+self.surface_vt)*0.79*0.2

        self.W_crewstat = 1426
        self.W_gond = 1.875*2*(55*10+55*10+10*10)
        self.W_eng = self.n_engines*4.848*(self.P_hp_per_engine)**0.7956
        self.W_eng_mount = 0.64*self.W_eng
        self.W_ec = 60.27*(150*self.n_engines/100)**0.724
        self.W_start = 98

        self.Kp = 31.92
        self.nblades = 3
        self.W_prop = self.Kp*self.n_engines*(self.nblades)**0.391*(self.D_prop*self.P_hp_per_engine/1000)**0.782
        self.W_fueltank = 2.49*(self.fuel_total/6)**0.6 *(2)**0.2 *self.n_engines**0.13
        self.W_pressuresys = 0.02*self.woe

        self.W_acls = 1.6*4057
        self.W_vms = 1493
        self.W_Elect = 33.73*(470+500)**0.51
        self.W_Msys = 0.05*self.woe
        self.W_crew = 1148
        self.W_fuel = 0.01*self.fuel_total
        self.W_margin = 0.06*self.woe
        self.W_both = self.W_fuel + self.W_margin
        self.woe2 = self.W_body+self.W_ball+self.W_tails +self.W_crewstat+self.W_gond+self.W_eng+self.W_eng_mount+self.W_ec+self.W_start+self.W_prop+self.W_fueltank+self.W_pressuresys+self.W_acls+self.W_vms+self.W_Elect+self.W_Msys+self.W_crew+self.W_both
        self.W_g2 = self.woe2+self.fuel_total+self.payload

        return
    def iterator(self,Volume):
        if Volume < abs(1e5):
            return 1e6
        self.volume = abs(Volume[0])
        self.geomertic_parameters()
        self.tailvolume()
        self.aerodynamic_properties()
        self.buoyant_lift()
        self.prelimanary_weight()
        self.fuel_calculations()
        self.weight_calculations()
        self.calculate_lift()
        self.engine()
        self.further_weight()
        self.ballonet()
        #print(self)
        return abs(self.wg - self.W_g2)
    def iterate_to_exact(self):
        #fsolve(self.iterator,2000000,xtol=1e-3)
        fmin(self.iterator,self.volume)
        return

    def __str__(self):
        return (f"Airship(volume={self.volume} ft³, wg ={self.wg} lb, wg2={self.W_g2} lb, difference= {self.wg - self.W_g2}")
               # f"length={self.length} m, dc={self.dc} m")

    