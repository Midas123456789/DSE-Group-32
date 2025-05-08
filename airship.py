import math

class Airship:
    def __init__(self, FR, volume, lobes, velocity, altitude):
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
        self.density_sl = 0.0765 #lb/ft^3 # density at SL
        self.mu = 3.66*10**-7               # find out
        self.n_engines = 4
        self.NL = 2.4               # find out
        self.gas_density =  0.0646 #lb/ft3 for helium at sea level
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
        self.ratio_dedc = -0.0178*self.n_lobes**2 + 0.361*self.n_lobes + 0.575

        # Calculate the diameter of the airship
        self.dc =self.de/self.ratio_dedc
        self.ht = self.dc                           # diameter of the lobes
        self.w = (1+self.n_lobes)*self.dc/2
        self.AR = 4*self.w**2/(math.pi*self.length*self.w)

        # Calculate wet surface area
        self.surface = math.pi*((self.length**self.p*self.w**self.p + self.length**self.p*self.ht**self.p + self.w**self.p*self.ht**self.p)/3)**(1/self.p)
        self.ratioperimenter = 1.122-0.1226*(1/self.n_lobes)
        self.wetsurface = self.ratioperimenter * self.surface

        return self.de, self.length, self.dc, self.w, self.ht, self.AR, self.surface, self.wetsurface
    
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
        self.CD0 = self.FFbody * self.Cf * self.wetsurface / self.reference_volume

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
    
    def buoyant_lift(self):
        """
        Calculates the buoyant lift of the airship.

        """
        self.buoyancy= self.gas_density*self.volume*self.density/self.density
    
    

    def print(self):
        return (f"Airship(FR={self.FR}, AR={self.AR}, volume={self.volume} m³, de={self.de} m, "
                f"length={self.length} m, dc={self.dc} m")
    
1
    