import math

class Airship:
    def __init__(self, FR, volume, lobes):
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
        #self.length = length
        #self.n_eng = n_eng
        #self.payload = payload
        #self.altitude = altitude
        #self.volume2_3 = self.volume ** (2 / 3)

        #Parameters for hybrid n_lobes >1
        
        self.de = (6*self.volume/(math.pi*self.FR))**(1/3)
        self.length = self.FR * self.de
        self.ratio_dedc = -0.0178*self.n_lobes**2 + 0.361*self.n_lobes + 0.575
        self.dc =self.de/self.ratio_dedc
        self.ht = self.dc
        self.w = (1+self.n_lobes)*self.dc/2
        self.AR = 4*self.w**2/(math.pi*self.length*self.w)
        self.surface = math.pi((self.length**self.p*self.w**p+self.length**self.p*self.ht**self.p+self.w**self.p+self.ht**self.p)/3)**(1/self.p)
        #self.Swet = 
    


    def print(self):
        return (f"Airship(FR={self.FR}, AR={self.AR}, volume={self.volume} m³, de={self.de} m, "
                f"length={self.length} m, dc={self.dc} m")
    
1
    