import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
import math


class Mass_wing(asb.Wing):
    def __init__(
        self,
        name=None,
        xsecs=None,
        symmetric=True,
        color=None,
        analysis_specific_options=None,
    ):
        super().__init__(name=name, xsecs=xsecs, symmetric=symmetric, color=color, analysis_specific_options=analysis_specific_options)
    
    def M_root(self, L, b):
        n = 100
        x = np.linspace(0, b, 100)
        return L * b / 3 / math.pi
    
    def root_R(self):
        return self.xsecs[0].chord * max(self.xsecs[0].airfoil.local_thickness()) / 2
    
    def tip_R(self):
        return self.xsecs[-1].chord * max(self.xsecs[-1].airfoil.local_thickness()) / 2
    
    def root_t(self, L, b, SF=3, max_stress = 100*10**6):
        M = self.M_root(L, b)
        root_R = self.root_R()
        return SF * M / math.pi / root_R ** 2 / max_stress

    def spar_mass(self, L, b, SF=3, max_stress=100*10**6, density=1600):
        root_R = self.root_R()
        tip_R = self.tip_R()

        R_list = []
        num_sec = len(self.xsecs)
        span = self.span()
        spar_mass = 0
        for i in range(num_sec):
            R = self.xsecs[i].chord * max(self.xsecs[i].airfoil.local_thickness()) / 2
            Lx = L / span
            d = span * (num_sec - i) / (2 * num_sec)
            M_local = Lx * d ** 2 / 2
            t = SF * M_local / (math.pi * R**2 * max_stress)
            sec_len = span / (2 * num_sec)
            section_mass = sec_len * (2 * math.pi * R * t) * density
            spar_mass += section_mass 
            #print(R)
        return spar_mass
 
    
