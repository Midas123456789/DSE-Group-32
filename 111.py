import numpy as np

class Wing:
    
    def __init__(self):
        
        self.S = S      # m
        self.taper = taper  # m
        self.b = b # m
    
    def root_chord(self):
        self.c_r = (2 * self.S) / ((1 + self.taper) * self.b)
        return self.c_r

    def tip_chord(self):
        self.c_t = 2 * self.taper
        return self.c_t