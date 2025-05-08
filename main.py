from airship import Airship

hybrid = Airship(3, 2000000, 3, 84.5, 4000)

hybrid.geomertic_parameters()
hybrid.tailvolume()
hybrid.aerodynamic_properties()
hybrid.buoyant_lift()
print(hybrid.aerodynamic_properties())