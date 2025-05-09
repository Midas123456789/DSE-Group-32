from airship import Airship

hybrid = Airship(3, 2000000, 3, 84.5, 60e3)
'''
hybrid.iterator(2e6)
print(hybrid)
hybrid.iterator(1.5e6)
print(hybrid)
hybrid.iterator(1.3e6)
print(hybrid)
hybrid.iterator(1e6)
print(hybrid)
'''
hybrid.iterate_to_exact()
print(hybrid)

'''
hybrid.geomertic_parameters()
hybrid.tailvolume()
hybrid.aerodynamic_properties()
hybrid.buoyant_lift()
hybrid.prelimanary_weight()
hybrid.fuel_calculations()
hybrid.weight_calculations()
hybrid.calculate_lift()
hybrid.engine()
hybrid.further_weight()
hybrid.ballonet()
'''