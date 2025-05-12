from airship import Airship
from visualise import visualise
from airship_simple import Airship as AS

hybrid = Airship(3, 2e8, 3, 84.5, 60e3,1000,1000)
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
#hybrid.iterate_to_exact()
#print(hybrid)

simplehybrid = AS(3, 2e8, 3, 84.5, 60e3,1000)
simplehybrid.complete()
graph = visualise.altitude_graph(simplehybrid)


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