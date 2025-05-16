from Aerosandbox_models.LH2_tandem import *
import aerosandbox as asb
plane = Tandem_LH()
plane.solve()
plane._get_solution_plane()
sol_plane = plane.airplane_sol
plane.print_solution()
sol_plane.draw()

