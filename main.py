from Aerosandbox_models.LH2_tandem import *
import aerosandbox as asb
plane = Tandem_LH()
plane.solve()
plane._get_solution_plane()
sol_plane = plane.airplane_sol
plane.print_solution()
#sol_plane.draw()

alpha = np.linspace(-20, 20, 1000)
aero = asb.AeroBuildup(
    airplane=sol_plane,
    op_point=asb.OperatingPoint(velocity=30, alpha=alpha, beta=0),
).run()

plt.plot(alpha, aero["CL"]/aero["CD"])
plt.show()