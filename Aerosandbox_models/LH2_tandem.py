import aerosandbox as asb
import aerosandbox.numpy as np
from mass_wing import Mass_wing
import matplotlib.pyplot as plt

class Tandem_LH:
    def __init__(self, N_cords = 5, wing_airfoil = asb.Airfoil("sd7037"), altitude = 18e3, mission_days = 7, W_payload = 1e3, P_payload = 5e3, tail_drag_factor = 1.05):
        #GEOMETRY
        self.N_cords = N_cords
        self.wing_airfoil = wing_airfoil
        
        #WEIGHT
        self.W_payload = W_payload
        
        #OPERATIONS
        self.altitude = altitude
        self.mission_days = mission_days
        self.second_in_day = 86400
        self.mission_seconds = mission_days * self.second_in_day
        
        #AERODYNAMICS
        self.tail_drag_factor = tail_drag_factor
        self.atm = asb.Atmosphere(altitude=altitude)
        
        #POWER
        self.P_payload = P_payload
        
        #CONSTANTS
        self.g = 9.81
        self.density_LH = 70.85 #kg/m3
        self.energy_density_LH = 142 * 10 ** 6 #J/kg
        self.power_density_PEMFCs = 40000 #W/kg
        self.PEMFCs_eff = 0.55
        self.motor_eff = 0.97 #Fundamentals of Aircraft and Airship Design
        self.propeller_eff = 0.85 #Fundamentals of Aircraft and Airship Design
        self.conversion_eff = self.PEMFCs_eff
        self.propulsive_eff = self.conversion_eff * self.motor_eff * self.propeller_eff
        self.fuel_tank_fuel_mass_fraction = 0.34  # from Brewer, Hydrogen Aircraft Technology pg. 29
        
        #OPTI
        self.opti = asb.Opti()
          
    def _variables_setup(self):
        """
        Setup variables
        """
        self.cords = self.opti.variable(init_guess=3 * np.ones(self.N_cords), n_vars=self.N_cords)
        self.b = self.opti.variable(init_guess=60, upper_bound=200, lower_bound=0)
        self.V = self.opti.variable(init_guess=30, upper_bound=200, lower_bound=1)
        self.alpha = self.opti.variable(init_guess=5, lower_bound=0, upper_bound=30)
        
    def _setup_op(self):
        self.op = asb.OperatingPoint(velocity=self.V, alpha=self.alpha, atmosphere=self.atm)
        
    def _plane_geom_setup(self):
        """
        Initialise parameterised mass wing. 
        """
        self.y_sections = np.linspace(0, self.b / 2, self.N_cords)
        self.wing1 = Mass_wing(
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[-0.25 * self.cords[i], self.y_sections[i], 0],
                    chord=self.cords[i],
                    airfoil=self.wing_airfoil
                )
                for i in range(self.N_cords)
            ]
        )
        
        self.wing2 = Mass_wing(
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[-0.25 * self.cords[i], self.y_sections[i], 0],
                    chord=self.cords[i],
                    airfoil=self.wing_airfoil
                )
                for i in range(self.N_cords)
            ]
        ).translate([9,0,0])
        self.airplane = asb.Airplane(wings=[self.wing1] + [self.wing2])
        self.wing_area = self.wing1.area() + self.wing2.area()
        
    def _setup_aero(self):
        self.aero = asb.AeroBuildup(airplane=self.airplane, op_point=self.op).run()
        self.L = 0.5 * self.atm.density() * self.V**2 * self.aero['CL'] * self.wing_area
        self.D = 0.5 * self.atm.density() * self.V**2 * self.aero['CD'] * self.wing_area * self.tail_drag_factor

    def _setup_constraints(self):
        self.opti.subject_to([
        self.cords > 0,
        np.diff(self.cords) <= 0,
        self.wing_area < 2000,
        self.L >= self.W_total
        ])
        
    def _get_spar(self):
        spar_mass, t_list = self.wing1.spar_mass(self.L, self.b)
        spar_mass *= 2
        self.W_spar = spar_mass * self.g
        self.t_list = t_list
        
    def compute_weights(self):
        # Required power
        self.P_required = self.V * self.D

        # Energy required for the mission
        self.E_LH = (1 / self.conversion_eff) * (self.P_required / self.propulsive_eff + self.P_payload) * self.mission_seconds

        # Hydrogen mass and weight
        self.M_LH = self.E_LH / self.energy_density_LH
        self.W_LH = self.M_LH * self.g

        # Fuel tank mass and weight
        self.M_tank = self.fuel_tank_fuel_mass_fraction * self.M_LH
        self.W_tank = self.M_tank * self.g

        # Fuel cell mass and weight
        self.M_fuel_cell = self.E_LH / self.mission_seconds / self.power_density_PEMFCs
        self.W_fuel_cell = self.M_fuel_cell * self.g

        # LH volume
        self.LH_volume = self.M_LH / self.density_LH

        # Skin weight
        skin_density = 0.25  # kg/m², adjust based on material
        self.W_skin = 2 * self.wing_area * skin_density * self.g

        # Spar weight
        self._get_spar()

        # Wing weight
        self.W_wing = self.W_spar + self.W_skin

        # Miscellaneous weight
        self.W_misc = self.W_wing

        # Total weight
        self.W_total = self.W_payload + self.W_wing + self.W_LH + self.W_misc + self.W_tank + self.W_fuel_cell


    def solve(self):
        self._variables_setup()
        self._setup_op()
        self._plane_geom_setup()
        self._setup_aero()
        self.compute_weights()
        self._setup_constraints()

        # Objective: minimize required power
        self.opti.minimize(self.P_required)

        # Solve the optimization problem
        sol = self.opti.solve(verbose=False)

        # Store the solution
        self.solution = {
            #GEOMETRY
            "cords": sol.value(self.cords),
            "b": sol.value(self.b),
            "wing_area": sol.value(self.wing_area),
            "t_list": sol.value(self.t_list),
            "Airfoil": self.wing_airfoil,
            "LH_volume": sol.value(self.LH_volume),
            #FLIGHT
            "V": sol.value(self.V),
            "alpha": sol.value(self.alpha),
            "CL": sol.value(self.aero['CL']),
            "CD": sol.value(self.aero['CD']) * self.tail_drag_factor,
            "L": sol.value(self.L),
            "D": sol.value(self.D),
            "P_required": sol.value(self.P_required),      
            "Altitude": self.altitude,
            "Atm": self.atm,
            #WEIGHT (self.W_payload + self.W_wing + self.W_LH + self.W_misc + self.W_tank + self.W_fuel_cell)
            "W_payload": sol.value(self.W_payload),
            "W_wing": sol.value(self.W_wing),
            "W_LH": sol.value(self.W_LH),
            "W_misc": sol.value(self.W_misc),
            "W_tank": sol.value(self.W_tank),
            "W_spar": sol.value(self.W_spar),
            "W_fuel_cell": sol.value(self.W_fuel_cell),
            "W_skin": sol.value(self.W_skin),
            "W_total": sol.value(self.W_total),
            #MASS
            "M_payload": sol.value(self.W_payload) / self.g,
            "M_wing": sol.value(self.W_wing) / self.g,
            "M_LH": sol.value(self.W_LH) / self.g,
            "M_misc": sol.value(self.W_misc) / self.g,
            "M_tank": sol.value(self.W_tank) / self.g,
            "M_spar": sol.value(self.W_spar) / self.g,
            "M_fuel_cell": sol.value(self.W_fuel_cell) / self.g,
            "M_skin": sol.value(self.W_skin) / self.g,
            "M_total": sol.value(self.W_total) / self.g
        }
        self._get_solution_plane()
        return self.solution
    
    def print_solution(self):
        sol = self.solution
        print("\n--- HYDROGEN TANDEM OPTIMIZATION RESULTS ---\n")

        print("GEOMETRY:")
        print(f"  Wing span (b):                 {sol['b']:.2f} m")
        print(f"  Wing area:                     {sol['wing_area']:.2f} m²")
        print(f"  Chord distribution (m):        {np.round(sol['cords'], 2)}")
        print(f"  Spar thicknesses (m):          {np.round(sol['t_list'], 4)}")
        print(f"  LH2 volume:                    {sol['LH_volume']:.2f} m³")
        print(f"  Airfoil:                       {sol['Airfoil'].name}\n")

        print("FLIGHT:")
        print(f"  Altitude:                      {sol['Altitude']} m")
        print(f"  Velocity:                      {sol['V']:.2f} m/s")
        print(f"  Angle of attack:               {sol['alpha']:.2f} deg")
        print(f"  Lift coefficient (CL):         {sol['CL']:.3f}")
        print(f"  Drag coefficient (CD):         {sol['CD']:.4f}")
        print(f"  Lift (L):                      {sol['L']:.2f} N")
        print(f"  Drag (D):                      {sol['D']:.2f} N")
        print(f"  Required power:                {sol['P_required'] / 1e3:.2f} kW\n")

        print("WEIGHTS (N):")
        print(f"  Payload:                       {sol['W_payload']:.2f}")
        print(f"  Wing total:                    {sol['W_wing']:.2f}")
        print(f"    └ Spar:                      {sol['W_spar']:.2f}")
        print(f"    └ Skin:                      {sol['W_skin']:.2f}")
        print(f"  LH2 fuel:                      {sol['W_LH']:.2f}")
        print(f"  Tank:                          {sol['W_tank']:.2f}")
        print(f"  Fuel cell:                     {sol['W_fuel_cell']:.2f}")
        print(f"  Miscellaneous:                 {sol['W_misc']:.2f}")
        print(f"  Total:                         {sol['W_total']:.2f}\n")

        print("MASSES (kg):")
        print(f"  Payload:                       {sol['M_payload']:.2f}")
        print(f"  Wing total:                    {sol['M_wing']:.2f}")
        print(f"    └ Spar:                      {sol['M_spar']:.2f}")
        print(f"    └ Skin:                      {sol['M_skin']:.2f}")
        print(f"  LH2 fuel:                      {sol['M_LH']:.2f}")
        print(f"  Tank:                          {sol['M_tank']:.2f}")
        print(f"  Fuel cell:                     {sol['M_fuel_cell']:.2f}")
        print(f"  Miscellaneous:                 {sol['M_misc']:.2f}")
        print(f"  Total:                         {sol['M_total']:.2f}")

    def _get_solution_plane(self):
        sol = self.solution
        b_sol = sol['b']
        cords_sol = sol['cords']
        wing_airfoil = self.wing_airfoil
        y_sections_sol = np.linspace(0, b_sol / 2, len(cords_sol))
        wing_sol = Mass_wing(
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[-0.25 * cords_sol[i], y_sections_sol[i], 0],
                    chord=cords_sol[i],
                    airfoil=wing_airfoil
                )
                for i in range(len(cords_sol))
            ]
        )

        self.airplane_sol = asb.Airplane(wings=[wing_sol])
        
    def draw(self):
        self.airplane_sol.draw()
        
    def plot_aero(self):
        airplane = self.airplane_sol
        
        alpha = np.linspace(-20, 20, 1000)
        aero = asb.AeroBuildup(
            airplane=airplane,
            op_point=asb.OperatingPoint(velocity=30, alpha=alpha, beta=0),
        ).run()

        plt.plot(alpha, aero["CL"])
        plt.show()
    
if __name__ == "__main__":
    LH_plane = Tandem_LH()
    LH_plane.solve()
    LH_plane.print_solution()
    LH_plane.draw()
    #LH_plane.plot_aero()