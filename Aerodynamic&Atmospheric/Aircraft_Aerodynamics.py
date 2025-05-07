class AircraftAerodynamic(Aerodynamic):
    def plot_feasible_S_V(self, CL):
        """
        Plot feasible wing area (S) vs. velocity (V) for different coefficients of lift.

        Parameters:
        - CL: float or list of floats (lift coefficients)
        """
        rho = self.altitude_data[self.altitude]["Density [kg/m³]"]
        V = np.linspace(10, 200, 500)

        if not isinstance(CL, list):
            CL = [CL]

        plt.figure(figsize=(10, 6))
        for cl in CL:
            S = (2 * self.weight) / (rho * V**2 * cl)
            plt.plot(V, S, label=f'CL = {cl}')
            plt.fill_between(V, 0, S, alpha=0.1)

        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Wing Area (m²)')
        plt.title(f'Feasible Wing Area (S) vs Velocity (V)\nAltitude: {self.altitude} m, Weight: {self.weight} N')
        plt.xlim(0, 100)
        plt.grid(True)
        plt.legend()
        plt.show()
        
    def estimate_drag(self, S, V, CD0, k):
        """
        Estimate drag using the drag equation.

        Parameters:
        - S: wing area in m²
        - V: velocity in m/s
        - CD0: zero-lift drag coefficient
        - k: induced drag factor

        Returns:
        - Drag force in Newtons
        """
        rho = self.rho
        q = 0.5 * rho * V**2
        CD = CD0 + k * (self.weight / (0.5 * rho * V**2 * S))**2









# Aircraft example
aircraft = AircraftAerodynamic(weight=4000, altitude=20000)
aircraft.plot_feasible_S_V([1, 1.5, 2, 2.5, 3])