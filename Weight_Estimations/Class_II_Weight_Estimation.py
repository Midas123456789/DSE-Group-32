class ClassIIWeightEstimation:
    """
    Class II Weight Estimation using component breakdown and iterative MTOW convergence.
    """
    def __init__(self, 
                payload_weight_kg, 
                battery_mass_kg,
                initial_mtow_guess_kg, 
                iteration_limit=100, 
                tolerance=0.1,
                g=9.80665,
                S=30,         # Wing area (m^2)
                A=25,         # Aspect ratio
                fuselage_length=10, 
                fuselage_diameter=1.5,
                W1_WTO=1, W2_W1=1, W3_W2=1, W4_W3=1, W5_W4=1, W6_W5=1, W7_W6=1, W8_W7=1, Wfinal_W8=1):

        self.payload_weight_kg = payload_weight_kg
        self.battery_mass_kg = battery_mass_kg
        self.MTOW_kg = initial_mtow_guess_kg
        self.iteration_limit = iteration_limit
        self.tolerance = tolerance
        self.g = g
        self.S = S
        self.A = A
        self.fuselage_length = fuselage_length
        self.fuselage_diameter = fuselage_diameter
        self.W1_WTO=W1_WTO
        self.W2_W1=W2_W1
        self.W3_W2=W3_W2
        self.W4_W3=W4_W3
        self.W5_W4=W5_W4
        self.W6_W5=W6_W5
        self.W7_W6=W7_W6
        self.W8_W7=W8_W7
        self.Wfinal_W8=Wfinal_W8
        self.results = {}

        self.converged = False
        self.Determine_MTOW()

    def Estimate_Component_Weights(self, MTOW):
        W_wing = 0.036 * MTOW**0.758 * self.S**0.6 * self.A**0.04  # [Raymer]
        W_fuselage = 0.052 * MTOW**0.689 * self.fuselage_length**0.924
        W_empennage = 0.019 * MTOW**0.889
        W_landing_gear = 0.04 * MTOW
        W_avionics = 0.01 * MTOW
        W_propulsion = 0.04 * MTOW  # motors, wiring, inverters

        total_structure_weight = W_wing + W_fuselage + W_empennage + W_landing_gear + W_avionics + W_propulsion
        return {
            "Wing": W_wing,
            "Fuselage": W_fuselage,
            "Empennage": W_empennage,
            "Landing Gear": W_landing_gear,
            "Avionics": W_avionics,
            "Propulsion": W_propulsion,
            "Total Structure": total_structure_weight
        }

    def Determine_MTOW(self):
        for _ in range(self.iteration_limit):
            components = self.Estimate_Component_Weights(self.MTOW_kg)
            OEW = components["Total Structure"]
            
            mission_fraction = (self.W1_WTO * self.W2_W1 * self.W3_W2 * self.W4_W3 * self.W5_W4 * self.W6_W5 * self.W7_W6 * self.W8_W7 * self.Wfinal_W8)
            fuel_fraction = 1 - mission_fraction
            fuel_weight = fuel_fraction * self.MTOW_kg
            new_MTOW = OEW + self.payload_weight_kg + self.battery_mass_kg + fuel_weight

            if abs(new_MTOW - self.MTOW_kg) < self.tolerance:
                self.converged = True
                break
            self.MTOW_kg = new_MTOW

        if self.converged:
            self.results["MTOW [kg]"] = round(self.MTOW_kg, 2)
            self.results["OEW [kg]"] = round(OEW, 2)
            self.results["Battery [kg]"] = round(self.battery_mass_kg, 2)
            self.results["Fuel [kg]"] = round(fuel_weight, 2)
            self.results["Payload [kg]"] = round(self.payload_weight_kg, 2)
            self.results.update({f"{k} [kg]": round(v, 2) for k, v in components.items()})
            
            self.estimated_MTOM = self.results["MTOW [kg]"]
            self.estimated_OEM = self.results["OEW [kg]"]
            self.estimated_OEM_fraction = self.estimated_OEM / self.estimated_MTOM
            self.estimated_battery_mass_kg = self.results["Battery [kg]"]
            self.estimated_fuel_weight_kg = self.results["Fuel [kg]"]

    def __str__(self):
        if not self.converged:
            return "MTOW did not converge."
        
        output = ["Class II Weight Estimation Results:\n"]
        max_key_length = max(len(k) for k in self.results)
        header = f"{'Component'.ljust(max_key_length)} | Value"
        output.append(header)
        output.append("-" * len(header))
        
        for key, val in self.results.items():
            output.append(f"{key.ljust(max_key_length)} | {val:>10,.2f}")
        return "\n".join(output)


# Example usage
if __name__ == "__main__":
    
    battery_mass = 300
    payload_mass = 100

    Class_II = ClassIIWeightEstimation(
        payload_weight_kg=payload_mass,
        battery_mass_kg=battery_mass,
        S=30,
        A=25,
        fuselage_length=10,
        fuselage_diameter=1.5,
        initial_mtow_guess_kg= 2000
    )

    print(Class_II)
