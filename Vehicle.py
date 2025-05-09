import numpy as np
import math
import pandas as pd


class Vehicle:
    def __init__(self, windproperties, radius, position_dd):
        """
        Initializes the vehicle with wind properties, radius, and position in DD.
        position_dd should be a list: [longitude_dd, latitude_dd]
        """
        self.windproperties = windproperties  # [u_wind (m/s), v_wind (m/s)]
        self.maxvelocity = 500  # m/s
        self.desiredvelocity = 200  # m/s
        self.lateralconvert = 30715 #m
        self.longitudinalconvert = 110600 #m


        self.radius = radius
        self.position_dd =  position_dd

    def calculate_distance(self, target_position_dd):
        """Calculates the distance in meters."""
        delta_lon = (target_position_dd[0] - self.position_dd[0])*self.longitudinalconvert
        delta_lat = (target_position_dd[1] - self.position_dd[1])*self.lateralconvert
        return delta_lon, delta_lat

    def calculate_bearing(self, target_position_dd):
        """Calculates the initial bearing to the target in decimal degrees."""
        delta_lon_rad = math.radians(target_position_dd[0] - self.position_dd[0])
        lat1_rad = math.radians(self.position_dd[1])
        lat2_rad = math.radians(target_position_dd[1])

        y = math.sin(delta_lon_rad) * math.cos(lat2_rad)
        x = math.cos(lat1_rad) * math.sin(lat2_rad) - math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(delta_lon_rad)
        bearing_rad = math.atan2(y, x)
        return math.degrees(bearing_rad + 360) % 360

    def trajectory(self, target_position_dd):
        """
        Determines the vehicle's movement towards the target position in DD,
        considering wind influence and velocity constraints.
        """
        #target_position_dd = np.array([dms_to_dd(*pos) for pos in target_position_dms]) #target position is already in dd
        delta_lon_m, delta_lat_m = self.calculate_distance(target_position_dd)
        distance_m = math.sqrt((delta_lon_m**2 + (delta_lat_m**2)))

        if distance_m < 100:  # Small tolerance for reaching target in meters
            self.position_dd = target_position_dd
            return self.position_dd

        bearing_to_target_rad = math.radians(self.calculate_bearing(target_position_dd))

        wind_u, wind_v = self.windproperties
        #wind_speed = math.sqrt(wind_u**2 + wind_v**2)

        # Calculate required ground velocity components in decimal degrees per second
        time_to_target = distance_m / self.desiredvelocity  # Approximate time

        u_flightspeed = -self.desiredvelocity*math.cos(bearing_to_target_rad) - wind_u
        v_flightspeed = self.desiredvelocity*math.sin(bearing_to_target_rad) - wind_v
        self.ground_speed = math.sqrt(u_flightspeed**2 + v_flightspeed**2) # Calculate ground speed

        if self.ground_speed <= self.maxvelocity:
            # We can achieve the desired ground velocity

            self.position_dd[0] += (u_flightspeed / self.longitudinalconvert)
            self.position_dd[1] += (v_flightspeed / self.lateralconvert)

        else:
            # Fly at max airspeed towards the target
            max_airspeed_u = self.maxvelocity * math.sin(bearing_to_target_rad)
            max_airspeed_v = self.maxvelocity * math.cos(bearing_to_target_rad)

            u_flightspeed = max_airspeed_u + wind_u
            v_flightspeed = max_airspeed_v + wind_v
            self.ground_speed = math.sqrt(u_flightspeed ** 2 + v_flightspeed ** 2)

            self.position_dd[0] += (u_flightspeed / self.longitudinalconvert)
            self.position_dd[1] += (v_flightspeed / self.lateralconvert)
        return  self.position_dd, self.ground_speed

    def vehicle_data(self, identifier, target_position):

        df = pd.DataFrame()
        df['identifier'] = [identifier]
        df['C. longitudinal'] = [self.position_dd[0]]
        df['C. lateral'] = [self.position_dd[1]]
        df['velocity'] = [self.ground_speed]
        df['T. longitudional'] = [target_position[0]]
        df['T. lateral'] = [target_position[1]]

        return df.set_index('identifier')




# Example usage with DD:
wind = [20, -13]  # Wind blowing at 5 m/s eastward
start_position_dd = [4.350, 52.01]  # Delft: 4.350, 52.01
target_position_dd = [4.9, 52.37] # Amsterdam:  4.9, 52.37
vehicle_dd = Vehicle(wind, 10, start_position_dd)

for _ in range(1000):
    current_position_dd = vehicle_dd.trajectory(target_position_dd)
    lon_dd, lat_dd = current_position_dd
    print(vehicle_dd.vehicle_data("01", target_position_dd))
    #print(f"Current Position (DD): {lon_dd:.6f} E, {lat_dd:.6f} N")
    if np.linalg.norm(vehicle_dd.position_dd - np.array(target_position_dd)) < 0.001:
        print("Reached target (DMS)!")
        break
