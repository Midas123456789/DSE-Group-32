import numpy as np
import math

from numpy.ma.core import arcsin


def dd_to_dms(dd):
    """Converts decimal degrees to degrees, minutes, seconds."""
    degrees = int(dd)
    minutes_decimal = (dd - degrees) * 60
    minutes = int(minutes_decimal)
    seconds = (minutes_decimal - minutes) * 60
    return degrees, minutes, seconds

def dms_to_dd(degrees, minutes, seconds):
    """Converts degrees, minutes, seconds to decimal degrees."""
    dd = float(degrees) + float(minutes) / 60 + float(seconds) / 3600
    return dd

class Vehicle:
    def __init__(self, windproperties, radius, position_dms):
        """
        Initializes the vehicle with wind properties, radius, and position in DMS.
        position_dms should be a list: [longitude_dms, latitude_dms]
        where each _dms is a tuple: (degrees, minutes, seconds).
        """
        self.windproperties = windproperties  # [u_wind (m/s), v_wind (m/s)]
        self.maxvelocity = 500  # m/s
        self.desiredvelocity = 200  # m/s
        self.lateralconvert = 30715 #m
        self.longitudinalconvert = 110600 #m


        self.radius = radius
        self.position_dms = position_dms
        self.position_dd = np.array([dms_to_dd(*pos) for pos in position_dms]) # Convert to DD for calculations

    def calculate_distance(self, target_position_dd):
        """Calculates the distance in decimal degrees."""
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

    def trajectory(self, target_position_dms):
        """
        Determines the vehicle's movement towards the target position in DMS,
        considering wind influence and velocity constraints.
        """
        target_position_dd = np.array([dms_to_dd(*pos) for pos in target_position_dms])
        delta_lon_dd, delta_lat_dd = self.calculate_distance(target_position_dd)
        distance_dd = math.sqrt((delta_lon_dd**2 + (delta_lat_dd**2)))

        if distance_dd < 100:  # Small tolerance for reaching target in DD
            self.position_dms = target_position_dms
            self.position_dd = target_position_dd
            return self.position_dms

        bearing_to_target_rad = math.radians(self.calculate_bearing(target_position_dd))

        wind_u, wind_v = self.windproperties
        #wind_speed = math.sqrt(wind_u**2 + wind_v**2)

        # Calculate required ground velocity components in decimal degrees per second
        time_to_target = distance_dd / self.desiredvelocity  # Approximate time based on DD distance

        u_flightspeed = self.desiredvelocity*math.sin(bearing_to_target_rad) + wind_u
        v_flightspeed = self.desiredvelocity*math.cos(bearing_to_target_rad) + wind_v
        v_flightspeed = math.sqrt(u_flightspeed**2 + v_flightspeed**2)

        if v_flightspeed <= self.maxvelocity:
            # We can achieve the desired ground velocity

            self.position_dd[0] += v_flightspeed/self.longitudinalconvert
            self.position_dd[1] += u_flightspeed/self.lateralconvert
            self.position_dms = [dd_to_dms(self.position_dd[0]), dd_to_dms(self.position_dd[1])]

        else:
            # Fly at max airspeed towards the target
            max_airspeed_u = self.maxvelocity * math.sin(bearing_to_target_rad)
            max_airspeed_v = self.maxvelocity * math.cos(bearing_to_target_rad)

            u_flightspeed = max_airspeed_u + wind_u
            v_flightspeed = max_airspeed_v + wind_v

            self.position_dd[0] += v_flightspeed / self.longitudinalconvert
            self.position_dd[1] += u_flightspeed / self.lateralconvert
            self.position_dms = [dd_to_dms(self.position_dd[0]), dd_to_dms(self.position_dd[1])]

        return self.position_dms

# Example usage with DMS:
wind = [55, -10]  # Wind blowing at 5 m/s eastward
start_position_dms = [(4, 21, 0.0), (52, 0, 36.0)]  # Delft: 4°21'00.0" E, 52°00'36.0" N
target_position_dms = [(4, 54, 0.0), (52, 22, 12.0)] # Amsterdam: 4°54'00.0" E, 52°22'12.0" N
vehicle_dms = Vehicle(wind, 10, start_position_dms)

for _ in range(1000):
    current_position_dms = vehicle_dms.trajectory(target_position_dms)
    lon_dms, lat_dms = current_position_dms
    print(f"Current Position (DMS): {lon_dms[0]}°{lon_dms[1]:02.0f}'{lon_dms[2]:04.1f}\" E, {lat_dms[0]}°{lat_dms[1]:02.0f}'{lat_dms[2]:04.1f}\" N")
    if np.linalg.norm(vehicle_dms.position_dd - np.array([dms_to_dd(*pos) for pos in target_position_dms])) < 0.0001:
        print("Reached target (DMS)!")
        break