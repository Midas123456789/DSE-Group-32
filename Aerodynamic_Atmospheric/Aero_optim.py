import numpy as np

def optimize_velocity_for_min_drag(aircraft, v_range):
    min_drag = float('inf')
    optimal_v = None
    for V in v_range:
        D = aircraft.Drag(V=V)
        if D < min_drag:
            min_drag = D
            optimal_v = V
    return optimal_v, min_drag

def optimize_altitude_for_min_drag(aircraft, h_range):
    min_drag = float('inf')
    optimal_h = None
    for h in h_range:
        aircraft.altitude = h
        V = aircraft.min_velocity()
        aircraft.V = V
        D = aircraft.Drag()
        if D < min_drag:
            min_drag = D
            optimal_h = h
        print(aircraft.altitude, aircraft.rho, aircraft.V, D)

    return optimal_h, min_drag

def determine_states_altitude(aircraft, h_range):
    """For range of altitudes, determine the velocity and drag"""
    results = []
    for h in h_range:
        aircraft.altitude = h
        rho = aircraft.rho
        V = aircraft.min_velocity()
        aircraft.V = V
        # print(aircraft.altitude, aircraft.rho, aircraft.V, aircraft.a, aircraft.M)
        print('h', h, 'rho', rho, 'V', V, 'a', aircraft.a, 'M', aircraft.M)
        D = aircraft.Drag()
        L = aircraft.Lift()
        M = aircraft.M
        results.append([h, rho, V, D, L, M])
    return results


def parameter_sweep(aircraft, CL_list, V_range, S_range):
    results = []
    for CL in CL_list:
        for V in V_range:
            for S in S_range:
                aircraft.CL = CL
                aircraft.V = V
                aircraft.S = S
                D = aircraft.Drag()
                L = aircraft.Lift()
                results.append((CL, V, S, D, L))
    return results



