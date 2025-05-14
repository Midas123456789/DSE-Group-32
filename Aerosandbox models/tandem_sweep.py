from full_model_hydrogen_tandem import *
W_pl = 1000
P_pl = 10e3
for altitude in np.linspace(12000, 18000, 6):
    print(1)
    plane = get_LH_conventional(altitude=altitude, P_payload=P_pl, W_payload=W_pl)
    