import numpy as np

Free_space_of_data_lists = 200
Num_observation = 5
converge_factor = 0.2e-8
JD_to_sec = 86400 # length of a day in sec
mu = np.float64(132712.0e+6) # sun's GM. (km^3/m^2)
#mu = np.float64(0.3986e+6) # earth's GM. (km^3/m^2)
ecliptic_obluquity_deg = 23.4392911
ecliptic_obluquity_rad = np.deg2rad(ecliptic_obluquity_deg)
au_to_km = np.float64(149597870.700)

conv_eq2ec_mat = np.array([[1,                             0,                               0],
                          [0, np.cos(ecliptic_obluquity_rad), np.sin(ecliptic_obluquity_rad)],
                          [0, -np.sin(ecliptic_obluquity_rad), np.cos(ecliptic_obluquity_rad)]])

conv_ec2eq_mat = np.array([[1,                              0,                              0],
                          [0, np.cos(ecliptic_obluquity_rad), -np.sin(ecliptic_obluquity_rad)],
                          [0, np.sin(ecliptic_obluquity_rad), np.cos(ecliptic_obluquity_rad)]])
