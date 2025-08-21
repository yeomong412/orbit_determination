import numpy as np
from constants import conv_eq2ec_mat

def coordinate_transform_ecliptic(RA,DEC,distance):
    RA_rad = np.deg2rad(RA)
    DEC_rad = np.deg2rad(DEC)

    x = distance * np.cos(DEC_rad) * np.cos(RA_rad)
    y = distance * np.cos(DEC_rad) * np.sin(RA_rad)
    z = distance * np.sin(DEC_rad)

    return np.dot(conv_eq2ec_mat, [x,y,z])
