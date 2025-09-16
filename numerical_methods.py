import numpy as np
from constants import *

#이곳에 수치해석에 필요한 함수들 저장

def stumpff_C(z):
    if np.abs(z) < 1e-8:  # z가 0에 가까운 경우
        return float(1/2)
    elif z > 0:
        return float((1 - np.cos(np.sqrt(z))) / z)
    else:  # z < 0
        return float((1 - np.cosh(np.sqrt(-z))) / z)

def stumpff_C_diff(z):
    if np.abs(z) < 1e-8: return float(-1/24)
    else: return float((1 / (2 * z)) * (1 - z * stumpff_S(z) - 2 * stumpff_C(z)))

def stumpff_S(z):
    if np.abs(z) < 1e-8:  # z가 0에 가까운 경우
        return float(1/6)
    elif z > 0:
        return float((np.sqrt(z) - np.sin(np.sqrt(z))) / (np.power(np.sqrt(z), 3)))
    else:  # z < 0
        return float((np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.power(np.sqrt(-z), 3)))
    
def stumpff_S_diff(z):
    if np.abs(z) < 1e-8: return float(-1/80)
    else: return float((1 / (2 * z)) * (stumpff_C(z) - 3 * stumpff_S(z)))

def lambert_y(z,r1,r2,A):
    return float(np.linalg.norm(r1) + np.linalg.norm(r2) + A * (z * stumpff_S(z) - 1) / np.sqrt(stumpff_C(z)))

def lambert_y_diff(z, r1, r2, A):
    C = stumpff_C(z)
    S = stumpff_S(z)
    Cp = stumpff_C_diff(z)
    Sp = stumpff_S_diff(z)
    sqrtC = np.sqrt(C)
    f = z * S - 1.0
    numerator = (S + z * Sp) * sqrtC - f * 0.5 * (Cp / sqrtC)
    return float(A * numerator / C)

def lambert_F(z,r1,r2,A,delta_t):
    return np.float64(np.power(float( lambert_y(z,r1,r2,A) / stumpff_C(z) ), 1.5) * stumpff_S(z) + A * np.sqrt(lambert_y(z,r1,r2,A)) - np.sqrt(mu) * delta_t)

def lambert_F_diff(z,r1,r2,A):
    y = lambert_y(z, r1, r2, A)
    C = stumpff_C(z)
    S = stumpff_S(z)
    Cp = stumpff_C_diff(z)
    Sp = stumpff_S_diff(z)
    ydz = lambert_y_diff(z, r1, r2, A)

    term1 = (3/2) * np.power((y / C), 0.5) * ((ydz * C - Cp * y) / (C**2)) * S
    term2 = np.power((y / C), 1.5) * Sp
    term3 = 0.5 * A * (ydz / np.sqrt(y))
    return term1 + term2 + term3

def lagrange_f(z, r1, r2, A):
    return 1 - (lambert_y(z,r1,r2,A) / np.linalg.norm(r1))

def lagrange_g(z,r1,r2, A):
    return A * np.sqrt(lambert_y(z,r1,r2,A) / mu)

def lagrange_f_diff(z,r1,r2,A):
    return (np.sqrt(mu) / (np.linalg.norm(r1) * np.linalg.norm(r2))) * np.sqrt(lambert_y(z,r1,r2,A) / stumpff_C(z)) * (z * stumpff_S(z) - 1)

def lagrange_g_diff(z,r1,r2,A):
    return 1 - (lambert_y(z,r1,r2,A) / np.linalg.norm(r2))
