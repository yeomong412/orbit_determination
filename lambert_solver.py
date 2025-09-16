import math
import numpy as np
from constants import *         # 상수값
from numerical_methods import * # 수치해석 함수들


def solve_lambert(r1, r2, delta_t):
    # 시간 단위 변환 (JD -> sec)
    delta_t *= JD_to_sec

    # 거리 단위 변환 (AU -> km)
    r1 = r1 * au_to_km
    r2 = r2 * au_to_km

    # 각도 및 A 계산
    cos_dth = np.dot(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2))
    cos_dth = np.clip(cos_dth, -1.0, 1.0)
    delta_theta = float(np.arccos(cos_dth))
    if delta_theta > np.pi:
        delta_theta = 2 * np.pi - delta_theta

    A = np.sin(delta_theta) * np.sqrt((np.linalg.norm(r1) * np.linalg.norm(r2)) / (1 - np.cos(delta_theta)))

    # 안전한 뉴턴법 설정
    z_cur = 0.0
    Z_MAX = 20.0
    MAX_IT = 60
    TOLF = 1e-8

    print('--------------------------')
    for _ in range(MAX_IT):
        try:
            F = lambert_F(z_cur, r1, r2, A, delta_t)
            Fp = lambert_F_diff(z_cur, r1, r2, A)
        except Exception:
            F, Fp = np.nan, np.nan

        if (not np.isfinite(F)) or (not np.isfinite(Fp)) or abs(Fp) < 1e-14:
            step = -np.sign(F) * 0.1 if np.isfinite(F) and F != 0 else 0.0
            z_try = float(np.clip(z_cur + step, -Z_MAX, Z_MAX))
        else:
            dz = -F / Fp
            alpha = 1.0
            z_try = float(np.clip(z_cur + alpha * dz, -Z_MAX, Z_MAX))
            while alpha > 1e-4:
                F_try = lambert_F(z_try, r1, r2, A, delta_t)
                if np.isfinite(F_try) and abs(F_try) < abs(F):
                    break
                alpha *= 0.5
                z_try = float(np.clip(z_cur + alpha * dz, -Z_MAX, Z_MAX))

        print(f"갱신된 z = {z_try}")
        if (abs(z_try - z_cur) < converge_factor) or (np.isfinite(F) and abs(F) < TOLF):
            z_cur = z_try
            print('성공!')
            break
        z_cur = z_try
    else:
        print('------------실패-----------')

    print('--------------------------')

    z_val = z_cur
    f = lagrange_f(z_val, r1, r2, A)
    g = lagrange_g(z_val, r1, r2, A)
    f_diff = lagrange_f_diff(z_val, r1, r2, A)
    g_diff = lagrange_g_diff(z_val, r1, r2, A)

    v1 = (1 / g) * (r2 - f * r1)
    v2 = (1 / g) * (g_diff * r2 - r1)

    return [v1, v2]

