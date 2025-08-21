import math
import numpy as np
from constants import *         # 상수값
from numerical_methods import * # 수치해석 함수들

# -----------------------------입력 예시 -------------------------------
# r1, r2 = np.array([-1.14933075,  1.18727464,  0.05321403]), np.array([-1.40015014,  0.90012888,  0.05335051])
# delta_t = 30.0 * 86400   #30일 (sec)
# ----------------------------------------------------------------------

def solve_lambert(r1, r2, delta_t):

    #시간단위 변환 (일 -> 초)
    delta_t *= JD_to_sec 

    #거리 단위 변환 (au -> km)
    r1, r2 = r1 * au_to_km, r2 * au_to_km

    # 델타세타 계산
    delta_theta = np.arccos(np.dot(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2)))
    if delta_theta > np.pi:
        delta_theta = 2 * np.pi - delta_theta

    # 상수 A 정의
    A = np.sin(delta_theta) * np.sqrt((np.linalg.norm(r1) * np.linalg.norm(r2)) / (1 - np.cos(delta_theta)))

    # z값 계산하기. 초기 z값은 적당히 유추하여 설정. 화성의 경우 4.7 ~ 5.6 을 초기값으로 잡음.
    z = [0] * 1
    z_count = 1
    # 뉴턴-랩슨 
    print('--------------------------')
    while(True):
        z.append(z[z_count-1] - (lambert_F(z[z_count-1], r1, r2, A, delta_t)/lambert_F_diff(z[z_count-1], r1, r2, A)))
        if np.abs(z[z_count] - z[z_count-1]) < converge_factor: break
        lambert_y(z[z_count-1],r1,r2,A)
        print(f"갱신된 z = {z[z_count]}") # 갱신된 z 값 출력
        if math.isnan(z[-1]) :            # z값이 nan인지 확인하고, nan이라면 False 반환.
            print(f"------------실패-----------")
            print(f"다른 z값 시도: {z[0] + 0.5}")
            z = [z[0] + 0.01]             # 0.01씩 초기 z값 조정
            z_count = 0
        z_count += 1
    print('--------------------------')

    z_val = z[-1]
    # 라그랑지 계수 계산
    f = lagrange_f(z_val, r1, r2, A)
    g = lagrange_g(z_val, r1, r2, A)
    f_diff = lagrange_f_diff(z_val, r1, r2, A)
    g_diff = lagrange_g_diff(z_val, r1, r2, A)

    #v1,v2 계산
    v1 = (1/g) * (r2 - f * r1)
    v2 = (1/g) * (g_diff * r2 - r1)

    return [v1,v2]