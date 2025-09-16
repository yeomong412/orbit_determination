import sys
import numpy as np
from scipy import stats
from constants import *                                         # 상수값
from coordinate_transform import coordinate_transform_ecliptic  # 적도좌표계 -> 황도좌표계 변환
from getcoordinate_earth import getcoordinate_earth             # NASA JPL Horizons에서 지구의 위치좌표 얻어오기 
from lambert_solver import solve_lambert                        # 램버트 문제의 해를 구함
asteriod_name = '2000AC6'

# ----------------------------------------------------------------------------------------
# 기준 평면은 지구 공전면이고, 황도 좌표계를 이용함.
# 관측 데이터 파일(input.txt) 에서 관측 정보를 읽어옴.
# 순서 및 단위: 적경(°), 적위(°), 거리(AU), 시간(JD)
# 참고: nasa_raw_input에 "https://ssd.jpl.nasa.gov/horizons/app.html#/" 에서 긁어온 데이터를 넣고
#       input_format.py를 실행하면 input.txt파일에 자동으로 관측정보가 입력됨.
# ----------------------------------------------------------------------------------------
line = [0] * Free_space_of_data_lists # 파일에서 읽어온 한 줄을 저장
RA, DEC, distance = [0] * Free_space_of_data_lists, [0] * Free_space_of_data_lists, [0] * Free_space_of_data_lists #위치 벡터
time = [0] * Free_space_of_data_lists # 시간 저장

with open('data/input.txt', 'r',encoding="utf-8") as input_file:
    for i in range(Num_observation):
        print(f"읽는 중...{i+1}")
        line[i] = input_file.readline() 
        # 공백 기준 분리(split) 이후 np.float64로 각각 저장.
        RA[i], DEC[i], distance[i], time[i] = map(np.float64, line[i].split())
# 시간은 getcoordinate_earth()의 입력 형식에 맞춰 문자열로 저장
time = list(map(str, time))


# ----------------------------------------------------------------------------------------
# 태양계 질량중심 좌표계로 천체의 위치벡터를 계산
# CM 기준 천체 벡터 = 지구 기준 천체 벡터 + 태양 기준 지구 벡터
# ----------------------------------------------------------------------------------------
coord_of_earth = list(map(getcoordinate_earth, time[:Num_observation])) # 관측 개수에 맞춰서, 각 시간의 지구 위치벡터를 저장
coord_of_body_from_earth = list( #RA, DEC, distance 데이터를 입력 -> 카르테시안(x,y,z)으로 변환
    map(coordinate_transform_ecliptic, RA[:Num_observation], DEC[:Num_observation], distance[:Num_observation]) 
)
coord_of_body_form_CM = [earth + body for earth, body in zip(coord_of_earth, coord_of_body_from_earth)]
# 계산 결과 출력
print(f"지구의 좌표(x[AU], y[AU], z[AU]):\n {coord_of_earth[:Num_observation]}")
print(f"지구 기준 천체 좌표(x[AU], y[AU], z[AU]):\n {coord_of_body_from_earth[:Num_observation]}")
print(f"CM 기준 천체 좌표(x[AU], y[AU], z[AU]):\n {coord_of_body_form_CM[:Num_observation]}")

# ----------------------------------------------------------------------------------------
# 두 관측 시점의 위치 벡터를 이용해 각 지점에서 천체 속도 벡터 계산
# velocity_of_body 리스트에 각 위치별 속도 벡터 저장
# ----------------------------------------------------------------------------------------
velocity_of_body = []
for i in range(Num_observation - 1):
    # vel = solve_lambert(coord_of_body_form_CM[i], coord_of_body_form_CM[i+1], float(time[i+1]) - float(time[i]))
    vel = solve_lambert(coord_of_body_form_CM[i], coord_of_body_form_CM[i+1], float(time[i+1]) - float(time[i]))
    if vel == False : sys.exit("초기 z값 조정 필요")
    print(
        f"천체의 속도벡터 {i+1}번:\n==== Date{i+1}: [{vel[0][0]} km/s {vel[0][1]} km/s {vel[0][2]} km/s] => {np.linalg.norm(vel[0])} km/s \n"
        f"==== Date{i+2}: [{vel[1][0]} km/s {vel[1][1]} km/s {vel[1][2]} km/s] => {np.linalg.norm(vel[1])} km/s"
    )
    velocity_of_body.append(vel)
    print('.')
print('--------------------------')

# ----------------------------------------------------------------------------------------
# 각 관측별로 궤도요소 6개 계산 후 저장.
# 순서대로 장반경(a), 이심률(e), 궤도경사(i), 승교점 경도(RAAN), 근일점 편각(omega), 진근점 이각(theta), 시간(JD).
#-----------------------------------------------------------------------------------------
orbital_elements = [[0] * 7 for _ in range(2 * (Num_observation - 1))]
for i in range(Num_observation - 1):
    for j in range(2):
        idx = i * 2 + j
        # r_vec = coord_of_body_form_CM[i+j] * au_to_km  # 위치 벡터
        r_vec = coord_of_body_form_CM[i+j] * au_to_km  # 위치 벡터
        r = np.linalg.norm(r_vec)  # 거리
        v_vec = velocity_of_body[i][j]  # 속도 벡터
        v = np.linalg.norm(v_vec)  # 속력
        epsilon = v * v / 2 - mu / r  # 총 역학적 에너지
        h = np.cross(r_vec, v_vec)  # 각운동량 벡터
        N = np.cross([0, 0, 1], h)  # 승교점 벡터

        # orbital_elements[i+j][0] = -mu / (2 * epsilon)  # 장반경
        # orbital_elements[i+j][1] = np.linalg.norm(np.cross(v_vec, h) / mu - r_vec / r)  # 이심률
        # orbital_elements[i+j][2] = np.arccos(np.dot(h, [0, 0, 1]) / np.linalg.norm(h))  # 궤도 경사각
        # orbital_elements[i+j][3] = np.arccos(np.dot(N, [1, 0, 0]) / np.linalg.norm(N))  # 승교점 경도
        # orbital_elements[i+j][4] = np.arccos(np.dot(N, np.cross(v_vec, h) / mu - r_vec / r) / (np.linalg.norm(N) * orbital_elements[i][1]))  # 근일점 편각

        # # ------------------------ 천체 위치 정보 계산--------------------------------
        # theta = np.arccos(np.dot(N, r_vec) / (np.linalg.norm(N) * r))
        # if np.dot(r_vec, v_vec) < 0:
        #     theta = 2 * np.pi - theta
        # orbital_elements[i+j][5] = theta  # 진근점 이각
        # orbital_elements[i+j][6] = float(time[i+j])  # 시간 #이전 코드

        # --- 공통 벡터 ---
        e_vec = np.cross(v_vec, h) / mu - r_vec / r
        e = np.linalg.norm(e_vec)

        # 1) 장반경 a, 이심률 e, 궤도경사 i
        orbital_elements[idx][0] = -mu / (2 * epsilon)
        orbital_elements[idx][1] = e
        orbital_elements[idx][2] = np.arccos(h[2] / np.linalg.norm(h))

        # 2) 노드벡터 N 및 RAAN Ω (atan2)
        N = np.cross([0.0, 0.0, 1.0], h)
        N_norm = np.linalg.norm(N)
        if N_norm < 1e-12:
            Omega = 0.0  # 적도궤도 특이치
        else:
            Omega = np.mod(np.arctan2(N[1], N[0]), 2*np.pi)
        orbital_elements[idx][3] = Omega

        # 3) 근일점 편각 ω (atan2, 사분면 안전)
        h_norm = np.linalg.norm(h)
        if e < 1e-12 or N_norm < 1e-12:
            omega = 0.0  # 원궤도/적도궤도 특이치
        else:
            cos_w = np.dot(N, e_vec) / (N_norm * e)
            sin_w = np.dot(np.cross(N, e_vec), h) / (N_norm * e * h_norm)
            omega = np.mod(np.arctan2(sin_w, cos_w), 2*np.pi)
        orbital_elements[idx][4] = omega

        # 4) 진근점이각 ν (atan2, 사분면 안전)
        if e < 1e-12:
            # 원궤도: ν 대신 u를 사용(참고용)
            cos_u = np.dot(N, r_vec) / (max(N_norm,1e-12) * r)
            sin_u = np.dot(np.cross(N, r_vec), h) / (max(N_norm,1e-12) * r * h_norm)
            nu = np.mod(np.arctan2(sin_u, cos_u), 2*np.pi)
        else:
            cos_nu = np.dot(e_vec, r_vec) / (e * r)
            sin_nu = np.dot(np.cross(e_vec, r_vec), h) / (e * r * h_norm)
            nu = np.mod(np.arctan2(sin_nu, cos_nu), 2*np.pi)
        orbital_elements[idx][5] = nu

        # 5) 시간
        orbital_elements[idx][6] = float(time[i+j])


        # 계산된 궤도 요소 출력
        print(
            f"계산된 {idx+1}차 궤도요소 : [{orbital_elements[idx][0]} km, {orbital_elements[idx][1]} # "
            f"{orbital_elements[idx][2]} rad, {orbital_elements[idx][3]} rad, {orbital_elements[idx][4]} rad, "
            f"{orbital_elements[idx][5]} rad, {orbital_elements[idx][6]} JD]"
        )

print('--------------------------')

# ---------------------------------------------------------
# 전체 관측에 대한 궤도요소의 최종 평균값 계산

# 1. 장반경 단위: km, 궤도 요소 각도: rad
# 2. 장반경 단위: km, 궤도 요소 각도: deg
# 3. 장반경 단위: au, 궤도 요소 각도: deg
# ---------------------------------------------------------
oe_km_rad = stats.trim_mean(orbital_elements, 0.1)
oe_km_deg = [oe_km_rad[0], oe_km_rad[1]] + list(map(np.rad2deg, oe_km_rad[2:6])) + [oe_km_rad[6]]
oe_au_deg = [oe_km_deg[0] / au_to_km, oe_km_deg[1]] + list(map(np.rad2deg, oe_km_rad[2:6])) + [oe_km_rad[6]]

print("궤도요소의 최종 절사 평균값:")
print(
    f"==== 장반경(km,rad): {oe_km_rad[0]} km, 이심률: {oe_km_rad[1]}, "
    f"궤도경사: {oe_km_rad[2]} rad, 승교점경도: {oe_km_rad[3]} rad, 근일점편각: {oe_km_rad[4]} rad"
)
print(
    f"==== 장반경(km,deg): {oe_km_deg[0]} km, 이심률: {oe_km_deg[1]}, "
    f"궤도경사: {oe_km_deg[2]} deg, 승교점경도: {oe_km_deg[3]} deg, 근일점편각: {oe_km_deg[4]} deg"
)
print(
    f"==== 장반경(au,deg): {oe_au_deg[0]} au, 이심률: {oe_au_deg[1]}, "
    f"궤도경사: {oe_au_deg[2]} deg, 승교점경도: {oe_au_deg[3]} deg, 근일점편각: {oe_au_deg[4]} deg"
)
el_input = open('data\el_input.txt','w')
el_input.write(f"{oe_au_deg[0]} {oe_au_deg[1]} {oe_au_deg[2]} {oe_au_deg[3]} {oe_au_deg[4]} {oe_au_deg[5]} {oe_au_deg[6]} {asteriod_name}")


print('-------------종료-------------')

