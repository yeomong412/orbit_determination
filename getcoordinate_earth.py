from astroquery.jplhorizons import Horizons
from astropy.time import Time
from datetime import datetime
import numpy as np

def getcoordinate_earth(date):

    # 원하는 Julian Date 값 (예: 2460000.5)
    jd_value = date
    # 단일 시점을 리스트로 전달하며, 'JD' 접두어 없이 숫자형 문자열로 지정
    epoch_str = str(jd_value)

    earth = Horizons(
        id='399',
        location='500@0',
        epochs=[epoch_str],  # 단일 시점 전달
        id_type=None
    )

    # refplane은 소문자 'ecliptic'으로 지정
    vectors = earth.vectors(refplane='ecliptic')

    # 결과 테이블에서 시간과 x, y, z 데이터를 추출
    data = vectors['datetime_str', 'x', 'y', 'z']
    #print("Horizons로부터 받아온 지구 좌표 데이터:")
    #print(data)

    # (선택사항) 시간 문자열을 datetime 객체 및 Julian Date로 변환
    time_str = data['datetime_str'][0]
    # 예: "A.D. 2025-Jan-01 00:00:00.0000" 형식 -> "A.D. " 제거 후 파싱
    time_dt = datetime.strptime(time_str.replace("A.D. ", ""), "%Y-%b-%d %H:%M:%S.%f")
    time_numeric = Time(time_dt).jd

    # x, y, z 값을 float64로 변환
    x_val = np.float64(data['x'][0])
    y_val = np.float64(data['y'][0])
    z_val = np.float64(data['z'][0])

    # 하나의 행 ([t, x, y, z])을 2차원 배열로 생성 (단일 행)
    matrix = np.array([x_val, y_val, z_val])
    #print("\n[t, x, y, z]로 구성된 지구 좌표 array:")
    #print(matrix)

    return matrix
