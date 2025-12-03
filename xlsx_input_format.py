from datetime import datetime, timezone
import math
import pandas as pd
from parallax_distance import *

def iso_to_julian(iso_string):
    s = str(iso_string).replace("Z", "+00:00")
    dt = datetime.fromisoformat(s)

    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)

    year, month = dt.year, dt.month
    day = dt.day + (dt.hour + (dt.minute + dt.second/60)/60)/24

    if month <= 2:
        year -= 1
        month += 12

    A = year // 100
    B = 2 - A + (A // 4)

    JD = (math.floor(365.25*(year+4716)) +
          math.floor(30.6001*(month+1)) +
          day + B - 1524.5)
    return JD


# ====== 데이터 불러오기 ======
excel_path = r"data\프기연통합_관측결과.xlsx"
df = pd.read_excel(excel_path)

df.columns = df.columns.str.strip()
df["관측자"] = df["관측자"].ffill().astype(str).str.strip()
df["시각"] = df["시각"].astype(str)

observer = "여은수"

df_obs = df[df["관측자"] == observer][["시각", "RA", "DEC"]].dropna()

groups = df_obs.groupby("시각")

lat1, lon1, h1 = 36.523, 127.248, 0.060
lat2, lon2, h2 = 34.5263, 127.4468, 0.060

results = []

for t, g in groups:
    if len(g) < 2:
        continue

    r1 = g.iloc[0]
    r2 = g.iloc[1]

    RA1 = float(r1["RA"])
    DEC1 = float(r1["DEC"])
    RA2 = float(r2["RA"])
    DEC2 = float(r2["DEC"])

    d_km, d_au = parallax_distance_from_absolute_positions(
        lat1, lon1, h1,
        lat2, lon2, h2,
        RA1, DEC1,
        RA2, DEC2
    )

    jd = iso_to_julian(t)

    results.append((jd, RA1, DEC1, RA2, DEC2, d_au))

# 시간순 정렬
results.sort(key=lambda x: x[0])

# ====== 숫자만 출력 ======
for jd, RA1, DEC1, RA2, DEC2, AU in results:
    print(f"{RA1:.8f} {DEC1:.8f} {AU:.10f} {jd:.8f}")
