import numpy as np

# ----- WGS84 -----
a = 6378.1370  # km
b = 6356.7523  # km
e2 = 1 - (b*b)/(a*a)


def ecef(lat_deg, lon_deg, h=0.0):
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    N = a / np.sqrt(1 - e2*np.sin(lat)**2)
    x = (N + h) * np.cos(lat) * np.cos(lon)
    y = (N + h) * np.cos(lat) * np.sin(lon)
    z = ((1 - e2) * N + h) * np.sin(lat)
    return np.array([x, y, z])


def chord_distance(lat1, lon1, h1, lat2, lon2, h2):
    r1 = ecef(lat1, lon1, h1)
    r2 = ecef(lat2, lon2, h2)
    return float(np.linalg.norm(r1 - r2))


# ===============================
#   ★★★ GUI 로직 기반 gp 계산 ★★★
# ===============================
def compute_gp_from_gui_logic(ra1_h, dec1_deg, ra2_h, dec2_deg, lat1_deg, lat2_deg):
    """
    Tkinter GUI의 compute_gp_from_values()와 동일한 계산을 수행.
    RA는 hour 단위 입력이다.
    """

    # RA hour → arcsec 변환 방식 (GUI 방식 그대로)
    ra1_arcsec = ra1_h * 3600 * 15 * np.cos(np.deg2rad(dec1_deg))
    ra2_arcsec = ra2_h * 3600 * 15 * np.cos(np.deg2rad(dec2_deg))

    # RA 차이
    ra_diff_arcsec = abs(
        (ra1_arcsec - ra2_arcsec + 12.0 * 15 * 3600) % (24.0 * 15 * 3600)
        - 12.0 * 15 * 3600
    )

    # DEC 차이
    dec_diff_arcsec = abs((dec2_deg - dec1_deg) * 3600.0)

    dec_diff_rad = np.deg2rad(dec_diff_arcsec / 3600.0)
    ra_ang_rad   = np.deg2rad(ra_diff_arcsec / 3600.0)

    # GUI 로직의 분모
    dec_mean_deg = 0.5 * (dec1_deg + dec2_deg)
    denom = np.cos(np.deg2rad((lat1_deg + lat2_deg)/2.0 - dec_mean_deg))

    # GUI 로직과 동일한 최종 각도
    val = np.arccos(np.cos(ra_ang_rad) * np.cos(dec_diff_rad)) / denom

    gp_arcsec = 3600.0 * np.rad2deg(val)
    return gp_arcsec


def distance_from_parallax(baseline_km, gp_arcsec):
    d_km = (baseline_km / gp_arcsec) * 206265.0
    d_AU = d_km / 1.496e8
    return d_km, d_AU


# ===========================================================
#          최종: GUI 계산식 100% 동일한 순수 함수
# ===========================================================
def parallax_distance_from_absolute_positions(
    lat1, lon1, h1,
    lat2, lon2, h2,
    RA1_h, DEC1_deg,
    RA2_h, DEC2_deg
):
    """
    이 함수는 Tkinter GUI 계산 방식과 완전히 동일하게 동작한다.
    RA는 hour 단위로 입력해야 한다.
    """

    # ① baseline 계산 (동일)
    baseline = chord_distance(lat1, lon1, h1, lat2, lon2, h2)

    # ② gp 계산 (GUI 로직 그대로)
    gp_arcsec = compute_gp_from_gui_logic(
        RA1_h/15, DEC1_deg,
        RA2_h/15, DEC2_deg,
        lat1, lat2
    )

    # ③ 거리 변환 (동일)
    d_km, d_au = distance_from_parallax(baseline, gp_arcsec)

    print(f"Baseline (km): {baseline:.6f}")
    print(f"Geometric parallax gp (arcsec): {gp_arcsec:.10f}")

    return d_km, d_au
