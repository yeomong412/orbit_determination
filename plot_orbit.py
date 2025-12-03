import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, atan2, pi
from datetime import datetime, timedelta, timezone

# ---------- Constants ----------
K_GAUSS = 0.01720209895  # Gaussian gravitational constant [rad/day]
twopi = 2*np.pi
deg = np.deg2rad

# Earth (approx) — for path + current point (J2000 elements, simple model)
EARTH = dict(
    a=1.00000011,              # AU
    e=0.0167086,
    i=0.0,                     # deg (ecliptic ≈ equatorial here)
    raan=0.0,                  # deg
    argp=102.9372,             # deg
    M0=357.52715,              # deg at J2000 (≈ L0 - ω)
    JD_epoch=2451545.0         # J2000.0
)

# ---------- Time helpers ----------
def jd_to_datetime_kst(jd):
    days = jd - 2440587.5
    dt_utc = datetime(1970,1,1,tzinfo=timezone.utc) + timedelta(days=days)
    return dt_utc.astimezone(timezone(timedelta(hours=9)))

# ---------- Angle helpers ----------
def wrap_pm_pi(x): return (x + np.pi) % (2*np.pi) - np.pi
def wrap_0_2pi(x): return x % (2*np.pi)

# ---------- Kepler / geometry ----------
def mean_motion(a_au):
    # n = k / a^{3/2}  [rad/day]
    return K_GAUSS / (a_au**1.5)

def solve_kepler_elliptic(M, e, tol=1e-12, maxit=50):
    M = wrap_pm_pi(M)
    E = M if e < 0.8 else np.pi
    for _ in range(maxit):
        f = E - e*np.sin(E) - M
        fp = 1 - e*np.cos(E)
        dE = -f/fp
        E += dE
        if abs(dE) < tol:
            break
    return wrap_0_2pi(E)

def true_anomaly_from_M(M, e):
    E = solve_kepler_elliptic(M, e)
    s = sqrt((1+e)/(1-e))
    nu = 2*atan2(sin(E/2)*s, cos(E/2))
    return wrap_0_2pi(nu), E

def eccentric_from_true(nu, e):
    # quadrant-safe: tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
    E = 2*atan2(sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2))
    return wrap_0_2pi(E)

def perifocal_position(a,e,nu):
    p = a*(1-e**2)
    r = p / (1 + e*np.cos(nu))
    return np.array([r*np.cos(nu), r*np.sin(nu), 0.0])

def R1(x):
    c,s = np.cos(x), np.sin(x)
    return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def R3(x):
    c,s = np.cos(x), np.sin(x)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def inertial_from_orbital(x_pf, i, raan, argp):
    return (R3(raan) @ R1(i) @ R3(argp)) @ x_pf

def orbit_curve(a,e,i,raan,argp,n=1200):
    f = np.linspace(0,2*np.pi,n)
    p = a*(1-e**2)
    r = p/(1+e*np.cos(f))
    xyz_pf = np.vstack([r*np.cos(f), r*np.sin(f), np.zeros_like(f)])
    R = R3(raan) @ R1(i) @ R3(argp)
    return R @ xyz_pf

# ---------- I/O ----------
def load_elements(filename):
    """
    입력 형식 (공백/콤마 구분):
    a  e  i  raan  argp  nu  JD_obs
    단위: AU, deg, JD
    """
    bodies = []
    with open(filename,'r',encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('//'):
                continue
            parts = [p for p in line.replace(',', ' ').split() if p]
            if len(parts) != 8:
                raise ValueError("각 줄은 a e i raan argp nu JD_obs obj_name (8개)여야 합니다.")
            a,e,i,raan,argp,nu,JD_obs = map(float, parts[:7])
            obj_name = parts[-1]
            bodies.append({
                'a':a, 'e':e,
                'i':deg(i), 'raan':deg(raan), 'argp':deg(argp),
                'nu0':deg(nu), 'JD_obs':JD_obs, 'obj_name':obj_name
            })
    return bodies

# ---------- Earth state at JD ----------
def earth_state_at(JD):
    a = EARTH['a']; e = EARTH['e']
    i = deg(EARTH['i']); O = deg(EARTH['raan']); w = deg(EARTH['argp'])
    M0 = deg(EARTH['M0']); JD0 = EARTH['JD_epoch']
    n = mean_motion(a)
    M_t = wrap_0_2pi(M0 + n*(JD - JD0))
    nu,_ = true_anomaly_from_M(M_t, e)
    r_pf = perifocal_position(a,e,nu)
    r = inertial_from_orbital(r_pf, i, O, w)
    return r, nu

# ---------- Main plotting ----------
def plot_from_file(filename, JD_ref=None, plot_earth=True, show_markers=True):
    """
    filename : 'data\\el_input.txt'
    JD_ref   : 기준 JD.
               None이면 각 줄의 'nu0 @ JD_obs'를 그대로 표시.
               값을 주면 nu0→E0→M0로 바꿔 n*Δt만큼 전파해 새 ν를 계산.
    """
    bodies = load_elements(filename)

    # --- 색상 설정 ---
    EARTH_COLOR = 'tab:blue'  # 지구: 궤도+점 모두 파란색
    AST_COLORS = ['tab:orange','tab:green','tab:red','tab:purple',
                  'tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']  # 파랑 제외

    # 궤도 곡선
    Xe = None
    if plot_earth:
        Xe = orbit_curve(EARTH['a'], EARTH['e'],
                         deg(EARTH['i']), deg(EARTH['raan']), deg(EARTH['argp']))

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    if Xe is not None:
        ax.plot(Xe[0], Xe[1], Xe[2], label='Earth orbit', color=EARTH_COLOR)

    # 각 소행성
    for idx, el in enumerate(bodies, start=1):
        color_i = AST_COLORS[(idx-1) % len(AST_COLORS)]
        X = orbit_curve(el['a'], el['e'], el['i'], el['raan'], el['argp'])
        ax.plot(X[0], X[1], X[2], label=f'{el['obj_name']}', color=color_i)

        # 기준 JD 결정
        JD_base = el['JD_obs'] if JD_ref is None else JD_ref

        # 관측시각의 nu0 -> E0 -> M0
        E0 = eccentric_from_true(el['nu0'], el['e'])
        M0 = wrap_0_2pi(E0 - el['e']*sin(E0))

        # 전파: M(t) = M0 + n * (JD_base - JD_obs)
        n = mean_motion(el['a'])
        M_t = wrap_0_2pi(M0 + n*(JD_base - el['JD_obs']))

        # 현재 진근점이각
        nu_t,_ = true_anomaly_from_M(M_t, el['e'])

        # 위치
        r_pf = perifocal_position(el['a'], el['e'], nu_t)
        r_inert = inertial_from_orbital(r_pf, el['i'], el['raan'], el['argp'])

        if show_markers:
            ax.scatter([r_inert[0]],[r_inert[1]],[r_inert[2]], s=35, marker='o',
                       label=f'_nolegend_', color=color_i)

    # 태양
    ax.scatter([0],[0],[0], s=40, color='red')

    # 지구 현재 위치 점
    if plot_earth:
        JD_base = (EARTH['JD_epoch'] if JD_ref is None else JD_ref)
        rE, nuE = earth_state_at(JD_base)
        ax.scatter([rE[0]],[rE[1]],[rE[2]], s=50, marker='o',
                   label=f'_nolegend_', color=EARTH_COLOR)

    ax.set_xlabel('X [AU]'); ax.set_ylabel('Y [AU]'); ax.set_zlabel('Z [AU]')

    # aspect 보정
    xs, ys, zs = [], [], []
    for line in ax.lines:
        x,y,z = line.get_data_3d()
        xs.extend(x); ys.extend(y); zs.extend(z)
    for p in ax.collections:
        Xo,Yo,Zo = p._offsets3d
        xs.extend(np.atleast_1d(Xo)); ys.extend(np.atleast_1d(Yo)); zs.extend(np.atleast_1d(Zo))
    xs, ys, zs = np.array(xs), np.array(ys), np.array(zs)
    r = np.ptp([xs.min(), xs.max(), ys.min(), ys.max(), zs.min(), zs.max()]) / 2
    cx, cy, cz = xs.mean(), ys.mean(), zs.mean()
    ax.set_xlim(cx-r, cx+r); ax.set_ylim(cy-r, cy+r); ax.set_zlim(cz-r, cz+r)

    # 제목(기준 시각)
    if JD_ref is None:
        title = "Keplerian Orbits @ each line's JD_obs"
    else:
        title = f"Keplerian Orbits @ JD {JD_ref:.5f} ({jd_to_datetime_kst(JD_ref).strftime('%Y-%m-%d %H:%M:%S KST')})"
    ax.set_title(title)
    ax.legend(loc='upper right')
    plt.show()

# -------- Usage --------
# el_input.txt: a e i raan argp nu JD_obs  (공백/콤마 구분, 단위: AU, deg, JD)
JD_ref =2460992.95986111        # None이면 각 라인의 JD_obs에서의 nu로 표시
plot_from_file('data\\el_input.txt', JD_ref=JD_ref, plot_earth=True, show_markers=True)
