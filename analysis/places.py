import argparse
import os
import re
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PLACES = [("the border", 20.49), ("Syabrubesi", 34.41), ("Betrawati", 64.68), ("Bidur", 73.68), ("Galchhi", 101.49)]
OFFS = (-0.2, -0.1, 0.0, 0.1, 0.2)
HALF = 400.0
RHO = (2700.0, 1000.0, 1080.0)

ap = argparse.ArgumentParser()
ap.add_argument("results"); ap.add_argument("prefix")
ap.add_argument("--tint", type=float, default=30.0, help="output interval of the run [s]")
ap.add_argument("--out", default="places"); ap.add_argument("--nproc", type=int, default=16)
ap.add_argument("--places", default="", help="other places as km:km:... (names become the km)")
ap.add_argument("--hmin", type=float, default=1.0,
                help="a phase thinner than this [m] carries no momentum in the sums (the paper: 1 m; 0 counts all flow)")
a = ap.parse_args()
if a.places:
    PLACES = [("%g km" % float(k), float(k)) for k in a.places.split(":")]
ASC = os.path.join(a.results, a.prefix + "_ascii")


def header(path):
    h = {}
    with open(path) as f:
        for _ in range(6):
            k, v = f.readline().split()
            h[k.lower()] = float(v)
    return h


def read(path):
    with open(path) as f:
        for _ in range(6):
            f.readline()
        z = np.array(f.read().split(), dtype=np.float64).reshape(NR, NC)
    return np.where(z == ND, 0.0, z)


g = header(os.path.join(ASC, a.prefix + "_elev.asc"))
NR, NC, CS = int(g["nrows"]), int(g["ncols"]), g["cellsize"]
ND = g.get("nodata_value", -9999.0)
X0 = g.get("xllcenter", g.get("xllcorner", 0.0) + CS / 2); Y0 = g.get("yllcenter", g.get("yllcorner", 0.0) + CS / 2)
XC = X0 + CS * np.arange(NC); YC = Y0 + CS * (NR - 1 - np.arange(NR))
elev = read(os.path.join(ASC, a.prefix + "_elev.asc"))
frames = sorted(int(m.group(1)) for f in os.listdir(ASC) if (m := re.fullmatch(a.prefix + r"_hflow1(\d{4})\.asc", f)))
TT = np.array(frames, float) * a.tint

P = np.loadtxt(os.path.join(HERE, "data", "path.csv"), delimiter=",", skiprows=1)
secs = []
for ip, (name, km) in enumerate(PLACES):
    for off in OFFS:
        i0 = int(np.argmin(np.abs(P[:, 2] - (km + off)))); px, py = P[i0, :2]
        t = P[min(i0 + 5, len(P) - 1), :2] - P[max(i0 - 5, 0), :2]; t /= np.linalg.norm(t)
        al = (XC[None, :] - px) * t[0] + (YC[:, None] - py) * t[1]
        ac = -(XC[None, :] - px) * t[1] + (YC[:, None] - py) * t[0]
        rr, cc = np.nonzero((np.abs(al) <= CS / 2) & (np.abs(ac) <= HALF))
        if rr.size == 0:
            raise SystemExit(f"section at {km + off:.2f} km lies outside the grid")
        secs.append((ip, rr, cc, t, float(elev[rr, cc].min())))
ns = len(secs)
from scipy.spatial import cKDTree
XX, YY = np.meshgrid(XC, YC)
dd, ii = cKDTree(P[:, :2]).query(np.c_[XX.ravel(), YY.ravel()])
KM = np.where(dd <= 2000.0, P[ii, 2], -1.0).reshape(NR, NC)


def one(q):
    f = lambda k: read(os.path.join(ASC, "%s_%s%04d.asc" % (a.prefix, k, q)))
    h = [f("hflow%d" % p) for p in (1, 2, 3)]
    me = [np.where(h[p - 1] > a.hmin, f("vflowy%d" % p), 0.0) for p in (1, 2, 3)]
    mn = [np.where(h[p - 1] > a.hmin, -f("vflowx%d" % p), 0.0) for p in (1, 2, 3)]
    bc = f("basechange")
    out = np.full((7, ns), np.nan)
    for s, (ip, rr, cc, t, bed0) in enumerate(secs):
        hk = [x[rr, cc] for x in h]; H = hk[0] + hk[1] + hk[2]
        out[0, s] = sum(float((me[p][rr, cc] * t[0] + mn[p][rr, cc] * t[1]).sum()) for p in range(3)) * CS
        out[1, s] = float(H.sum()) * CS
        out[2, s] = sum(RHO[p] * float(hk[p].sum()) for p in range(3)) / float(H.sum()) if H.sum() > 0 else 0.0
        w = H > 1.0
        if w.any():
            out[3, s] = float(np.median((elev[rr, cc] + bc[rr, cc] + H)[w]))
        out[4, s] = float(np.percentile(hk[0], 90))
        out[5, s] = float(np.percentile(H, 90))
    m = (h[0] + h[1]) > 1.0
    front = float(KM[m].max()) if m.any() else 0.0
    return out, front


nf = len(frames)
Q, AA, RH, SURF, ROCK, THK = (np.zeros((nf, ns)) for _ in range(6)); FRONT = np.zeros(nf)
with Pool(a.nproc) as pool:
    for j, (o, fr) in enumerate(pool.imap(one, frames, chunksize=2)):
        Q[j], AA[j], RH[j], SURF[j], ROCK[j], THK[j] = o[0], o[1], o[2], o[3], o[4], o[5]; FRONT[j] = fr
        if j % 80 == 0:
            print("output", j, "of", nf, flush=True)

V = np.where(AA > 0, Q / np.where(AA > 0, AA, 1), 0.0)
PR = RH * V ** 2
BED = np.array([s[4] for s in secs])
DEP = SURF - BED[None, :]
RISE = SURF - np.where(np.isfinite(SURF[0]), SURF[0], BED)[None, :]
dt = np.gradient(TT) / 60.0


def first(mask):
    k = np.flatnonzero(mask)
    return TT[k[0]] / 60 if k.size else np.nan


def vol80(q):
    c = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(TT))])
    if c[-1] <= 0:
        return np.nan, 0.0
    return (np.interp(0.9 * c[-1], c, TT) - np.interp(0.1 * c[-1], c, TT)) / 60, float(c[-1])


rows = []
for ip, (name, km) in enumerate(PLACES):
    per = []
    for s in [s for s in range(ns) if secs[s][0] == ip]:
        qs = np.convolve(Q[:, s], np.ones(3) / 3, "same")
        pr = PR[:, s]; prs = np.median(np.c_[np.r_[pr[:1], pr[:-1]], pr, np.r_[pr[1:], pr[-1:]]], axis=1)
        rise = np.nan_to_num(RISE[:, s], nan=0.0)
        fl = (THK[:, s] - THK[0, s]) > 1.0
        ta = first(rise > 1.0)
        win = (TT / 60 >= ta - 0.5) & fl if np.isfinite(ta) else np.zeros_like(fl)
        if not win.any():
            continue
        iq = int(np.argmax(np.where(win, qs, -np.inf))); pk = qs[iq]
        jp = int(np.argmax(np.where(win, prs, -np.inf)))
        qf = np.where(win, np.clip(Q[:, s], 0, None), 0.0)
        d80, vol = vol80(qf)
        per.append(dict(arrival=ta, peak_Q=pk, t_peak_Q=TT[iq] / 60,
                        depth=float(np.nanmax(np.where(win, DEP[:, s], np.nan))),
                        p=prs[jp], rho_p=RH[jp, s], V_p=abs(V[jp, s]), dur_flood=float(dt[fl].sum()),
                        dur_Q50=float(dt[(qs >= 0.5 * pk) & win].sum()), dur_Q80=float(dt[(qs >= 0.8 * pk) & win].sum()),
                        dur_vol80=d80, volume=vol, end=int(iq >= nf - 3)))
    if not per:
        rows.append((name, km, None)); continue
    rows.append((name, km, {k: float(np.nanmedian([d[k] for d in per])) for k in per[0]}))

keys = ["arrival", "peak_Q", "t_peak_Q", "depth", "p", "rho_p", "V_p", "dur_flood", "dur_Q50", "dur_Q80", "dur_vol80",
        "volume", "end"]
with open(a.out + ".csv", "w") as o:
    o.write("place,km," + ",".join(keys) + "\n")
    for name, km, m in rows:
        o.write("%s,%.2f," % (name, km) + (",".join("%.4g" % m[k] for k in keys) if m else ",".join("" for _ in keys)) + "\n")
np.savez_compressed(a.out + "_series.npz", t=TT, front_km=FRONT, places=np.array([p[0] for p in PLACES]),
                    **{f"{v}_{ip}": np.nanmedian(A[:, [s for s in range(ns) if secs[s][0] == ip]], axis=1)
                       for v, A in (("Q", Q), ("depth", DEP), ("p", PR), ("rho", RH), ("V", V), ("rise", RISE), ("thk", THK))
                       for ip in range(len(PLACES))})
print("\nplace          km    arrival  peak Q      at     depth  p      dur>80%Q  middle 80% volume  flood present")
print("                    [min]    [m3 s-1]    [min]  [m]    [MPa]  [min]     [min]              [min]")
for name, km, m in rows:
    if m is None:
        print(f"{name:13s} {km:6.2f}  not reached"); continue
    print(f"{name:13s} {km:6.2f}  {m['arrival']:6.1f}  {m['peak_Q']:9.0f}  {m['t_peak_Q']:6.1f}  {m['depth']:5.1f}  "
          f"{m['p'] / 1e6:5.2f}  {m['dur_Q80']:8.1f}  {m['dur_vol80']:17.1f}  {m['dur_flood']:13.1f}"
          + ("   (peak in the last outputs: run too short)" if m['end'] else ""))
print("\nfront position [km]: " + "  ".join(f"{TT[j] / 60:.0f} min {FRONT[j]:.1f}" for j in range(0, nf, 20)))
