"""Scenarios A and B of Section 5.6: corner peaks, numerical predictions and certified decisions.
Run from a directory that contains the solver modules (scripts/reproduce_results.py arranges this)."""
import json, sys
import e4, e3, solver
a = 0.2/(1200*1800); Lm = 0.004; k = 0.2; T0, Tinf = 25.0, 100.0
SC = {"A": dict(h=(60, 15), ts=30, tp=180, lim=65), "B": dict(h=(500, 5), ts=10, tp=180, lim=60)}
def cls(pU, pL, L): return 'SAFE_ALL' if pU < L else 'EXCEEDS_ALL' if pL >= L else 'MIXED'
out = {}
for name, s in SC.items():
    B1, B2 = s['h'][0]*Lm/k, s['h'][1]*Lm/k; Fs, Fp = a*s['ts']/Lm**2, a*s['tp']/Lm**2
    L = (s['lim']-T0)/(Tinf-T0); box = ((0.8*B1, 1.2*B1), (0.8*B2, 1.2*B2))
    corners = {'lower': (box[0][0], box[1][0]), 'upper': (box[0][1], box[1][1])}
    ev = {c: e4.endpoints_float(*b, Fs, Fp) for c, b in corners.items()}           # (theta_s(ts), theta_s(tp))
    g20 = {c: e4.full_history_peak(*b, Fs, Fp, K=300, nt=20) for c, b in corners.items()}
    dec, info = solver.decide(box, Fs, Fp, L)
    enc = {}
    for c, b in corners.items():
        (lo, hi), ths, thp, M = e3.theta_peak(*b, Fs, Fp, 1e-10); S = 2.0**96
        enc[c] = (float(hi-lo)/S)
    T = lambda x: T0 + (Tinf-T0)*x
    res = dict(B1=B1, B2=B2, Fs=Fs, Fp=Fp, theta_lim=L,
        switch_C={c: T(v[0]) for c, v in ev.items()}, end_C={c: T(v[1]) for c, v in ev.items()},
        peak_C={c: T(max(v)) for c, v in ev.items()},
        predictions={'end-time only': cls(ev['upper'][1], ev['lower'][1], L), '20-point grid': cls(g20['upper'], g20['lower'], L),
                     'float endpoint': cls(max(ev['upper']), max(ev['lower']), L)},
        certified=dec, decision_width_theta=info['width'], enclosure_width_tol1e10=enc)
    out[name] = res
    print(name, json.dumps(res, indent=1, default=float))
json.dump(out, open('scenarios.json', 'w'), indent=1, default=float)
