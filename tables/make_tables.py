"""Recompute Tables 5 and 9 of the paper from the archived per-case results. These are tables of the ORIGINAL runs: the
decisions and timings are read from benchmarks/cases/e4_G*.jsonl (Table 5; the baselines ran two processes at once) and from
benchmarks/results/e5final_main.json (Table 9; single process, cold start). Nothing is recomputed here; new calculations are in reproduce_results.py."""
import json, glob, statistics as st, pathlib
R = pathlib.Path(__file__).resolve().parent.parent
C = [json.loads(l) for f in sorted(glob.glob(str(R/"benchmarks/cases/e4_G*.jsonl"))) for l in open(f)]
print("Table 5 (agreement with the 60-digit reference; archived times, two-process runs)")
for key, name in (("B1", "finite difference, backward Euler"), ("B2", "full-history floating-point spectral, 200 points"), ("B3", "floating-point endpoint solver")):
    ok = sum(c[key]["dec"] == c["truth"] for c in C); bad = len(C)-ok
    un = sum(c[key]["dec"] == "SAFE_ALL" and c["truth"] != "SAFE_ALL" for c in C)
    print(f"  {name:50s} consistent {ok:5d}  disagreement {bad:4d}  incorrect universal-safety {un:4d}  median {st.median(c[key]['t'] for c in C)*1e3:6.2f} ms")
E = json.load(open(R/"benchmarks/results/e5final_main.json"))
print("Table 9 (certified implementations; archived single-process cold-start times)")
import numpy as np
for key, name in (("fr", "original 256-bit"), ("v7", "optimized 96-bit")):
    t = [r[key]["t"] for r in E]
    ok = sum(r[key]["dec"] == r["truth"] for r in E); inc = sum(r[key]["dec"] == "INCONCLUSIVE" for r in E)
    print(f"  {name:20s} consistent {ok:5d}  disagreement {len(E)-ok-inc}  inconclusive {inc}  median {np.median(t)*1e3:6.2f} ms  P90 {np.percentile(t,90)*1e3:6.2f} ms")
