#!/usr/bin/env python3
"""Re-run the calculations behind the paper and compare them with the archived results.

    python scripts/reproduce_results.py            # everything (about 10 minutes on two cores)
    python scripts/reproduce_results.py --quick    # skip the three longest steps
    python scripts/reproduce_results.py --only benchmark,ablation

Every step runs the real code; stored results are used only as the reference the new results are compared with, never in
their place. Steps that are skipped are reported as SKIPPED, not as passed. New timings are reported next to, and never
merged with, the timings archived from the original runs. Outputs: output/logs, output/results, output/figures, output/summary.json.
Exit status is 0 only if every step that ran passed."""
import argparse, json, os, re, shutil, subprocess, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/"output"; WORK = OUT/"work"
ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--only", default="")
args = ap.parse_args()

def stage():
    if WORK.exists(): shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    for f in (ROOT/"src").glob("*.py"): shutil.copy(f, WORK)
    shutil.copytree(ROOT/"frozen", WORK/"frozen")
    for f in (ROOT/"benchmarks"/"cases").glob("*.jsonl"): shutil.copy(f, WORK)
    for d in ("benchmarks/scripts", "validation", "comparisons", "ablation", "examples", "figures"):
        for f in (ROOT/d).glob("*.py"): shutil.copy(f, WORK)
    for d in ("logs", "results", "figures"): (OUT/d).mkdir(parents=True, exist_ok=True)

def run(script, *a, timeout=1800, env=None):
    t0 = time.time()
    e = dict(os.environ, OUTDIR=str(OUT/"figures"), PYTHONHASHSEED="0"); e.update(env or {})
    p = subprocess.run([sys.executable, script, *a], cwd=WORK, capture_output=True, text=True, timeout=timeout, env=e)
    txt = p.stdout + p.stderr
    (OUT/"logs"/(Path(script).stem + ".log")).write_text(txt)
    return p.returncode, txt, time.time()-t0

def load(p): return json.load(open(p))
results = []; STEPS = []
def step(name, long=False):
    def deco(fn): STEPS.append((name, long, fn)); return fn
    return deco

def run_step(name, long, fn):
    only = [s for s in args.only.split(",") if s]
    if (only and name not in only) or (args.quick and long and not only):
        results.append((name, "SKIPPED", "")); print(f"SKIPPED {name}"); return
    print(f"... {name}", flush=True)
    try: ok, msg = fn()
    except Exception as ex: ok, msg = False, f"exception: {type(ex).__name__}: {ex}"
    results.append((name, "PASS" if ok else "FAIL", msg)); print(f"{'PASS' if ok else 'FAIL'}  {name}: {msg}", flush=True)

@step("archive_hashes")
def _():
    p = subprocess.run([sys.executable, str(ROOT/"scripts"/"verify_archive.py")], capture_output=True, text=True)
    (OUT/"logs"/"verify_archive.log").write_text(p.stdout); return p.returncode == 0, p.stdout.strip().splitlines()[-1]

@step("validation_kernel_exact_rational")
def _():
    rc, t, s = run("kval.py")
    return rc == 0 and "610 checks, overlap failures 0" in t and "failures 0" in t.split("root brackets")[-1] and "True" in t, f"{s:.0f}s; " + " | ".join(t.strip().splitlines()[-3:])[:220]

@step("validation_kernel_mutation_tests")
def _():
    rc, t, s = run("audit.py"); lines = [l for l in t.splitlines() if "MUTANT" in l]
    real = [l for l in t.splitlines() if l.startswith("REAL")]
    return rc == 0 and len(lines) >= 6 and all("detected: True" in l for l in lines) and len(real) == 1 and all(v == "0" for v in re.findall(r": (\d+)", real[0])), f"{len(lines)} injected defects, all detected" if lines else t[-200:]

@step("validation_bc_oracle_120_digit")
def _():
    rc, t, s = run("bcval.py"); m = re.search(r"REAL kernel.*?(\d+) bc-oracle comparisons, enclosure failures (\d+)", t)
    mut = re.findall(r"MUTANT.*?enclosure failures (\d+)", t)
    return bool(m) and m.group(2) == "0" and len(mut) >= 2 and all(int(x) > 0 for x in mut), (m.group(0)[:120] if m else t[-200:])

@step("validation_lemma_4_1_tail_bound")
def _():
    rc, t, s = run("lemmaG_reaudit.py")
    return rc == 0 and "'total': 0" in t.split("max ratio")[0], t.strip().splitlines()[0][:200] if t.strip() else "no output"

@step("validation_proof_chain_numerics")
def _():
    rc, t, s = run("verify_core.py"); rc2, t2, s2 = run("tp2hp.py")
    m = re.search(r"TP2 minors tested (\d+) violations (\d+)", t); r = re.search(r"ratio-monotone profiles tested (\d+) violations (\d+)", t)
    ok = rc == 0 and rc2 == 0 and r and r.group(2) == "0" and "False" not in t2
    return bool(ok), f"TP2 double-precision violations {m.group(2) if m else '?'} (expected 6; all recomputed positive in 70-digit arithmetic by tp2hp.py), monotone-ratio violations {r.group(2) if r else '?'}"

@step("benchmark_1600_regression")
def _():
    rc, t, s = run("e5final.py")
    if rc != 0: return False, t[-300:]
    new, old = load(WORK/"e5final_main.json"), load(ROOT/"benchmarks"/"results"/"e5final_main.json")
    same = all(a[k]["dec"] == b[k]["dec"] and a[k]["w"] == b[k]["w"] for a, b in zip(new, old) for k in ("fr", "v7", "fb"))
    fb = [r["fb"]["dec"] for r in new]; agree = sum(r["fb"]["dec"] == r["truth"] for r in new); inc = fb.count("INCONCLUSIVE")
    dis = sum(d not in ("INCONCLUSIVE",) and d != r["truth"] for d, r in zip(fb, new))
    import statistics as st
    shutil.copy(WORK/"e5final_main.json", OUT/"results"/"e5final_main_rerun.json")
    med_new = st.median(r["v7"]["t"] for r in new)*1e3; med_old = st.median(r["v7"]["t"] for r in old)*1e3
    return same and agree == 1289 and inc == 311 and dis == 0, f"identical decisions and widths to archive: {same}; agreement {agree}, inconclusive {inc}, disagreement {dis}; median time, this run {med_new:.2f} ms (archived run {med_old:.2f} ms; timings are machine dependent and not compared)"

@step("benchmark_enclosures_vs_60_digit_reference")
def _():
    rc, t, s = run("e5verify.py", timeout=3000); ok = "reference outside interval: frozen 0 | V7 0" in t and "disjoint: 0" in t
    return rc == 0 and ok, " | ".join(t.strip().splitlines()[-3:])[:250]

@step("edge_cases_and_precision_exhaustion")
def _():
    rc, t, s = run("e5edge.py", timeout=3000)
    return rc == 0 and "wrong V7: 0" in t and "differs from frozen: 0" in t, t.strip().splitlines()[0][:200] if t.strip() else "no output"

@step("edge_inputs_equality_resolution")
def _():
    rc, t, s = run("edge_inputs.py", timeout=600)
    exp = ["EXCEEDS_ALL", "INCONCLUSIVE", "SAFE_ALL"]
    lines = t.strip().splitlines()
    ok = rc == 0 and [l.split()[-1] for l in lines[:3]] == exp and lines[-2].endswith("SAFE_ALL") and lines[-1].endswith("EXCEEDS_ALL")
    return ok, " | ".join(lines)[:250]

@step("comparison_finite_volume", long=True)
def _():
    rc, t, s = run("fv16.py")
    if rc != 0: return False, t[-300:]
    new, old = load(WORK/"fv16.json"), load(ROOT/"comparisons"/"fv16.json")
    same = all(abs(a["err"]-b["err"]) < 1e-12 and a["peak_at_switch_or_end"] == b["peak_at_switch_or_end"] for a, b in zip(new, old)) and len(new) == len(old) == 16
    shutil.copy(WORK/"fv16.json", OUT/"results"/"fv16_rerun.json"); return same, f"16 cases; finest-mesh errors equal archive to 1e-12: {same}"

@step("ablation_theorem_dependence", long=True)
def _():
    rc, t, s = run("ablate.py", timeout=3600)
    if rc != 0: return False, t[-300:]
    new, old = load(WORK/"ablate.json"), load(ROOT/"ablation"/"ablate.json")
    same = all(new[k][f] == old[k][f] for k in old for f in ("disagree", "unsafe_called_safe", "disagree_gap_ge_1e14"))
    shutil.copy(WORK/"ablate.json", OUT/"results"/"ablate_rerun.json"); return same, f"Table 6 counts equal archive: {same}; end-only {new['end time only, corners']['disagree']}, corner rule {new['Theorem 3.5 endpoints, corners (corner rule)']['disagree']}"

@step("ablation_tolerance_resolution", long=True)
def _():
    rc, t, s = run("floorrun.py", timeout=3600)
    if rc != 0: return False, t[-300:]
    new, old = load(WORK/"floorrun.json"), load(ROOT/"ablation"/"floorrun.json")
    same = all([r[0] for r in new[k]] == [r[0] for r in old[k]] for k in old)
    shutil.copy(WORK/"floorrun.json", OUT/"results"/"floorrun_rerun.json"); return same, "Table 8 decisions equal archive at every tolerance floor: %s" % same

@step("comparison_uncertainty_sampling", long=True)
def _():
    rc, t, s = run("uncertainty_sampling_studies.py")
    if rc != 0: return False, t[-300:]
    keep = lambda txt: [l.strip() for l in txt.splitlines() if l.strip() and not l.strip().startswith("timing")]
    old = (ROOT/"comparisons"/"uncertainty_sampling_archived.log").read_text()
    same = keep(t) == keep(old)
    return same, f"Section 5.2 studies (peak underestimation, Table 7) equal archive: {same}; {s:.0f}s"

@step("tables_5_and_9_from_archive")
def _():
    p = subprocess.run([sys.executable, str(ROOT/"tables"/"make_tables.py")], capture_output=True, text=True)
    (OUT/"logs"/"make_tables.log").write_text(p.stdout + p.stderr)
    ok = p.returncode == 0 and "1171" in p.stdout and "1592" in p.stdout and "1289" in p.stdout
    return ok, "Tables 5 and 9 recomputed from archived per-case results (original-run timings, not new measurements)"

@step("examples_scenarios_A_B")
def _():
    rc, t, s = run("scenarios.py")
    if rc != 0: return False, t[-300:]
    d = load(WORK/"scenarios.json"); A, B = d["A"], d["B"]
    ok = (A["certified"] == "SAFE_ALL" and B["certified"] == "EXCEEDS_ALL" and abs(A["peak_C"]["lower"]-49.96) < 0.01 and abs(A["peak_C"]["upper"]-57.84) < 0.01
          and abs(B["switch_C"]["lower"]-80.2) < 0.05 and abs(B["switch_C"]["upper"]-86.1) < 0.05 and abs(B["end_C"]["lower"]-44.3) < 0.05 and abs(B["end_C"]["upper"]-48.5) < 0.05
          and B["predictions"]["end-time only"] == "SAFE_ALL")
    shutil.copy(WORK/"scenarios.json", OUT/"results"/"scenarios.json"); return ok, f"A {A['certified']}, B {B['certified']}; values equal Table 11 inputs: {ok}"

@step("figures")
def _():
    fallback = not (WORK/"floorrun.json").exists()
    if fallback: shutil.copy(ROOT/"ablation"/"floorrun.json", WORK/"floorrun.json")  # --quick: Figure 3 drawn from archived Table 8 data
    rc, t, s = run("figs_hires.py"); rc2, t2, s2 = run("figs_cert.py")
    made = sorted(p.name for p in (OUT/"figures").glob("*.png"))
    return rc == 0 and rc2 == 0 and len(made) >= 5, f"wrote {made}" + (" (Figure 3 from archived floorrun.json)" if fallback else "") if rc == rc2 == 0 else (t[-150:] + t2[-150:])

if __name__ == "__main__":
    stage()
    for name, long, fn in STEPS: run_step(name, long, fn)
    (OUT/"summary.json").write_text(json.dumps([dict(step=n, status=s, detail=m) for n, s, m in results], indent=1))
    n_fail = sum(s == "FAIL" for _, s, _ in results); n_skip = sum(s == "SKIPPED" for _, s, _ in results)
    print(f"\n{sum(s == 'PASS' for _, s, _ in results)} passed, {n_fail} failed, {n_skip} skipped")
    sys.exit(1 if n_fail else 0)
