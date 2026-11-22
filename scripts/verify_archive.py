#!/usr/bin/env python3
"""Check that the archived files are the ones used in the paper: SHA-256 of the frozen original implementation and of the
optimized implementation, and the structure of the 1,600-case benchmark. Exit status 0 only if every check passes."""
import hashlib, json, glob, sys, collections
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_sums(p):
    d = {}
    for line in Path(p).read_text().splitlines():
        parts = line.split()
        if len(parts) == 2: d[parts[1]] = parts[0]
    return d
ok = True
def check(label, cond, extra=""):
    global ok; ok &= bool(cond); print(("PASS  " if cond else "FAIL  ") + label + (("  " + extra) if extra else ""))
# 1. frozen original implementation (e4_freeze.sha256 names the files; e4ref.py and closed.py live in src/ and are shared)
for name, h in read_sums(ROOT/"frozen"/"e4_freeze.sha256").items():
    f = ROOT/"frozen"/name
    if not f.exists(): f = ROOT/"src"/name
    check(f"frozen design file {name}", f.exists() and sha(f) == h)
# 2. optimized implementation (default_V7.sha256 lists ivk.py, e3.py, e4.py of the final default configuration)
for name, h in read_sums(ROOT/"src"/"default_V7.sha256").items():
    check(f"optimized implementation file {name}", sha(ROOT/"src"/name) == h)
# 3. benchmark structure
cases = [json.loads(l) for f in sorted(glob.glob(str(ROOT/"benchmarks"/"cases"/"e4_G*.jsonl"))) for l in open(f)]
cnt = collections.Counter(c['group'] for c in cases)
check("1,600 benchmark cases", len(cases) == 1600, str(len(cases)))
check("group sizes G1 600, G2 400, G3 400, G4 200", dict(cnt) == {'G1': 600, 'G2': 400, 'G3': 400, 'G4': 200}, str(dict(cnt)))
check("seeds unique", len({c['seed'] for c in cases}) == 1600)
print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED")
sys.exit(0 if ok else 1)
