# Certified peak-temperature computation for a heat equation with switching Robin boundary conditions

Code and data behind the paper *Structure-Exploiting Certified Peak-Temperature Computation for a Heat Equation with Switching
Robin Boundary Conditions* (S. S. Avro, submitted to Mathematics and Computers in Simulation; citation and DOI below once released).
The repository reproduces the numerical evidence in the paper. The mathematical proofs are in the manuscript and are not repeated here.

## 1. Problem and assumptions
Slab `u_t = u_xx`, `u_x(0,t)=0`, `u_x(1,t)+B(t)u(1,t)=0`, `u(x,0)=1`, with `θ = 1-u`, and a Robin coefficient that switches once:
`B = B1` on `(0, t_s)`, `B = B2` on `(t_s, t_p]`, with `B1, B2 >= 0` (dimensionless Biot numbers) and `0 < t_s < t_p` (Fourier numbers).
The paper proves that the maximum surface temperature over `[0, t_p]` equals `max{θ_s(t_s), θ_s(t_p)}` and that, over a box of
coefficients, the extremes occur at the lower–lower and upper–upper corners. See `docs/mathematical_assumptions.md`.

## 2. What the software decides
`solver.decide(box, Fs, Fp, L)` takes a box `((B1lo,B1hi),(B2lo,B2hi))`, the times `Fs=t_s`, `Fp=t_p` and a dimensionless threshold `L`, and returns one of
* `SAFE_ALL`: the peak is strictly below `L` for every coefficient pair in the box,
* `EXCEEDS_ALL`: the peak is at or above `L` for every pair,
* `MIXED`: the upper corner is certified at or above `L` and the lower corner certified below it,
* `INCONCLUSIVE`: the enclosures cannot decide at the requested tolerance (never a guess).

Invalid input raises an error. The guarantee is conditional on the proofs in the paper, the truncation bound of Lemma 4.1 and the correctness of the interval kernel; the tests in this repository support the implementation and do not prove it.

## 3. Requirements
Python 3.13, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2 (`requirements.txt`, `environment.yml`) and the `bc` calculator for the 120-digit oracle test. These are the versions used for the reported results (two virtual cores of an Intel Xeon 2.10 GHz, about 8 GB). The solver itself uses only Python integers, NumPy and SciPy for the floating-point baselines.

## 4. Quick start
```
pip install -r requirements.txt
python scripts/verify_archive.py            # SHA-256 of the frozen and optimized implementations, 1,600 cases
python scripts/reproduce_results.py --quick # about 5 minutes; skips the three longest steps
python scripts/reproduce_results.py         # everything, about 10 minutes
```
Reproduction runs the real code and compares it with the archived results (`benchmarks/results`, `ablation`, `comparisons`); it never substitutes stored results, and reports skipped steps as SKIPPED. Results go to `output/`.

## 5. Example
```
cd src
python -c "import solver; print(solver.decide(((3.0,3.2),(1.0,1.1)), 0.05, 0.2, 0.5))"
```
Returns the decision and a dictionary with the decision-sufficient enclosure width and whether the 256-bit repeat was used.

## 6. Contents
`src/` optimized (96-bit) implementation and fallback; `frozen/` original 256-bit implementation with design and SHA-256; `benchmarks/` the 1,600 cases, reference classifications and per-method results; `validation/` kernel, root, oracle, mutation, tail-bound and edge-case tests; `comparisons/` finite-volume study (the finite-difference, full-history and floating-point baselines are in `src/e4.py`); `ablation/` theorem-dependence and tolerance-resolution experiments; `examples/` Scenarios A and B; `figures/` figure code and final images; `docs/` assumptions, protocol and verification scope. The mapping from every table and figure of the paper to code and data is in `REPRODUCIBILITY.md`.

## 7. Limits
No independent Arb/FLINT cross-check has been run, the theorem is two-stage only, and the certified solver is slower than the floating-point endpoint solver. Timings depend on the machine; the paper's timings are archived separately from any new measurement.

## 8. Citation and archive
See `CITATION.cff`. Zenodo DOI: [to be inserted at release]. License: MIT (`LICENSE`).
