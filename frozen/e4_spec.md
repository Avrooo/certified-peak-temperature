# E4 pre-registered design (frozen before any E4 result was seen)

Question: does the structure of two-stage convective heating let us certify peak-temperature decisions more efficiently, accurately, or simply than existing numerical approaches? Three structural elements are measured separately.

## Decision task
Box of uncertain (B1,B2), fixed Fs<Fp, threshold L (a binary float). Truth class from the independent reference:
SAFE_ALL (peak at upper corner < L), EXCEEDS_ALL (peak at lower corner >= L), MIXED (otherwise).
Peak at a corner = max{theta_s(Fs), theta_s(Fp)} (theorem), evaluated by the reference at 60 digits (Decimal; own Newton roots; own Taylor sin/cos).

## Methods
- B1 finite-difference (backward Euler, N=100 nodes, 1200 steps total, peak over all steps), run at both corners.
- B2 full-history float spectral (K=300 modes, uniform 200-point time grid on (0,Fp], NOT containing Fs), both corners.
- B3 theorem-aware float spectral: only theta_s(Fs), theta_s(Fp), K from the Lemma-G float estimate (<=1e-14), both corners.
- B4 fixed-mode interval evaluator (M=64), tail bound included, no refinement; inconclusive if the interval straddles L.
- B5 adaptive E3 (tolerance schedule 1e-3..1e-13, both endpoints at each level, corner order U then Lc).
- B6 decision-first A2: float screening (B3 values) chooses corner order and the starting tolerance (|screen-L|/4 clipped to [1e-13,1e-3]); evaluates only the endpoints needed (one endpoint suffices to prove "above"; both needed for "below"); same final certificate rule as B5.
B1,B2 use the corner rule (otherwise they would need parameter sampling); this credits the corner reduction to them.

## Timing
Per case and method, wall time incl. root finding / preprocessing; caches cleared before each method.

## Groups (cases)
G1 ordinary 600; G2 extreme 400 (Fs 3e-4..3e-3, B up to 500, B1~B2 to 1e-9, B=0); G3 near-threshold 400 (margins 1e-8..1e-13); G4 independent 200 (different generator: wider B, Fs, long pulses; seed set).
Class mix per group: 25% safe, 25% exceed, 30% mixed, 20% near-equality (|L-peak|<=3e-15).
Threshold: L = float(reference peak +/- margin); the ACTUAL separation (L - reference, 60 digits) is stored and truth is taken from it.

## Reported
Correct / wrong / inconclusive per method; false-safe (SAFE_ALL declared when truth is not SAFE_ALL); time median/P90/max; certificate width; A0 (B4) vs A1 (B5) vs A2 (B6) ablation; false-safe by float methods vs near-threshold margin.

## Not done / limits (declared in advance)
No Arb/FLINT cross-check (not installable here). No general validated PDE integrator baseline (no implementation for switched Robin available). Reference shares the closed-form series and Lemma-G truncation with E3 (independent arithmetic, roots, evaluation order).
