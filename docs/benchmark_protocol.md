# Benchmark protocol
* 1,600 synthetic cases: G1 standard (600), G2 stress-test (400), G3 threshold-proximity (400), G4 separately generated (200). Each record stores the seed, the box, `Fs`, `Fp`, the threshold `L`, the group, the reference classification and the separation between `L` and the nearest corner peak.
* The generator is `make_case`/`gen_box` in `src/e4.py` (seeds 1000+, 2000+, 3000+, 9000+). The archived files in `benchmarks/cases/` are the cases used in the paper; do not regenerate them.
* Reference classification: 60-digit decimal evaluator (`src/e4ref.py`) with its own root finder and Taylor sums, sharing the closed-form series and the Lemma 4.1 truncation with the certified solver. It is a numerical reference, not an independent proof.
* Scoring: reference-consistent, reference disagreement, inconclusive. Invalid input raises an error and is never counted as inconclusive.
* Tolerance schedule 1e-3, 1e-5, ..., 1e-13; the optimized solver uses 96 fixed-point bits and 24 Taylor terms, and a 256-bit repeat on unexpected failure (never triggered in the final regression).
* The benchmark design and fingerprints were written before results were seen (`frozen/e4_spec.md`, `frozen/e4_freeze.sha256`).
