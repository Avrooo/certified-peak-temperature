# Reproducibility

## Source-to-table and source-to-figure map
Numbers in the paper come from the following code and data. "Re-run" means `scripts/reproduce_results.py` recomputes it and compares with the archive.

| Item in the paper | Code | Archived data | Status |
|:--|:--|:--|:--|
| Table 5, baselines vs reference (1,600 cases) | `src/e4.py` (`fd_peak`, `full_history_peak`, `endpoints_float`); `tables/make_tables.py` | `benchmarks/cases/e4_G*.jsonl` | recomputed from archived original-run results (two-process timings) |
| Table 6, theorem-dependence ablation | `ablation/ablate.py` | `ablation/ablate.json` | re-run |
| Table 7, parameter-sampling study; peak underestimation study (Section 5.2) | `comparisons/uncertainty_sampling_studies.py` | `comparisons/uncertainty_sampling_archived.log` | re-run (log compared except timing line) |
| Table 8 and Figure 3, tolerance floor | `ablation/floorrun.py`, `ablation/e4_floor.py`, `ablation/figs_cert.py` | `ablation/floorrun.json` | re-run |
| Table 9, certified implementations | `benchmarks/scripts/e5final.py`; `tables/make_tables.py` | `benchmarks/results/e5final_main.json` | decisions and widths re-run and identical; timings archived, re-measured separately |
| Figure 2, 1,600-case decisions by gap | `figures/figs_hires.py` | `benchmarks/results/e5final_main.json` | re-run on the new benchmark output |
| Table 10, finite-volume cross-check | `comparisons/fv16.py` | `comparisons/fv16.json` | re-run |
| Table 11, Figures 4 and 5, Scenarios A and B | `examples/scenarios.py`, `figures/figs_hires.py` | (computed) | re-run |
| Figure 1 | `figures/figs_hires.py` | (computed) | re-run |
| Section 5.3: comparison of enclosures with the 60-digit reference | `benchmarks/scripts/e5verify.py` | `benchmarks/results/e5verify.log` | re-run |
| Section 5.3: zero modes, degenerate input, reduced precision | `validation/e5edge.py` | `validation/logs/e5edge.log` | re-run (the extreme-input outcomes of the equality/resolution cases are in `validation/edge_inputs.py` output, see below) |
| Appendix A: 610 transcendental checks, 42 root brackets | `validation/kval.py` | printed counts | re-run |
| Appendix A: mutation tests | `validation/audit.py`, `validation/bcval.py` | printed counts | re-run |
| Appendix A: 120-digit `bc` oracle, 247 comparisons | `validation/bcval.py` | printed counts | re-run |
| Remark 4.1: Lemma 4.1 numerical check (1,500 sets) | `validation/lemmaG_reaudit.py` | `validation/logs/lemmaG_reaudit.log` | re-run (the script's internal labels G1..G5 are the paper's (L1)..(L5)) |
| Section 3 numerical checks of the proof chain | `validation/verify_core.py`, `validation/tp2hp.py` | printed counts | re-run |
| Table 7 of the optimization ladder (first certified vs optimized) timing | `benchmarks/scripts/e5.py`, `benchmarks/scripts/e5run.sh` | `benchmarks/results/e5_V*.json` | scripts and archived results included; timing-based, so not part of the automatic comparison |

## Original run versus re-run
* The decisions and enclosure widths of the 1,600-case regression are deterministic and reproduce bit for bit.
* Timings are machine dependent. The paper's timings were measured on two virtual cores of an Intel Xeon 2.10 GHz (Python 3.13, NumPy 2.5.3, SciPy 1.18.1). The re-run reports its own medians next to the archived ones and does not compare them.
* The baselines' times in Table 5 and the float row of Table 9 come from the original runs, which executed two processes at once; the certified implementations' times in Table 9 are single-process, cold start.

## What the tests establish
They support the implementation. They do not prove Theorem 3.5, Lemma 4.1 or the correctness of the interval kernel; those rest on the proofs in the paper and on independent mathematical reading.

## Files kept unchanged
`frozen/` holds the original 256-bit implementation, the written benchmark design and its SHA-256 sums (`frozen/e4_freeze.sha256`); `src/default_V7.sha256` lists the optimized files. `scripts/verify_archive.py` checks both. `ablation/e4_floor.py` is a modified copy of `src/e4.py` in which the 1e-13 tolerance floor is a variable; it is not a fingerprinted file.

## Not included
Working drafts, abandoned algorithms and one-off diagnostics from development.
