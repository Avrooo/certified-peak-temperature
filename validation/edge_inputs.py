"""Outcomes for threshold equal to the peak or within the arithmetic resolution (Section 5.3, degenerate and extreme input)."""
import solver, e4
for n, b, L in (("B1=B2=0, L=0", ((0.0, 0.0), (0.0, 0.0)), 0.0), ("B1=B2=0, L=1e-30", ((0.0, 0.0), (0.0, 0.0)), 1e-30), ("B1=B2=0, L=0.5", ((0.0, 0.0), (0.0, 0.0)), 0.5)):
    print(n, solver.decide(b, 0.05, 0.2, L)[0])
b = ((3.0, 3.0), (1.0, 1.0)); pk = max(e4.endpoints_float(3.0, 1.0, 0.05, 0.2))
for off in (0.0, 1e-15, -1e-15, 1e-10, -1e-10):
    print("point box (3,1), threshold = float peak %+.0e ->" % off, solver.decide(b, 0.05, 0.2, pk+off)[0])
