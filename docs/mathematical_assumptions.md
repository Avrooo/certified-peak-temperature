# Mathematical assumptions
* Homogeneous slab of unit dimensionless thickness, constant properties, constant ambient and initial temperature.
* `u = 1 - θ` solves `u_t = u_xx`, `u_x(0,t) = 0`, `u_x(1,t) + B(t) u(1,t) = 0`, `u(x,0) = 1`.
* `B(t) = B1` on `(0,t_s)`, `B2` on `(t_s,t_p]`, with `B1, B2 >= 0` and `0 < t_s < t_p`. In the code `t_s`, `t_p` are `Fs`, `Fp`.
* Only `B1` and `B2` may be uncertain (rectangular box). Uncertain switch time, end time, diffusivity, thickness, ambient or initial temperature are not covered.
* Exactly one switch. Three or more stages are not covered by the peak theorem.
* Decision convention: `SAFE_ALL` = peak strictly below the threshold; `EXCEEDS_ALL` = peak at or above it; a decision is returned only if the enclosures establish it.
* The theorems (peak location and corner reduction) and the truncation bound are proved in the paper; the software certificate is conditional on them.
