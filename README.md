# C45 Romanoff upper-density certificate

This repository contains reproducibility code for the C45 Romanoff upper-density certificate

```text
overline d(A) < 0.490249407811155.
```

Here

```text
A = {n odd : n = p + 2^k for some prime p and k >= 0}.
```

The verifier is split into two stages.

1. `generate_seed_histogram.cpp` recomputes the 13-prime seed histogram from the residue conditions. The seed histogram is not hard-coded.
2. `verify_c45_romanoff_upper_bound.py` reads the generated seed histogram, applies the exact coprime-order cluster updates for the remaining primes, and checks the final rational upper-density inequality.

Run everything with:

```bash
./verify_all.sh
```

This command compiles `generate_seed_histogram.cpp`, writes
`seed_histogram_13_prime.txt`, runs exact C++ and Python self-tests, and then
runs the Python verifier. It also writes `certificate_summary.json`, a
machine-readable summary containing the exact rational upper sum numerator and
denominator plus the generated seed-histogram SHA-256 digest. These generated
files are intentionally ignored by git; rerun the command to reproduce them.

The scripts use exact integer and rational arithmetic for the proof. Decimal arithmetic is only used for display.

The self-tests are intentionally small and independent of the final
36-prime computation:

- the C++ seed-histogram dynamic program is checked against direct enumeration
  on a small prime set;
- the Python cluster-polynomial update is checked against brute-force residue
  enumeration for small clusters;
- the Python verifier tests exercise histogram parsing and exact helper
  invariants.

Expected final output includes:

```text
36-prime Romanoff upper-density certificate
number of primes = 36
M = 77647987881031766653638954957516345525837812179564689471444926820062741622939917219985
phi(M) = 20439232269111000863616155005739376004784767346781800757207096895187762282496000000000
T = 116570053844283433719485160
seed histogram sha256 = f35f2f3d0a14c04ecc0a8f65a0b2e8b469f358a344b3e0038bee83e1e65c7d08
seed histogram entries = 2104
added multiplier entries = 302400
rational upper sum = 0.4902494078111542556868594292410419912779819272368602249728959220606257217212393446088018129903300560
proved: upper density < 0.490249407811155
improvement over current displayed upper bound = 0.000091681047089
```

## Mathematical outline

Let `P` be the 36-prime obstruction set. Put

```text
M = product_{q in P} q,
T = lcm_{q in P} ord_q(2).
```

For each residue class `a mod M`, define

```text
F(a) = {r mod T : gcd(a - 2^r, M) = 1},
nu(a) = |F(a)|,
delta(nu) = #{a mod M : nu(a) = nu}.
```

The obstruction argument gives

```text
overline d(A) <= sum_nu delta(nu) min(1/(2M), nu/(T phi(M) log 2)).
```

The 13-prime seed histogram is computed from first principles. The added prime clusters have orders coprime to the base period and to each other, so their cluster polynomials multiply the seed histogram exactly by the Chinese remainder theorem.

The final inequality uses the first 40 terms of

```text
log(2) = 2 sum_{j >= 0} 1 / ((2j+1) 3^(2j+1))
```

as a rational lower bound for `log(2)`. Its reciprocal is therefore a rigorous rational upper bound for `1/log(2)`.
