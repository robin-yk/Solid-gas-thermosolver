# Pilot: apparent TOF range in place of the fixed TOF

SI Note 2a divides every sample's initial CO rate by one fixed site count,
2.31 umol/g (0.17 ML of 13.6 umol/g bridging O). This pilot replaces that
single number with the range of apparent TOF allowed by the published site
energies, the measured vacancy inventory, charge compensation, bulk vacancy
pairing, the surface-reconstruction limit and vacancy transport.

## Result (`outputs/sample_tof_range.csv`)

Core = 40 cases per sample: 900 nm sphere, five sourced energy maps, NEUTRAL
and LOCAL charge closures, four reactive-site definitions.

| Sample | Fixed TOF, SI 2a (s⁻¹) | Core min | Core median | Core max | With all sensitivities |
|---|---:|---:|---:|---:|---:|
| A600 | 0.978 | 0.981 | 2.07 | 5.75 | 0.327–7.35 |
| R500 | 1.63 | 1.64 | 3.41 | 7.26 | 0.545–7.60 |
| R600 | 1.70 | 1.70 | 3.03 | 7.55 | 0.567–7.55 |
| R800 | 0.870 | 0.872 | 1.55 | 3.87 | 0.291–3.87 |
| R1000 | 0.0197 | 0.0197 | 0.0350 | 0.0875 | 0.0066–0.0875 |
| R1100 | 0.0105 | not calculated | | | |

- In every core case for R600, R800 and R1000, and in 36 of 40 for R500 and
  32 of 40 for A600, the equilibrium bridging coverage exceeds 17%. An
  unreconstructed (110) surface cannot hold that much (BIR2024; Note 2a), so
  the reactive count is capped at 0.17 ML.
  - With all bridging vacancies counted, the cap gives the SI 2a denominator:
    2.305 against 2.31 umol/g. The difference is the bridging capacity,
    13.558 here against 13.6 in the manuscript.
  - The SI 2a value is therefore the lower edge of each core range.
- The width of the core range comes from the reactive-site definition. Isolated
  vacancies with z = 2, 4, 8 at the cap give 1.45, 2.1 and 4.4 times the SI 2a
  value.
- Particle size is the largest lever in the other direction. At 300 nm every
  TOF is 3 times lower (−0.48 decades). No sample size distribution was
  measured.
- R1000 values are lower bounds. The sample was treated at 1000 °C, and Yuan
  2024 finds Ti2O3-(1×2) above 900 °C, which removes bridging rows.
- R1100 has no inventory with a source; the v12 workbook marks it BLOCKED.

`outputs/sensitivity_effects.csv` gives each change's effect in decades,
matched case by case to its core counterpart:

| Change | Median | Largest |
|---|---:|---:|
| 600 nm / 300 nm sphere | −0.18 / −0.48 | −0.52 |
| GLOBAL closure, eps 64 / 107 | 0 | +0.26 |
| ZHA2017 bulk pairs, 600 °C or frozen at treatment | 0 | +0.23 |
| Note 2b decay length 0.25 / 1 nm | 0 | ±0.05 |
| Basal (in-plane) vacancies also reactive | 0 | −0.13 |
| R600 inventory 94.6 instead of 94.0 | 0 | 0 |

Most medians are zero because the capped cases do not move. The largest
shifts occur where a change pulls the coverage below 17%.

## What is solved

One free energy over a finite state table, with the measured inventory fixed:

    F = sum x_a E_a + kT sum x_a ln(x_a / (C_d g_a)) + 1/2 q^T H q
    sum v_a x_a = N (measured)

- **NEUTRAL**: electrons implicit.
- **LOCAL**: each trilayer (discrete maps) or shell (continuum) is neutral,
  with its own Ti³⁺ pool and ideal Ti entropy. This is the manuscript Note 2b
  convention.
- **GLOBAL**: only the particle is neutral. Every atomic plane is a charged
  shell: vacancies (+2) sit on O planes and Ti³⁺ (−1) on Ti planes, and the
  shells interact through H_lm = K / max(r_l, r_m).
  - The inverse of H is tridiagonal: it is the capacitance matrix of nested
    spherical capacitors. Each Newton step therefore costs O(shells).
  - The permittivity is the rutile solid's own value from Parker 1961, Fig. 1
    at 873 K: 64 along a, which is the direction of a field normal to (110),
    and 107 along c.
  - The measured bed ε′ cannot stand in for it, for two reasons.
    - A bed that is 23% solid by volume saturates at 1.89 (Maxwell Garnett)
      or 3.21 (Bruggeman) for any solid permittivity. The measured 3.2–4.0
      therefore does not fix the solid value.
    - The measured ε′ already contains the vacancy electrons that GLOBAL
      places explicitly.
  - Surface-trilayer Ti³⁺ sits 0.2 eV above the subsurface optimum (RET2018).
- **Bulk pairs**: vacancy pairs below 1.30 nm use the ZHA2017 bulk potential
  E(r) = A/r − B/r² − C/r⁶. It is taken at all 267 rutile O–O neighbours
  within its 1 nm range, each with the Mayer weight (n_k/2)(e^(−E_k/kT) − 1),
  which is exact to second order in fugacity.
  - Pairs hold more than half of the eligible bulk vacancies in most cases.
    Those cases are BOUNDARY: larger clusters would lower the surface
    population further, so the TOF is a lower bound.
- **Transport**: every link carries J = Λ(F, R)(μ_i − μ_j)/kT.
  - F and R are the forward and backward site-exclusion hop fluxes, with
    barriers from `specification/transport_edges.csv` and prefactor kT/h.
  - For NEUTRAL this is exact hopping. For LOCAL it gives the ambipolar
    factor 3.
  - The free energy falls monotonically, so the only rest point is the
    equilibrium solver's answer.
- **Reactive sites**: bridging coverage capped at 0.17, then
  - C_bri θ (all bridging vacancies), or
  - C_bri θ(1−θ)^z for isolated vacancies, z = 2, 4, 8.
  - With no pair interaction on the surface this count is exact.

The measured CO rate enters only at the end: TOF = r_CO / N_react.

## Transport gate (`outputs/transport_gate.csv`)

Each case starts from equilibrium at the treatment temperature and is held at
600 °C. The surface coverage reaches within 1% of its 600 °C value in:

| Sample | Iddir 2007, 1.10 eV (lowest bulk barrier) | Uberuaga 2011, 1.50 eV (highest) |
|---|---:|---:|
| R500 | 0.010 s | 2.0 s |
| R800 | 0.0015 s | 0.30 s |
| R1000 | 0.0003 s | 0.054 s |

The shortest observation time on the v12 grid is 1 s (the actual time was not
recorded). Only R500 at 1.50 eV lags at 1 s, by up to 4% in coverage. That
coverage is above the cap, so the TOF does not change. The 600 °C equilibrium is
therefore the state at the rate measurement for every hop-connected
population.

## Checks (`python3 -m pytest tests -q`, 27 tests)

- Every O and Ti site is counted once.
- NEUTRAL agrees with an independent bisection. LOCAL satisfies its
  stationarity condition at every site.
- Manuscript Note 2b is reproduced digit for digit: x(R), the R600 interior,
  the top-2-nm share and its ξ scan, and the 0.86–1.04 ratio.
- The coupled_v1 neutral R600 equilibrium is reproduced to 1e-9. Its 4-site
  open chain overcounts isolated vacancies by (2−θ)/(2(1−θ)).
- Capacitance inverse and potentials are exact against dense algebra.
- GLOBAL converges to Poisson and neutrality, and tends to LOCAL as eps → 0.
- The screening length equals Debye–Hückel from the bulk site statistics
  (0.4315 vs 0.4311 nm).
- Pairs only lower the surface population. Their weights are the Mayer
  second virial. Freezing them at the same temperature reproduces
  equilibrium.
- Transport ends at the equilibrium solver (1e-9) and conserves mass.
- Its slowest radial mode matches exact sphere diffusion, D k₁², including
  the LOCAL ambipolar factor. The free energy falls along every trajectory.
- Regenerating the outputs is byte-identical.

## Not in this model

- Surface vacancy–vacancy energies: none with a source (ZHA2017 is bulk;
  BIR2024 gives pattern statistics, not energies).
- Reconstruction energies: the 0.17 ML cap and the Yuan 900 °C rule are
  limits, not free energies. Reconstructed area is taken as zero, so upper
  TOF values would rise if it is not.
- Cooling rate and exact observation time: not recorded. The v12 grid of
  1–600 s is used.
- Particle size distribution and facet fractions: not measured.
- Kinetics of pair formation: not sourced. Pairs are either equilibrated at
  600 °C or frozen at the treatment temperature.
- Registry rows kept as evidence but unused: Matsunaga, V8 (−0.665,
  M10 code only), the 0.271 eV pairing (Note 2c, bulk, no geometry), aggregate
  capacity fractions (M10 convention), cooling-rate grid, unmatched BET areas.

```
python3 run.py      # about 1 minute
```
