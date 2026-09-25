# Pilot: apparent TOF range in place of the fixed TOF

SI Note 2a divides every sample's initial CO rate by one fixed site count,
2.31 umol/g (0.17 ML of 13.6 umol/g bridging O). This pilot replaces that
single number with the range of apparent TOF allowed by:

- the published site energies and the measured vacancy inventory;
- charge compensation;
- bulk aggregates of every size;
- explicit surface reconstruction;
- particle size;
- vacancy transport.

## Result (`outputs/sample_tof_range.csv`)

Core = 40 cases per sample: 900 nm sphere, five sourced energy maps, NEUTRAL
and LOCAL charge closures, four reactive-site definitions.

| Sample | Fixed TOF, SI 2a (s⁻¹) | Core min | Core median | Core max |
|---|---:|---:|---:|---:|
| A600 | 0.978 | 0.981 | 2.07 | 5.75 |
| R500 | 1.63 | 1.64 | 3.41 | 7.26 |
| R600 | 1.70 | 1.70 | 3.03 | 7.55 |
| R800 | 0.870 | 0.872 | 1.55 | 3.87 |
| R1000 | 0.0197 | 0.0197 | 0.0350 | 0.0875 |
| R1100 | 0.0105 | not calculated | | |

- In every core case for R600, R800 and R1000, and in 36 of 40 for R500 and
  32 of 40 for A600, the equilibrium bridging coverage exceeds 17%. An
  unreconstructed (110) surface cannot hold that much (BIR2024; Note 2a), so
  the core caps the reactive coverage at 0.17 ML.
  - This cap is a rule. It does not compute a reconstructed population; the
    two reconstruction families below do.
  - With all bridging vacancies counted, the cap gives the SI 2a denominator:
    2.305 against 2.31 umol/g. The difference is the bridging capacity,
    13.558 here against 13.6 in the manuscript.
  - The SI 2a value is therefore the lower edge of each core range.
- The width of the core range comes from the reactive-site definition. Isolated
  vacancies with z = 2, 4, 8 at the cap give 1.45, 2.1 and 4.4 times the SI 2a
  value.
- R1000 values are lower bounds. The sample was treated at 1000 °C, and Yuan
  2024 finds Ti2O3-(1×2) above 900 °C, which removes bridging rows.
- R1100 has no inventory with a source; the v12 workbook marks it BLOCKED.

### Every family (`outputs/cases.csv`)

TOF range over the valid cases (s⁻¹), with valid/total cases. A case is
INVALID when fewer than 0.01 umol/g, or fewer than 1% of the bridging
capacity, are reactive; its TOF is then undefined.

| Family | A600 | R500 | R600 | R800 | R1000 |
|---|---|---|---|---|---|
| Core | 0.981–5.75 (40/40) | 1.64–7.26 (40/40) | 1.70–7.55 (40/40) | 0.872–3.87 (40/40) | 0.0197–0.0875 (40/40) |
| Diameter 1250, 1600 nm | 1.36–10.1 | 2.27–12.9 | 2.36–13.4 | 1.21–6.88 | 0.0274–0.155 |
| Equal mass, 900–1600 nm | 1.31–7.65 | 2.19–9.74 | 2.28–10.1 | 1.17–5.19 | 0.0264–0.117 |
| Reconstruction, fixed area 0.25/0.5/0.75 | 0.223–16.0 (71/120) | 0.371–27.2 (56/120) | 0.386–24.4 (48/120) | 0.198–11.2 (43/120) | 0.00446–0.266 (40/120) |
| Reconstruction, (1×2) state, dG −0.4…+0.4 eV | 0.171–12.6 (102/160) | 0.278–24.1 (78/160) | 0.289–27.9 (64/160) | 0.148–14.7 (54/160) | 0.00335–0.264 (44/160) |
| Aggregates, MC, cutoff 1 or 0.6 nm | none (0/80) | none (0/80) | none (0/80) | none (0/80) | none (0/80) |
| GLOBAL, eps 64 / 107 | 0.981–5.37 | 1.64–7.26 | 1.70–7.55 | 0.872–3.87 | 0.0197–0.0875 |
| Bulk pairs (dilute), 600 °C or frozen | 0.981–7.35 | 1.64–7.60 | 1.70–7.55 | 0.872–3.87 | 0.0197–0.0875 |

- **Size.** The samples span 900–1600 nm. Every TOF rises with diameter,
  because bridging capacity scales as 1/D: 1.39× at 1250 nm, 1.78× at
  1600 nm, 1.34× for equal mass over the range.
- **Reconstruction, explicit.** No 17% cap is applied here.
  - Without the cap, saturated surfaces leave almost no isolated vacancies.
    Most INVALID cases are isolated-site definitions on such surfaces.
  - The computed reconstructed fraction at R600 (LOCAL, discrete maps) is
    0.02–0.27 at dG = 0 and 0.60–0.73 at dG = −0.4 eV. At R1000 it is
    0.12–0.55 at dG = 0.
- **Aggregates.** With pairwise-additive ZHA2017 energies, a vacancy in a
  dense bulk aggregate has a free energy of −4.3 eV (cutoff 1 nm) or −2.3 eV
  (0.6 nm), relative to an isolated bulk vacancy. Surface sites sit at −0.8
  to −1.3 eV. The aggregates therefore take nearly every vacancy: bridging
  coverage falls below 0.0014, and all 400 aggregate cases are INVALID.
  - Under this assumption the reactive-site count, and so the CO rate, would
    be about zero. The measured rates are not zero.

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
- **Bulk aggregates**: bulk O below 1.30 nm is grouped into boxes of 700
  sites.
  - A box holding N vacancies has free energy −kT ln Q(N). Q(N) comes from
    Monte Carlo with Widom insertion summed over every site, and the chain
    grows one vacancy at a time.
  - The box is then one state table in the same free energy, so clusters of
    every size are in the equilibrium.
- **Reconstruction**:
  - Fixed area: a fraction f of the surface is Ti2O3-(1×2). That area holds
    its vacancies and its own Ti³⁺, and has no reactive sites.
  - Explicit state: each 1×2 cell is either two bridging sites or
    reconstructed. The reconstructed fraction is then computed.
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

## Checks (`python3 -m pytest tests -q`, 33 tests)

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
- Aggregate Monte Carlo equals exact enumeration for one and two vacancies,
  and two seeds agree within 0.05 eV per vacancy.
- Fixed reconstruction conserves the inventory. The reconstructed fraction
  falls monotonically with dG and vanishes as dG → ∞, recovering the
  unreconstructed result.
- The size mixture equals the mass average of its diameters.
- Regenerating the outputs is byte-identical.

## Assumptions that stand in for missing values

Each is an explicit input in `specification/parameter_registry.csv`.

| Item | Assumption | Why |
|---|---|---|
| Aggregates of 3+ | Pair energies add (no many-body term); cutoff 1 nm (ZHA2017 range) or 0.6 nm | only pair energies are sourced |
| Aggregate Ti3+ | two Ti3+ per vacancy on its own three Ti; one electron per Ti | caps local vacancy fraction at 1/4 |
| Aggregate box | 700 O sites, boxes independent | smallest box wider than twice the cutoff |
| Reconstruction, fixed | area fraction 0.25, 0.5, 0.75; 0.5 vacancy per 1×1 cell; no reactive sites; own Ti3+ | Ti2O3 row stoichiometry |
| Reconstruction, state | 1×2 cell energy eps_BRI + dG, dG from −0.4 to +0.4 eV | no sourced energy; the fraction is an output |
| Size | 900–1600 nm (user); equal mass at 8 diameters; same inventory per gram | distribution shape not given |
| Surface pairs | none | no sourced surface energy (ZHA2017 is bulk) |
| Observation time | 1–600 s grid | not recorded |

Registry rows kept as evidence but unused: Matsunaga, V8 (−0.665, M10 code
only), the 0.271 eV pairing (Note 2c, bulk, no geometry), aggregate capacity
fractions (M10 convention), cooling-rate grid, unmatched BET areas.

```
python3 run.py      # about 3 minutes; aggregate Monte Carlo is cached in outputs/mc
```
