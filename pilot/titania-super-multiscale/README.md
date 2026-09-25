# Pilot: apparent TOF range in place of the fixed TOF

SI Note 2a divides every sample's initial CO rate by one fixed site count,
2.31 umol/g (0.17 ML of 13.6 umol/g bridging O). This pilot replaces that
single number with the range of apparent TOF allowed by:

- the published site energies and the measured vacancy inventory;
- charge compensation;
- bulk aggregates (pairwise-additive pair energies, Monte Carlo over a periodic box);
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
| R1100 | 0.0105 | 0.0106 | 0.0188 | 0.0469 |

- In every core case for R600, R800, R1000 and R1100, and in 36 of 40 for R500 and
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
- R1000 and R1100 values are lower bounds. The samples were treated at 1000
  and 1100 °C, and Yuan 2024 finds Ti2O3-(1×2) above 900 °C, which removes
  bridging rows.
- R1100 inputs: inventory 1343 umol/g and rate 0.02433 umol/g/s (v12 record,
  Origin extraction). Manuscript 092426 Fig. 2b plots R1100 at 1336 umol/g
  read off the figure; the known points read back within 3%. Treatment 5% H2,
  1100 °C, 1 h (Fig. 2 and Fig. S5 captions).

### Every family (`outputs/cases.csv`)

TOF range (s⁻¹) over cases at or above the site threshold, with their count
where some fall below it. The threshold is a reporting rule of this pilot:
fewer than 0.01 umol/g, or fewer than 1% of the bridging capacity, reactive.
Cases below it (class BELOW_SITE_THRESHOLD) keep their TOF = rate / sites and
Q in `cases.csv`, and their range is reported apart in
`sample_tof_range.csv`. They stay out of the ranges below because a near-zero
site count sets a TOF of up to 1e48 s⁻¹. Such a value means the assumed
distribution needs a very large rate per site to give the measured rate.

| Family | A600 | R500 | R600 | R800 | R1000 | R1100 |
|---|---|---|---|---|---|---|
| Core | 0.981–5.75 | 1.64–7.26 | 1.70–7.55 | 0.872–3.87 | 0.0197–0.0875 | 0.0106–0.0469 |
| Diameter 1250, 1600 nm | 1.36–10.1 | 2.27–12.9 | 2.36–13.4 | 1.21–6.88 | 0.0274–0.155 | 0.0147–0.0833 |
| Equal mass, 900–1600 nm | 1.31–7.65 | 2.19–9.74 | 2.28–10.1 | 1.17–5.19 | 0.0264–0.117 | 0.0142–0.0628 |
| Reconstruction, fixed area 0.25/0.5/0.75 | 0.223–16.4 (70/120) | 0.371–26.1 (58/120) | 0.386–24.1 (50/120) | 0.198–11.2 (44/120) | 0.00446–0.266 (41/120) | 0.00239–0.143 (40/120) |
| Reconstruction, (1×2) state, dG −0.4…+0.4 eV | 0.171–12.6 (102/160) | 0.278–25.9 (79/160) | 0.289–27.7 (66/160) | 0.148–14.4 (60/160) | 0.00335–0.321 (47/160) | 0.00179–0.127 (44/160) |
| Aggregates, MC, cutoff 0.28 nm | 0.981–8.45 | 1.64–8.40 | 1.70–7.55 | 0.872–3.87 | 0.0197–0.0875 | 0.0106–0.0469 |
| Aggregates, MC, cutoff 0.34 nm | 0.981–6.53 (36/40) | 1.64–17.2 | 1.70–10.9 | 0.872–4.16 | 0.0197–0.0875 | 0.0106–0.0469 |
| Aggregates, MC, cutoff 0.40 nm | 0.981–11.6 (36/40) | 1.64–13.2 (36/40) | 1.70–19.5 | 0.872–6.58 | 0.0197–0.0875 | 0.0106–0.0469 |
| Aggregates, MC, cutoff 0.45 nm | none (0/40) | none (0/40) | 23.0–25.5 (8/40) | 7.25–12.5 (12/40) | 0.0619–0.187 (16/40) | 0.0220–0.0731 (16/40) |
| Aggregates, MC, cutoff 0.6 or 1 nm | none (0/80) | none (0/80) | none (0/80) | none (0/80) | none (0/80) | none (0/80) |
| Combined: aggregates 0.34/0.40 nm + (1×2) state, dG −0.2/0/+0.2 eV + equal mass 900–1600 nm | 0.226–15.3 (115/192) | 0.375–24.1 (127/192) | 0.390–26.2 (139/192) | 0.200–14.0 (132/192) | 0.00451–0.237 (110/192) | 0.00242–0.138 (101/192) |
| (110) share 0.75 / 0.5 | 1.31–11.3 | 2.18–14.5 | 2.27–15.1 | 1.16–7.74 | 0.0263–0.175 | 0.0141–0.0937 |
| Li maps, Matsunaga basal +0.11 eV | 0.273–4.35 | 0.146–7.26 | 0.133–7.55 | 0.0684–3.87 | 0.00154–0.0875 | 0.000827–0.0469 |
| GLOBAL, eps 64 / 107 | 0.981–5.37 | 1.64–7.26 | 1.70–7.55 | 0.872–3.87 | 0.0197–0.0875 | 0.0106–0.0469 |
| Bulk pairs (dilute), 600 °C or frozen | 0.981–7.35 | 1.64–7.60 | 1.70–7.55 | 0.872–3.87 | 0.0197–0.0875 | 0.0106–0.0469 |

- **Size.** The samples span 900–1600 nm. Every TOF rises with diameter,
  because bridging capacity scales as 1/D: 1.39× at 1250 nm, 1.78× at
  1600 nm, 1.34× for equal mass over the range.
- **Reconstruction, explicit.** No 17% cap is applied here.
  - Without the cap, saturated surfaces leave almost no isolated vacancies.
    Most cases below the site threshold are isolated-site definitions on
    such surfaces.
  - The computed reconstructed fraction at R600 (LOCAL, discrete maps) is
    0.012–0.088 at dG = 0 and 0.66–0.76 at dG = −0.4 eV. At R1000 it is
    0.024–0.13 at dG = 0, and at R1100 0.014–0.15.
  - In LOCAL each 1×2 cell carries its own four layer-1 Ti. The reconstructed
    row keeps two of them for its own Ti³⁺, so no Ti is counted in two pools.
    The fixed-area branch removes the same Ti from the free layer-1 pool.
- **Combined.** Aggregates, the (1×2) state and the 900–1600 nm mixture are
  switched on together, so their interaction is computed. The mixture
  averages reactive sites, not TOF, over the 8 diameters. The combined ranges
  lie inside the union of the single-change ranges for every sample.
- **Aggregates.** These are results of one stated model: the ZHA2017 energy
  of two vacancies, summed over every pair within a cutoff. ZHA2017 computed
  two vacancies only; whether the sum holds at high density is not
  established. In this model the free energy of a vacancy in a dense bulk
  aggregate, relative to an isolated bulk vacancy, depends strongly on the
  cutoff:

  | Cutoff (nm) | 0.28 | 0.34 | 0.40 | 0.45 | 0.6 | 1.0 |
  |---|---:|---:|---:|---:|---:|---:|
  | eV per vacancy | −0.54 | −0.75 | −0.86 | −1.46 | −2.3 | −4.3 |

  - The surface bridging energies are −0.81 (PAB), −1.21 (HAM) and −1.31 eV
    (Li). At 0.28 and 0.34 nm the aggregate value lies above all three; at
    0.40 nm (−0.86 eV) it lies below PAB. Even so, at 0.40 nm 36 to 40 of the
    40 cases per sample stay above the site threshold. The split is set by
    each state's full partition function and capacity, not by the mean
    energy per vacancy.
  - From 0.45 nm the aggregates take nearly every vacancy, and most cases
    fall below the site threshold. Dividing the measured rate by that small
    count needs a rate per site many decades above the core values.
  - Going from a 0.6 to a 1.0 nm cutoff moves the value by 2 eV (26.6 kT at
    600 °C). The −4.3 eV figure is therefore a property of the 1 nm pairwise
    sum, not a verified rutile aggregate energy.
  - Q(N) gives the box free energy for each vacancy count N. It does not
    report how those vacancies split into clusters.
- **Q (Note 11.2).** Every case also carries Q = TOF / TOF(R600) in the same
  scenario (`Q_vs_R600`).
  - Core Q for R1000 is 0.0116, the rate ratio itself, because every core
    case sits at the 17% cap.
  - Core Q for R1100 is 0.0062.
  - Across all families Q(R1000) spans 0.0019–0.089 and Q(R1100)
    0.00068–0.051. Q values whose R600 reference is below the site threshold
    are marked (`Q_ref_below_threshold`) and left out of these ranges.

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
  - A box holding N vacancies has free energy −kT ln Q(N). At each N,
    canonical Monte Carlo moves vacancies between sites. Every 50 moves the
    Widom insertion weight is summed over every allowed empty site, which
    gives Q(N+1)/Q(N). The chain then grows by one vacancy.
  - The box is then one state table in the same free energy, so clusters of
    every size that fit in the box are in the equilibrium.
  - The Ti rule (two Ti³⁺ per vacancy on its own Ti, one per Ti) only admits
    or rejects a configuration. The energy and number of the possible Ti³⁺
    arrangements are not summed.
- **Reconstruction**:
  - Fixed area: a fraction f of the surface is Ti2O3-(1×2). That area holds
    its vacancies and its own Ti³⁺, and has no reactive sites. Those Ti leave
    the free layer-1 Ti pool.
  - Explicit state: each 1×2 cell is either two bridging sites or
    reconstructed. The reconstructed fraction is then computed. In LOCAL the
    cell also holds its four layer-1 Ti; the reconstructed state keeps two.
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
| R1100 | 0.0001 s | 0.021 s |

The shortest observation time on the v12 grid is 1 s (the actual time was not
recorded). Only R500 at 1.50 eV lags at 1 s, by up to 4% in coverage. That
coverage is above the cap, so the TOF does not change. The 600 °C equilibrium is
therefore the state at the rate measurement for every hop-connected
population.

## Checks (`python3 -m pytest tests -q`, 36 tests)

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
- Aggregate Monte Carlo equals exact enumeration for one and two vacancies.
  At N = 170, two seeds agree within 0.05 eV per vacancy. Other N are not
  compared across seeds.
- Fixed reconstruction conserves the inventory. The reconstructed fraction
  falls monotonically with dG and vanishes as dG → ∞, recovering the
  unreconstructed result to 1e-9. No layer-1 Ti holds two electrons or sits
  in two pools, in either reconstruction branch.
- The size mixture equals the mass average of its diameters.
- A (110) share of 1 reproduces the base model. Q equals the TOF ratio to
  R600 in the same scenario.
- Regenerating the outputs is byte-identical.

## Assumptions that stand in for missing values

Each is an explicit input in `specification/parameter_registry.csv`.

| Item | Assumption | Why |
|---|---|---|
| Aggregates of 3+ | Pair energies add (no many-body term); cutoff scanned 0.28–1.0 nm (1.0 = ZHA2017 range) | only pair energies are sourced |
| Aggregate Ti3+ | two Ti3+ per vacancy on its own three Ti; one electron per Ti | caps local vacancy fraction at 1/4 |
| Aggregate box | 700 O sites, boxes independent | smallest box wider than twice the cutoff |
| Reconstruction, fixed | area fraction 0.25, 0.5, 0.75; 0.5 vacancy per 1×1 cell; no reactive sites; own Ti3+ | Ti2O3 row stoichiometry |
| Reconstruction, state | 1×2 cell energy eps_BRI + dG, dG from −0.4 to +0.4 eV; row keeps 2 of the cell's 4 layer-1 Ti | no sourced energy; the fraction is an output |
| Combined family | cutoff 0.34, 0.40 nm; dG −0.2, 0, +0.2 eV; discrete maps | cutoffs where the surface keeps vacancies; state needs explicit sites |
| Site threshold | 0.01 umol/g and 1% of bridging capacity | reporting rule; cases below are kept and reported apart |
| Size | 900–1600 nm (user); equal mass at 8 diameters; same inventory per gram | distribution shape not given |
| Surface pairs | none | no sourced surface energy; v12 forbids reusing bulk values |
| Facets | (110) share 0.75 or 0.5; the rest is bulk-like (energy 0) and has no reactive sites | no energies or reactivity for other facets |
| Observation time | 1–600 s grid | not recorded |
| R1100 rate | 0.02433 umol/g/s (v12) | not in manuscript 092426 Fig. 2a or Table S1 |

Registry rows kept as evidence but unused: Matsunaga sensitivity branches (M10), V8 (−0.665, M10 code
only), the 0.271 eV pairing (Note 2c, bulk, no geometry), aggregate capacity
fractions (M10 convention), cooling-rate grid, unmatched BET areas.

```
python3 run.py      # about 16 minutes; aggregate Monte Carlo is cached in outputs/mc
```
