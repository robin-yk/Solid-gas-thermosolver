# Pilot: apparent TOF range from one model

SI Note 2a divides every sample's initial CO rate by one fixed site count,
2.31 umol/g. This pilot replaces that number with a range. One model
distributes each sample's measured vacancy inventory over every site a
vacancy can occupy; the unknowns of that model are its parameters, and the
model is solved at every combination of them. The range of rate / reactive
sites over those solutions is the apparent TOF range.

## Result (`outputs/sample_tof_range.csv`)

720 parameter points x 5 reactive-site definitions = 3600 cases per sample.
The range is over the cases whose reactive-site count is countable (at least
0.01 umol/g and 1% of the bridging capacity).

| Sample | Inventory (umol/g) | Fixed TOF, SI 2a (s⁻¹) | TOF min | TOF median | TOF max | Countable |
|---|---:|---:|---:|---:|---:|---:|
| A600 | 13.4 | 0.978 | 1.2 | 5.71 | 33.3 | 1670 / 3600 |
| R500 | 36.9 | 1.63 | 1.37 | 8.58 | 54.3 | 1744 / 3600 |
| R600 | 94 | 1.7 | 1 | 8.15 | 53.1 | 1793 / 3600 |
| R800 | 193 | 0.87 | 0.382 | 4.31 | 29.4 | 1798 / 3600 |
| R1000 | 781 | 0.0197 | 0.00518 | 0.101 | 0.662 | 1900 / 3600 |
| R1100 | 1.34e+03 | 0.0105 | 0.00234 | 0.0489 | 0.35 | 1766 / 3600 |

| Sample | Q min | Q max | Bridging vacancies, median (umol/g) | Bulk vacancies, median (umol/g) |
|---|---:|---:|---:|---:|
| A600 | 0.132 | 1.23 | 0.0613 | 13.4 |
| R500 | 0.321 | 1.4 | 0.0776 | 36.7 |
| R600 | 1 | 1 | 0.0991 | 93.8 |
| R800 | 0.356 | 1.57 | 0.107 | 193 |
| R1000 | 0.00296 | 0.0748 | 0.158 | 781 |
| R1100 | 0.00117 | 0.0534 | 0.145 | 1.34e+03 |

Q = TOF / TOF(R600) at the same parameter point and definition.

- The bulk holds almost the whole inventory: median share 99.4% (A600),
  99.8% (R600), 99.9% (R1100) of the vacancies. The surface holds about
  0.06 to 0.16 umol/g of bridging vacancies at the median.
- The fixed SI 2a value lies inside the range for R500 to R1100, in its
  lower part. For A600 it is below the range: 0.978 against a minimum of 1.20 s⁻¹.
- Cases below the site threshold are kept in `outputs/model_cases.csv` and
  reported apart. They come from the aggregate cutoff. Each cutoff value has 3600
  cases over the six samples; below threshold are 3600 at 0.6 nm, 3600 at
  1.0 nm, 2480 at 0.45 nm, and 287 to 620 at 0.28 to 0.40 nm.
- R1000 and R1100 were treated at 1000 and 1100 °C and are flagged
  YUAN_STRONG_REDUCTION: Yuan 2024 finds Ti2O3-(1×2) above 900 °C.
- R1100 inputs: inventory 1343 umol/g and rate 0.02433 umol/g/s (v12 record,
  Origin extraction). Manuscript 092426 Fig. 2b plots R1100 at 1336 umol/g
  read off the figure; the known points read back within 3%. Treatment 5% H2,
  1100 °C, 1 h (Fig. 2 and Fig. S5 captions).

## The model (`tofrange/model.py`)

Every part below is on at once, at equilibrium at 873.15 K with the
measured inventory fixed.

- **Surface.** A share f110 of the surface is (110), as four explicit
  O–Ti2O2–O trilayers (1.30 nm) with bridging (BRI), in-plane (IPL) and
  sub-bridging (SBR) O. The rest is bulk-like and has no reactive sites.
- **Reconstruction.** Bridging rows form (1×2) cells that own their four
  layer-1 Ti. A cell is intact, has one or both bridging O missing with 0 to
  4 Ti³⁺, or is reconstructed (energy E + dG, one vacancy, two internal Ti³⁺
  and up to two free).
- **Bulk.** Below 1.30 nm, O sites are grouped into Monte Carlo boxes of 700
  sites. A box holding N vacancies has free energy −kT ln Q(N) with the
  ZHA2017 pair energy E(r) = A/r − B/r² − C/r⁶ up to the cutoff; the chain is
  grown one vacancy at a time with Widom insertion. The cache is in
  `outputs/mc` (not committed).
- **Charge.** Each vacancy leaves two electrons as Ti³⁺. Every atomic plane
  is a charged shell and the shells interact through the Poisson equation
  with the rutile permittivity (Parker 1961 at 873 K: 64 along a, 107 along
  c). Layer-1 Ti³⁺ sits 0.2 eV above the subsurface optimum (RET2018). Only
  the particle is neutral.
- **Size.** Equal-mass mixture of 900 to 1600 nm at four diameters (eight
  change the size factor by 1%).
- **Reactive sites.** C θ (every bridging vacancy), C θ(1−θ)^z with z = 2, 4, 8
  (isolated bridging vacancies), or bridging plus layer-1 in-plane
  vacancies. There is no coverage cap; reconstruction is in the model.

The measured CO rate enters only at the end: TOF = r_CO / N_react.

## Parameters solved over

| Unknown | Values |
|---|---|
| Vacancy energy map | PAB, HAM, LI_SBR1, LI_L2 (`specification/energy_sets.csv`) |
| Aggregate cutoff | 0.28, 0.34, 0.40, 0.45, 0.6, 1.0 nm (1.0 = ZHA2017 range) |
| Reconstruction dG | −0.4, −0.2, 0, +0.2, +0.4 eV |
| (110) share | 1.0, 0.75, 0.5 |
| Permittivity | 64 (a), 107 (c) |

## Checks (`python3 -m pytest tests -q`, 13 tests)

- Every O and Ti site is counted once.
- Aggregate Monte Carlo equals exact enumeration for one and two vacancies;
  at N = 170 two seeds agree within 0.05 eV per vacancy.
- The capacitance inverse is exact against dense algebra.
- The split-shell gradient and Hessian agree with finite differences.
- Solutions conserve the inventory and charge and satisfy Poisson; the pools
  add up to the inventory.
- At large permittivity the solution tends to the scalar μ, μ_e limit as
  1/ε (below 1e-7 at ε = 1e14).
- The screening length equals Debye–Hückel from the site statistics.
- Reconstruction keeps its own Ti and vanishes at dG = 5 eV.
- Q equals the TOF ratio to R600 at the same point.
- Saved rows regenerate for two parameter points.

## Assumptions that stand in for missing values

Each is an input in `specification/parameter_registry.csv`.

| Item | Assumption | Why |
|---|---|---|
| Aggregates of 3+ | pair energies add; cutoff scanned | only pair energies are sourced |
| Aggregate box | 700 O sites, boxes independent | smallest box wider than twice the cutoff |
| Reconstruction energy | dG from −0.4 to +0.4 eV | no sourced energy; the fraction is an output |
| Site threshold | 0.01 umol/g and 1% of bridging capacity | reporting rule; cases below are kept |
| Size | 900–1600 nm, equal mass | distribution shape not given |
| Other facets | bulk-like energy, no reactive sites | no energies or reactivity for them |
| Transport | complete at 600 °C | observation time not recorded |
| R1100 rate | 0.02433 umol/g/s (v12) | not in manuscript 092426 Fig. 2a or Table S1 |

```
python3 run.py      # about 44 minutes on four cores
```
