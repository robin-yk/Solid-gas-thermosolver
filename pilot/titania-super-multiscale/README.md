# Pilot: apparent TOF range in place of the fixed TOF

SI Note 2a divides every sample's initial CO rate by one fixed site count,
2.31 umol/g. This pilot replaces that single number with the range of apparent
TOF that the published site energies, the measured vacancy inventory and the
particle size allow.

It is an equilibrium model only: fixed measured inventory, common temperature
600 °C, no transport, no thermal history.

## What is solved

A particle is a set of domains. Each domain is a group of identical cells, and
each cell can be in one of a small number of states. A state carries an energy
E, a vacancy count v and a Ti3+ count e. The populations x minimise one free
energy:

    F = sum x E + kT sum x ln(x / (C g))
    subject to   sum v x = N                  (measured inventory)
                 sum (2v - e) x = 0 per region  (LOCAL closure only)

The minimum has Gibbs form. The only unknowns are therefore the vacancy
chemical potential mu and one charge multiplier per region. Region charge
rises monotonically with its multiplier, and inventory rises monotonically
with mu once every region is neutral. Both are solved as bracketed scalar
roots, so the solver cannot diverge.

The measured CO rate enters only at the end:

    TOF = r_CO / N_react

## Scenarios (168 per sample)

| Dimension | Values | Source |
|---|---|---|
| Particle diameter | 900, 600, 300 nm | 900: manuscript; 600, 300: v12 Q1 sensitivity |
| Site energies | Pabisiak L1–4, Hameeuw L1–3, Li on L1 SBR, Li on all of L2, Li continuum (xi 0.25 / 0.5 / 1 nm) | `specification/energy_sets.csv`, manuscript Note 2b |
| Charge closure | NEUTRAL; LOCAL (each trilayer, or each continuum shell, neutral with Ti3+ entropy) | manuscript Note 2b for LOCAL |
| Reactive sites | all bridging vacancies; isolated bridging vacancies with z = 2, 4, 8 | manuscript Note 2b |

Site capacities come from rutile(110): per 1×1 cell and per trilayer, 1
bridging, 2 in-plane and 1 sub-bridging O plus 2 Ti. The bulk takes the
remainder, so every O and Ti of the particle is counted once.

Isolated vacancies use C θ(1−θ)^z. With no pair interaction in the energies,
occupancies are independent, and this count is exact for an infinite row (a
test enumerates a 12-site ring to show it).

A case whose denominator is below 0.01 umol/g or below 1 % of the bridging
capacity is marked INVALID (v12 boundary rule 2) and left out of the range.

## Result (`outputs/sample_tof_range.csv`)

| Sample | Fixed TOF, SI 2a (s⁻¹) | TOF min | median | max | valid / 168 |
|---|---:|---:|---:|---:|---:|
| A600 | 0.978 | 0.0680 | 1.33 | 6.72 | 133 |
| R500 | 1.63 | 0.0940 | 1.39 | 9.04 | 101 |
| R600 | 1.70 | 0.0964 | 1.18 | 27.1 | 90 |
| R800 | 0.870 | 0.0494 | 0.614 | 14.3 | 87 |
| R1000 | 0.0197 | 0.00112 | 0.00795 | 0.245 | 75 |
| R1100 | 0.0105 | not calculated | | | |

R1100 has no inventory with a recorded source; the v12 workbook marks it
BLOCKED.

`outputs/scenario_tof.csv` holds every case with its bridging coverage, its
denominator and its solver residuals.

## Checks (`python3 -m pytest tests -q`)

- Every O and Ti site is counted once (relative error below 1e-13).
- NEUTRAL solutions agree with an independent bisection on mu.
- LOCAL solutions satisfy eps + kT ln θ/(1−θ) + 2kT ln y/(1−y) = mu at every site.
- Manuscript Note 2b is reproduced digit for digit:
  - x(R) = 10.13, 17.59, 22.49, 24.1, 24.9 %;
  - R600 interior 0.00338 and top-2-nm share 0.111;
  - top-2-nm share 0.064 and 0.193 at xi = 0.25 and 1 nm;
  - R1000/R600 isolated ratio 0.86–1.04.
- The coupled_v1 neutral R600 equilibrium is reproduced to 1e-9.
- Regenerating the outputs is byte-identical.

## Not in this model

- Global charge separation: there are no sourced charged-state energies and
  no permittivity.
- Surface aggregates and reconstruction: there are no sourced energies. The
  0.271 eV pairing in Note 2c is a bulk value.
- Transport and thermal history.

The engine accepts any finite state table, so these can be added as states
once their energies have a source.

```
python3 run.py
```
