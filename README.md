# solid-gas-thermosolver

Equilibrium solver for a gas phase in contact with a reducible oxide, built to
answer one question: **can a reverse water-gas-shift feed reduce rutile TiO2 to
the Magneli phase Ti4O7?**

The answer is no, and the margin is not close. Getting there needs methane.

## What it does

Total Gibbs energy is minimised over the mole vector of

| phase | species |
|---|---|
| gas (ideal, mixing) | CO2, H2, CO, H2O, CH4 |
| condensed (pure, unit activity) | TiO2 (rutile), Ti4O7 (Magneli) |

subject to elemental balance on C, H, O and Ti. SLSQP is run from several
starting points at each temperature, because the reduced-solid branch and the
untouched-rutile branch are separate local minima and a single start lands in
whichever one it was seeded near.

Thermodynamic data is NIST-JANAF (4th ed.) Shomate polynomials throughout.

## Results

Conditions in both runs: 100 mg TiO2, 25 mL, 1 atm, 500-1500 C.

### RWGS feed, CO2/H2 = 1:1

Ti(3+) is **0.0000% at every temperature**. Rutile is untouched across the whole
range. The gas equilibrates to ordinary RWGS behaviour — CO2 conversion climbs
from 25% at 500 C to 66% at 1500 C, with CH4 significant only below ~600 C.

```
  T(C)   CO2%    H2%    CO%   H2O%    CH4%  conv%   Ti3+%     Kp      Q
   500  44.39  15.92   5.61  24.59   9.492   25.4   0.000  0.195  0.195
   804  25.43  25.40  24.57  24.59   0.009   49.2   0.000  0.936  0.936
  1006  21.78  21.78  28.22  28.22   0.000   56.4   0.000  1.680  1.680
  1500  17.10  17.10  32.90  32.90   0.000   65.8   0.000  3.703  3.703
```

`Kp` is computed from standard potentials and `Q` from the converged mole
fractions. They agree to three decimals at every temperature, which is the
check that the minimiser actually found equilibrium rather than a nearby point.

### Reducing feeds

| feed | Ti(3+) at 900 C | Ti(3+) at 1500 C |
|---|---|---|
| CO2/H2 = 1:1 | 0.000% | 0.000% |
| pure H2 | 0.386% | 8.01% |
| pure CH4 | 41.2% | 46.2% |

Methane opens the Magneli branch roughly 500 C below hydrogen. The reason is
visible in the product gas: CH4 cracks and the carbon leaves as **CO**, not CO2,
and removing oxygen as CO becomes more favourable the hotter it gets, while
removing it as H2O does not. Under CH4 the equilibrium gas at 900 C is 65% H2 /
33% CO — the reduction is being driven by carbon, with hydrogen along for the
ride.

The CH4 curve is not monotonic: Ti(3+) peaks at 41.2% near 905 C, dips to 39.2%
by 1108 C, then climbs again to 46.2% at 1500 C. The dip is where CH4 runs out
(0.92% remaining at 905 C, 0.06% at 1006 C); past that point the reduction is
limited by the CO/CO2 ratio the residual gas can sustain rather than by
available carbon.

## Layout

```
solidgas/
  shomate.py      NIST-JANAF coefficients; H, S, Cp, mu0, Kp
  equilibrium.py  species and element matrix, G_total, constrained solve
run_Ti4O7_RWGS.py       the CO2/H2 sweep, writes Ti4O7_RWGS.png
run_Ti4O7_reducing.py   the pure-H2 and pure-CH4 sweeps
tests/                  data checks and solver checks
```

## Running

```bash
pip install -r requirements.txt
python3 run_Ti4O7_RWGS.py        # ~40 s
python3 run_Ti4O7_reducing.py    # ~7 min
python3 -m pytest tests/ -q
```

The reducing sweep is the slow one: 80 temperatures x 4 starting points x 2
feeds, each an SLSQP solve.

## Assumptions, and where they bite

- **Condensed phases are pure and immiscible.** TiO2 and Ti4O7 are treated as
  separate stoichiometric solids with unit activity. Real Magneli formation goes
  through a homologous series (Ti_nO_2n-1, n = 4..10) with mutual solubility, so
  the true reduction onset is gradual where this model puts a sharp branch.
- **Only Ti4O7 represents the reduced phase.** Ti3O5, Ti2O3 and the other
  Magneli members are absent from the species list. Their inclusion would lower
  the onset temperature further.
- **Equilibrium only.** Nothing here says how fast any of it happens. The CH4
  results in particular assume carbon deposition equilibrates rather than
  passivating the surface, which is not what happens in a real reactor.
- **The Shomate branch switches at 1000 K for every species**, while the JANAF
  breakpoints differ per species (CO2 at 1200 K, CH4 and CO at 1300 K). The
  resulting discontinuity in mu0 is at most 0.027 kJ/mol, on CH4 — three orders
  of magnitude below RT at that temperature, so it does not move any result.
- **Range.** The Ti4O7 fit is valid to 1950 K (1677 C), so the 1500 C ceiling is
  inside it. The H2O entry uses the 500-1700 K fit and should not be evaluated
  at room temperature.
- **Data note.** The rutile fit reproduces Cp(298) = 55.20 J/mol/K against a
  literature 55.06, but its S(298) comes out 49.86 against a literature
  ~50.3-50.6. That offset is at the very bottom of the fitted range and well
  below the 773 K floor used here.
