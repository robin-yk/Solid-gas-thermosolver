# Numerical verification

## Thermodynamic oracle

`scripts/oracle_tio.py` is an independent mpmath implementation with 90 working digits and 80 reported digits. It shares thermodynamic coefficient tables with production but does not import the production Gibbs solver. Each candidate active set is solved as an element-potential root problem and classified by KKT conditions.

The committed reference is `data/reference_results_high_precision.json`. It covers six gas cases, six temperatures, and boundary probes. Production agrees with the reference to the following measured worst deviations:

| Quantity | Worst observed deviation |
| --- | ---: |
| Gas mole fraction | 2e-13 |
| Inactive-phase reduced cost | 8e-12 kJ |
| Browser active-set result | 6e-15 |

Trace solid amounts down to about 1e-41 mol are compared in logarithmic form.

## Defect oracle

`scripts/oracle_statmech.py` runs at 50 digits and shares only `data/rutile_dft.json` with production. It certifies:

- analytic three-class site-exclusion matching, phase-aware, with the closed-form boundary ladder
- the declared-shell-thickness ladder (one to four oxygen layers), which fixes the shell/bulk split in the three-class engine
- exact finite-lattice canonical ensembles from generating polynomials
- exact enumeration of small interacting lattices
- the chemical potential estimated by Widom insertion

The committed reference is `data/reference_statmech.json`. Same-seed Python and JavaScript Monte Carlo trajectories are compared as integer states. Class occupancies and histograms must match exactly; the vacancy chemical potential uses a 1e-9 tolerance.

## Release gates

The test suite checks:

- elemental balance and KKT feasibility
- phase assemblage and inactive-phase stability
- Python-to-JavaScript parity
- high-precision oracle agreement
- deterministic Monte Carlo trajectories
- data export freshness
- generated-page freshness
- absence of external runtime requests
- required model disclosures and validity warnings

Run the complete gate with:

```bash
python3 -m pytest tests/ -q
```

Generated reference files should be updated only by their commands:

```bash
python3 scripts/oracle_tio.py
python3 scripts/oracle_statmech.py
python3 scripts/export_activeset.py
python3 scripts/export_thermo.py
python3 scripts/export_oxide.py
```
