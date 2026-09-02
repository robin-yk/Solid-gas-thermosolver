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

## The population model

The browser mirror `web/population.js` integrates the same scalar equation as
`solidgas/vacancy_population.py` with a different integrator (Strang splitting
against Radau), and `tests/test_population_port.py` holds the two to a relative
1e-6 across the reduction series and the parameter search.

## Release gates

The test suite checks:

- elemental balance and KKT feasibility
- phase assemblage and inactive-phase stability
- high-precision oracle agreement
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
python3 scripts/reproduce_paper.py
python3 scripts/build_site.py
```
