# Results

## Conditions

Unless stated otherwise, calculations use 100 mg TiO2, a 25 mL gas volume, and 1 atm total pressure. Temperatures span 500 to 1500 °C.

## Ti3+ fraction

| Feed | 500 °C | 700 °C | 900 °C | 1100 °C | 1300 °C | 1500 °C |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CO2/H2 = 1:1 | 0, exact | 0, exact | 0, exact | 0, exact | 0, exact | 0, exact |
| Pure H2 | 0.0005% | 0.0168% | 0.1702% | 0.8694% | 2.82% | 6.37% |
| Sealed N2 | 1.5e-36% | 1.7e-26% | 7.2e-20% | 3.9e-15% | 1.4e-11% | 8.7e-9% |

The active-set solver returns exact zero for every inactive reduced phase under the RWGS feed. The sealed-N2 trace is a real two-phase-buffer solution rather than an optimizer floor.

## Stable reduced phases

| Feed | Reduced phases present |
| --- | --- |
| CO2/H2 = 1:1 | None; rutile throughout |
| Pure H2 | Ti10O19 |
| Sealed N2 | Ti10O19 at trace amount |

Pure H2 remains on the TiO2/Ti10O19 buffer across the temperature range. It converts more of the solid at higher temperature but does not pass the first reduction step under the stated finite gas charge.

## RWGS stability margin

Positive values mean the reduced phase remains unstable.

| Temperature | 500 °C | 700 °C | 900 °C | 1100 °C | 1300 °C | 1500 °C |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Nearest reduced-phase KKT cost, kJ mol-O⁻¹ | 70.40 | 62.79 | 54.74 | 46.17 | 37.03 | 27.30 |

The nearest implemented phase is Ti10O19. Higher-`n` Magneli members narrow the formal line-compound margin, but their limit is the unimplemented rutile TiO2-x solution. The result supported by the current registry is that the RWGS feed remains on the rutile side throughout 500 to 1500 °C.

## Sources

Production values are stored in `data/reference_results_high_precision.json`, and the table above is checked against it by the test suite. Legacy regression values are stored in `data/results.json`; that run carries the incomplete CH4 case and reports an oxygen-potential ladder margin rather than a KKT reduced cost, so its numbers are not interchangeable with the table above and are not quoted.
