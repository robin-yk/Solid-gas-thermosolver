# Thermodynamic data

## Sources

Gas standard states use NIST-JANAF 4th edition Shomate coefficients. Titanium-bearing condensed phases use Waldner and Eriksson, Calphad 23 (1999) 189-218, including the Magneli series.

The two condensed-phase assessments available in the package are kept explicit. They differ by several kilojoules per mole for rutile and reduced oxides, which is comparable to some reported stability margins. A correction derived on one basis must not be applied to a margin calculated on the other.

## Implemented titanium oxide model

The production active-set solver treats condensed phases as pure line compounds. Its registry contains rutile and reduced titanium oxide phases through Ti2O3. The reported RWGS margin is measured against the easiest available reduction step from rutile.

The shallow Magneli members move this boundary relative to a Ti4O7-only comparison. At high `n`, however, the line-compound series approaches oxygen-deficient rutile. A rutile TiO2-x solution is needed to represent that limit.

## Known data issue

The printed Ti7O13 formation enthalpy in the source assessment duplicates the Ti4O7 value. The package retains a series-based repair for internal ordering checks but does not use Ti7O13 to establish the quoted RWGS margins.

## Carbon scope

The legacy gas registry includes CH4 but does not include graphite. Carbon therefore has no condensed sink in that calculation. CH4 results remain in regression data and must not be quoted as complete carbon equilibrium.

## Files

- `solidgas/shomate.py`: NIST-JANAF gas and legacy condensed-phase coefficients
- `solidgas/waldner.py`: Ti-O assessment and integrated solid chemical potentials
- `data/thermo_data.json`: browser export for the legacy dashboard
- `data/activeset_data.json`: browser export for the production solver
- `data/oxide_data.json`: browser export for the legacy oxide tool
- `data/rutile_dft.json`: defect-model inputs and provenance

The JSON exports are generated from the Python package. Tests reject stale exports.
