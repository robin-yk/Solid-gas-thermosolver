# Vacancy distribution at fixed inventory

The distribution workspace implements Supplementary Note 2b supplied on 2026-09-19. It accepts a measured oxygen removal in micromoles of oxygen per gram of initial TiO2, particle radius, and a common distribution temperature. Sample presets change the inventory only. The gas–solid equilibrium and kinetic population models remain separate calculations.

For oxygen-site vacancy fraction x, local charge compensation gives a Ti3+ fraction y = 4x. The mean oxygen-site vacancy fraction is inventory × 79.866 × 10^-6 / 2. This x is half the stoichiometric oxygen deficit δ in TiO2−δ.

The neutral vacancy energy is −A exp[−(R−r)/ξ]. Ideal mixing of oxygen vacancies and compensating Ti3+ gives the stationarity equation

    G + kBT [ln(x/(1−x)) + 2 ln(4x/(1−4x))] = λ

The multiplier λ enforces the volume-integrated inventory. The local free energy is strictly convex for 0 < x < 0.25, so the local occupation and inventory-constrained solution are unique. The browser uses nested bisection; Python uses a bracketed scalar root with local bisection. Spherical shell weights integrate volume exactly; occupations are sampled at shell midpoints. A refined surface mesh resolves the exponential energy, with a cell boundary at the 2 nm integration depth.

The default radius is 450 nm and temperature is 873.15 K. The assumed decay length is 0.5 nm. A = 1.31 eV is the supplied manuscript's difference between bulk and surface bridging-vacancy energies, 5.70 and 4.39 eV. These zero-temperature energies are used at finite temperature. Dense surface occupations extrapolate dilute-defect energies; vacancy interactions and structural reconstruction are omitted.

## Validation

`tests/test_vacancy_distribution.py` compares Python and browser outputs for all five measured inventories, changed radius and temperature, both sensitivity decay lengths, zero stabilization, and a particle wholly within the 2 nm shell. Invalid inputs are rejected. An independent 80-digit local calculation and Gauss–Legendre integration check the default solution against the continuum inventory constraint. The quadrature acceptance is 10^-6 relative inventory error; the discrete inventory closure acceptance is 10^-8.

The R600 default returns 22.4856% vacant oxygen sites at the surface, 0.338014% in the interior, and 11.1152% of all vacancies in the outermost 2 nm. These are model outputs calculated at the supplied measured inventory of 94 micromoles O per gram, not measurements of surface coverage.

The kinetic tab retains the original `#population` URL. The new tab uses `#distribution`; previous `#ws-thermo` and `#ws-population` links are also accepted.

### Repository check on 2026-09-19

The full suite initially returned 236 passed and 4 failed. One failure was the browser-file inventory assertion, which required the two new distribution files. The other three arose from regenerated paper tables differing from the committed tables and the corresponding embedded page data. The regeneration changed 278 numeric JSON entries and no keys, list lengths, or text. The largest relative change was 8.20 × 10^-7 in the existing population-model fitted prefactor; its objective changed by 2.74 × 10^-12 relative. The equilibrium conversion differences were approximately 10^-12 percentage points. This check does not establish byte-identical reproduction across environments.

The original committed paper outputs were restored for the preview, and the page was rebuilt from those outputs. The targeted distribution and paper-output tests then passed: 29 passed, with the byte-for-byte regeneration test excluded. The browser-file assertion now includes the distribution engine and interface.

For release, the tables and embedded page data were regenerated together using Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1 and mpmath 1.4.1. This records the small numerical differences above; no gas–solid or kinetic model equations were changed. The original byte-for-byte reproducibility gate remains unchanged.

The release run passed all 240 tests with no skips or expected failures. One existing intermediate overflow warning was emitted by the gas–solid solver in the full-series test; that test passed. Browser checks confirmed all three workspace schematics loaded, the R600 distribution result, the source-code-only footer, and no horizontal overflow at 390 px viewport width.
