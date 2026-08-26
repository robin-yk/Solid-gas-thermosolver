# Changelog

## Unreleased

- Added a manuscript figure-plate pipeline: a registry of what each figure claims, a frozen extract of every value the plates carry, a print drawing kit at journal column widths, and a self-contained proof sheet with per-plate SVG and 600 dpi export.
- Corrected the 500 and 700 C entries of the RWGS stability-margin table, which quoted the legacy oxygen-potential ladder margin from the CH4-containing run instead of the certified KKT reduced cost, and added a gate that holds the table and the README to the reference.
- Exposed the declared subsurface shell thickness (1 to 4 oxygen layers) as a workspace control, certified the 600 C / 95 umol-O/g ladder against the oracle, and documented the shell/bulk sensitivity it carries.
- Redrew the magnified slab so one dot is one vacancy on one oxygen site: the oxygen sublattice is shown at four O per (1x1) cell per 0.325 nm layer, on an isotropic scale, with the shell bracketed at its real thickness.
- Replaced the reconstruction onset warning with a phase-aware partition: the (1x2) added-row Ti2O3 reconstruction enters as a first-order surface phase (lever rule at pinned mu_V, areal deficiency 0.5 bridging-ML), and beyond the point-defect solubility mu_V pins at the estimated rutile/CS coexistence with the excess reported as extended defects. Closed-form boundary ladder (2.31 / 6.78 / ~160 umol-O/g at 600 C) certified by the 50-digit oracle.
- Decoupled the layer partition (analytic, formation energies only) from the surface-ordering Monte Carlo (effective pair interactions constrained to the Birschitzky statistics; arrangement exhibits only), restoring exact oracle certification of the shipped default.
- Added a coarse-grained kinetic Monte Carlo reoxidation-transport layer after Katsoulakis & Vlachos (2003): detailed-balance exchange on radial coarse cells, expm-oracle anchored, with an effective-barrier fit from measured recovery.
- Marked every CO2-recoverable quantity illustrative until the refill kinetics are fitted.
- Added a two-workspace browser portal and full-viewport light interface.
- Redesigned the equilibrium and defect workspaces with Helvetica-family sans-serif typography and larger scientific figures.
- Replaced the defect particle graphic with spherical site-class shells tied to calculated occupancies.
- Corrected the rutile (110) reconstruction threshold from a hard 20% cap to an onset warning.
- Added thermodynamic validity bridging between the gas-solid and defect workspaces.
- Reorganized browser sources, generated pages, data, scripts, and technical documentation into dedicated directories.
- Split the former monolithic README into focused method, verification, data, results, legacy, and roadmap documents.
