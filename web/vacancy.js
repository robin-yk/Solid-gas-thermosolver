/* The surface-capacity bound, in the browser.

   A line-for-line mirror of solidgas/vacancy_inventory.py. It is a dozen
   lines of arithmetic on the rutile lattice constants, which is why it is
   worth mirroring at all: the alternative is a precomputed grid that
   cannot answer a BET area the user actually measured. A parity gate
   holds this against the Python module to the last place. */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.Vacancy = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  var N_AVO = 6.02214076e23;

  function cellArea(g) {                       /* cm^2, the (110) 1x1 cell */
    var a = g.lattice_constants_A.a * 1e-8;
    var c = g.lattice_constants_A.c * 1e-8;
    return a * Math.sqrt(2) * c;
  }
  function arealDensity(g) {
    return g.surface_cell.bridging_O_per_1x1_cell / cellArea(g);
  }
  function density(g) {                        /* g/cm^3, from the cell */
    var a = g.lattice_constants_A.a * 1e-8;
    var c = g.lattice_constants_A.c * 1e-8;
    return g.formula_units_per_cell * g.molar_mass_g_mol / (N_AVO * a * a * c);
  }
  /* The grouping below is the Python module's, parenthesis for
     parenthesis. Floating-point multiplication does not associate, and
     reordering these chains moves the answer by a couple of places - which
     is enough to turn a bit-identical gate into a tolerance. */
  function area(g, d_um, bet) {                /* m^2/g */
    if (bet != null) return bet;
    if (d_um == null) throw new Error('give a BET area or a particle diameter');
    var d_cm = d_um * 1e-4;
    return 6.0 / (density(g) * d_cm) * 1e-4;
  }
  function capacity(g, d_um, bet) {            /* umol-O/g, whole row */
    var area_cm2_g = area(g, d_um, bet) * 1e4;
    return arealDensity(g) * area_cm2_g / N_AVO * 1e6;
  }
  function deficiency(g, d_um, bet) {          /* umol-O/g, (1x2) */
    return g.reconstruction.areal_O_deficiency_ML * capacity(g, d_um, bet);
  }
  function siteDensity(g) {
    return g.oxygen_per_formula_unit / g.molar_mass_g_mol * 1e6;
  }
  function toX(g, umol) {
    return umol / siteDensity(g) * g.oxygen_per_formula_unit;
  }
  function toTi3(g, umol) {
    return Math.min(2.0 * umol / (1.0 / g.molar_mass_g_mol * 1e6), 1.0);
  }

  function bound(g, umol, d_um, bet) {
    if (!(umol > 0)) throw new Error('the measured inventory must be positive');
    var full = capacity(g, d_um, bet), recon = deficiency(g, d_um, bet);
    return {
      measured_umol_o_g: umol,
      area_m2_g: area(g, d_um, bet),
      density_g_cm3: density(g),
      bridging_capacity_umol_o_g: full,
      reconstruction_deficiency_umol_o_g: recon,
      inventory_over_capacity: umol / full,
      surface_fraction_max: Math.min(recon / umol, 1),
      below_surface_fraction_min: Math.max(1 - recon / umol, 0),
      surface_fraction_max_full_row: Math.min(full / umol, 1),
      below_surface_fraction_min_full_row: Math.max(1 - full / umol, 0),
      x_in_TiO2mx: toX(g, umol),
      Ti3_fraction: toTi3(g, umol)
    };
  }

  return { cellArea: cellArea, arealDensity: arealDensity, density: density,
           area: area, capacity: capacity, deficiency: deficiency,
           siteDensity: siteDensity, toX: toX, toTi3: toTi3, bound: bound };
});
