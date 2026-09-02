/* Run web/vacancy.js over the grid the parity gate asks for. */
'use strict';
const fs = require('fs');
const ROOT = process.argv[2];
const spec = JSON.parse(process.argv[3]);
const V = require(ROOT + '/web/vacancy.js');
const D = JSON.parse(fs.readFileSync(ROOT + '/paper_outputs/site_data.json', 'utf8'));
const g = D.surface.geometry_card;

const out = {
  cell_area_cm2: V.cellArea(g),
  areal_cm2: V.arealDensity(g),
  density_g_cm3: V.density(g),
  site_density: V.siteDensity(g),
  areas: spec.diameters.map(d => ({ d_um: d, area: V.area(g, d, null) })),
  capacities: [],
  bounds: [],
};
for (const d of spec.diameters)
  out.capacities.push({ d_um: d, bet: null,
    full: V.capacity(g, d, null), recon: V.deficiency(g, d, null) });
for (const bet of spec.bets)
  out.capacities.push({ d_um: null, bet: bet,
    full: V.capacity(g, null, bet), recon: V.deficiency(g, null, bet) });
for (const inv of spec.inventories) {
  for (const d of spec.diameters)
    out.bounds.push({ d_um: d, bet: null, inv: inv, out: V.bound(g, inv, d, null) });
  for (const bet of spec.bets)
    out.bounds.push({ d_um: null, bet: bet, inv: inv, out: V.bound(g, inv, null, bet) });
}
const threw = f => { try { f(); return false; } catch (e) { return true; } };
out.rejects_zero_inventory = threw(() => V.bound(g, 0, 0.9, null));
out.rejects_missing_area = threw(() => V.bound(g, 95, null, null));
console.log(JSON.stringify(out));
