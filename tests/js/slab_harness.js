/* Run web/slab_canvas.js's site fractions and vacancy placement over the
   stored cases the parity gate asks for. */
'use strict';
const fs = require('fs');
const ROOT = process.argv[2];
const points = JSON.parse(process.argv[3]);
const S = require(ROOT + '/web/slab_canvas.js');
const D = JSON.parse(fs.readFileSync(ROOT + '/paper_outputs/tof_range.json', 'utf8'));
const at = JSON.parse(fs.readFileSync(ROOT + '/web/render/slab_atoms.json', 'utf8'));
const ne = D.axes.eps.length, nf = D.axes.f110.length;
const out = {};
for (const smp of D.samples) for (const p of points) {
  const f110 = D.axes.f110[Math.floor(p / ne) % nf];
  const fr = S.fractions(D, smp, p, f110);
  out[smp.sample + ':' + p] = { fractions: fr, vacant: S.vacancies(at, fr) };
}
console.log(JSON.stringify(out));
