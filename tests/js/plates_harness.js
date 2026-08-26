/* Render every implemented plate under node and print the SVG sources.
   Usage: node plates_harness.js <repo root>  ->  {id: svg, ...} */
const path = require('path');
const fs = require('fs');

const root = process.argv[2];
const PL = require(path.join(root, 'web', 'plates.js'));
const D = JSON.parse(
  fs.readFileSync(path.join(root, 'data', 'figure_data.json'), 'utf8'));

const out = { implemented: [], svg: {}, errors: {} };
Object.keys(PL.PLATES).sort().forEach(id => {
  try {
    out.svg[id] = PL.render(id, D);
    out.implemented.push(id);
  } catch (e) {
    out.errors[id] = e.message;
  }
});
process.stdout.write(JSON.stringify(out));
