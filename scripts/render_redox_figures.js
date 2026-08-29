/* Draw the redox figures from the exported rows. No physics here. */
const fs = require('fs');
const path = require('path');
const root = process.argv[2];
const data = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const out = process.argv[4];
const K = require(path.join(root, 'web', 'figkit.js'));
global.self = { FigKit: K };
const F = require(path.join(root, 'web', 'figures_redox.js'));
const made = [];
Object.keys(F.FIGURES).forEach(id => {
  const svg = F.FIGURES[id](data);
  fs.writeFileSync(path.join(out, 'redox_' + id + '.svg'), svg);
  made.push('redox_' + id + '.svg');
});
console.log(made.join(', ') + ' written');
