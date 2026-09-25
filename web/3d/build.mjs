/* Bundle the 3D figures into web/vacancy3d.js. The first line of the
   bundle records a hash of every source it was built from, so the test
   suite can tell a stale bundle without running node. */
import { build } from 'esbuild';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';

export const SOURCES = ['common.js', 'main.js', 'particle.js', 'slab.js', 'assets/studio.jpg', 'package-lock.json'];
const h = createHash('sha256');
for (const f of SOURCES) h.update(f + '\0').update(readFileSync(new URL(f, import.meta.url)));
const tag = h.digest('hex');

await build({
  entryPoints: ['main.js'], bundle: true, minify: true, format: 'iife', globalName: 'Vacancy3D',
  target: 'es2019', legalComments: 'eof', loader: { '.jpg': 'dataurl' },
  outfile: '../vacancy3d.js',
  banner: { js: '/* vacancy3d sources ' + tag + ' - built from web/3d by `npm run build`; three.js (MIT) inside */' }
});
console.log('web/vacancy3d.js', tag.slice(0, 12));
