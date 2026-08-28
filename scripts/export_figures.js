/* Write every plate to docs/figures as print SVG and 600 dpi PNG.

   The plate coordinate system is points and the SVG declares its width in
   millimetres, so a browser lays it out at true print size: 89 mm is
   3.504 in, which is 336.4 CSS px at the nominal 96 dpi. A device scale
   factor of 600/96 therefore rasterises at exactly 600 dpi.

     node scripts/export_figures.js [--svg-only]
*/
const fs = require('fs');
const path = require('path');

const ROOT = path.resolve(__dirname, '..');
const OUT = path.join(ROOT, 'docs', 'figures');
const DPI = 600;
const SCALE = DPI / 96;

const PL = require(path.join(ROOT, 'web', 'plates.js'));
const D = JSON.parse(fs.readFileSync(
  path.join(ROOT, 'data', 'figure_data.json'), 'utf8'));

function main() {
  fs.mkdirSync(OUT, { recursive: true });
  const ids = Object.keys(PL.PLATES).sort();
  const made = [];
  ids.forEach(id => {
    const svg = PL.render(id, D);
    const f = path.join(OUT, id + '.svg');
    fs.writeFileSync(f, svg + '\n');
    made.push({ id: id, svg: f });
  });
  console.log(`${made.length} SVG written to docs/figures`);
  if (process.argv.includes('--svg-only')) return Promise.resolve();

  const { chromium } = require('playwright');
  return chromium.launch().then(async browser => {
    const page = await browser.newPage({ deviceScaleFactor: SCALE });
    for (const m of made) {
      const svg = fs.readFileSync(m.svg, 'utf8');
      await page.setContent(
        '<style>html,body{margin:0;padding:0;background:#fff}'
        + 'svg{display:block}</style>' + svg,
        { waitUntil: 'load' });
      const el = await page.$('svg');
      const png = path.join(OUT, m.id + '.png');
      await el.screenshot({ path: png, scale: 'device' });
      const { width, height } = await el.boundingBox();
      console.log(`  ${m.id}: ${Math.round(width * SCALE)} x `
        + `${Math.round(height * SCALE)} px at ${DPI} dpi`);
    }
    await browser.close();
  });
}

Promise.resolve(main()).catch(e => {
  console.error(e);
  process.exit(1);
});
