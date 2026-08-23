/* Runs the browser engine under node and prints its answers as JSON, so the
   Python test can compare them against the reference without a browser. */
const fs = require('fs'); global.fs = fs;
const ROOT = process.argv[2]; global.ROOT = ROOT;
global.document = {
  getElementById: () => ({ textContent: fs.readFileSync(ROOT + '/data/oxide_data.json', 'utf8') }),
};
const harness = `
const ref = JSON.parse(fs.readFileSync(ROOT + '/data/oxide_reference.json','utf8'));
const GAS = ['CO2','H2','CO','H2O','CH4','O2'];
const out = [];
for (const name of Object.keys(Oxide.D.systems)) {
  const sy = Oxide.D.systems[name];
  const S = new Oxide.GibbsSystem(GAS, sy.solids, ['CO2','H2','H2O', sy.host]);
  const nMetal = 0.100 / sy.mw;
  for (const row of ref.rows.filter(r => r.system === name)) {
    const T = row.T_C + 273.15, ng = 25.0 / (Oxide.D.R_ATM * T), n0 = ng / 2;
    const b = S.bFrom({CO2:n0, H2:n0, [sy.host]: nMetal});
    const r = Oxide.minimise(S, b, T, Oxide.buildSeeds(S, ng, nMetal), 1.0);
    const y = S.gasFractions(r.n);
    const gas = {}; GAS.forEach(g => { gas[g] = y[g] * 100; });
    const bal = {};
    S.elements.forEach((e, j) => {
      let v = 0;
      for (let i = 0; i < S.species.length; i++) v += S.E[i][j] * r.n[i];
      bal[e] = v - b[j];
    });
    out.push({
      system: name, T_C: row.T_C, method: r.method, formulation: r.formulation,
      reduced_pct: S.reducedPercent(r.n, nMetal),
      phase_split_pct: S.phaseSplit(r.n, nMetal),
      gas_pct: gas,
      conversion_CO2_pct: (n0 - r.n[S.species.indexOf('CO2')]) / n0 * 100,
      Q_rwgs: (y.CO * y.H2O) / (y.CO2 * y.H2),
      G_rel_kJ: r.G_rel_kJ,
      balance_residual_mol: r.residual,
      balance_by_element: bal,
      elements: S.elements,
    });
  }
}
console.log(JSON.stringify(out));
`;
(0, eval)(fs.readFileSync(ROOT + '/web/oxide.js', 'utf8') + harness);
