/* Worker shell for the defect stat-mech engine. The page assembles this
   worker from a Blob whose source is particle.js, then statmech.js, then
   this file, so the StatMech global is in scope and so is the Particle
   global its free energy comes from; the dataset arrives in the message,
   so the worker makes no requests of its own. */

'use strict';

self.onmessage = function (ev) {
  var msg = ev.data;
  if (!msg || msg.cmd !== 'run') return;
  try {
    var pts = StatMech.isothermScan(msg.data, msg.T_C, {
      seed: msg.seed,
      quality: msg.quality,
      eps: msg.eps,
      ordering: msg.ordering,
      mc: msg.mc,
      progress: function (i, n) {
        self.postMessage({ kind: 'progress', done: i, total: n,
                           token: msg.token });
      }
    });
    self.postMessage({ kind: 'done', isotherm: pts, token: msg.token });
  } catch (e) {
    self.postMessage({ kind: 'error',
                       message: String(e && e.message ? e.message : e),
                       token: msg.token });
  }
};
