/* Redox cycling with memory, in the browser.

   Line-for-line port of solidgas/cycling.py. The state is two numbers, not
   one:

       A   accessible   inventory a CO2 exposure can still reach
       L   locked       inventory that survived an exposure and stopped
                        responding to them

   A reduction adds its oxygen deficit to A. An exposure recovers a
   fraction f of A and locks a fraction `lock` of what is left. At
   lock = 0 nothing is ever locked and the model is the memoryless one
   exactly, so the memory is a switch and not an assumption.

   Two cycles constrain two parameters and no more. `lock` is a survival
   probability, not a mechanism: migration out of reach, aggregation into
   an ordered defect and a local relaxation all produce it, and cycle
   integrals alone cannot tell them apart. */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.Cycling = api;
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  function initialState(accessible, locked) {
    return { accessible: accessible || 0.0, locked: locked || 0.0 };
  }
  function total(st) { return st.accessible + st.locked; }

  function reduceStep(st, created) {
    if (created < 0.0) {
      throw new Error('a reduction step cannot create negative inventory');
    }
    return { accessible: st.accessible + created, locked: st.locked };
  }

  function reoxidiseStep(st, fRec, lock) {
    lock = lock == null ? 1.0 : lock;
    if (!(fRec >= 0.0 && fRec <= 1.0)) {
      throw new Error('f_rec is a fraction of the accessible pool');
    }
    if (!(lock >= 0.0 && lock <= 1.0)) {
      throw new Error('lock is a probability');
    }
    var a = st.accessible;
    var recovered = fRec * a;
    var left = a - recovered;
    return { state: { accessible: left * (1.0 - lock),
                      locked: st.locked + left * lock },
             recovered: recovered };
  }

  function run(schedule, fRec, lock, state) {
    lock = lock == null ? 1.0 : lock;
    var st = state || initialState();
    var rows = [], i;
    for (i = 0; i < schedule.length; i++) {
      var beforeTotal = total(st);
      st = reduceStep(st, schedule[i]);
      var pre = { accessible: st.accessible, locked: st.locked };
      var f = typeof fRec === 'function' ? fRec(i, st) : fRec;
      var step = reoxidiseStep(st, f, lock);
      st = step.state;
      rows.push({ cycle: i + 1, created_umol_g: schedule[i],
                  carried_in_umol_g: beforeTotal,
                  pre_total_umol_g: total(pre),
                  pre_accessible_umol_g: pre.accessible,
                  pre_locked_umol_g: pre.locked,
                  f_rec_used: f,
                  recovered_umol_g: step.recovered,
                  residual_umol_g: total(st),
                  accessible_umol_g: st.accessible,
                  locked_umol_g: st.locked,
                  f_vs_created: schedule[i]
                    ? step.recovered / schedule[i] : 0.0,
                  f_vs_pre_total: total(pre)
                    ? step.recovered / total(pre) : 0.0 });
    }
    return rows;
  }

  function fit(observations, lockMax) {
    lockMax = lockMax == null ? 1.0 : lockMax;
    var obs = observations.slice();
    if (obs.length < 2) {
      throw new Error('two cycles are the minimum that constrains lock');
    }
    var first = obs[0];
    if (!(first.created > 0)) {
      throw new Error('the first cycle must create some inventory');
    }
    var f = first.recovered / first.created;
    if (!(f > 0.0 && f <= 1.0)) {
      throw new Error('the first cycle recovers more than it created');
    }

    /* Only the part of the ACCESSIBLE pool that an exposure failed to
       recover is available to lock; inventory locked in an earlier cycle
       is already gone from A and must not be counted again. So the
       denominator is A_n - R_n and not the cumulative residual - those
       coincide for two cycles and diverge from the third onwards. */
    var locks = [], raw = [], notes = [];
    var left = first.created - first.recovered;
    var n, row, need, carried, lk, over;
    for (n = 1; n < obs.length; n++) {
      row = obs[n];
      need = row.recovered / f;
      carried = need - row.created;
      lk = left <= 0.0 ? lockMax : 1.0 - carried / left;
      raw.push(lk);
      if (lk > lockMax) {
        over = f * row.created - row.recovered;
        notes.push('cycle ' + (n + 1) + ': complete locking is the least '
          + 'recovery this model can predict and it still over-predicts by '
          + over.toFixed(2) + ' umol-O/g ('
          + (100.0 * over / row.recovered).toFixed(1) + '%), so f is '
          + 'falling as well as the residual surviving');
        lk = lockMax;
      }
      if (lk < 0.0) {
        notes.push('cycle ' + (n + 1) + ': implies negative locking, i.e. '
          + 'the residual is more accessible than fresh inventory');
      }
      locks.push(lk);
      left = need - row.recovered;
    }

    var lock = locks.reduce(function (a, v) { return a + v; }, 0)
      / locks.length;
    lock = Math.min(lockMax, Math.max(0.0, lock));
    var pred = run(obs.map(function (o) { return o.created; }), f, lock);
    var resid = obs.map(function (o, i) {
      return pred[i].recovered_umol_g - o.recovered;
    });
    return { f_rec: f, lock: lock, lock_per_cycle: locks,
             lock_unclamped: raw,
             lock_at_bound: raw.some(function (x) { return x >= lockMax; }),
             predicted: pred,
             residual_umol_g: resid,
             rms_umol_g: Math.sqrt(resid.reduce(function (a, r) {
               return a + r * r; }, 0) / resid.length),
             notes: notes,
             f_per_cycle: obs.map(function (o) {
               return o.created ? o.recovered / o.created : 0.0; }) };
  }

  /* The experiment that separates the two models: the same total V_O made
     in one long reduction and in several short cycles. A memoryless
     inventory says the two must recover the same, because the total is its
     whole state; a locked reservoir says the cycled sample recovers less,
     by exactly what has stopped responding. */
  function twoSamplePrediction(totalVo, fRec, lock, split) {
    lock = lock == null ? 1.0 : lock;
    var recFresh = reoxidiseStep(initialState(totalVo), fRec, lock).recovered;
    var n = split == null ? 2 : Math.round(split);
    if (n < 1) throw new Error('the cycled sample needs at least one cycle');
    var st = initialState();
    var per = totalVo / (n + 1), i;
    for (i = 0; i < n; i++) {
      st = reduceStep(st, per);
      st = reoxidiseStep(st, fRec, lock).state;
    }
    st = reduceStep(st, totalVo - total(st));
    var cycled = { accessible: st.accessible, locked: st.locked };
    var recCycled = reoxidiseStep(cycled, fRec, lock).recovered;
    return { total_umol_g: totalVo, f_rec: fRec, lock: lock,
             cycles_to_build: n,
             fresh: { accessible: totalVo, locked: 0.0,
                      recovered_umol_g: recFresh },
             cycled: { accessible: cycled.accessible,
                       locked: cycled.locked,
                       recovered_umol_g: recCycled },
             recovery_gap_umol_g: recFresh - recCycled,
             recovery_ratio: recFresh ? recCycled / recFresh : NaN };
  }

  return { initialState: initialState, total: total,
           reduceStep: reduceStep, reoxidiseStep: reoxidiseStep,
           run: run, fit: fit,
           twoSamplePrediction: twoSamplePrediction };
}));
