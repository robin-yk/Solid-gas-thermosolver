/* The two figures the particle engine exists to produce.

   depth   theta against depth below the bridging plane: the layer stack
           the engine actually sums, the bulk plateau it decays to, and
           the Ti2O3 ceiling the cation sublattice imposes.
   where   how much of the oxygen sits inside each depth, against how
           much of the volume that depth is - the gap between the two
           curves is the enrichment, and it is the reason a shell too
           thin to see is worth drawing at all.

   Neither computes a site fraction. web/particle.js does that and these
   two only reshape what it returned, so a figure here cannot report a
   profile the engine never produced.

   One thing the depth figure has to be careful about and says on its
   face: the bridging row and the interior do not share a denominator.
   The row is one oxygen per (1x1) cell and an interior layer is four,
   so theta on the row is a fraction of a quarter-density sublattice.
   Putting them on one axis without saying so would suggest the surface
   is four times more reduced than it is. */

(function (root, factory) {
  'use strict';
  var kit = (typeof module !== 'undefined' && module.exports)
    ? require('./figkit.js') : root.FigKit;
  if (!kit) throw new Error('figures_particle.js needs figkit.js');
  var api = factory(kit);
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.ParticleFigures = api;
})(typeof self !== 'undefined' ? self : this, function (K) {
  'use strict';

  var C = K.C, LW = K.LW, T = K.FIGTYPE, MARK = K.MARK;
  var lin = K.lin, lg = K.lg, frame = K.frame, axisX = K.axisX;
  var axisY = K.axisY, series = K.series, tint = K.tint, decades = K.decades;

  /* --------------------------------------------------------- payloads */

  function depthData(res, part, layers) {
    var seg = part.segregation;
    var pts = layers.map(function (q) {
      return { z: q.z_nm, theta: q.theta, k: q.k };
    });
    return {
      layers: pts,
      surface_theta: res.surface.theta,
      surface_phase: res.surface.phase,
      bulk: res.theta_bulk_plateau,
      cap: 1.0 / res.polaron_ratio,
      xi_nm: seg.xi_nm,
      d110_nm: seg.d110_nm,
      dE_seg_eV: seg.dE_seg_eV,
      loading: res.loading_umol_g,
      T_C: res.T_C,
      limited_by: res.limited_by
    };
  }

  function whereData(rows, res) {
    return { rows: rows.slice(), loading: res.loading_umol_g, T_C: res.T_C };
  }

  /* ---------------------------------------------------------- figures */

  function depth(d) {
    var f = K.square();
    f.pane = { x0: 46, y0: 36, x1: 244, y1: 206 };
    var p = f.pane;

    /* The bridging row gets its own column.  It is not a deeper layer of
       the same lattice gas - past the transition it is the added-row
       (1x2) phase, whose oxygen deficiency is counted per bridging site
       and whose stoichiometry is Ti2O3 outright.  Sharing one x axis
       with the interior would invite reading it as a continuation of the
       decay, and sharing the interior's Ti2O3 ceiling would make a phase
       that IS Ti2O3 look like a violation. */
    var STRIP = 24, GAP = 9;
    var sx = p.x0 + 0.5 * STRIP;
    var ix0 = p.x0 + STRIP + GAP;

    var zmax = Math.max(6.0, 10.0 * d.xi_nm);
    var lo = Math.max(1e-6, Math.pow(10, Math.floor(
      Math.log(d.bulk / 3.0) / Math.LN10)));
    var X = lin(0, zmax, ix0, p.x1);
    var Y = lg(lo, 1.0, p.y1, p.y0);

    f.rect(p.x0, p.y0, STRIP, p.y1 - p.y0, tint(C.surface, 0.10));

    /* the cation sublattice is full at the cap - for the interior gas */
    f.rect(ix0, p.y0, p.x1 - ix0, Y(d.cap) - p.y0, tint(C.extended, 0.14));
    f.line(ix0, Y(d.cap), p.x1, Y(d.cap), C.extended, LW.guide, '3 2');
    f.text(p.x1 - 5, Y(d.cap) - 4, 'Ti₂O₃ — cation sublattice full',
           { size: T.small, anchor: 'end', fill: C.extended });

    f.line(ix0, Y(d.bulk), p.x1, Y(d.bulk), C.bulk, LW.guide, '2 2');
    f.text(p.x1 - 3, Y(d.bulk) + 9, 'bulk plateau',
           { size: T.small, anchor: 'end', fill: C.bulk });

    f.line(X(d.xi_nm), Y(d.cap), X(d.xi_nm), p.y1, C.structure,
           LW.hair, '1 2');
    f.text(X(d.xi_nm) + 3, Y(d.cap) + 11, 'ξ',
           { size: T.body, fill: C.structure });

    frame(f, [p.x0, p.x1], Y);
    f.line(ix0 - 0.5 * GAP, p.y0, ix0 - 0.5 * GAP, p.y1, C.ink, LW.hair, '2 3');
    axisY(f, Y, p.x0, decades(lo, 1.0), 'oxygen missing, as a fraction',
          K.expLabel);
    K.niceTicks(0, zmax, 6).forEach(function (t) {
      if (t > zmax) return;
      f.line(X(t), p.y1, X(t), p.y1 - K.TICK.major, C.ink, LW.axis);
      f.line(X(t), p.y0, X(t), p.y0 + K.TICK.major, C.ink, LW.axis);
      f.text(X(t), p.y1 + T.tick + 3, String(t),
             { size: T.tick, anchor: 'middle' });
    });
    f.text(sx, p.y1 + T.tick + 3, 'row',
           { size: T.tick, anchor: 'middle', fill: C.surface });
    f.text(0.5 * (ix0 + p.x1), p.y1 + T.tick + T.body + 7.5,
           'depth below the bridging plane (nm)',
           { size: T.body, anchor: 'middle' });

    var pts = d.layers.filter(function (q) { return q.z <= zmax; });
    series(f, pts.map(function (q) { return [q.z, q.theta]; }), X, Y,
           C.subsurface, LW.curve);
    pts.forEach(function (q) {
      f.dot(X(q.z), Y(q.theta), MARK * 0.7, C.subsurface);
    });

    f.dot(sx, Y(d.surface_theta), MARK * 1.4, C.surface);
    f.text(sx, Y(d.surface_theta) - 7, d.surface_phase,
           { size: T.small, anchor: 'middle', fill: C.surface });

    f.text(p.x0, 13, d.loading.toPrecision(3) + ' μmol-O/g at ' + d.T_C
           + ' °C', { size: T.body });
    f.text(p.x0, 23, 'one dot is one d₁₁₀ oxygen layer',
           { size: T.small, fill: C.structure });
    return f.done();
  }

  function where(d) {
    var f = K.square();
    f.pane = { x0: 46, y0: 20, x1: 244, y1: 206 };
    var p = f.pane;

    var rows = d.rows.filter(function (q) { return q.depth_nm > 0; });
    var z0 = rows[0].depth_nm, z1 = rows[rows.length - 1].depth_nm;
    var X = lg(z0, z1, p.x0, p.x1);
    var Y = lin(0, 1, p.y1, p.y0);

    frame(f, X, Y);
    axisX(f, X, p.y1, decades(z0, z1), 'within this depth of the surface (nm)',
          K.expLabel);
    axisY(f, Y, p.x0, [0, 0.25, 0.5, 0.75, 1.0], 'cumulative fraction');

    /* the gap between the curves is the whole point, so fill it */
    var poly = rows.map(function (q) {
      return X(q.depth_nm) + ',' + Y(q.inventory_fraction);
    }).concat(rows.slice().reverse().map(function (q) {
      return X(q.depth_nm) + ',' + Y(q.volume_fraction);
    }));
    f.path('M' + poly.join('L') + 'Z', 'none', 0, null,
                 tint(C.subsurface, 0.16));

    series(f, rows.map(function (q) {
      return [q.depth_nm, q.volume_fraction]; }), X, Y,
           C.structure, LW.curve, '3 2');
    series(f, rows.map(function (q) {
      return [q.depth_nm, q.inventory_fraction]; }), X, Y,
           C.subsurface, LW.curve);

    var mid = rows[Math.floor(rows.length * 0.5)];
    f.text(X(mid.depth_nm) - 6, Y(mid.inventory_fraction) - 6, 'oxygen',
           { size: T.body, anchor: 'end', fill: C.subsurface });
    var vy = Math.min(p.y1 - 6, Y(mid.volume_fraction) + 12);
    f.text(X(mid.depth_nm) + 6, vy, 'volume',
           { size: T.body, fill: C.structure });

    var one = rows.filter(function (q) { return q.depth_nm >= 1.0; })[0];
    if (one) {
      f.text(p.x0, 13,
                   'at 1 nm: ' + (100 * one.inventory_fraction).toFixed(1)
                   + '% of the oxygen in '
                   + (100 * one.volume_fraction).toFixed(2) + '% of the volume',
                   { size: T.small });
    }
    return f.done();
  }

  return { depthData: depthData, whereData: whereData,
           depth: depth, where: where };
});
