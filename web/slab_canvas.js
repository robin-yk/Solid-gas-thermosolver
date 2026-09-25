/* The (110) surface panel of the vacancy-distribution workspace, drawn live
   on a canvas. Atoms come from web/render/slab_atoms.json (rutile CIF,
   unrelaxed, four trilayers of 8 x 12 cells). A site of class c is vacant
   when its fixed rank is below the calculated fraction of c for the
   parameter point shown; scripts/slab_atoms.py holds the same rule and a
   parity test runs both.

   Motion is for display: a slow sway of the view and atomic vibration of
   arbitrary amplitude. Vacancy positions do not move; they change only
   when the calculated fractions change. */
(function (root) {
  'use strict';

  /* Site fractions of one stored case: theta for BRI, pool / capacity for
     the other classes (capacities at the case's (110) share). */
  function fractions(D, smp, p, f110) {
    var cap = D.capacity_umol_g[String(f110)], pools = {};
    D.pools.forEach(function (q, k) { pools[q] = smp.cases.pools[p][k]; });
    return { BRI: smp.cases.theta[p],
             IPL: pools.basal_L1 / cap.basal_L1,
             SBR: pools.L1_subbridging / cap.L1_subbridging,
             L24: pools.subsurface_L2_4 / cap.subsurface_L2_4 };
  }
  function vacancies(doc, frac) {
    var out = [];
    for (var i = 0; i < doc.site.length; i++) {
      var f = frac[doc.site[i]];
      if (f !== undefined && doc.rank[i] < f) out.push(i);
    }
    return out;
  }
  function counts(doc, vac) {
    var n = { BRI: 0, IPL: 0, SBR: 0, L24: 0 };
    vac.forEach(function (i) { n[doc.site[i]] += 1; });
    return n;
  }

  var COL = { O: '#ece6dc', Ti: '#5d6b78', Ti3: '#CC79A7' };
  var RING = { surface: '#0072B2', subsurface: '#D55E00' };
  var RAD = { O: 1.15, Ti: 0.72 };
  var TILT = 1.0, SWAY = 0.32, PERIOD = 16, DIST = 150, VIB = 0.09;

  function hexRGB(h) { return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)]; }
  function sprite(hex) {
    var S = 96, c = document.createElement('canvas');
    c.width = c.height = S;
    var g = c.getContext('2d'), rgb = hexRGB(hex), r = S / 2;
    var lit = rgb.map(function (v) { return Math.min(255, v + 70); });
    var dark = rgb.map(function (v) { return (v * 0.55) | 0; });
    var gr = g.createRadialGradient(r * 0.62, r * 0.58, r * 0.06, r, r, r);
    gr.addColorStop(0, 'rgb(' + lit + ')');
    gr.addColorStop(0.4, 'rgb(' + rgb + ')');
    gr.addColorStop(1, 'rgb(' + dark + ')');
    g.fillStyle = gr;
    g.beginPath(); g.arc(r, r, r - 0.5, 0, 2 * Math.PI); g.fill();
    return c;
  }

  function disc(hex) {
    var S = 96, c = document.createElement('canvas');
    c.width = c.height = S;
    var g = c.getContext('2d');
    g.fillStyle = hex;
    g.beginPath(); g.arc(S / 2, S / 2, S / 2 - 0.5, 0, 2 * Math.PI); g.fill();
    return c;
  }

  /* layer groups the particle panel can point at */
  function inGroup(doc, i, name) {
    var L = doc.layer[i], s = doc.site[i];
    if (name === 'surface') return L === 1 && (s === 'BRI' || s === 'IPL' || s === 'Ti');
    if (name === 'subsurface') return (L === 1 && s === 'SBR') || L > 1;
    return false;
  }

  function mount(canvas, doc, opts) {
    opts = opts || {};
    var n = doc.el.length, g = canvas.getContext('2d');
    var reduced = !!(root.matchMedia && root.matchMedia('(prefers-reduced-motion: reduce)').matches);
    var SPR = { O: sprite(COL.O), Ti: sprite(COL.Ti), Ti3: sprite(COL.Ti3), fog: disc('#ffffff') };
    var zc = 0, i;
    for (i = 0; i < n; i++) zc += doc.xyz[i][2];
    zc /= n;
    var ztop = -Infinity;
    for (i = 0; i < n; i++) ztop = Math.max(ztop, doc.xyz[i][2]);
    var ph = [], seed = 7;
    function rnd() { seed = (seed * 16807) % 2147483647; return seed / 2147483647; }
    for (i = 0; i < n; i++) ph.push([6.283 * rnd(), 6.283 * rnd(), 6.283 * rnd()]);
    var vac = new Float32Array(n), vacT = new Uint8Array(n), ti3 = new Float32Array(n), ti3T = new Uint8Array(n);
    var focus = null, dim = new Float32Array(n).fill(1), vacList = [];
    var running = false, paused = reduced, raf = 0, last = 0, t = 0, visible = true;
    var order = new Array(n);
    var P = new Float32Array(4 * n);

    function set(frac) {
      vacList = vacancies(doc, frac);
      vacT.fill(0); ti3T.fill(0);
      vacList.forEach(function (k) {
        vacT[k] = 1;
        doc.ti3[k].forEach(function (j) { ti3T[j] = 1; });
      });
      if (reduced || !running) { for (i = 0; i < n; i++) { vac[i] = vacT[i]; ti3[i] = ti3T[i]; } }
      kick();
      return counts(doc, vacList);
    }
    function highlight(name) { focus = name || null; kick(); }

    function size() {
      var dpr = Math.min(2, root.devicePixelRatio || 1), w = canvas.clientWidth, h = canvas.clientHeight;
      if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
        canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
      }
      return dpr;
    }

    function draw() {
      var dpr = size(), W = canvas.width, H = canvas.height;
      var yaw = reduced ? 0 : SWAY * Math.sin(2 * Math.PI * t / PERIOD);
      var cy = Math.cos(yaw), sy = Math.sin(yaw), ct = Math.cos(TILT), st = Math.sin(TILT);
      var sc = W / 58, ox = W / 2, oy = H * 0.56, vib = reduced ? 0 : VIB;
      g.clearRect(0, 0, W, H);
      for (i = 0; i < n; i++) {
        var a = doc.xyz[i], x = a[0], y = a[1], z = a[2] - zc;
        if (vib) {
          x += vib * Math.sin(t * 9.1 + ph[i][0]);
          y += vib * Math.sin(t * 10.3 + ph[i][1]);
          z += vib * Math.sin(t * 8.2 + ph[i][2]);
        }
        var x1 = x * cy - y * sy, y1 = x * sy + y * cy;
        var d = -y1 * st + z * ct, u = y1 * ct + z * st, f = DIST / (DIST - d);
        P[4 * i] = ox + x1 * sc * f; P[4 * i + 1] = oy - u * sc * f; P[4 * i + 2] = d; P[4 * i + 3] = f;
        order[i] = i;
      }
      order.sort(function (p, q) { return P[4 * p + 2] - P[4 * q + 2]; });
      for (var k = 0; k < n; k++) {
        i = order[k];
        var X = P[4 * i], Y = P[4 * i + 1], F = P[4 * i + 3], isO = doc.el[i] === 'O';
        var r = (isO ? RAD.O : RAD.Ti) * sc * F, base = dim[i];
        var alpha = base * (isO ? 1 - vac[i] : 1);
        if (alpha < 0.01) continue;
        g.globalAlpha = alpha;
        var spr = isO ? SPR.O : SPR.Ti;
        g.drawImage(spr, X - r, Y - r, 2 * r, 2 * r);
        if (!isO && ti3[i] > 0.01) {
          g.globalAlpha = alpha * ti3[i];
          g.drawImage(SPR.Ti3, X - r, Y - r, 2 * r, 2 * r);
        }
        /* depth shading: deeper atoms fade toward the background */
        var fog = Math.min(0.6, -(doc.xyz[i][2] - ztop) / 11 * 0.6);
        if (fog > 0.02) {
          g.globalAlpha = alpha * fog;
          g.drawImage(SPR.fog, X - r, Y - r, 2 * r, 2 * r);
        }
      }
      /* vacancy markers on top, so a site under its neighbours stays visible */
      for (var m = 0; m < n; m++) {
        i = order[m];
        if (doc.el[i] !== 'O' || vac[i] < 0.01) continue;
        var surf = doc.layer[i] === 1 && (doc.site[i] === 'BRI' || doc.site[i] === 'IPL');
        var R = (surf ? 1.55 : 1.2) * RAD.O * sc * P[4 * i + 3], cx = P[4 * i], cyy = P[4 * i + 1];
        var pulse = reduced ? 1 : 0.75 + 0.25 * Math.sin(t * 2.4 + ph[i][0]);
        var col = surf ? RING.surface : RING.subsurface, w = vac[i] * dim[i] * (surf ? 1 : 0.75);
        var gl = g.createRadialGradient(cx, cyy, 0, cx, cyy, R * 1.25);
        gl.addColorStop(0, 'rgba(' + hexRGB(col) + ',0.28)');
        gl.addColorStop(1, 'rgba(' + hexRGB(col) + ',0)');
        g.globalAlpha = w * pulse;
        g.fillStyle = gl;
        g.beginPath(); g.ellipse(cx, cyy, R * 1.25, R * 1.25 * ct, 0, 0, 2 * Math.PI); g.fill();
        g.setLineDash([5 * dpr, 3.5 * dpr]);
        g.lineDashOffset = reduced ? 0 : -t * 6 * dpr;
        g.lineWidth = (surf ? 2.6 : 2) * dpr;
        g.strokeStyle = col;
        g.beginPath(); g.ellipse(cx, cyy, R, R * ct, 0, 0, 2 * Math.PI); g.stroke();
        g.setLineDash([]);
      }
      g.globalAlpha = 1;
    }

    function step(dt) {
      var s = Math.min(1, dt / 0.6), moving = false;
      for (i = 0; i < n; i++) {
        var want = focus === null ? 1 : (inGroup(doc, i, focus) ? 1 : 0.22);
        if (Math.abs(dim[i] - want) > 1e-3) { dim[i] += (want - dim[i]) * Math.min(1, dt / 0.25); moving = true; } else dim[i] = want;
        if (vac[i] !== vacT[i]) { vac[i] += Math.sign(vacT[i] - vac[i]) * Math.min(Math.abs(vacT[i] - vac[i]), s); moving = true; }
        if (ti3[i] !== ti3T[i]) { ti3[i] += Math.sign(ti3T[i] - ti3[i]) * Math.min(Math.abs(ti3T[i] - ti3[i]), s); moving = true; }
      }
      return moving;
    }
    function frame(now) {
      raf = 0;
      var dt = last ? Math.min(0.1, (now - last) / 1000) : 0;
      last = now;
      var moving = step(dt);
      if (!paused) t += dt;
      draw();
      if (visible && (!paused || moving)) raf = root.requestAnimationFrame(frame);
      else last = 0;
    }
    function kick() {
      if (!running) { draw(); return; }
      if (!raf && visible) raf = root.requestAnimationFrame(frame);
    }
    function start() {
      running = true;
      if (root.IntersectionObserver) {
        new root.IntersectionObserver(function (es) {
          visible = es[0].isIntersecting; if (visible) kick();
        }).observe(canvas);
      }
      root.addEventListener('resize', kick);
      kick();
    }
    function pause(v) { paused = !!v; kick(); return paused; }

    start();
    return { set: set, highlight: highlight, pause: pause,
             paused: function () { return paused; }, reduced: reduced,
             time: function () { return t; }, vacant: function () { return vacList.slice(); } };
  }

  var api = { fractions: fractions, vacancies: vacancies, counts: counts, mount: mount };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.SlabCanvas = api;
})(typeof window !== 'undefined' ? window : this);
