/* Figure 2: the rutile (110) surface the reactive-site count refers to.
   Atoms come from the bulk structure (a = 4.594 A, c = 2.959 A,
   u = 0.305) cut into four O-Ti2O2-O trilayers; nothing is relaxed.
   Every site carries a fixed random rank and is drawn vacant when its
   rank is below the site fraction the model gives, so a change of
   scenario removes or restores atoms without reshuffling the rest.
   The (1x2) added rows are schematic: drawn where the model puts its
   reconstructed cells, not at relaxed positions. */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { ROLE, rng, stage, orbit, sway, REDUCED_MOTION } from './common.js';

const A = 4.594, C = 2.959, U = 0.305, S2 = Math.SQRT2;
const PX = A * S2, D110 = A / S2;            /* 6.497 A cell, 3.248 A trilayer */
export const NX = 8, NY = 16, LAYERS = 4;
const R_O = 0.62, R_TI = 0.44, BOND = 2.1;

const BASIS = [['Ti', 0, 0, 0], ['Ti', 0.5, 0.5, 0.5], ['O', U, U, 0], ['O', 1 - U, 1 - U, 0],
  ['O', 0.5 + U, 0.5 - U, 0.5], ['O', 0.5 - U, 0.5 + U, 0.5]];

/* Atoms of the slab in surface coordinates: X along [1-10], Y along
   [001], Z along [110] with the top Ti2O2 plane at Z = 0. */
function atoms() {
  const out = [], seen = new Set();
  for (let i = -40; i <= 40; i++) for (let j = -40; j <= 40; j++) for (let k = 0; k < NY; k++) {
    for (const [el, x, y, z] of BASIS) {
      const X = ((x + i) - (y + j)) * A / S2, Z = ((x + i) + (y + j)) * A / S2, Y = (z + k) * C;
      const home = Math.round(Z / D110);
      if (home > 0 || home < 1 - LAYERS) continue;
      if (X < -1e-6 || X >= NX * PX - 1e-6) continue;
      const key = el + X.toFixed(2) + ',' + Y.toFixed(2) + ',' + Z.toFixed(2);
      if (seen.has(key)) continue;
      seen.add(key);
      const dz = Z - home * D110;
      let site = el;
      if (el === 'O') site = home < 0 ? 'L24' : dz > 0.5 ? 'BRI' : dz < -0.5 ? 'SBR' : 'IPL';
      out.push({ el, X, Y, Z, layer: 1 - home, site, col: Math.floor(X / PX + 1e-6) });
    }
  }
  return out;
}

function glowTexture() {
  const c = document.createElement('canvas');
  c.width = c.height = 128;
  const g = c.getContext('2d');
  const grd = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  grd.addColorStop(0, 'rgba(255,255,255,1)');
  grd.addColorStop(0.28, 'rgba(255,255,255,0.85)');
  grd.addColorStop(0.5, 'rgba(255,255,255,0.28)');
  grd.addColorStop(1, 'rgba(255,255,255,0)');
  g.fillStyle = grd; g.fillRect(0, 0, 128, 128);
  g.strokeStyle = 'rgba(255,255,255,1)'; g.lineWidth = 7;
  g.beginPath(); g.arc(64, 64, 34, 0, 2 * Math.PI); g.stroke();
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

export function createSlab(host) {
  const S = stage(host, 26);
  const scene = S.scene;
  const cx = NX * PX / 2, cy = NY * C / 2;
  const W = function (a, dy) { return new THREE.Vector3(a.X - cx, a.Z + (dy || 0), a.Y - cy); };

  scene.add(new THREE.HemisphereLight(0xffffff, 0xd6dde3, 0.75));
  const key = new THREE.DirectionalLight(0xffffff, 1.6);
  key.position.set(-30, 70, 25);
  key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048);
  Object.assign(key.shadow.camera, { left: -45, right: 45, top: 45, bottom: -45, near: 20, far: 160 });
  key.shadow.bias = -0.0004;
  key.shadow.normalBias = 0.05;
  scene.add(key);
  S.renderer.shadowMap.enabled = true;
  S.renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  S.camera.far = 1000;
  S.camera.near = 1;
  S.camera.updateProjectionMatrix();
  S.camera.position.set(-60, 50, 80);
  const target = new THREE.Vector3(0, -7, 0);
  S.camera.lookAt(target);
  const ctl = orbit(OrbitControls, S, target);
  ctl.minDistance = 45; ctl.maxDistance = 190;
  const root = new THREE.Group();
  scene.add(root);

  const all = atoms();
  const r = rng(2959);
  all.forEach(function (a) { a.rank = r(); });
  const O = all.filter(function (a) { return a.el === 'O'; });
  const Ti = all.filter(function (a) { return a.el === 'Ti'; });

  /* depth cue: each trilayer down is a little darker */
  const shade = function (base, layer) { return new THREE.Color(base).multiplyScalar(1 - 0.11 * (layer - 1)); };
  function instanced(list, radius, mat, base) {
    const m = new THREE.InstancedMesh(new THREE.SphereGeometry(radius, 28, 18), mat, list.length);
    const M = new THREE.Matrix4();
    list.forEach(function (a, i) {
      M.makeTranslation(W(a));
      m.setMatrixAt(i, M);
      m.setColorAt(i, shade(base, a.layer));
    });
    m.instanceMatrix.needsUpdate = true;
    m.castShadow = m.receiveShadow = true;
    return m;
  }
  const oMat = new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.32, clearcoat: 0.6, clearcoatRoughness: 0.2 });
  const tMat = new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.28, metalness: 0.25, clearcoat: 0.4 });
  const oMesh = instanced(O, R_O, oMat, 0xe6eaee);
  const tMesh = instanced(Ti, R_TI, tMat, 0x75828d);
  root.add(oMesh, tMesh);

  /* Ti-O bonds, as thin cylinders; a bond disappears with its O */
  const bonds = [];
  Ti.forEach(function (t) {
    O.forEach(function (o, j) {
      const d = Math.hypot(t.X - o.X, t.Y - o.Y, t.Z - o.Z);
      if (d < BOND) bonds.push([t, o, j]);
    });
  });
  const bMesh = new THREE.InstancedMesh(new THREE.CylinderGeometry(0.13, 0.13, 1, 10, 1, true),
    new THREE.MeshStandardMaterial({ color: 0xaab4bc, roughness: 0.5 }), bonds.length);
  bMesh.castShadow = bMesh.receiveShadow = true;
  const up = new THREE.Vector3(0, 1, 0);
  const bondMatrix = function (t, o) {
    const a = W(t), b = W(o), mid = a.clone().add(b).multiplyScalar(0.5), dir = b.clone().sub(a);
    const len = dir.length();
    const q = new THREE.Quaternion().setFromUnitVectors(up, dir.normalize());
    return new THREE.Matrix4().compose(mid, q, new THREE.Vector3(1, len, 1));
  };
  const bondM = bonds.map(function (b) { return bondMatrix(b[0], b[1]); });
  root.add(bMesh);

  /* (1x2) added rows: two Ti and three O over each reconstructed cell,
     a cell being two bridging rows wide and one c long */
  const addTi = [], addO = [];
  for (let p = 0; p < NX / 2; p++) for (let j = 0; j < NY; j++) {
    const X0 = (2 * p + 1) * PX;          /* over the fivefold Ti row between the two bridging rows */
    addTi.push({ X: X0 - 1.6, Y: j * C + C / 2, Z: 2.05, layer: 1, cell: p * NY + j });
    addTi.push({ X: X0 + 1.6, Y: j * C + C / 2, Z: 2.05, layer: 1, cell: p * NY + j });
    addO.push({ X: X0, Y: j * C, Z: 3.0, layer: 1, cell: p * NY + j });
    addO.push({ X: X0 - 2.7, Y: j * C, Z: 1.9, layer: 1, cell: p * NY + j });
    addO.push({ X: X0 + 2.7, Y: j * C, Z: 1.9, layer: 1, cell: p * NY + j });
  }
  const recTint = new THREE.Color(ROLE.surface).lerp(new THREE.Color(0xffffff), 0.38);
  const aoMesh = instanced(addO, R_O, oMat.clone(), recTint.getHex());
  const atMesh = instanced(addTi, R_TI, tMat.clone(), new THREE.Color(ROLE.surface).lerp(new THREE.Color(0x4a5560), 0.4).getHex());
  root.add(aoMesh, atMesh);

  /* cells take the reconstruction a whole [001] row pair at a time, in a
     fixed order, so the added rows form domains as on the real surface */
  const pairOrder = [2, 0, 3, 1];
  const recRank = new Float32Array(NX / 2 * NY);
  pairOrder.forEach(function (p, k) {
    for (let j = 0; j < NY; j++) recRank[p * NY + j] = (k * NY + j + 0.5) / (NX / 2 * NY);
  });
  const cellOf = function (o) { return Math.floor(o.col / 2) * NY + (Math.round(o.Y / C) % NY); };

  /* vacancy markers */
  const tex = glowTexture();
  const markers = new THREE.Group();
  markers.renderOrder = 10;
  root.add(markers);
  const pool = [];
  function marker(i) {
    if (!pool[i]) {
      const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false, depthWrite: false }));
      markers.add(m);
      pool[i] = m;
    }
    return pool[i];
  }

  const HIDE = new THREE.Matrix4().makeScale(0, 0, 0);
  const api = { stage: S, counts: {} };
  let glow = [];
  const neighbours = { ISO_z2: 1, ISO_z4: 2, ISO_z8: 4 };

  api.set = function (v) {
    const rec = new Uint8Array(NX / 2 * NY);
    for (let c = 0; c < rec.length; c++) rec[c] = recRank[c] < v.f_rec ? 1 : 0;
    const frac = { BRI: v.theta, IPL: v.x_basal, SBR: v.x_sbr, L24: v.x_l24 };
    const vac = new Uint8Array(O.length);
    const M = new THREE.Matrix4();
    const bri = {};                          /* vacant bridging sites by (col,row) */
    const counts = { BRI: 0, IPL: 0, SBR: 0, L24: 0, cells: 0, reactive: 0 };
    O.forEach(function (o, i) {
      const cut = o.site === 'BRI' && rec[cellOf(o)];
      const gone = !cut && o.layer === 1 ? o.rank < frac[o.site] : o.layer > 1 && o.rank < frac.L24;
      vac[i] = gone ? 1 : 0;
      if (gone) counts[o.site] += 1;
      if (gone && o.site === 'BRI') bri[o.col + ',' + Math.round(o.Y / C) % NY] = 1;
      oMesh.setMatrixAt(i, cut || gone ? HIDE : M.makeTranslation(W(o)));
    });
    oMesh.instanceMatrix.needsUpdate = true;
    bonds.forEach(function (b, k) { bMesh.setMatrixAt(k, vac[b[2]] || (b[1].site === 'BRI' && rec[cellOf(b[1])]) ? HIDE : bondM[k]); });
    bMesh.instanceMatrix.needsUpdate = true;
    [[aoMesh, addO], [atMesh, addTi]].forEach(function (z) {
      z[1].forEach(function (a, i) { z[0].setMatrixAt(i, rec[a.cell] ? M.makeTranslation(W(a)) : HIDE); });
      z[0].instanceMatrix.needsUpdate = true;
    });
    for (let c = 0; c < rec.length; c++) counts.cells += rec[c];

    /* markers: bridging and in-plane vacancies in the surface hue, the
       sub-bridging and deeper ones in the subsurface hue; the ones that
       count as reactive sites under the chosen definition pulse */
    const z = neighbours[v.reactive] || 0;
    let n = 0;
    glow = [];
    O.forEach(function (o, i) {
      if (!vac[i]) return;
      let reactive = false;
      if (o.site === 'BRI') {
        reactive = true;
        if (z) {
          const row = Math.round(o.Y / C) % NY;
          for (let d = 1; d <= z && reactive; d++) {
            if (bri[o.col + ',' + (row + d) % NY] || bri[o.col + ',' + (row - d + NY) % NY]) reactive = false;
          }
        }
      } else if (o.site === 'IPL' && o.layer === 1 && v.reactive === 'BRI+BASAL') reactive = true;
      const m = marker(n++);
      const surf = o.layer === 1 && (o.site === 'BRI' || o.site === 'IPL');
      m.material.color.set(surf ? ROLE.surface : ROLE.subsurface);
      m.position.copy(W(o));
      m.visible = true;
      m.userData = { reactive, base: surf ? 2.6 : 2.0 };
      m.scale.setScalar(m.userData.base);
      m.material.opacity = reactive ? 1 : 0.55;
      if (reactive) { glow.push(m); counts.reactive += 1; }
    });
    for (let k = n; k < pool.length; k++) pool[k].visible = false;
    api.counts = counts;
    S.dirty = true;
  };

  S.tick = function (t, still) {
    let moving = false;
    const a = still ? 0 : sway(S, t, 0.35);
    if (a !== null) { root.rotation.y = a; moving = true; }
    if (!REDUCED_MOTION) {
      const s = 1 + 0.22 * (0.5 + 0.5 * Math.sin(t * 3.1));
      glow.forEach(function (m) { m.scale.setScalar(m.userData.base * s); });
      if (glow.length) moving = true;
    }
    if (!still && ctl.update()) moving = true;
    return moving;
  };
  api.sites = { O: O.length, Ti: Ti.length, BRI: O.filter(function (o) { return o.site === 'BRI'; }).length };
  return api;
}
