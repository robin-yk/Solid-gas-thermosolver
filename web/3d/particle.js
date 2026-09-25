/* Figure 1: one rutile particle with an octant cut away. Every bead is a
   fixed amount of vacancies per gram, so the bead count is the pool size:
   bulk beads fill the volume, surface and subsurface beads sit in a band
   under the (110) tiles only, because the explicit trilayers exist only
   there. The band is drawn thicker than scale; the caption says by how
   much. */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { ROLE, rng, stage, ease, contactShadow, beadMaterial, orbit, sway } from './common.js';

export const DOT_UMOL = 0.04;              /* umol O per g of oxide per bead */
const CAP = { bulk: 34000, subsurface: 6000, surface: 3000 };
const BAND = { surface: [0.978, 0.995], subsurface: [0.952, 0.975], bulk: [0, 0.948] };
const TILES = 96;

function inOctant(x, y, z) { return x > 0 && y > 0 && z > 0; }

function fibonacci(n) {
  const out = [], g = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < n; i++) {
    const y = 1 - 2 * (i + 0.5) / n, r = Math.sqrt(1 - y * y);
    out.push(new THREE.Vector3(Math.cos(g * i) * r, y, Math.sin(g * i) * r));
  }
  return out;
}

function tileOf(seeds, v) {
  let best = 0, bd = -2;
  for (let i = 0; i < seeds.length; i++) {
    const d = seeds[i].dot(v);
    if (d > bd) { bd = d; best = i; }
  }
  return best;
}

/* The body: a sphere minus the +x+y+z octant, and the three quarter
   discs that close the cut. */
function body() {
  const g = new THREE.Group();
  const glass = new THREE.MeshPhysicalMaterial({
    color: 0xffffff, transmission: 1, thickness: 0.35, ior: 1.4, roughness: 0.04,
    attenuationColor: new THREE.Color(0xe9eef5), attenuationDistance: 4,
    clearcoat: 1, clearcoatRoughness: 0.05, envMapIntensity: 1.2, side: THREE.FrontSide
  });
  /* the inside wall, seen through the cut: a faint matte shell */
  const wall = new THREE.MeshStandardMaterial({
    color: 0xeef2f6, roughness: 0.9, transparent: true, opacity: 0.45,
    depthWrite: false, side: THREE.BackSide
  });
  const a = function () { return new THREE.SphereGeometry(1, 160, 110, Math.PI, 1.5 * Math.PI, 0, Math.PI); };
  const b = function () { return new THREE.SphereGeometry(1, 56, 56, 0.5 * Math.PI, 0.5 * Math.PI, 0.5 * Math.PI, 0.5 * Math.PI); };
  const inner = [new THREE.Mesh(a(), wall), new THREE.Mesh(b(), wall)];
  inner.forEach(function (m) { m.renderOrder = 1; });
  g.add(new THREE.Mesh(a(), glass), new THREE.Mesh(b(), glass), inner[0], inner[1]);
  /* the cut edges, drawn as structure */
  const lines = new THREE.LineBasicMaterial({ color: 0x7f8c96 });
  for (let k = 0; k < 3; k++) {
    const pts = [];
    for (let i = 0; i <= 96; i++) {
      const t = 0.5 * Math.PI * i / 96, c = Math.cos(t), s = Math.sin(t);
      pts.push(new THREE.Vector3(...[[c, s, 0], [0, c, s], [s, 0, c]][k]));
    }
    g.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), lines));
  }
  const dashes = new THREE.LineDashedMaterial({ color: 0x9aa6af, dashSize: 0.04, gapSize: 0.03 });
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]].forEach(function (v) {
    const l = new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3(...v)]), dashes);
    l.computeLineDistances();
    g.add(l);
  });
  return g;
}

/* The (110) tiles: a Voronoi patchwork on the sphere, drawn per pixel
   so the boundaries are smooth. The lit share is f110; tiles light in a
   fixed order, so changing f110 only adds or removes tiles. */
function tiles(seeds, order) {
  const lit = new Float32Array(seeds.length);
  const mat = new THREE.ShaderMaterial({
    transparent: true, depthWrite: false,
    uniforms: {
      uSeeds: { value: seeds },
      uLit: { value: lit },
      uColor: { value: new THREE.Color(ROLE.surface) }
    },
    vertexShader: `
      varying vec3 vP; varying vec3 vN; varying vec3 vV;
      void main() {
        vP = position;
        vec4 mv = modelViewMatrix * vec4(position, 1.0);
        vN = normalize(normalMatrix * normal);
        vV = normalize(-mv.xyz);
        gl_Position = projectionMatrix * mv;
      }`,
    fragmentShader: `
      #define N ${seeds.length}
      uniform vec3 uSeeds[N]; uniform float uLit[N]; uniform vec3 uColor;
      varying vec3 vP; varying vec3 vN; varying vec3 vV;
      void main() {
        vec3 p = normalize(vP);
        if (p.x > 0.0 && p.y > 0.0 && p.z > 0.0) discard;
        float d1 = -2.0, d2 = -2.0; int b = 0, c = 0;
        for (int i = 0; i < N; i++) {
          float d = dot(p, uSeeds[i]);
          if (d > d1) { d2 = d1; c = b; d1 = d; b = i; }
          else if (d > d2) { d2 = d; c = i; }
        }
        float on = uLit[b], other = uLit[c];
        float edge = 1.0 - smoothstep(0.0012, 0.0045, d1 - d2);
        float fres = pow(1.0 - abs(dot(vN, vV)), 2.0);
        float a = on * (0.13 + 0.22 * fres) + edge * max(on, other) * 0.75;
        if (a < 0.004) discard;
        gl_FragColor = vec4(uColor, a);
        #include <colorspace_fragment>
      }`
  });
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(1.003, 160, 110), mat);
  mesh.renderOrder = 3;
  const rank = new Int32Array(seeds.length);
  order.forEach(function (t, i) { rank[t] = i; });
  function set(f110) {
    const n = Math.round(f110 * seeds.length);
    for (let i = 0; i < seeds.length; i++) lit[i] = rank[i] < n ? 1 : 0;
  }
  return { mesh, set };
}

function cloud(n, band, accept, seed) {
  const r = rng(seed), pos = new Float32Array(3 * n), s = new Float32Array(n);
  const v = new THREE.Vector3();
  let i = 0, guard = 0;
  while (i < n && guard++ < 200 * n) {
    const u = 2 * r() - 1, ph = 2 * Math.PI * r(), q = Math.sqrt(1 - u * u);
    v.set(q * Math.cos(ph), u, q * Math.sin(ph));
    const rad = band[0] === 0 ? band[1] * Math.cbrt(r())
      : Math.cbrt(band[0] ** 3 + (band[1] ** 3 - band[0] ** 3) * r());
    v.multiplyScalar(rad);
    if (inOctant(v.x, v.y, v.z) || !accept(v)) continue;
    pos[3 * i] = v.x; pos[3 * i + 1] = v.y; pos[3 * i + 2] = v.z; s[i] = r();
    i++;
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('aSeed', new THREE.BufferAttribute(s, 1));
  g.setDrawRange(0, 0);
  return g;
}

export function createParticle(host) {
  const S = stage(host, 30);
  const scene = S.scene;
  const root = new THREE.Group();
  scene.add(root);
  S.camera.position.set(3.05, 2.35, 3.05);
  const target = new THREE.Vector3(0, -0.05, 0);
  S.camera.lookAt(target);
  const ctl = orbit(OrbitControls, S, target);
  ctl.minDistance = 2.6; ctl.maxDistance = 8;

  scene.add(new THREE.HemisphereLight(0xffffff, 0xdfe6ec, 0.5));
  const key = new THREE.DirectionalLight(0xffffff, 1.3);
  key.position.set(-3, 5, 2);
  scene.add(key);
  const shadow = contactShadow(1.5, 0.28);
  shadow.position.y = -1.22;
  scene.add(shadow);

  root.add(body());
  /* the opening: the octant that is cut away starts in place and slides
     out along (1,1,1) as the figure first comes into view */
  const lid = new THREE.Group();
  const lidMat = new THREE.MeshPhysicalMaterial({
    color: 0xffffff, transmission: 1, thickness: 0.35, ior: 1.4, roughness: 0.04,
    clearcoat: 1, clearcoatRoughness: 0.05, transparent: true, opacity: 1, side: THREE.DoubleSide
  });
  lid.add(new THREE.Mesh(new THREE.SphereGeometry(1, 48, 48, 0.5 * Math.PI, 0.5 * Math.PI, 0, 0.5 * Math.PI), lidMat));
  const q = function () { return new THREE.CircleGeometry(1, 48, 0, 0.5 * Math.PI); };
  const lz = new THREE.Mesh(q(), lidMat), lx = new THREE.Mesh(q(), lidMat), ly = new THREE.Mesh(q(), lidMat);
  lx.rotation.y = -0.5 * Math.PI; ly.rotation.x = 0.5 * Math.PI;
  lid.add(lz, lx, ly);
  root.add(lid);
  const OPEN = new THREE.Vector3(1, 1, 1).normalize();
  let opened = -1;
  const seeds = fibonacci(TILES);
  const order = seeds.map(function (_, i) { return i; });
  const rr = rng(110);
  for (let i = order.length - 1; i > 0; i--) {
    const j = Math.floor(rr() * (i + 1)); const t = order[i]; order[i] = order[j]; order[j] = t;
  }
  const T = tiles(seeds, order);
  root.add(T.mesh);
  const rank = new Int32Array(TILES);
  order.forEach(function (t, i) { rank[t] = i; });

  /* surface and subsurface beads are generated per tile rank, so a bead
     belongs to the first lit share that includes its tile */
  const tileRank = function (v) { return rank[tileOf(seeds, v.clone().normalize())] / TILES; };
  function rankedCloud(n, band, seed) {
    const g = cloud(n, band, function () { return true; }, seed);
    const p = g.attributes.position, key = new Float32Array(n), idx = [];
    const v = new THREE.Vector3();
    for (let i = 0; i < n; i++) { v.fromBufferAttribute(p, i); key[i] = tileRank(v); idx.push(i); }
    return { g, key };
  }
  const bulkG = cloud(CAP.bulk, BAND.bulk, function () { return true; }, 7);
  const sub = rankedCloud(CAP.subsurface, BAND.subsurface, 11);
  const sur = rankedCloud(CAP.surface, BAND.surface, 13);

  const mats = {
    bulk: beadMaterial(ROLE.bulk, 0.019),
    subsurface: beadMaterial(ROLE.subsurface, 0.05),
    surface: beadMaterial(ROLE.surface, 0.06)
  };
  mats.bulk.uniforms.uJitter.value = 0.0025;
  const pts = {
    bulk: new THREE.Points(bulkG, mats.bulk),
    subsurface: new THREE.Points(sub.g, mats.subsurface),
    surface: new THREE.Points(sur.g, mats.surface)
  };
  Object.keys(pts).forEach(function (k) { pts[k].frustumCulled = false; root.add(pts[k]); });

  /* A ranked cloud draws its first n beads whose tile is lit. Sorting
     the beads once by (lit-by, index) makes that a draw range. */
  function sortBy(c) {
    const n = c.key.length, ix = Array.from({ length: n }, function (_, i) { return i; });
    ix.sort(function (a, b) { return c.key[a] - c.key[b] || a - b; });
    const p = c.g.attributes.position, s = c.g.attributes.aSeed;
    const P = new Float32Array(3 * n), Q = new Float32Array(n), K = new Float32Array(n);
    ix.forEach(function (j, i) {
      P[3 * i] = p.getX(j); P[3 * i + 1] = p.getY(j); P[3 * i + 2] = p.getZ(j);
      Q[i] = s.getX(j); K[i] = c.key[j];
    });
    c.g.setAttribute('position', new THREE.BufferAttribute(P, 3));
    c.g.setAttribute('aSeed', new THREE.BufferAttribute(Q, 1));
    c.key = K;
  }
  sortBy(sub); sortBy(sur);

  /* Beads of a ranked cloud are scattered over the lit tiles only: the
     lit tiles are the first `lit` ranks, so take beads from the sorted
     list in a stride that spreads them over all lit tiles evenly. */
  function litCount(c, f110) {
    let n = 0;
    const lim = Math.round(f110 * TILES) / TILES - 1e-9;
    while (n < c.key.length && c.key[n] < lim) n++;
    return n;
  }
  function spread(c, avail, want) {
    /* reorder the first `avail` beads so that any prefix is spread over
       all lit tiles: a fixed pseudo-random permutation of that prefix */
    const g = c.g;
    if (c._avail === avail) return;
    c._avail = avail;
    const p = g.attributes.position.array, s = g.attributes.aSeed.array;
    const base = c._base || (c._base = { p: p.slice(), s: s.slice() });
    const r = rng(avail * 7919 + 3);
    const ix = Array.from({ length: avail }, function (_, i) { return i; });
    for (let i = avail - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); const t = ix[i]; ix[i] = ix[j]; ix[j] = t; }
    for (let i = 0; i < avail; i++) {
      const j = ix[i];
      p[3 * i] = base.p[3 * j]; p[3 * i + 1] = base.p[3 * j + 1]; p[3 * i + 2] = base.p[3 * j + 2]; s[i] = base.s[j];
    }
    g.attributes.position.needsUpdate = true;
    g.attributes.aSeed.needsUpdate = true;
  }

  const now = { bulk: 0, subsurface: 0, surface: 0 };
  let from = Object.assign({}, now), to = Object.assign({}, now), t0 = -1, f110 = -1;
  const api = { dotUmol: DOT_UMOL, shown: now, stage: S };

  api.set = function (view) {
    if (view.f110 !== f110) {
      f110 = view.f110;
      T.set(f110);
      spread(sub, litCount(sub, f110));
      spread(sur, litCount(sur, f110));
    }
    const want = {
      bulk: view.pools.bulk / DOT_UMOL,
      subsurface: (view.pools.L1_subbridging + view.pools.subsurface_L2_4) / DOT_UMOL,
      surface: (view.pools.bridging + view.pools.reconstructed_row + view.pools.basal_L1) / DOT_UMOL
    };
    api.want = want;
    from = Object.assign({}, now);
    to = {
      bulk: Math.min(CAP.bulk, Math.round(want.bulk)),
      subsurface: Math.min(sub._avail, Math.round(want.subsurface)),
      surface: Math.min(sur._avail, Math.round(want.surface))
    };
    t0 = view.animate === false ? -2 : -1;
    S.dirty = true;
  };

  const scaleOf = function () {
    return S.renderer.getPixelRatio() * S.size / (2 * Math.tan(THREE.MathUtils.degToRad(S.camera.fov / 2)));
  };
  S.tick = function (t, still) {
    let moving = false;
    if (t0 === -1) t0 = t;
    if (t0 === -2 || still) t0 = t - 10;
    const k = ease((t - t0) / 1.2);
    Object.keys(pts).forEach(function (n) {
      now[n] = Math.round(from[n] + (to[n] - from[n]) * k);
      pts[n].geometry.setDrawRange(0, now[n]);
      mats[n].uniforms.uScale.value = scaleOf();
      mats[n].uniforms.uTime.value = t;
    });
    if (k < 1) moving = true;
    if (lid.visible) {
      if (opened < 0) opened = still ? t - 10 : t + 0.4;
      const o = ease((t - opened) / 1.6);
      lid.position.copy(OPEN).multiplyScalar(1.4 * o);
      lidMat.opacity = 1 - o;
      if (o >= 1) lid.visible = false;
      moving = true;
    }
    const a = still ? 0 : sway(S, t, 0.42);
    if (a !== null) { root.rotation.y = a; moving = true; }
    if (!still && ctl.update()) moving = true;
    return moving || mats.bulk.uniforms.uJitter.value > 0;
  };
  return api;
}
