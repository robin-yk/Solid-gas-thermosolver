/* Shared pieces of the two 3D figures: renderer, studio light, seeded
   random numbers, and a render loop that sleeps while the figure is off
   screen or its workspace is hidden. Nothing here knows any physics. */
import * as THREE from 'three';
import studio from './assets/studio.jpg';

export const ROLE = {                 /* the figure kit's role hues */
  surface: 0x0072B2, subsurface: 0xD55E00, bulk: 0x6A51A3,
  extended: 0x777777, structure: 0xAAAAAA
};

export function rng(seed) {           /* mulberry32 */
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function webglAvailable() {
  try {
    const c = document.createElement('canvas');
    return !!(window.WebGLRenderingContext &&
      (c.getContext('webgl2') || c.getContext('webgl')));
  } catch (e) { return false; }
}

export const REDUCED_MOTION = typeof matchMedia === 'function' &&
  matchMedia('(prefers-reduced-motion: reduce)').matches;

/* A square stage inside `host`: the canvas takes the host's width. */
export function stage(host, fov) {
  const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.0;
  renderer.domElement.className = 'v3d-canvas';
  host.appendChild(renderer.domElement);

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0xffffff);
  const camera = new THREE.PerspectiveCamera(fov, 1, 0.01, 100);

  const pmrem = new THREE.PMREMGenerator(renderer);
  new THREE.TextureLoader().load(studio, function (tex) {
    tex.mapping = THREE.EquirectangularReflectionMapping;
    tex.colorSpace = THREE.SRGBColorSpace;
    scene.environment = pmrem.fromEquirectangular(tex).texture;
    tex.dispose();
    pmrem.dispose();
    S.dirty = true;
  });

  const S = {
    renderer, scene, camera, host, dirty: true, tick: null, visible: true,
    size: 0, t0: performance.now()
  };

  function fit() {
    const w = Math.round(host.clientWidth);
    if (!w || w === S.size) return;
    S.size = w;
    renderer.setSize(w, w, true);
    camera.aspect = 1;
    camera.updateProjectionMatrix();
    S.dirty = true;
  }
  if (typeof ResizeObserver === 'function') new ResizeObserver(fit).observe(host);
  if (typeof IntersectionObserver === 'function') {
    new IntersectionObserver(function (e) { S.visible = e[0].isIntersecting; })
      .observe(host);
  }

  function frame() {
    requestAnimationFrame(frame);
    if (!host.clientWidth || !S.visible || document.hidden) return;
    fit();
    const t = (performance.now() - S.t0) / 1000;
    /* a still page (the browser gate sets this) renders only on change,
       at the rest pose, with every transition already finished */
    const still = !!window.V3D_STILL;
    if (still && !S.dirty) return;
    const moving = S.tick ? S.tick(still ? 1e4 : t, still) : false;
    if (moving || S.dirty) {
      renderer.render(scene, camera);
      S.dirty = false;
      S.frames = (S.frames || 0) + 1;
    }
  }
  requestAnimationFrame(frame);
  return S;
}

export function ease(t) { return t < 0 ? 0 : t > 1 ? 1 : 1 - Math.pow(1 - t, 3); }

/* A soft round contact shadow under an object, as a textured plane. */
export function contactShadow(radius, strength) {
  const c = document.createElement('canvas');
  c.width = c.height = 128;
  const g = c.getContext('2d');
  const grd = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  grd.addColorStop(0, 'rgba(24,55,53,' + strength + ')');
  grd.addColorStop(0.55, 'rgba(24,55,53,' + strength * 0.35 + ')');
  grd.addColorStop(1, 'rgba(24,55,53,0)');
  g.fillStyle = grd;
  g.fillRect(0, 0, 128, 128);
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  const m = new THREE.Mesh(new THREE.PlaneGeometry(radius * 2, radius * 2),
    new THREE.MeshBasicMaterial({ map: tex, transparent: true, depthWrite: false }));
  m.rotation.x = -Math.PI / 2;
  return m;
}

/* Round, lit beads for point clouds: each sprite shades itself as a
   small sphere, so thousands of points read as objects, not pixels. */
export function beadMaterial(colour, size) {
  return new THREE.ShaderMaterial({
    uniforms: {
      uColor: { value: new THREE.Color(colour) },
      uSize: { value: size },
      uScale: { value: 1 },
      uTime: { value: 0 },
      uJitter: { value: 0 }
    },
    vertexShader: `
      attribute float aSeed;
      uniform float uSize, uScale, uTime, uJitter;
      varying float vFog;
      void main() {
        vec3 p = position;
        float s = aSeed * 6.2831853;
        p += uJitter * vec3(sin(uTime * 1.3 + s), sin(uTime * 1.7 + 2.0 * s), sin(uTime * 1.1 + 3.0 * s));
        vec4 mv = modelViewMatrix * vec4(p, 1.0);
        gl_Position = projectionMatrix * mv;
        gl_PointSize = uSize * uScale / -mv.z;
        vFog = clamp((-mv.z - 2.4) / 3.0, 0.0, 1.0);
      }`,
    fragmentShader: `
      uniform vec3 uColor;
      varying float vFog;
      void main() {
        vec2 q = gl_PointCoord * 2.0 - 1.0;
        float r2 = dot(q, q);
        if (r2 > 1.0) discard;
        vec3 n = vec3(q.x, -q.y, sqrt(1.0 - r2));
        vec3 l = normalize(vec3(-0.45, 0.6, 0.66));
        float diff = 0.55 + 0.45 * max(dot(n, l), 0.0);
        float spec = pow(max(dot(reflect(-l, n), vec3(0.0, 0.0, 1.0)), 0.0), 24.0);
        vec3 c = uColor * diff + vec3(0.9) * spec * 0.55;
        c = mix(c, vec3(1.0), vFog * 0.35);
        gl_FragColor = vec4(c, 1.0);
        #include <colorspace_fragment>
      }`
  });
}

export function orbit(OrbitControls, S, target) {
  const c = new OrbitControls(S.camera, S.renderer.domElement);
  c.target.copy(target);
  c.enableDamping = true;
  c.dampingFactor = 0.08;
  c.enablePan = false;
  c.addEventListener('change', function () { S.dirty = true; });
  c.addEventListener('start', function () { S.interact = performance.now(); S.grabbed = true; });
  c.addEventListener('end', function () { S.interact = performance.now(); S.grabbed = false; });
  return c;
}

/* The idle sway: a slow swing about the vertical, paused while the user
   holds the figure and for a few seconds after. */
export function sway(S, t, amp) {
  if (REDUCED_MOTION || S.grabbed) return null;
  if (S.interact && performance.now() - S.interact < 6000) return null;
  return amp * Math.sin(t * 0.22);
}
