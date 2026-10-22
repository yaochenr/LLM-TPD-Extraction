// White Molecular Universe — background scene for the AutoTPD+ page.
// Molecule geometry is generated offline by scripts/build_molecule_universe.py.
import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js';

const canvas = document.getElementById('universe');
const dialog = document.getElementById('download');
const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');
const small = window.matchMedia('(max-width: 640px)').matches;

const FOV = 45;
const DEPTH = 120;          // how far ahead molecules live
const BEHIND = 8;           // how far behind the camera before wrapping
const CLEAR_RADIUS = 4.2;   // keep the axis clear so near molecules never cover the title
const COUNT = small ? 60 : 100;

// Without WebGL this throws and the page simply stays plain white.
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, powerPreference: 'high-performance' });
renderer.setClearColor(0xffffff, 1);
let pixelRatio = Math.min(window.devicePixelRatio || 1, 2);
renderer.setPixelRatio(pixelRatio);

const scene = new THREE.Scene();
scene.fog = new THREE.Fog(0xffffff, DEPTH * 0.25, DEPTH * 0.98);
// Soft studio light: shading on the atom spheres is what makes the molecules read as 3D.
scene.add(new THREE.HemisphereLight(0xffffff, 0xa8a8a8, 1.7));
const sun = new THREE.DirectionalLight(0xffffff, 1.6);
sun.position.set(-0.6, 1, 0.8);
scene.add(sun);
const camera = new THREE.PerspectiveCamera(FOV, 1, 0.1, DEPTH + 20);

// Half-extents of the box molecules are scattered in, sized to the frustum at mid depth.
const box = { x: 1, y: 1 };
function resize() {
  const w = window.innerWidth, h = window.innerHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  const halfH = Math.tan(THREE.MathUtils.degToRad(FOV / 2)) * DEPTH * 0.45;
  box.y = halfH;
  box.x = halfH * Math.min(Math.max(camera.aspect, 0.6), 2.2);
}

const rand = (a, b) => a + Math.random() * (b - a);

// Monochrome ball-and-stick: one shared low-poly sphere, shaded grey by element class
// (C, N, O, S/P, halogen, other) and sized slightly larger for heteroatoms.
const atomGeometry = new THREE.IcosahedronGeometry(1, 1);
const ELEMENT_GREY = [0.82, 0.5, 0.36, 0.28, 0.32, 0.42];
const ELEMENT_RADIUS = [0.17, 0.21, 0.21, 0.24, 0.23, 0.22];
const _m = new THREE.Matrix4(), _p = new THREE.Vector3(), _s = new THREE.Vector3(), _q = new THREE.Quaternion();
const _c = new THREE.Color();

function buildAtoms(mol, scale) {
  const n = mol.e.length;
  const material = new THREE.MeshLambertMaterial({ transparent: true, opacity: 0 });
  const mesh = new THREE.InstancedMesh(atomGeometry, material, n);
  for (let i = 0; i < n; i++) {
    const e = mol.e[i];
    _p.set(mol.a[i * 3] / scale, mol.a[i * 3 + 1] / scale, mol.a[i * 3 + 2] / scale);
    _s.setScalar(ELEMENT_RADIUS[e]);
    mesh.setMatrixAt(i, _m.compose(_p, _q, _s));
    mesh.setColorAt(i, _c.setScalar(ELEMENT_GREY[e]));
  }
  return mesh;
}

function buildGeometry(coords, scale) {
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(coords, (v) => v / scale), 3));
  g.computeBoundingSphere();
  return g;
}

const items = [];
let camZ = 0;

function place(item, depth) {
  // Uniform in the box, but outside a cylinder around the view axis.
  const minR = CLEAR_RADIUS + item.radius;
  let x, y;
  do {
    x = rand(-box.x, box.x);
    y = rand(-box.y, box.y);
  } while (Math.hypot(x, y * 1.4) < minR);
  item.base.set(x, y, camZ - depth);
}

function build(data) {
  const mols = data.molecules.slice().sort(() => Math.random() - 0.5).slice(0, COUNT);
  mols.forEach((mol, i) => {
    const geometry = buildGeometry(mol.s, data.scale);
    const size = rand(0.3, 0.52);
    const color = new THREE.Color().setScalar(rand(0.08, 0.3));
    const material = new THREE.LineBasicMaterial({ color, transparent: true, depthWrite: false, opacity: 0 });
    const obj = new THREE.LineSegments(geometry, material);
    obj.scale.setScalar(size);
    const atoms = buildAtoms(mol, data.scale);
    obj.add(atoms);
    scene.add(obj);

    const item = {
      obj, material, atoms,
      radius: geometry.boundingSphere.radius * size,
      peak: rand(0.65, 0.95),
      base: new THREE.Vector3(),
      driftZ: rand(-0.12, 0.12),
      // Bounded wander: two slow, incommensurate sines per axis.
      wander: { a: rand(0.5, 1.4), w: rand(0.03, 0.07), p: rand(0, Math.PI * 2), q: rand(0, Math.PI * 2) },
      rot0: new THREE.Euler(rand(0, Math.PI * 2), rand(0, Math.PI * 2), rand(0, Math.PI * 2)),
      spin: new THREE.Vector3(rand(-1, 1), rand(-1, 1), rand(-1, 1)).multiplyScalar(0.12),
    };
    // Spread the initial cloud through the whole depth range.
    place(item, (i + Math.random()) / mols.length * DEPTH);
    items.push(item);
  });
}

const smooth = (e0, e1, x) => {
  const t = Math.min(Math.max((x - e0) / (e1 - e0), 0), 1);
  return t * t * (3 - 2 * t);
};

const ndc = new THREE.Vector3();

function update(t, dt) {
  // Camera: wheel-driven travel plus subtle pointer parallax.
  velocity *= Math.exp(-dt * 2.4);
  camZ -= velocity * dt;
  parallax.x += (pointer.x * 1.4 - parallax.x) * Math.min(dt * 2.5, 1);
  parallax.y += (pointer.y * 0.9 - parallax.y) * Math.min(dt * 2.5, 1);
  camera.position.set(parallax.x, parallax.y, camZ);
  camera.lookAt(parallax.x * 0.3, parallax.y * 0.3, camZ - 60);

  for (const it of items) {
    // Wrap molecules through the visible depth range so space feels endless.
    it.base.z += it.driftZ * dt;
    let depth = camZ - it.base.z;
    if (depth < -BEHIND) { place(it, depth + DEPTH + BEHIND); depth = camZ - it.base.z; }
    else if (depth > DEPTH) { place(it, depth - DEPTH - BEHIND); depth = camZ - it.base.z; }

    const w = it.wander;
    it.obj.position.set(
      it.base.x + (Math.sin(t * w.w + w.p) + 0.5 * Math.sin(t * w.w * 2.3 + w.q)) * w.a,
      it.base.y + (Math.sin(t * w.w * 0.8 + w.q) + 0.5 * Math.sin(t * w.w * 1.9 + w.p)) * w.a,
      it.base.z + Math.sin(t * w.w * 1.3 + w.p) * w.a * 0.5,
    );
    it.obj.rotation.set(it.rot0.x + t * it.spin.x, it.rot0.y + t * it.spin.y, it.rot0.z + t * it.spin.z);

    // Fade in from the far fog and out just before passing the camera (soft depth of field),
    // and recede behind the title block in the middle of the screen.
    const d = camZ - it.obj.position.z;
    ndc.copy(it.obj.position).project(camera);
    const centre = Math.hypot(ndc.x / 0.55, ndc.y / 0.42);
    it.material.opacity = it.peak * smooth(DEPTH, DEPTH * 0.6, d) * smooth(1.5, 10, d)
      * (0.18 + 0.82 * smooth(0.45, 1.05, centre));
    it.atoms.material.opacity = it.material.opacity;
    it.obj.visible = it.material.opacity > 0.003;
  }
}

// ---- input ----
let velocity = 0;
const pointer = { x: 0, y: 0 };
const parallax = { x: 0, y: 0 };
const MAX_SPEED = 90;
const push = (amount) => {
  velocity = Math.max(-MAX_SPEED, Math.min(MAX_SPEED, velocity + amount));
};

window.addEventListener('wheel', (e) => {
  if (dialog && dialog.open) return;
  const unit = e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? window.innerHeight : 1;
  push(e.deltaY * unit * 0.16);
}, { passive: true });

let touchY = null;
window.addEventListener('touchstart', (e) => { touchY = e.touches[0].clientY; }, { passive: true });
window.addEventListener('touchmove', (e) => {
  if (touchY === null || (dialog && dialog.open)) return;
  const y = e.touches[0].clientY;
  push((touchY - y) * 0.12);
  touchY = y;
}, { passive: true });
window.addEventListener('touchend', () => { touchY = null; }, { passive: true });

window.addEventListener('keydown', (e) => {
  if (dialog && dialog.open) return;
  if (e.key === 'ArrowDown' || e.key === 'PageDown') push(14);
  if (e.key === 'ArrowUp' || e.key === 'PageUp') push(-14);
});

window.addEventListener('pointermove', (e) => {
  if (e.pointerType !== 'mouse') return;
  pointer.x = (e.clientX / window.innerWidth) * 2 - 1;
  pointer.y = -((e.clientY / window.innerHeight) * 2 - 1);
});

// ---- loop with performance-based quality adjustment ----
const clock = new THREE.Clock();
let time = rand(0, 200), running = false, frames = 0, slowFrames = 0;

function frame() {
  if (!running) return;
  const dt = Math.min(clock.getDelta(), 0.1);
  time += dt;
  update(time, dt);
  renderer.render(scene, camera);

  frames++;
  if (dt > 1 / 40) slowFrames++;
  if (frames === 120) {
    if (slowFrames > 60) degrade();
    frames = slowFrames = 0;
  }
  requestAnimationFrame(frame);
}

function degrade() {
  if (pixelRatio > 1) {
    pixelRatio = 1;
    renderer.setPixelRatio(1);
    resize();
  } else if (items.length > 50) {
    // Drop every other molecule, but keep the scene populated.
    for (let i = items.length - 1; i >= 0; i -= 2) {
      scene.remove(items[i].obj);
      items[i].obj.geometry.dispose();
      items[i].material.dispose();
      items[i].atoms.material.dispose();
      items[i].atoms.dispose();
      items.splice(i, 1);
    }
  }
}

function start() {
  update(time, 0);
  renderer.render(scene, camera);
  if (reducedMotion.matches) {
    running = false;
    return;
  }
  if (running) return;
  running = true;
  clock.getDelta();
  requestAnimationFrame(frame);
}

window.addEventListener('resize', () => { resize(); if (!running) renderer.render(scene, camera); });
reducedMotion.addEventListener('change', start);

resize();
fetch(new URL('molecules.json?v=3d-2', import.meta.url))
  .then((r) => r.json())
  .then((data) => {
    build(data);
    start();
    canvas.classList.add('ready');
  });
