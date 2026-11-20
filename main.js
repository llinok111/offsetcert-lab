// Renderer + interaction only. All geometry, volume, certificates and the
// containment oracle come from src/core.mjs (exact BigInt rationals).
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import * as core from './src/core.mjs';

const MM = core.q(10); // samples are authored at unit scale; show millimetres
const shapes = core.shapes().map(p => core.scalePoly(p, MM));

const view = document.getElementById('view');
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(view.clientWidth, view.clientHeight);
view.appendChild(renderer.domElement);

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0d0f14);
scene.fog = new THREE.Fog(0x0d0f14, 55, 130);

const camera = new THREE.PerspectiveCamera(42, view.clientWidth / view.clientHeight, 0.1, 500);
camera.position.set(26, 21, 30);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
controls.target.set(0, 0, 0);

scene.add(new THREE.HemisphereLight(0xbfd4ff, 0x141a22, 1.15));
const key = new THREE.DirectionalLight(0xffffff, 2.1); key.position.set(18, 30, 22); scene.add(key);
const rim = new THREE.DirectionalLight(0xff9f6a, 1.0); rim.position.set(-24, 10, -18); scene.add(rim);

const grid = new THREE.GridHelper(80, 40, 0x2a3342, 0x1a212c);
grid.position.y = -18;
scene.add(grid);

const group = new THREE.Group();
scene.add(group);

const material = new THREE.MeshStandardMaterial({
  color: 0xff9f43, metalness: 0.18, roughness: 0.52, side: THREE.DoubleSide,
  polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1,
});
const edgeMaterial = new THREE.LineBasicMaterial({ color: 0x141a22, transparent: true, opacity: 0.85 });
const wireMaterial = new THREE.LineBasicMaterial({ color: 0x79c0ff, transparent: true, opacity: 0.22 });

let mesh = null, edges = null, wire = null, law = null;

function triangulate(poly) {
  const pos = [];
  for (const face of poly.faces) {
    for (let i = 1; i + 1 < face.length; i++) {
      for (const vi of [face[0], face[i], face[i + 1]]) {
        const v = poly.verts[vi];
        pos.push(core.qnum(v[0]), core.qnum(v[1]), core.qnum(v[2]));
      }
    }
  }
  return new Float32Array(pos);
}

function edgeSegments(poly) {
  const seen = new Set(), pos = [];
  for (const face of poly.faces) {
    for (let i = 0; i < face.length; i++) {
      const a = face[i], b = face[(i + 1) % face.length];
      const k = a < b ? `${a}_${b}` : `${b}_${a}`;
      if (seen.has(k)) continue;
      seen.add(k);
      const va = poly.verts[a], vb = poly.verts[b];
      pos.push(core.qnum(va[0]), core.qnum(va[1]), core.qnum(va[2]));
      pos.push(core.qnum(vb[0]), core.qnum(vb[1]), core.qnum(vb[2]));
    }
  }
  return new Float32Array(pos);
}

function disposeAll() {
  for (const o of [mesh, edges, wire]) {
    if (!o) continue;
    group.remove(o);
    o.geometry.dispose();
  }
  mesh = edges = wire = null;
}

function rebuild(shapeIdx, d) {
  const base = shapes[shapeIdx];
  const poly = core.offsetPoly(base, d);

  disposeAll();
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(triangulate(poly), 3));
  geo.computeVertexNormals();
  mesh = new THREE.Mesh(geo, material);
  group.add(mesh);

  const eg = new THREE.BufferGeometry();
  eg.setAttribute('position', new THREE.BufferAttribute(edgeSegments(poly), 3));
  edges = new THREE.LineSegments(eg, edgeMaterial);
  group.add(edges);

  const wg = new THREE.BufferGeometry();
  wg.setAttribute('position', new THREE.BufferAttribute(triangulate(poly), 3));
  wire = new THREE.LineSegments(new THREE.WireframeGeometry(geo), wireMaterial);
  wg.dispose();
  group.add(wire);

  const wt = core.watertight(poly);
  const cvx = core.convexCertificate(poly);
  const vol = core.volume(poly);
  const slope = core.evalPoly([law.coeffs[1], core.qmul(core.q(2), law.coeffs[2]), core.qmul(core.q(3), law.coeffs[3])], d);

  set('vol', core.qnum(vol).toFixed(4) + ' mm³', 'big');
  set('frac', core.qstr(vol));
  set('slope', core.qnum(slope).toFixed(4) + ' mm³/mm');
  set('coeffs', law.coeffs.map(c => core.qstr(c)).join(' , '));
  set('wt', wt.ok ? `true (${wt.edges} edges)` : `false b=${wt.boundary} nm=${wt.nonManifold}`, wt.ok ? 'ok' : 'bad');
  set('cvx', cvx.ok ? `true (violations 0, min tight ${cvx.minTight})` : `false (${cvx.violations})`, cvx.ok ? 'ok' : 'bad');
  const tris = poly.faces.reduce((n, f) => n + f.length - 2, 0);
  set('topo', `${wt.edges} / ${tris}`);

  // Independent oracle over deterministic lattice points, in the browser.
  const pts = core.makePoints(20261121, 2000, 20);
  let dis = 0;
  for (const p of pts) {
    const h = core.containsHalfspace(poly, p), t = core.tetraOracle(poly, p);
    if (h.boundary !== t.boundary || (!h.boundary && h.inside !== t.inside)) dis++;
  }
  set('oracle', `${dis} disagreements / 2000`, dis === 0 ? 'ok' : 'bad');

  // Exact scale equivariance at the current d.
  let exact = 0;
  const scales = ['1/1000', '1/100', '1', '100', '1000'];
  for (const ss of scales) {
    const s = core.q(ss);
    const a = core.offsetPoly(core.scalePoly(base, s), core.qmul(s, d));
    const b = core.scalePoly(poly, s);
    if (core.polyEqual(a, b)) exact++;
  }
  set('scale', `${exact}/${scales.length} exact`, exact === scales.length ? 'ok' : 'bad');
}

function set(id, text, cls = '') {
  const el = document.getElementById(id);
  el.textContent = text;
  el.className = cls;
}

const shapeSel = document.getElementById('shape');
const dRange = document.getElementById('d');

function update(rebuildLaw = false) {
  const idx = Number(shapeSel.value);
  const d = core.q(dRange.value);
  document.getElementById('dval').textContent = Number(dRange.value).toFixed(2) + ' mm';
  if (rebuildLaw) {
    law = core.cubicLaw(shapes[idx]);
    frameCamera(idx);
  }
  rebuild(idx, d);
}

function frameCamera(idx) {
  const p = shapes[idx];
  let r = 1;
  for (const v of p.verts) for (let k = 0; k < 3; k++) r = Math.max(r, Math.abs(core.qnum(v[k])));
  const dist = r * 3.1;
  camera.position.set(dist * 0.72, dist * 0.58, dist * 0.82);
  camera.near = 0.05; camera.far = dist * 20;
  camera.updateProjectionMatrix();
  controls.target.set(0, 0, 0);
  controls.update();
  grid.position.y = -r * 1.55;
}

shapeSel.addEventListener('change', () => update(true));
dRange.addEventListener('input', () => update(false));
addEventListener('resize', () => {
  camera.aspect = view.clientWidth / view.clientHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(view.clientWidth, view.clientHeight);
});

update(true);
renderer.setAnimationLoop(() => {
  controls.update();
  renderer.render(scene, camera);
});