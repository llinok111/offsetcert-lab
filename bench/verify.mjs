// Assertion suite. Every assertion prints the value it actually measured.
// Exit 0 only if all assertions pass.
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import * as core from '../src/core.mjs';
import * as naive from '../baseline/naive.mjs';
import { buildPayload, serialize, fnv1a } from '../tools/gen-sample.mjs';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
let pass = 0, fail = 0;
const failed = [];
function ok(cond, label, detail) {
  if (cond) { pass++; console.log(`PASS  ${label}  ::  ${detail}`); }
  else { fail++; failed.push(label); console.log(`FAIL  ${label}  ::  ${detail}`); }
}
const eq = (a, b) => a === b;
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

const shapes = core.shapes();
const OFFSETS = ['0.1', '0.25', '0.5', '0.75', '1', '2'];
const SCALES = ['1/1000000', '1/1000', '1', '1000', '1000000'];

// ---- A1: inputs ----------------------------------------------------------
ok(shapes.length === 3 && shapes.every(p => p.verts.length > 0), 'A1 3 sample shapes', `shapes=${shapes.map(p => p.name).join(',')}`);
const cube = shapes[0], pyramid = shapes[1], prism = shapes[2];
for (const p of shapes) {
  const w = core.watertight(p), c = core.convexCertificate(p);
  ok(w.ok && c.ok, `A1 watertight+convex input ${p.name}`, `edges=${w.edges} boundary=${w.boundary} violations=${c.violations} minTight=${c.minTight}`);
}

// ---- A2: known volumes ---------------------------------------------------
ok(core.qeq(core.volume(cube), core.q(8)), 'A2 cube V(0)=8', core.qstr(core.volume(cube)));
ok(core.qeq(core.volume(pyramid), core.q(1)), 'A2 pyramid V(0)=1', core.qstr(core.volume(pyramid)));
ok(core.qeq(core.volume(prism), core.q(95, 12)), 'A2 prism V(0)=95/12', core.qstr(core.volume(prism)));

// ---- A3/A4: every offset body isValid -----------------------------------
let offsetsOk = 0, minEdges = Infinity;
for (const p of shapes) {
  for (const ds of OFFSETS) {
    const off = core.offsetPoly(p, core.q(ds));
    const w = core.watertight(off), c = core.convexCertificate(off);
    if (w.ok && c.ok && core.qsign(core.volume(off)) > 0) offsetsOk++;
    minEdges = Math.min(minEdges, w.edges);
    if (!w.ok || !c.ok) console.log(`      offset detail ${p.name} d=${ds} wt=${JSON.stringify(w)} cvx=${JSON.stringify(c)}`);
  }
}
ok(offsetsOk === 18, 'A3 18/18 offsets watertight+convex', `ok=${offsetsOk}/18 minEdges=${minEdges}`);
ok(minEdges >= 8, 'A3 topology preserved', `min edges across 18 offsets=${minEdges}`);

// ---- A5/A6: cubic volume law --------------------------------------------
const closed = { cube: ['8', '24', '24', '8'], pyramid: ['1', '9', '27', '27'] };
const laws = [];
for (const p of shapes) {
  const law = core.cubicLaw(p);
  laws.push(law);
  ok(law.mismatches === 0, `A5 cubic law ${p.name}`, `coeffs=[${law.coeffs.map(core.qstr).join(', ')}] recheck=${law.points} mismatches=${law.mismatches} maxDiff=${law.maxDiffNum}`);
  if (closed[p.name]) {
    ok(law.coeffs.map(core.qstr).join(',') === closed[p.name].join(','), `A5 closed form ${p.name}`, `got=[${law.coeffs.map(core.qstr)}] want=[${closed[p.name]}]`);
  }
}
const prismLaw = laws[2];
ok(core.qeq(prismLaw.coeffs[0], core.q(95, 12)), 'A5 prism a0 exact', core.qstr(prismLaw.coeffs[0]));
ok(laws.every(l => l.maxDiffNum === 0), 'A6 0 volume error on 63 recheck points', `points=${laws.map(l => l.points).join('+')} maxDiff=${laws.map(l => l.maxDiffNum).join(',')}`);

// ---- A7: monotonicity ----------------------------------------------------
for (const l of laws) {
  ok(l.monotone && core.qsign(l.minSlope) > 0, 'A7 V\'(d)>0', `minSlope=${core.qnum(l.minSlope).toFixed(6)} (exact ${core.qstr(l.minSlope)})`);
}

// ---- A8: exact scale equivariance ---------------------------------------
let scaleOk = 0, scaleTotal = 0, scaleMinVol = Infinity;
for (const p of shapes) {
  for (const ss of SCALES) {
    const s = core.q(ss);
    const d = core.q(1, 10);
    const a = core.offsetPoly(core.scalePoly(p, s), core.qmul(s, d));
    const b = core.scalePoly(core.offsetPoly(p, d), s);
    scaleTotal++;
    if (core.polyEqual(a, b) && core.watertight(a).ok) scaleOk++;
    scaleMinVol = Math.min(scaleMinVol, Math.abs(core.qnum(core.volume(a))));
  }
}
ok(scaleOk === scaleTotal, 'A8 scale equivariance 5 scales x 3 shapes', `exact=${scaleOk}/${scaleTotal} smallestVolume=${scaleMinVol.toExponential(3)}`);

// ---- A9: independent tetra oracle ---------------------------------------
let oraclePoints = 0, disagreements = 0, boundary = 0, insideN = 0;
for (const p of shapes) {
  const pts = core.makePoints(20261121, 20000);
  for (const pt of pts) {
    oraclePoints++;
    const h = core.containsHalfspace(p, pt);
    const t = core.tetraOracle(p, pt);
    if (h.boundary !== t.boundary) { disagreements++; continue; }
    if (h.boundary) { boundary++; continue; }
    if (h.inside !== t.inside) { disagreements++; continue; }
    if (h.inside) insideN++;
  }
}
ok(disagreements === 0, 'A9 halfspace vs tetra oracle', `points=${oraclePoints} disagreements=${disagreements} boundary=${boundary} inside=${insideN} outside=${oraclePoints - boundary - insideN}`);
ok(boundary > 0, 'A10 boundary points classified separately', `boundary=${boundary} (${(100 * boundary / oraclePoints).toFixed(2)}%)`);

// exact boundary witnesses: vertices and edge midpoints
let bWitness = 0, bTotal = 0;
for (const p of shapes) {
  for (const pt of core.boundaryPoints(p)) {
    bTotal++;
    const h = core.containsHalfspace(p, pt);
    const t = core.tetraOracle(p, pt);
    if (h.boundary && t.boundary) bWitness++;
  }
}
ok(bWitness === bTotal, 'A10 exact boundary witnesses', `boundary=${bWitness}/${bTotal}`);

// ---- A11: baseline scale window -----------------------------------------
const floatShapes = shapes.map(naive.toFloatPoly);
let baseSmallDegenerate = 0, exactSmallOk = 0;
for (const p of shapes) {
  const small = core.scalePoly(p, core.q(1, 1000000));
  const fp = naive.toFloatPoly(small);
  const st = naive.faceStats(fp);
  if (st.degenerate === st.total) baseSmallDegenerate++;
  const off = core.offsetPoly(small, core.q(1, 10000000));
  if (core.watertight(off).ok) exactSmallOk++;
}
ok(baseSmallDegenerate === 3, 'A11 baseline rejects 1e-6 geometry', `shapes with all faces < AREA_EPS=${baseSmallDegenerate}/3 (AREA_EPS=${naive.AREA_EPS})`);
ok(exactSmallOk === 3, 'A11 exact core certifies 1e-6 geometry', `watertight=${exactSmallOk}/3, zero thresholds`);

// ---- A12: baseline is accurate where it is usable -----------------------
const unitNaive = naive.analyzeNaive(naive.toFloatPoly(cube), 0.1);
const exactV = core.qnum(core.volume(core.offsetPoly(cube, core.q(1, 10))));
const vErr = unitNaive.usable ? rel(unitNaive.volume, exactV) : NaN;
ok(unitNaive.usable && vErr < 1e-9, 'A12 baseline @unit scale agrees', `naive=${unitNaive.volume} exact=${exactV} relErr=${vErr.toExponential(2)}`);
const naiveCubeCoeffs = naive.cubicNaive(naive.toFloatPoly(cube));
const coeffErr = Math.max(...naiveCubeCoeffs.map((c, i) => rel(c, core.qnum(laws[0].coeffs[i]))));
ok(coeffErr < 1e-9, 'A12 baseline cubic @unit scale agrees', `coeffs=[${naiveCubeCoeffs.map(x => x.toFixed(6))}] maxRelErr=${coeffErr.toExponential(2)}`);
const smallCube = naive.analyzeNaive(naive.toFloatPoly(core.scalePoly(cube, core.q(1, 1000))), 0.0001);
ok(smallCube.usable, 'A12 baseline usable at 1e-3', `faces=${smallCube.totalFaces} minArea=${smallCube.minArea.toExponential(3)} volume=${smallCube.volume.toExponential(3)}`);
const tinyCube = naive.analyzeNaive(naive.toFloatPoly(core.scalePoly(cube, core.q(1, 1000000))), 0.0000001);
ok(!tinyCube.usable && tinyCube.degenerateFaces === tinyCube.totalFaces, 'A12 baseline refuses 1e-6 cube', `${tinyCube.reason} (exact volume=${core.qstr(core.volume(core.offsetPoly(core.scalePoly(cube, core.q(1, 1000000)), core.q(1, 10000000))))})`);

// ---- A13: generator determinism -----------------------------------------
const p1 = buildPayload(), p2 = buildPayload();
ok(fnv1a(JSON.stringify(p1)) === fnv1a(JSON.stringify(p2)) && p1.sha256 === p2.sha256, 'A13 generator deterministic', `hash=${p1.sha256} twice`);
const onDisk = readFileSync(join(ROOT, 'data', 'shapes.json'), 'utf8');
ok(onDisk === serialize(p1), 'A13 data file matches generator', `bytes=${onDisk.length} hash=${p1.sha256}`);

// ---- A14: core purity ----------------------------------------------------
const importOut = execFileSync(process.execPath, ['-e', "import('./src/core.mjs').then(m=>console.log(Object.keys(m).length))"], { cwd: ROOT, encoding: 'utf8' }).trim();
ok(Number(importOut) > 20, 'A14 core imports under plain node', `exports=${importOut} no DOM/three`);

// ---- A15: exact arithmetic identities -----------------------------------
const a = core.q(7, 3), b = core.q(-5, 11), c = core.q(13, 17);
ok(core.qeq(core.qadd(a, b), core.q(62, 33)), 'A15 qadd normalization', `${core.qstr(a)}+${core.qstr(b)}=${core.qstr(core.qadd(a, b))}`);
ok(core.qeq(core.qmul(a, core.qadd(b, c)), core.qadd(core.qmul(a, b), core.qmul(a, c))), 'A15 distributivity', `a*(b+c)=${core.qstr(core.qmul(a, core.qadd(b, c)))}`);
ok(core.qeq(core.qsub(a, a), core.q(0)), 'A15 cancellation', `a-a=${core.qstr(core.qsub(a, a))}`);
ok(core.qsign(core.qdiv(core.q(1), core.q(3))) > 0 && core.qeq(core.qdiv(core.q(2, 6), core.q(1, 3)), core.q(1)), 'A15 division', `(2/6)/(1/3)=${core.qstr(core.qdiv(core.q(2, 6), core.q(1, 3)))}`);

// ---- A16: endpoint volumes match closed forms ---------------------------
const d2 = core.q(2);
ok(core.qeq(core.volume(core.offsetPoly(cube, d2)), core.q(216)), 'A16 cube V(2)=216', core.qstr(core.volume(core.offsetPoly(cube, d2))));
ok(core.qeq(core.volume(core.offsetPoly(pyramid, d2)), core.q(343)), 'A16 pyramid V(2)=343', core.qstr(core.volume(core.offsetPoly(pyramid, d2))));

console.log(`\n${pass} passed, ${fail} failed, ${pass + fail} assertions`);
if (fail) { console.log('FAILED: ' + failed.join(' | ')); process.exit(1); }
process.exit(0);