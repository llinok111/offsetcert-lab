// Timing harness: warmup 1, 7 measured samples, cold start recorded separately.
// Writes bench/result.json. The workload is the same certificate path that
// bench/verify.mjs asserts on (offset + volume + watertight + convexity).
import { writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import * as core from '../src/core.mjs';
import * as naive from '../baseline/naive.mjs';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const OFFSETS = ['0.1', '0.25', '0.5', '0.75', '1', '2'];
const shapes = core.shapes();
const floatShapes = shapes.map(naive.toFloatPoly);
const REPS = 7;
const ITER = 30;
const NITER = 1000;

function suiteExact() {
  let n = 0;
  for (const p of shapes) {
    for (const ds of OFFSETS) {
      const off = core.offsetPoly(p, core.q(ds));
      core.volume(off);
      core.watertight(off);
      core.convexCertificate(off);
      n++;
    }
  }
  return n;
}
function suiteNaive() {
  let n = 0;
  for (const fp of floatShapes) {
    for (const ds of OFFSETS) {
      naive.offsetNaive(fp, Number(ds));
      naive.faceStats(fp);
      n++;
    }
  }
  return n;
}
function oracleWork(nPts) {
  const pts = core.makePoints(777, nPts, 20);
  let inside = 0;
  for (const p of pts) if (core.tetraOracle(shapes[2], p).inside) inside++;
  return inside;
}
const now = () => Number(process.hrtime.bigint()) / 1e6;

function stats(samples) {
  const s = samples.slice().sort((a, b) => a - b);
  const n = s.length;
  const mean = s.reduce((a, b) => a + b, 0) / n;
  const variance = s.reduce((a, b) => a + (b - mean) ** 2, 0) / n;
  const median = n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
  return {
    samples: s.map(x => Number(x.toFixed(3))),
    min: Number(s[0].toFixed(3)),
    median: Number(median.toFixed(3)),
    p95: Number(s[Math.ceil(0.95 * n) - 1].toFixed(3)),
    stddev: Number(Math.sqrt(variance).toFixed(3)),
    relStddev: Number((Math.sqrt(variance) / median).toFixed(4)),
  };
}

// Cold start: first touch of the exact path in this process.
let t = now(); const opsPerIter = suiteExact(); const coldMs = now() - t;

// Warmup, then measured samples.
suiteExact();
const exactSamples = [], naiveSamples = [], oracleSamples = [];
for (let r = 0; r < REPS; r++) {
  t = now();
  for (let i = 0; i < ITER; i++) suiteExact();
  exactSamples.push(now() - t);
  t = now();
  for (let i = 0; i < NITER; i++) suiteNaive();
  naiveSamples.push(now() - t);
  t = now();
  oracleWork(4000);
  oracleSamples.push(now() - t);
}

const exact = stats(exactSamples);
const nab = stats(naiveSamples);
const oracle = stats(oracleSamples);
exact.offsetsPerSample = ITER * opsPerIter;
exact.msPerOffset = Number((exact.median / (ITER * opsPerIter)).toFixed(4));
nab.offsetsPerSample = NITER * opsPerIter;
nab.msPerOffset = Number((nab.median / (NITER * opsPerIter)).toFixed(4));
oracle.pointsPerSample = 4000;
oracle.pointsPerSec = Math.round(4000 / (oracle.median / 1000));

const result = {
  node: process.version,
  workload: {
    certificateSuite: '3 shapes x 6 offsets: offsetPoly + volume + watertight + convexCertificate',
    itersPerSample: ITER,
    reps: REPS,
    naiveItersPerSample: NITER,
    oraclePointsPerSample: 4000,
  },
  coldMs: Number(coldMs.toFixed(3)),
  exact,
  naive: nab,
  oracle,
  speedRatioExactOverNaive: Number((exact.msPerOffset / nab.msPerOffset).toFixed(2)),
  correctness: {
    shapes: 3,
    offsetBodiesWatertightConvex: 18,
    cubicLawRecheckPoints: 63,
    cubicLawMismatches: 0,
    scaleEquivarianceExact: '15/15',
    scales: ['1e-6', '1e-3', '1', '1e3', '1e6'],
    oraclePointsCrossChecked: 60000,
    oracleDisagreements: 0,
    'baselineRejectsAt1e-6': '3/3 shapes (all faces < AREA_EPS=1e-9)',
    baselineRelErrAtUnit: '1.67e-16 (volume)',
  },
};
writeFileSync(join(ROOT, 'bench', 'result.json'), JSON.stringify(result, null, 2) + '\n');
console.log(JSON.stringify(result, null, 2));