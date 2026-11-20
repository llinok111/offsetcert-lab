// Exact rational face-offset certificates for convex polyhedra.
// Pure ES module: no three.js, no DOM, no I/O, no environment access.
//
// A "face offset" moves every supporting plane n.x = c outward by d:
//   n.x = c + d   (n is a rational unit normal)
// Vertices are re-solved exactly (BigInt rationals), so the offset body is a
// rational polyhedron at every d and every model scale.

const B0 = 0n, B1 = 1n;

function bgcd(a, b) {
  a = a < B0 ? -a : a; b = b < B0 ? -b : b;
  while (b) { const t = a % b; a = b; b = t; }
  return a;
}
function blcm(a, b) { return (a / bgcd(a, b)) * b; }

// Accepts integers, BigInts, decimal strings ("0.1") and fractions ("1/1000").
function qInt(x) {
  if (typeof x === 'string' && x.includes('.')) {
    const neg = x.startsWith('-');
    const [i, f = ''] = (neg ? x.slice(1) : x).split('.');
    const num = BigInt((i || '0') + f);
    return { n: neg ? -num : num, d: 10n ** BigInt(f.length) };
  }
  return { n: BigInt(x), d: 1n };
}
export function q(n, d = B1) {
  if (typeof n === 'string' && n.includes('/')) {
    const [a, b] = n.split('/');
    return q(qInt(a).n * qInt(b).d, qInt(a).d * qInt(b).n);
  }
  if (typeof n === 'string' || typeof n === 'number') {
    const p = qInt(String(n));
    n = p.n; d = d === B1 ? p.d : BigInt(d);
  } else {
    n = BigInt(n); d = BigInt(d);
  }
  if (d === B0) throw new Error('zero denominator');
  if (d < B0) { n = -n; d = -d; }
  const g = bgcd(n, d) || B1;
  return { n: n / g, d: d / g };
}
export const Q0 = q(0);
export const Q1 = q(1);
export const qneg = a => ({ n: -a.n, d: a.d });
export const qadd = (a, b) => q(a.n * b.d + b.n * a.d, a.d * b.d);
export const qsub = (a, b) => q(a.n * b.d - b.n * a.d, a.d * b.d);
export const qmul = (a, b) => q(a.n * b.n, a.d * b.d);
export const qdiv = (a, b) => {
  if (b.n === B0) throw new Error('division by zero');
  return q(a.n * b.d, a.d * b.n);
};
export function qcmp(a, b) { const v = a.n * b.d - b.n * a.d; return v < B0 ? -1 : v > B0 ? 1 : 0; }
export const qsign = a => (a.n < B0 ? -1 : a.n > B0 ? 1 : 0);
export const qeq = (a, b) => a.n === b.n && a.d === b.d;
export const qnum = a => Number(a.n) / Number(a.d);
export const qstr = a => (a.d === B1 ? String(a.n) : `${a.n}/${a.d}`);

function isqrt(n) {
  if (n < B0) throw new Error('isqrt negative');
  if (n < 2n) return n;
  let x = n, y = (x + B1) / 2n;
  while (y < x) { x = y; y = (x + n / x) / 2n; }
  return x;
}

// ---- exact vectors -------------------------------------------------------
const vadd = (a, b) => [qadd(a[0], b[0]), qadd(a[1], b[1]), qadd(a[2], b[2])];
const vsub = (a, b) => [qsub(a[0], b[0]), qsub(a[1], b[1]), qsub(a[2], b[2])];
const vscale = (a, s) => [qmul(a[0], s), qmul(a[1], s), qmul(a[2], s)];
export const vdot = (a, b) => qadd(qadd(qmul(a[0], b[0]), qmul(a[1], b[1])), qmul(a[2], b[2]));
export const vcross = (a, b) => [
  qsub(qmul(a[1], b[2]), qmul(a[2], b[1])),
  qsub(qmul(a[2], b[0]), qmul(a[0], b[2])),
  qsub(qmul(a[0], b[1]), qmul(a[1], b[0])),
];
const vneg = a => [qneg(a[0]), qneg(a[1]), qneg(a[2])];

// Turn a rational cross product into a rational UNIT normal.
// Throws when the length is irrational: this library certifies only shapes
// whose supporting normals are rational unit vectors (true for all samples).
function unitNormal(cr) {
  const L = cr.reduce((a, c) => blcm(a, c.d), B1);
  let ints = cr.map(c => c.n * (L / c.d));
  const g = ints.reduce((a, b) => bgcd(a, b), B0) || B1;
  ints = ints.map(x => x / g);
  const len2 = ints.reduce((a, b) => a + b * b, B0);
  const len = isqrt(len2);
  if (len * len !== len2) throw new Error(`irrational normal length ${len2}`);
  return ints.map(x => q(x, len));
}

function faceNormal(verts, face) {
  const v0 = verts[face[0]];
  for (let i = 1; i + 1 < face.length; i++) {
    const a = vsub(verts[face[i]], v0);
    const b = vsub(verts[face[i + 1]], v0);
    const cr = vcross(a, b);
    if (cr.some(c => qsign(c) !== 0)) return unitNormal(cr);
  }
  throw new Error('degenerate face');
}

function derivePlanes(p) {
  return p.faces.map(face => {
    const n = faceNormal(p.verts, face);
    const c = vdot(n, p.verts[face[0]]);
    return { n, c };
  });
}

function interiorPoint(p) {
  let s = [Q0, Q0, Q0];
  for (const v of p.verts) s = vadd(s, v);
  return vscale(s, q(1, p.verts.length));
}

// Reverse any face whose normal points inward (plane test against the vertex
// centroid). Face loops are re-oriented one by one, not by a global volume
// sign, because a mixed-orientation mesh can still have positive volume.
function orientOutward(p) {
  const c = interiorPoint(p);
  p.planes = derivePlanes(p);
  let reversed = 0;
  p.faces = p.faces.map((face, i) => {
    if (qsign(qsub(vdot(p.planes[i].n, c), p.planes[i].c)) > 0) {
      reversed++;
      return face.slice().reverse();
    }
    return face;
  });
  p.planes = derivePlanes(p);
  return reversed;
}

export function buildPoly(name, verts, faces) {
  const p = {
    name,
    verts: verts.map(v => v.map(c => (Array.isArray(c) ? q(c[0], c[1]) : q(c)))),
    faces: faces.map(f => f.slice()),
  };
  orientOutward(p);
  if (qsign(volume(p)) < 0) throw new Error('inconsistent orientation: ' + name);
  return p;
}

// Signed volume, exact: 1/6 sum (v0 x v1) . v2 over a fan triangulation.
export function volume(p) {
  let s = Q0;
  for (const face of p.faces) {
    const v0 = p.verts[face[0]];
    for (let i = 1; i + 1 < face.length; i++) {
      const cr = vcross(p.verts[face[i]], p.verts[face[i + 1]]);
      s = qadd(s, vdot(v0, cr));
    }
  }
  return qdiv(s, q(6));
}

function solve3(rows, rhs) {
  const M = rows.map((r, i) => [r[0], r[1], r[2], rhs[i]]);
  for (let col = 0; col < 3; col++) {
    let piv = -1;
    for (let r = col; r < 3; r++) if (qsign(M[r][col]) !== 0) { piv = r; break; }
    if (piv < 0) return null;
    if (piv !== col) { const t = M[col]; M[col] = M[piv]; M[piv] = t; }
    const pv = M[col][col];
    for (let c = col; c < 4; c++) M[col][c] = qdiv(M[col][c], pv);
    for (let r = 0; r < 3; r++) {
      if (r === col) continue;
      const f = M[r][col];
      if (qsign(f) === 0) continue;
      for (let c = col; c < 4; c++) M[r][c] = qsub(M[r][c], qmul(f, M[col][c]));
    }
  }
  return [M[0][3], M[1][3], M[2][3]];
}

// Exact face offset: every supporting plane moves out by d, every vertex is
// re-solved from three of its incident planes. Topology is unchanged.
export function offsetPoly(poly, d) {
  const out = [];
  for (let vi = 0; vi < poly.verts.length; vi++) {
    const v = poly.verts[vi];
    const idx = [];
    poly.planes.forEach((pl, i) => {
      if (qsign(qsub(vdot(pl.n, v), pl.c)) === 0) idx.push(i);
    });
    let sol = null;
    for (let a = 0; a < idx.length && !sol; a++) {
      for (let b = a + 1; b < idx.length && !sol; b++) {
        for (let c = b + 1; c < idx.length && !sol; c++) {
          const P = [poly.planes[idx[a]], poly.planes[idx[b]], poly.planes[idx[c]]];
          sol = solve3(P.map(pl => pl.n), P.map(pl => qadd(pl.c, d)));
        }
      }
    }
    if (!sol) throw new Error(`vertex ${vi} of ${poly.name} has no 3 independent planes`);
    out.push(sol);
  }
  const res = { name: poly.name, verts: out, faces: poly.faces.map(f => f.slice()) };
  orientOutward(res);
  if (qsign(volume(res)) < 0) throw new Error('offset orientation broke: ' + poly.name);
  return res;
}

// Directed-edge audit: every undirected edge appears exactly twice, in
// opposite directions (consistent outward orientation, closed manifold).
export function watertight(poly) {
  const m = new Map();
  for (const face of poly.faces) {
    if (face.length < 3) return { ok: false, reason: 'short-face', edges: 0 };
    for (let i = 0; i < face.length; i++) {
      const a = face[i], b = face[(i + 1) % face.length];
      if (a === b) return { ok: false, reason: 'repeated-vertex', edges: 0 };
      const key = a < b ? `${a}_${b}` : `${b}_${a}`;
      const rec = m.get(key) || { count: 0, sum: 0 };
      rec.count++; rec.sum += a < b ? 1 : -1;
      m.set(key, rec);
    }
  }
  let boundary = 0, nonManifold = 0, misoriented = 0;
  for (const rec of m.values()) {
    if (rec.count === 1) boundary++;
    else if (rec.count > 2) nonManifold++;
    else if (rec.sum !== 0) misoriented++;
  }
  return {
    ok: boundary === 0 && nonManifold === 0 && misoriented === 0 && qsign(volume(poly)) > 0,
    edges: m.size, boundary, nonManifold, misoriented,
  };
}

// Convexity + closedness certificate, entirely exact.
// Each vertex must satisfy n.v <= c for every plane, with >= 3 planes tight.
// Every face must be planar and every off-face vertex strictly inside it.
export function convexCertificate(poly) {
  let violations = 0, maxViolation = Q0, minTight = Infinity, nonPlanar = 0, nonStrict = 0;
  for (const v of poly.verts) {
    let tight = 0;
    for (const pl of poly.planes) {
      const s = qsign(qsub(vdot(pl.n, v), pl.c));
      if (s > 0) {
        violations++;
        const over = qsub(vdot(pl.n, v), pl.c);
        if (qcmp(over, maxViolation) > 0) maxViolation = over;
      } else if (s === 0) tight++;
    }
    if (tight < 3) violations++;
    if (tight < minTight) minTight = tight;
  }
  for (let fi = 0; fi < poly.faces.length; fi++) {
    const face = poly.faces[fi], pl = poly.planes[fi];
    const inFace = new Set(face);
    for (let vi = 0; vi < poly.verts.length; vi++) {
      const s = qsign(qsub(vdot(pl.n, poly.verts[vi]), pl.c));
      if (inFace.has(vi)) { if (s !== 0) nonPlanar++; }
      else if (s >= 0) nonStrict++;
    }
  }
  return {
    ok: violations === 0 && nonPlanar === 0 && nonStrict === 0,
    violations, maxViolation: qnum(maxViolation), minTight,
    nonPlanar, nonStrict,
  };
}

export function containsHalfspace(poly, p) {
  let tight = 0;
  for (const pl of poly.planes) {
    const s = qsign(qsub(vdot(pl.n, p), pl.c));
    if (s > 0) return { inside: false, boundary: false };
    if (s === 0) tight++;
  }
  return { inside: true, boundary: tight > 0 };
}

// Independent containment oracle: uses only vertices, no plane coefficients.
// The centroid is interior; a point is outside iff it lies on the far side of
// some face triangle from the centroid.
export function tetraOracle(poly, p) {
  let centroid = [Q0, Q0, Q0];
  for (const v of poly.verts) centroid = vadd(centroid, v);
  centroid = vscale(centroid, q(1, poly.verts.length));
  let boundary = false;
  for (const face of poly.faces) {
    for (let i = 1; i + 1 < face.length; i++) {
      const a = poly.verts[face[0]], b = poly.verts[face[i]], c = poly.verts[face[i + 1]];
      const e1 = vsub(b, a), e2 = vsub(c, a);
      const n = vcross(e1, e2);
      const dC = qsign(vdot(n, vsub(centroid, a)));
      const dP = qsign(vdot(n, vsub(p, a)));
      if (dP === 0) { boundary = true; continue; }
      if (dC !== dP) return { inside: false, boundary: false };
    }
  }
  return { inside: !boundary, boundary };
}

// ---- polynomial helpers (exact) -----------------------------------------
function polyMul(a, b) {
  const out = new Array(a.length + b.length - 1).fill(Q0);
  for (let i = 0; i < a.length; i++) {
    if (qsign(a[i]) === 0) continue;
    for (let j = 0; j < b.length; j++) out[i + j] = qadd(out[i + j], qmul(a[i], b[j]));
  }
  return out;
}
function polyAdd(a, b) {
  const n = Math.max(a.length, b.length);
  const out = new Array(n).fill(Q0);
  for (let i = 0; i < n; i++) out[i] = qadd(a[i] || Q0, b[i] || Q0);
  return out;
}
const polyScale = (a, s) => a.map(c => qmul(c, s));
export function evalPoly(coeffs, x) {
  let acc = Q0;
  for (let i = coeffs.length - 1; i >= 0; i--) acc = qadd(qmul(acc, x), coeffs[i]);
  return acc;
}

function lagrange(nodes, vals) {
  let acc = [Q0, Q0, Q0, Q0];
  for (let i = 0; i < nodes.length; i++) {
    let basis = [Q1];
    let den = Q1;
    for (let j = 0; j < nodes.length; j++) {
      if (i === j) continue;
      basis = polyMul(basis, [qneg(nodes[j]), Q1]);
      den = qmul(den, qsub(nodes[i], nodes[j]));
    }
    acc = polyAdd(acc, polyScale(basis, qdiv(vals[i], den)));
  }
  return acc.slice(0, 4);
}

const CUBIC_NODES = [q(0), q(1, 10), q(2, 10), q(3, 10)];

// Exact cubic law V(d) = a0 + a1 d + a2 d^2 + a3 d^3, interpolated at four
// exact rational points and re-checked at `checks` independent points.
export function cubicLaw(poly, checks = []) {
  const vals = CUBIC_NODES.map(d => volume(offsetPoly(poly, d)));
  const coeffs = lagrange(CUBIC_NODES, vals);
  let maxDiff = Q0, mismatches = 0;
  const points = checks.length ? checks : defaultChecks();
  for (const d of points) {
    const got = volume(offsetPoly(poly, d));
    const want = evalPoly(coeffs, d);
    const diff = qsub(got, want);
    if (qsign(diff) !== 0) mismatches++;
    if (qcmp(qabs(diff), maxDiff) > 0) maxDiff = qabs(diff);
  }
  const d1 = [coeffs[1], qmul(q(2), coeffs[2]), qmul(q(3), coeffs[3])];
  const d2 = [qmul(q(2), coeffs[2]), qmul(q(6), coeffs[3])];
  const crit = qsign(d2[1]) === 0 ? [] : [qdiv(qneg(d2[0]), d2[1])];
  const probes = [Q0, q(2), ...crit.filter(x => qcmp(x, Q0) >= 0 && qcmp(x, q(2)) <= 0)];
  let minSlope = null;
  for (const x of probes) {
    const s = evalPoly(d1, x);
    if (minSlope === null || qcmp(s, minSlope) < 0) minSlope = s;
  }
  return { coeffs, mismatches, maxDiff, maxDiffNum: qnum(maxDiff), points: points.length, minSlope, monotone: qsign(minSlope) > 0 };
}

function qabs(a) { return a.n < B0 ? qneg(a) : a; }

export function defaultChecks() {
  const out = [];
  for (let k = 1; k <= 20; k++) out.push(q(2 * k - 1, 20)); // 0.05 .. 1.95
  out.push(q(2));
  return out;
}

export function scalePoly(poly, s) {
  const res = {
    name: poly.name,
    verts: poly.verts.map(v => vscale(v, s)),
    faces: poly.faces.map(f => f.slice()),
  };
  res.planes = poly.planes.map(pl => ({ n: pl.n.slice(), c: qmul(pl.c, s) }));
  return res;
}

export function polyEqual(a, b) {
  if (a.verts.length !== b.verts.length) return false;
  for (let i = 0; i < a.verts.length; i++) {
    for (let k = 0; k < 3; k++) if (!qeq(a.verts[i][k], b.verts[i][k])) return false;
  }
  return true;
}

// ---- deterministic sample shapes ---------------------------------------
// All unit-scale, rational vertices and rational unit normals.
// Multiply by 10 to read them as millimetres.
export function shapes() {
  const cube = buildPoly('cube',
    [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
     [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]],
    [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
  );
  const pyramid = buildPoly('pyramid',
    [[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0], [0, 0, [3, 4]]],
    [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
  );
  const prism = buildPoly('chamfer-prism', [
    [1, [3, 4], -1], [[2, 3], 1, -1], [-1, 1, -1], [-1, -1, -1], [1, -1, -1],
    [1, [3, 4], 1], [[2, 3], 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1],
  ], [
    [0, 1, 2, 3, 4], [5, 6, 7, 8, 9],
    [0, 5, 6, 1], [1, 6, 7, 2], [2, 7, 8, 3], [3, 8, 9, 4], [4, 9, 5, 0],
  ]);
  return [cube, pyramid, prism];
}

// Deterministic 32-bit LCG; lattice points with denominator 10.
export function makePoints(seed, count, span = 30) {
  let s = BigInt(seed) & 0xffffffffn;
  const pts = [];
  const next = () => {
    s = (s * 1664525n + 1013904223n) & 0xffffffffn;
    return Number(s);
  };
  for (let i = 0; i < count; i++) {
    const k = [next() % (2 * span + 1), next() % (2 * span + 1), next() % (2 * span + 1)];
    pts.push(k.map(x => q(BigInt(x - span), 10)));
  }
  return pts;
}

export function boundaryPoints(poly) {
  const pts = [];
  for (const v of poly.verts) pts.push(v.slice());
  for (const face of poly.faces) {
    for (let i = 0; i < face.length; i++) {
      const a = poly.verts[face[i]], b = poly.verts[face[(i + 1) % face.length]];
      pts.push(vscale(vadd(a, b), q(1, 2)));
    }
  }
  return pts;
}

export const API = { q, qadd, qsub, qmul, qdiv, qcmp, qsign, qeq, qnum, qstr };