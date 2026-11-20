// Float64 baseline for convex face offset, with a fixed absolute epsilon.
// This is the "naive but plausible" pipeline the exact core is compared to:
//   * plane membership and determinant tests use EPS = 1e-9 (absolute)
//   * face validity uses area >= AREA_EPS = 1e-9 (absolute)
// Both constants are absolute, so the pipeline has a hidden scale window.
// It refuses to invent geometry: a rejected face is reported, never patched.

export const EPS = 1e-9;
export const AREA_EPS = 1e-9;

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];

function solve3(rows, rhs) {
  const M = rows.map((r, i) => [r[0], r[1], r[2], rhs[i]]);
  for (let col = 0; col < 3; col++) {
    let piv = col;
    for (let r = col + 1; r < 3; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
    if (Math.abs(M[piv][col]) < EPS) return null;
    if (piv !== col) { const t = M[col]; M[col] = M[piv]; M[piv] = t; }
    const pv = M[col][col];
    for (let c = col; c < 4; c++) M[col][c] /= pv;
    for (let r = 0; r < 3; r++) {
      if (r === col) continue;
      const f = M[r][col];
      if (Math.abs(f) < 1e-300) continue;
      for (let c = col; c < 4; c++) M[r][c] -= f * M[col][c];
    }
  }
  return [M[0][3], M[1][3], M[2][3]];
}

export function toFloatPoly(poly) {
  return {
    name: poly.name,
    verts: poly.verts.map(v => v.map(c => Number(c.n) / Number(c.d))),
    faces: poly.faces.map(f => f.slice()),
    planes: poly.planes.map(p => ({
      n: p.n.map(c => Number(c.n) / Number(c.d)),
      c: Number(p.c.n) / Number(p.c.d),
    })),
  };
}

export function faceStats(poly) {
  const areas = [];
  for (const face of poly.faces) {
    let a = 0;
    const v0 = poly.verts[face[0]];
    for (let i = 1; i + 1 < face.length; i++) {
      const cr = cross(sub(poly.verts[face[i]], v0), sub(poly.verts[face[i + 1]], v0));
      a += 0.5 * Math.sqrt(dot(cr, cr));
    }
    areas.push(a);
  }
  return { areas, degenerate: areas.filter(a => a < AREA_EPS).length, total: areas.length };
}

export function offsetNaive(poly, d) {
  const verts = [];
  let degenerateVerts = 0;
  for (let vi = 0; vi < poly.verts.length; vi++) {
    const v = poly.verts[vi];
    const idx = [];
    poly.planes.forEach((p, i) => {
      if (Math.abs(dot(p.n, v) - p.c) <= EPS) idx.push(i);
    });
    let sol = null;
    for (let a = 0; a < idx.length && !sol; a++) {
      for (let b = a + 1; b < idx.length && !sol; b++) {
        for (let c = b + 1; c < idx.length && !sol; c++) {
          const P = [poly.planes[idx[a]], poly.planes[idx[b]], poly.planes[idx[c]]];
          sol = solve3(P.map(p => p.n), P.map(p => p.c + d));
        }
      }
    }
    if (!sol) { degenerateVerts++; sol = [0, 0, 0]; }
    verts.push(sol);
  }
  return { name: poly.name, verts, faces: poly.faces.map(f => f.slice()) };
}

export function volumeNaive(poly) {
  let s = 0;
  for (const face of poly.faces) {
    const v0 = poly.verts[face[0]];
    for (let i = 1; i + 1 < face.length; i++) {
      s += dot(v0, cross(poly.verts[face[i]], poly.verts[face[i + 1]]));
    }
  }
  return s / 6;
}

// Fit V(d) = a0 + a1 d + a2 d^2 + a3 d^3 through four samples, in Float64.
export function cubicNaive(poly, nodes = [0, 0.1, 0.2, 0.3]) {
  const vals = nodes.map(d => volumeNaive(offsetNaive(poly, d)));
  const a = new Array(4).fill(0);
  nodes.forEach((d, i) => {
    let p = [1];
    for (let j = 0; j < 4; j++) {
      if (j === i) continue;
      const q = [-nodes[j], 1];
      const out = new Array(p.length + 1).fill(0);
      for (let u = 0; u < p.length; u++) for (let w = 0; w < 2; w++) out[u + w] += p[u] * q[w];
      p = out;
    }
    let den = 1;
    for (let j = 0; j < 4; j++) if (j !== i) den *= nodes[i] - nodes[j];
    for (let k = 0; k < 4; k++) a[k] += vals[i] * p[k] / den;
  });
  return a;
}

// Full naive analysis: rejects bodies whose faces fall below AREA_EPS.
export function analyzeNaive(floatPoly, d) {
  const st = faceStats(floatPoly);
  if (st.degenerate > 0) {
    return {
      usable: false, reason: `degenerate-faces ${st.degenerate}/${st.total} (area < ${AREA_EPS})`,
      degenerateFaces: st.degenerate, totalFaces: st.total,
      minArea: Math.min(...st.areas), volume: null, coeffs: null,
    };
  }
  const off = offsetNaive(floatPoly, d);
  return {
    usable: true, reason: 'ok', degenerateFaces: 0, totalFaces: st.total,
    minArea: Math.min(...st.areas), volume: volumeNaive(off), coeffs: null,
  };
}