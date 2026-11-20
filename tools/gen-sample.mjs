// Deterministic sample generator. Writes data/shapes.json only.
// Same command twice -> byte-identical output (no timestamps, no randomness
// beyond the fixed-seed LCG).
import { writeFileSync, mkdirSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import * as core from '../src/core.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const ROOT = join(HERE, '..');

const OFFSETS = ['0.1', '0.25', '0.5', '0.75', '1', '2'];
const SCALES = ['1/1000000', '1/1000', '1', '1000', '1000000'];

export function fnv1a(str) {
  let h = 0x811c9dc5;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 0x01000193) >>> 0;
  }
  return h.toString(16).padStart(8, '0');
}

export function buildPayload() {
  const shapes = core.shapes();
  const out = { format: 'offsetcert-shapes-v1', offsets: OFFSETS, scales: SCALES, shapes: [] };
  for (const p of shapes) {
    const law = core.cubicLaw(p);
    const offsets = OFFSETS.map(ds => {
      const d = core.q(ds);
      const off = core.offsetPoly(p, d);
      const wt = core.watertight(off);
      const cv = core.convexCertificate(off);
      return {
        d: ds,
        volume: core.qstr(core.volume(off)),
        watertight: wt.ok,
        convex: cv.ok,
        edges: wt.edges,
      };
    });
    out.shapes.push({
      name: p.name,
      verts: p.verts.map(v => v.map(core.qstr)),
      faces: p.faces,
      planes: p.planes.map(pl => ({ n: pl.n.map(core.qstr), c: core.qstr(pl.c) })),
      volume0: core.qstr(core.volume(p)),
      cubic: law.coeffs.map(core.qstr),
      offsets,
    });
  }
  out.sha256 = fnv1a(JSON.stringify(out));
  return out;
}

export function serialize(payload) {
  return JSON.stringify(payload, null, 2) + '\n';
}

if (import.meta.url === `file://${process.argv[1]}`) {
  mkdirSync(join(ROOT, 'data'), { recursive: true });
  const payload = buildPayload();
  writeFileSync(join(ROOT, 'data', 'shapes.json'), serialize(payload));
  console.log(`wrote data/shapes.json  ${serialize(payload).length} bytes  hash=${payload.sha256}`);
}