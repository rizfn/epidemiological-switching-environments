// Trajectories of the SIQS model under a periodically switching environment.
//
//   dS/dt = omega (1 - S - I) - v S - beta S I + gamma I
//   dI/dt = I (beta S - gamma)
//
// The environment alternates every `dwellSec` seconds. Click anywhere to release a cloud
// of initial conditions; each one is integrated forever and drags a comet trail behind it,
// coloured by whichever environment it is travelling through. The trails fade away, the
// particles themselves never do.
//
// Turning the switching all the way up lands on dwellSec = 0, which is the time-averaged
// environment: the alpha -> infinity limit, drawn in purple.
//
// Two stacked canvases: `trails` fades towards transparent every frame, `heads` is cleared
// every frame for the axes, the fixed points and the moving dots. All text is HTML.

// ---------------------------------------------------------------- parameters
const GAMMA = 1.0;
const ENVS = [
  { beta: 1.0, omega: 10.0, v: 0.0, gamma: GAMMA, name: 'Environment 1' },
  { beta: 1.5, omega: 1.0, v: 1.0, gamma: GAMMA, name: 'Environment 2' },
];

const COLORS = {
  env: ['#1e3a8a', '#991b1b'],     // deep blue, deep red
  avg: '#6b21a8',                  // deep purple
  ink: '#000000',
};

const AVG = {
  beta: (ENVS[0].beta + ENVS[1].beta) / 2,
  omega: (ENVS[0].omega + ENVS[1].omega) / 2,
  v: (ENVS[0].v + ENVS[1].v) / 2,
  gamma: (ENVS[0].gamma + ENVS[1].gamma) / 2,
};

for (const p of [...ENVS, AVG]) {
  p.mu = p.omega + p.v;
  p.Sstar = p.omega / p.mu;                    // disease-free susceptible fraction
  p.R = p.beta * p.Sstar / p.gamma;
  p.Sc = p.gamma / p.beta;                     // herd-immunity threshold
  p.Iend = 1 - p.Sc * (1 + p.v / p.omega);     // endemic level
}

// ---------------------------------------------------------------- view
const xMin = 0.0, xMax = 1.0;
const yMin = -0.020, yMax = 0.50;
const X_TICKS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
const Y_TICKS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
const MARGIN = { left: 84, right: 66, top: 20, bottom: 66 };

const TIME_SCALE = 6.0;          // model time units per wall-clock second
const SUBSTEPS = 10;
const TRAIL_FADE = 0.06;         // fraction of the trail canvas erased each frame
const BATCH = 15;                // trajectories released per click
const SPREAD = 0.010;            // radius of the released cloud, in phase-space units
const MAX_TRAJ = 900;

// the dwell time is stepped multiplicatively, so a handful of presses spans the whole
// useful range; dropping below DWELL_MIN snaps to 0, the averaged environment
const DWELL_FACTOR = 1.15, DWELL_MIN = 0.03, DWELL_MAX = 1.0;

const trailCanvas = document.getElementById('trails');
const headCanvas = document.getElementById('heads');
const tctx = trailCanvas.getContext('2d');
const hctx = headCanvas.getContext('2d');
const overlay = document.getElementById('overlay');
const envText = document.getElementById('env-text');
const paramText = document.getElementById('param-text');
const verdictText = document.getElementById('verdict-text');

let width, height;               // CSS pixels
let activeEnv = 0;
let lastSwitchTime = Date.now();
let dwellSec = 0.12;             // wall-clock seconds spent in each environment
let isAverageMode = false;
let trajectories = [];

function toScreen(x, y) {
  const w = width - MARGIN.left - MARGIN.right;
  const h = height - MARGIN.top - MARGIN.bottom;
  return {
    px: MARGIN.left + ((x - xMin) / (xMax - xMin)) * w,
    py: MARGIN.top + h - ((y - yMin) / (yMax - yMin)) * h,
  };
}

function fromScreen(px, py) {
  const w = width - MARGIN.left - MARGIN.right;
  const h = height - MARGIN.top - MARGIN.bottom;
  return {
    x: xMin + ((px - MARGIN.left) / w) * (xMax - xMin),
    y: yMin + ((MARGIN.top + h - py) / h) * (yMax - yMin),
  };
}

function sizeCanvas(cv, c) {
  const dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(width * dpr);
  cv.height = Math.round(height * dpr);
  cv.style.width = width + 'px';
  cv.style.height = height + 'px';
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function resize() {
  width = window.innerWidth;
  height = window.innerHeight;
  sizeCanvas(trailCanvas, tctx);
  sizeCanvas(headCanvas, hctx);
  layoutLabels();
}
window.addEventListener('resize', resize);

// ---------------------------------------------------------------- dynamics
const k1 = [0, 0], k2 = [0, 0], k3 = [0, 0], k4 = [0, 0];

function deriv(S, I, p, out) {
  out[0] = p.omega * (1 - S - I) - p.v * S - p.beta * S * I + p.gamma * I;
  out[1] = I * (p.beta * S - p.gamma);
}

function rk4(tr, h, p) {
  const S = tr.x, I = tr.y;
  deriv(S, I, p, k1);
  deriv(S + 0.5 * h * k1[0], I + 0.5 * h * k1[1], p, k2);
  deriv(S + 0.5 * h * k2[0], I + 0.5 * h * k2[1], p, k3);
  deriv(S + h * k3[0], I + h * k3[1], p, k4);
  tr.x = S + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
  tr.y = I + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
  if (!(tr.y > 1e-12)) tr.y = 1e-12;            // also catches NaN
  if (!isFinite(tr.x)) tr.dead = true;
}

// ---------------------------------------------------------------- seeding
function inRegion(x, y) {
  return x >= 0 && x <= 1 && y >= 0 && x + y <= 1;
}

function seedBatch(cx, cy) {
  if (!inRegion(cx, cy)) return;
  for (let i = 0; i < BATCH; i++) {
    let x, y, tries = 0;
    do {
      const r = SPREAD * Math.sqrt(Math.random());
      const th = Math.random() * 2 * Math.PI;
      x = cx + r * Math.cos(th);
      y = cy + r * Math.sin(th);
    } while (!inRegion(x, y) && ++tries < 12);
    if (!inRegion(x, y)) { x = cx; y = cy; }
    trajectories.push({ x, y: Math.max(y, 1e-6), dead: false });
  }
  if (trajectories.length > MAX_TRAJ) {
    trajectories.splice(0, trajectories.length - MAX_TRAJ);
  }
}

function seedRandom() {
  // somewhere in the accessible region that is actually on screen, and high enough up
  // that the trajectory has something to show
  for (let attempt = 0; attempt < 50; attempt++) {
    const x = Math.random();
    const hi = Math.min(1 - x, yMax);
    if (hi < 0.03) continue;
    seedBatch(x, 0.02 + Math.random() * (hi - 0.02));
    return;
  }
}

function reset() {
  trajectories = [];
  tctx.clearRect(0, 0, width, height);
}

// ---------------------------------------------------------------- critical switching rate
function floquet(alpha) {
  const [p1, p2] = ENVS;
  const E1 = Math.exp(-p1.mu / alpha);
  const E2 = Math.exp(-p2.mu / alpha);
  const swing = (p1.Sstar - p2.Sstar) * (1 - E1) * (1 - E2) / (1 - E1 * E2);
  const a1 = p1.beta * p1.Sstar - p1.gamma;
  const a2 = p2.beta * p2.Sstar - p2.gamma;
  return (a1 + a2) / 2 - (alpha / 2) * (p1.beta / p1.mu - p2.beta / p2.mu) * swing;
}

function criticalAlpha() {
  let lo = 1e-3, hi = 1e3;
  if (floquet(hi) <= 0) return Infinity;
  if (floquet(lo) > 0) return 0;
  for (let i = 0; i < 80; i++) {
    const mid = Math.sqrt(lo * hi);
    if (floquet(mid) > 0) hi = mid; else lo = mid;
  }
  return Math.sqrt(lo * hi);
}
const ALPHA_C = criticalAlpha();

// ---------------------------------------------------------------- html labels
function makeLabel(cls, text) {
  const el = document.createElement('div');
  el.className = cls;
  el.textContent = text;
  overlay.appendChild(el);
  return el;
}

const LABELS = {
  xTicks: X_TICKS.map((v) => makeLabel('tick tick-x', v.toFixed(1))),
  yTicks: Y_TICKS.map((v) => makeLabel('tick tick-y', v.toFixed(1))),
  xTitle: makeLabel('axis-title title-x', 'S  (susceptible)'),
  yTitle: makeLabel('axis-title', 'I  (infected)'),
  cut: makeLabel('cut-label', 'Q = 0'),
};

function place(el, px, py) {
  el.style.left = px + 'px';
  el.style.top = py + 'px';
}

function cutLine() {
  const s0 = Math.max(xMin, 1 - yMax);
  const s1 = Math.min(xMax, 1);
  if (s1 <= s0) return null;
  return { a: toScreen(s1, 1 - s1), b: toScreen(s0, 1 - s0) };
}

function layoutLabels() {
  X_TICKS.forEach((v, i) => {
    const t = toScreen(v, 0);
    place(LABELS.xTicks[i], t.px, t.py + 12);
  });
  Y_TICKS.forEach((v, i) => {
    const t = toScreen(xMin, v);
    place(LABELS.yTicks[i], t.px - 12, t.py);
  });

  const o = toScreen(xMin, 0), xEnd = toScreen(xMax, 0), yEnd = toScreen(xMin, yMax);
  place(LABELS.xTitle, (o.px + xEnd.px) / 2, height - 12);
  LABELS.yTitle.style.left = '22px';
  LABELS.yTitle.style.top = ((o.py + yEnd.py) / 2) + 'px';
  LABELS.yTitle.style.transform = 'translate(-50%, -50%) rotate(-90deg)';

  const c = cutLine();
  if (c) place(LABELS.cut, c.b.px + 10, c.b.py + 4);
  // the fixed points are marked but deliberately left unlabelled
}

// ---------------------------------------------------------------- drawing
function plotRect() {
  return [MARGIN.left, MARGIN.top,
          width - MARGIN.left - MARGIN.right,
          height - MARGIN.top - MARGIN.bottom];
}

function drawAxes(c) {
  c.strokeStyle = 'rgba(0, 0, 0, 0.07)';
  c.lineWidth = 1;
  c.beginPath();
  for (const x of X_TICKS) {
    const a = toScreen(x, 0), b = toScreen(x, yMax);
    c.moveTo(a.px, a.py); c.lineTo(b.px, b.py);
  }
  for (const y of Y_TICKS) {
    const a = toScreen(xMin, y), b = toScreen(xMax, y);
    c.moveTo(a.px, a.py); c.lineTo(b.px, b.py);
  }
  c.stroke();

  const cut = cutLine();
  if (cut) {
    c.save();
    c.setLineDash([7, 6]);
    c.strokeStyle = 'rgba(0,0,0,0.42)';
    c.lineWidth = 1.6;
    c.beginPath();
    c.moveTo(cut.a.px, cut.a.py); c.lineTo(cut.b.px, cut.b.py);
    c.stroke();
    c.restore();
  }

  const o = toScreen(xMin, 0), xEnd = toScreen(xMax, 0), yEnd = toScreen(xMin, yMax);
  c.strokeStyle = COLORS.ink;
  c.lineWidth = 2;
  c.beginPath();
  c.moveTo(o.px, o.py); c.lineTo(xEnd.px, o.py);
  c.moveTo(o.px, o.py); c.lineTo(o.px, yEnd.py);
  for (const x of X_TICKS) {
    const t = toScreen(x, 0);
    c.moveTo(t.px, t.py); c.lineTo(t.px, t.py + 7);
  }
  for (const y of Y_TICKS) {
    const t = toScreen(xMin, y);
    c.moveTo(t.px, t.py); c.lineTo(t.px - 7, t.py);
  }
  c.stroke();
}

function marker(c, x, y, color, style) {
  const p = toScreen(x, y);
  c.save();
  if (style === 'active') {
    const throb = 1 + 0.15 * Math.sin(Date.now() / 200);
    c.globalAlpha = 0.3;
    c.beginPath(); c.arc(p.px, p.py, 13 * throb, 0, Math.PI * 2);
    c.fillStyle = color; c.fill();
    c.globalAlpha = 1;
    c.beginPath(); c.arc(p.px, p.py, 6.5, 0, Math.PI * 2);
    c.fillStyle = color; c.fill();
  } else {
    c.beginPath(); c.arc(p.px, p.py, 5.5, 0, Math.PI * 2);
    c.fillStyle = '#ffffff'; c.fill();
    c.lineWidth = 2.5; c.strokeStyle = color; c.stroke();
  }
  c.restore();
}

// ---------------------------------------------------------------- loop
let lastFrameTime = performance.now();

function animate() {
  const now = performance.now();
  let dt = (now - lastFrameTime) / 1000;
  if (dt > 0.1) dt = 0.1;
  lastFrameTime = now;

  if (!isAverageMode && Date.now() - lastSwitchTime >= dwellSec * 1000) {
    activeEnv = 1 - activeEnv;
    lastSwitchTime = Date.now();
  }

  const p = isAverageMode ? AVG : ENVS[activeEnv];
  const color = isAverageMode ? COLORS.avg : COLORS.env[activeEnv];
  const [rx, ry, rw, rh] = plotRect();

  // --- trails: fade what is already there towards TRANSPARENT, then add this frame's
  // segments. The page background shows through, and the particles are untouched.
  tctx.globalCompositeOperation = 'destination-out';
  tctx.fillStyle = `rgba(0, 0, 0, ${TRAIL_FADE})`;
  tctx.fillRect(0, 0, width, height);
  tctx.globalCompositeOperation = 'source-over';

  if (trajectories.length) {
    tctx.save();
    tctx.beginPath(); tctx.rect(rx, ry, rw, rh); tctx.clip();
    tctx.strokeStyle = color;
    tctx.globalAlpha = 0.85;
    tctx.lineWidth = 1.5;
    tctx.lineCap = 'round';
    tctx.lineJoin = 'round';
    tctx.beginPath();

    const h = dt * TIME_SCALE / SUBSTEPS;
    for (const tr of trajectories) {
      if (tr.dead) continue;
      let pos = toScreen(tr.x, tr.y);
      tctx.moveTo(pos.px, pos.py);
      for (let k = 0; k < SUBSTEPS; k++) {
        rk4(tr, h, p);
        if (tr.dead) break;
        pos = toScreen(tr.x, tr.y);
        tctx.lineTo(pos.px, pos.py);
      }
    }
    tctx.stroke();
    tctx.restore();
  }

  // --- axes, fixed points and the particles themselves, redrawn crisp every frame
  hctx.clearRect(0, 0, width, height);
  drawAxes(hctx);
  // the two disease-free points are plain static rings: which environment is live is
  // already obvious from the colour of the particles
  marker(hctx, ENVS[0].Sstar, 0, COLORS.env[0], 'idle');
  marker(hctx, ENVS[1].Sstar, 0, COLORS.env[1], 'idle');
  marker(hctx, AVG.Sc, AVG.Iend, COLORS.avg, isAverageMode ? 'active' : 'idle');

  hctx.save();
  hctx.beginPath(); hctx.rect(rx, ry, rw, rh); hctx.clip();
  hctx.fillStyle = color;
  for (const tr of trajectories) {
    if (tr.dead) continue;
    const pos = toScreen(tr.x, tr.y);
    hctx.beginPath();
    hctx.arc(pos.px, pos.py, 2.6, 0, Math.PI * 2);
    hctx.fill();
  }
  hctx.restore();

  requestAnimationFrame(animate);
}

// ---------------------------------------------------------------- ui
function fmt(x) { return Number.isInteger(x) ? x.toFixed(0) : x.toString(); }

// The header deliberately does NOT name the instantaneous environment: at these dwell
// times it would flicker several times a second. The colour of the particles and of the
// throbbing fixed point carries that information instead.
function updateText() {
  isAverageMode = (dwellSec === 0);

  if (isAverageMode) {
    envText.innerText = 'Averaged environment';
    envText.style.color = COLORS.avg;
    paramText.innerText =
      `β = ${fmt(AVG.beta)},  γ = ${fmt(AVG.gamma)},  ω = ${fmt(AVG.omega)}, ` +
      ` v = ${fmt(AVG.v)}    →    S* = ${AVG.Sstar.toFixed(3)}`;
    verdictText.innerText =
      `α → ∞   (α_c = ${ALPHA_C.toFixed(2)})   →   disease persists`;
    verdictText.style.color = COLORS.avg;
    return;
  }

  const alpha = 1 / (dwellSec * TIME_SCALE);
  const alive = alpha > ALPHA_C;
  envText.innerText = `Switching every ${dwellSec.toFixed(2)} s`;
  envText.style.color = COLORS.ink;
  paramText.innerText =
    `env 1: β=${fmt(ENVS[0].beta)} ω=${fmt(ENVS[0].omega)} v=${fmt(ENVS[0].v)}   ` +
    `env 2: β=${fmt(ENVS[1].beta)} ω=${fmt(ENVS[1].omega)} v=${fmt(ENVS[1].v)}   ` +
    `γ=${fmt(GAMMA)}`;
  verdictText.innerText =
    `α = ${alpha.toFixed(2)}   (α_c = ${ALPHA_C.toFixed(2)})   →   ` +
    (alive ? 'disease persists' : 'disease dies out');
  verdictText.style.color = alive ? COLORS.avg : COLORS.ink;
}

// dir = +1 slower, -1 faster. Faster than DWELL_MIN is the averaged environment.
function changeDwell(dir) {
  if (dwellSec === 0) {
    if (dir > 0) dwellSec = DWELL_MIN;           // step back out of the averaged case
  } else {
    const d = dwellSec * (dir > 0 ? DWELL_FACTOR : 1 / DWELL_FACTOR);
    dwellSec = d < DWELL_MIN * 0.999 ? 0 : Math.min(DWELL_MAX, d);
    if (dwellSec) lastSwitchTime = Date.now();
  }
  updateText();
}

headCanvas.addEventListener('mousedown', (e) => {
  const r = headCanvas.getBoundingClientRect();
  const q = fromScreen(e.clientX - r.left, e.clientY - r.top);
  seedBatch(q.x, q.y);
});

document.addEventListener('keydown', (e) => {
  if (e.code === 'Space') { e.preventDefault(); seedRandom(); }
  else if (e.key === 'r' || e.key === 'R') { e.preventDefault(); reset(); }
  else if (e.code === 'ArrowUp' || e.key === 'w' || e.key === 'W') { e.preventDefault(); changeDwell(+1); }
  else if (e.code === 'ArrowDown' || e.key === 's' || e.key === 'S') { e.preventDefault(); changeDwell(-1); }
});

for (const [id, fn] of [['btn-seed', seedRandom],
                        ['btn-reset', reset],
                        ['btn-t-plus', () => changeDwell(+1)],
                        ['btn-t-minus', () => changeDwell(-1)]]) {
  const el = document.getElementById(id);
  if (el) el.addEventListener('click', (e) => { e.preventDefault(); fn(); el.blur(); });
}

resize();
updateText();
animate();
