// Flow field of the SIQS model under a switching environment.
//
//   dS/dt = omega (1 - S - I) - v S - beta S I + gamma I
//   dI/dt = I (beta S - gamma)
//
// with S + I + Q = 1, so the physical region is the triangle S, I >= 0, S + I <= 1.
// The disease dies in each fixed environment, but the time-averaged environment has a
// stable endemic state, so switching fast enough keeps it alive.
//
// The canvas draws only flow streaks, the axes and the fixed-point markers; every piece
// of text is an absolutely positioned HTML element, so it never gets rasterised.

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

// average environment: the vector field is linear in the parameters, so averaging the
// fields is the same as averaging beta, omega, v and gamma
const AVG = {
  beta: (ENVS[0].beta + ENVS[1].beta) / 2,
  omega: (ENVS[0].omega + ENVS[1].omega) / 2,
  v: (ENVS[0].v + ENVS[1].v) / 2,
  gamma: (ENVS[0].gamma + ENVS[1].gamma) / 2,
  name: 'Averaged environment',
};

for (const p of [...ENVS, AVG]) {
  p.mu = p.omega + p.v;
  p.Sstar = p.omega / p.mu;                    // disease-free susceptible fraction
  p.R = p.beta * p.Sstar / p.gamma;            // growth of a small outbreak there
  p.Sc = p.gamma / p.beta;                     // herd-immunity threshold
  p.Iend = 1 - p.Sc * (1 + p.v / p.omega);     // endemic level (negative => unreachable)
}

// ---------------------------------------------------------------- view
const xMin = 0.0, xMax = 1.0;
const yMin = -0.020, yMax = 0.50;
const X_TICKS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
const Y_TICKS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
// the right margin has to clear the fixed point sitting exactly at S = 1, and its label
const MARGIN = { left: 84, right: 66, top: 20, bottom: 66 };

const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
const overlay = document.getElementById('overlay');
const envText = document.getElementById('env-text');
const paramText = document.getElementById('param-text');
const verdictText = document.getElementById('verdict-text');
const timeControls = document.getElementById('time-controls');

let width, height;               // CSS pixels
let currentEnv = 0;              // 0, 1 = fixed environments, 2 = switching
let activeEnv = 0;
let lastSwitchTime = Date.now();
// 0 means the averaged environment (infinitely fast switching), which is where the third
// panel starts; pressing "slower" walks it down to finite switching rates
let dwellSec = 0;                // wall-clock seconds spent in each environment
let isAverageMode = false;

// Environment 1 relaxes at rate mu = 10, so a particle can cross the whole plot in a
// couple of frames. Slow the clock a little, and integrate in substeps so that each
// frame draws a short curved STREAK rather than a scatter of disconnected dots.
const TIME_SCALE = 0.7;          // model time units per wall-clock second
const NUM_PARTICLES = 3000;
const SUBSTEPS = 6;
let particles = [];

function toScreen(x, y) {
  const w = width - MARGIN.left - MARGIN.right;
  const h = height - MARGIN.top - MARGIN.bottom;
  return {
    px: MARGIN.left + ((x - xMin) / (xMax - xMin)) * w,
    py: MARGIN.top + h - ((y - yMin) / (yMax - yMin)) * h,
  };
}

function resize() {
  // back the canvas with real device pixels, then work in CSS pixels everywhere else --
  // without this the lines and markers are drawn at 1x and upscaled, which looks blurry
  const dpr = window.devicePixelRatio || 1;
  width = window.innerWidth;
  height = window.innerHeight;
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);
  canvas.style.width = width + 'px';
  canvas.style.height = height + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  initParticles();
  layoutLabels();
}
window.addEventListener('resize', resize);

// ---------------------------------------------------------------- dynamics
const _d = [0, 0];               // scratch, so the hot loop allocates nothing
function deriv(S, I, p, out) {
  out[0] = p.omega * (1 - S - I) - p.v * S - p.beta * S * I + p.gamma * I;
  out[1] = I * (p.beta * S - p.gamma);
}

function activeParams() {
  return isAverageMode ? AVG : ENVS[activeEnv];
}

class Particle {
  constructor() { this.reset(true); }

  reset(randomizeAge = false) {
    this.x = xMin + Math.random() * (xMax - xMin);
    const top = Math.min(1 - this.x, yMax);          // stay inside S + I <= 1
    this.y = top > 0 ? Math.random() * top : 0;
    this.age = randomizeAge ? Math.random() * 200 : 0;
    this.maxAge = 150 + Math.random() * 200;
  }

  outside() {
    return this.x < xMin || this.x > xMax || this.y < 0 || this.y > yMax ||
           this.x + this.y > 1.002;
  }

  // Advance one frame, appending the path travelled to the current canvas path.
  // Every particle contributes to a single batched stroke, so this stays cheap.
  advance(dt, p) {
    const h = dt * TIME_SCALE / SUBSTEPS;
    let pos = toScreen(this.x, this.y);
    ctx.moveTo(pos.px, pos.py);

    for (let k = 0; k < SUBSTEPS; k++) {
      deriv(this.x, this.y, p, _d);
      this.x += _d[0] * h;
      this.y += _d[1] * h;
      if (this.outside()) { this.reset(); return; }   // end the streak at the boundary
      pos = toScreen(this.x, this.y);
      ctx.lineTo(pos.px, pos.py);
    }

    this.age += dt * 60;
    if (this.age > this.maxAge) this.reset();
  }
}

function initParticles() {
  particles = [];
  for (let i = 0; i < NUM_PARTICLES; i++) particles.push(new Particle());
}

// ---------------------------------------------------------------- critical switching rate
// Floquet exponent on the I = 0 boundary, where dS/dt = -mu (S - S*) is linear:
//   lambda(alpha) = (a1 + a2)/2 - (alpha/2) (b1/mu1 - b2/mu2) (Sa - Sb)
//   Sa - Sb       = (S1* - S2*) (1 - E1)(1 - E2) / (1 - E1 E2),   Ei = exp(-mu_i/alpha)
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
  if (floquet(hi) <= 0) return Infinity;         // fast switching cannot rescue
  if (floquet(lo) > 0) return 0;                 // always rescued
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
  yTitle: makeLabel('axis-title title-y', 'I  (infected)'),
  cut: makeLabel('cut-label', 'Q = 0'),
  fp: [makeLabel('fp-label', ''), makeLabel('fp-label', '')],
};

function place(el, px, py) {
  el.style.left = px + 'px';
  el.style.top = py + 'px';
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
  const yc = (o.py + yEnd.py) / 2;
  LABELS.yTitle.style.left = '22px';
  LABELS.yTitle.style.top = yc + 'px';
  LABELS.yTitle.style.transform = 'translate(-50%, -50%) rotate(-90deg)';

  const c = cutLine();
  if (c) place(LABELS.cut, c.b.px + 10, c.b.py + 4);

  labelKey = '';                                 // force the fixed-point labels to move too
}

// the edge of the physical region, I = 1 - S, clipped to the view
function cutLine() {
  const s0 = Math.max(xMin, 1 - yMax);
  const s1 = Math.min(xMax, 1);
  if (s1 <= s0) return null;
  return { a: toScreen(s1, 1 - s1), b: toScreen(s0, 1 - s0) };
}

let labelKey = '';
function updateFixedPointLabels() {
  const p = activeParams();
  const color = isAverageMode ? COLORS.avg : COLORS.env[activeEnv];
  const key = `${p.name}|${color}`;
  if (key === labelKey) return;                  // only touch the DOM when it changes
  labelKey = key;

  const items = [{ x: p.Sstar, y: 0, text: `S* = ${p.Sstar.toFixed(2)}` }];
  if (p.Iend > 1e-4) {
    items.push({ x: p.Sc, y: p.Iend,
                 text: `(${p.Sc.toFixed(2)}, ${p.Iend.toFixed(2)})` });
  }
  LABELS.fp.forEach((el, i) => {
    if (i < items.length) {
      const s = toScreen(items[i].x, items[i].y);
      el.textContent = items[i].text;
      el.style.color = color;
      place(el, s.px, s.py - 18);
      el.style.display = 'block';
    } else {
      el.style.display = 'none';
    }
  });
}

// ---------------------------------------------------------------- canvas drawing
function drawAxes() {
  const o = toScreen(xMin, 0);
  const xEnd = toScreen(xMax, 0);
  const yEnd = toScreen(xMin, yMax);

  ctx.strokeStyle = 'rgba(0, 0, 0, 0.07)';
  ctx.lineWidth = 1;
  ctx.beginPath();
  for (const x of X_TICKS) {
    const a = toScreen(x, 0), b = toScreen(x, yMax);
    ctx.moveTo(a.px, a.py); ctx.lineTo(b.px, b.py);
  }
  for (const y of Y_TICKS) {
    const a = toScreen(xMin, y), b = toScreen(xMax, y);
    ctx.moveTo(a.px, a.py); ctx.lineTo(b.px, b.py);
  }
  ctx.stroke();

  const c = cutLine();
  if (c) {
    ctx.save();
    ctx.setLineDash([7, 6]);
    ctx.strokeStyle = 'rgba(0,0,0,0.42)';
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    ctx.moveTo(c.a.px, c.a.py); ctx.lineTo(c.b.px, c.b.py);
    ctx.stroke();
    ctx.restore();
  }

  ctx.strokeStyle = COLORS.ink;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(o.px, o.py); ctx.lineTo(xEnd.px, o.py);
  ctx.moveTo(o.px, o.py); ctx.lineTo(o.px, yEnd.py);
  for (const x of X_TICKS) {
    const t = toScreen(x, 0);
    ctx.moveTo(t.px, t.py); ctx.lineTo(t.px, t.py + 7);
  }
  for (const y of Y_TICKS) {
    const t = toScreen(xMin, y);
    ctx.moveTo(t.px, t.py); ctx.lineTo(t.px - 7, t.py);
  }
  ctx.stroke();
}

function drawFixedPoint(x, y, color, isSink) {
  const p = toScreen(x, y);
  ctx.save();
  if (isSink) {
    const throb = 1 + 0.15 * Math.sin(Date.now() / 200);
    ctx.globalAlpha = 0.35;
    ctx.beginPath();
    ctx.arc(p.px, p.py, 13 * throb, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.globalAlpha = 1;
    ctx.beginPath();
    ctx.arc(p.px, p.py, 6.5, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
  } else {
    ctx.beginPath();
    ctx.arc(p.px, p.py, 6.5, 0, Math.PI * 2);
    ctx.fillStyle = '#ffffff';
    ctx.fill();
    ctx.lineWidth = 2.5;
    ctx.strokeStyle = color;
    ctx.stroke();
  }
  ctx.restore();
}

function drawFixedPoints() {
  const p = activeParams();
  const color = isAverageMode ? COLORS.avg : COLORS.env[activeEnv];
  // disease-free point: a sink while a small outbreak decays there, a saddle once it grows
  // (at R = 1 exactly it still attracts, just algebraically slowly, so count it as a sink)
  drawFixedPoint(p.Sstar, 0, color, p.R <= 1 + 1e-9);
  if (p.Iend > 1e-4) drawFixedPoint(p.Sc, p.Iend, color, true);
}

// ---------------------------------------------------------------- loop
let lastFrameTime = performance.now();

function animate() {
  const now = performance.now();
  let dt = (now - lastFrameTime) / 1000;
  if (dt > 0.1) dt = 0.1;
  lastFrameTime = now;

  if (currentEnv === 2 && !isAverageMode) {
    if (Date.now() - lastSwitchTime >= dwellSec * 1000) {
      activeEnv = 1 - activeEnv;
      lastSwitchTime = Date.now();
    }
  } else if (currentEnv !== 2) {
    activeEnv = currentEnv;
    lastSwitchTime = Date.now();
  }

  // fade the previous frame towards TRANSPARENT, so the page background shows through
  ctx.globalCompositeOperation = 'destination-out';
  ctx.fillStyle = 'rgba(0, 0, 0, 0.07)';
  ctx.fillRect(0, 0, width, height);
  ctx.globalCompositeOperation = 'source-over';

  drawAxes();

  const p = activeParams();
  ctx.strokeStyle = isAverageMode ? COLORS.avg : COLORS.env[activeEnv];
  ctx.lineWidth = 1.5;
  ctx.lineCap = 'round';
  ctx.lineJoin = 'round';
  ctx.globalAlpha = 0.7;
  ctx.beginPath();
  for (const particle of particles) particle.advance(dt, p);
  ctx.stroke();
  ctx.globalAlpha = 1;

  drawFixedPoints();
  updateFixedPointLabels();
  requestAnimationFrame(animate);
}

// ---------------------------------------------------------------- ui
function fmt(x) { return Number.isInteger(x) ? x.toFixed(0) : x.toString(); }

function updateText() {
  isAverageMode = (currentEnv === 2 && dwellSec === 0);
  const p = isAverageMode ? AVG : (currentEnv === 2 ? null : ENVS[currentEnv]);

  if (currentEnv === 2 && !isAverageMode) {
    const alpha = 1 / (dwellSec * TIME_SCALE);       // dwell in model time units
    envText.innerText = `Switching every ${dwellSec.toFixed(1)} s`;
    envText.style.color = COLORS.ink;
    paramText.innerText =
      `env 1: β=${fmt(ENVS[0].beta)} ω=${fmt(ENVS[0].omega)} v=${fmt(ENVS[0].v)}   ` +
      `env 2: β=${fmt(ENVS[1].beta)} ω=${fmt(ENVS[1].omega)} v=${fmt(ENVS[1].v)}`;
    const alive = alpha > ALPHA_C;
    verdictText.innerText =
      `α = ${alpha.toFixed(2)}   (α_c = ${ALPHA_C.toFixed(2)})   →   ` +
      (alive ? 'disease persists' : 'disease dies out');
    verdictText.style.color = alive ? COLORS.avg : COLORS.ink;
    timeControls.style.display = 'block';
    return;
  }

  envText.innerText = p.name;
  envText.style.color = isAverageMode ? COLORS.avg : COLORS.env[currentEnv];
  paramText.innerText =
    `β = ${fmt(p.beta)},  γ = ${fmt(p.gamma)},  ω = ${fmt(p.omega)},  v = ${fmt(p.v)}` +
    `    →    S* = ${p.Sstar.toFixed(3)}`;
  verdictText.innerText =
    `β S* / γ = ${p.R.toFixed(2)}   →   ` +
    (p.R > 1 ? 'disease persists' : p.R < 1 ? 'disease dies out' : 'exactly at threshold');
  verdictText.style.color = isAverageMode ? COLORS.avg : COLORS.env[currentEnv];
  timeControls.style.display = currentEnv === 2 ? 'block' : 'none';
}

function toggleEnvironment() {
  currentEnv = (currentEnv + 1) % 3;
  if (currentEnv === 2) {
    activeEnv = 0;
    lastSwitchTime = Date.now();
  }
  updateText();
  ctx.clearRect(0, 0, width, height);
}

function changeDwell(delta) {
  if (currentEnv !== 2) return;
  dwellSec = Math.min(10, Math.max(0, dwellSec + delta));
  if (dwellSec < 0.05) dwellSec = 0;
  updateText();
}

document.addEventListener('keydown', (e) => {
  if (e.code === 'Space') { e.preventDefault(); toggleEnvironment(); }
  else if (e.code === 'ArrowUp' || e.key === 'w' || e.key === 'W') { e.preventDefault(); changeDwell(+0.1); }
  else if (e.code === 'ArrowDown' || e.key === 's' || e.key === 'S') { e.preventDefault(); changeDwell(-0.1); }
});

for (const [id, fn] of [['toggle-btn', toggleEnvironment],
                        ['btn-t-plus', () => changeDwell(+0.1)],
                        ['btn-t-minus', () => changeDwell(-0.1)]]) {
  const el = document.getElementById(id);
  if (el) el.addEventListener('click', (e) => { e.preventDefault(); fn(); el.blur(); });
}

resize();
updateText();
animate();
