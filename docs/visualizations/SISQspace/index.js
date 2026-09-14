// SIQS in two neighbourhoods coupled by commuting, drawn as individual people.
//
//   dS_i/dt = w_i Q_i - v_i S_i - b_i S_i I_i + g I_i + D (S_j - S_i)
//   dQ_i/dt = v_i S_i - w_i Q_i                        <- no D: the protected stay put
//   dI_i/dt = b_i S_i I_i - g I_i       + D (I_j - I_i)
//
// The equations are integrated with RK4 and are the ground truth. The dots are a faithful
// rendering of that solution, not an independent simulation: each frame the six densities
// are rounded to whole people (largest-remainder, so the head count is exactly conserved)
// and the dots are reassigned to match. Blue = susceptible, red = infected, purple =
// protected. A dot therefore crosses the gap only when the equations call for a NET
// migration -- there is no cosmetic two-way traffic.
//
// Environment 1 has v = 0: nobody there is ever protected, so it holds no purple at all.
// The recovery rate is the same in both, so infinitely fast commuting is guaranteed to
// kill the disease (../../../src/SIQS_space/math/criticalDiffusion.md). The disease dies
// at D = 0 too. It survives only in the window 0 < D < D_c.

// ---------------------------------------------------------------- model
const GAMMA = 1.0;
const ENVS = [
  { beta: 0.9, omega: 1.0, v: 0.0, name: 'Environment 1' },   // nobody protected
  { beta: 2.5, omega: 1.0, v: 2.0, name: 'Environment 2' },   // protected, but infectious
];

const COLORS = {
  S: '#486ccf',                  // deep blue
  I: '#d23f3f',                  // deep red
  Q: '#7f1e91',                  // deep purple
  ink: '#000000',
};

for (const p of ENVS) {
  p.mu = p.omega + p.v;
  p.Sstar = p.omega / p.mu;                     // susceptible fraction when isolated
}
const SBAR = 2 / (1 / ENVS[0].Sstar + 1 / ENVS[1].Sstar);   // harmonic mean = 0.5
const A = ENVS.map((p) => p.beta * SBAR - GAMMA);           // = (-0.55, +0.25)
const D_C = A[0] * A[1] / (A[0] + A[1]);                    // 1/D_c = 1/a_1 + 1/a_2

const S = 0, I = 1, Q = 2;
const DOT_R = 5;                 // radius of one person, in pixels
const N0 = 300;                  // dots per unit of density
const NTOT = 2 * N0;             // total head count, conserved
const TIME_SCALE = 4.0;          // model time units per wall-clock second
const SUBSTEPS = 8;
const SEED = 0.06;               // density of infection added by a re-seed

const D_FACTOR = 1.4, D_MIN = 0.02, D_MAX = 20;
let D = 0.10;

// state vector: [S1, S2, Q1, Q2, I1, I2]
let y = new Float64Array(6);

function initialState() {
  const z = new Float64Array(6);
  for (let i = 0; i < 2; i++) {
    z[i] = ENVS[i].Sstar - SEED;               // take the seed out of the susceptibles,
    z[2 + i] = 1 - ENVS[i].Sstar;              // so each neighbourhood still holds 1 unit
    z[4 + i] = SEED;
  }
  return z;
}

function deriv(z, out) {
  for (let i = 0; i < 2; i++) {
    const j = 1 - i, p = ENVS[i];
    const Si = z[i], Qi = z[2 + i], Ii = z[4 + i];
    out[i] = p.omega * Qi - p.v * Si - p.beta * Si * Ii + GAMMA * Ii + D * (z[j] - Si);
    out[2 + i] = p.v * Si - p.omega * Qi;
    out[4 + i] = p.beta * Si * Ii - GAMMA * Ii + D * (z[4 + j] - Ii);
  }
}

const _a = new Float64Array(6), _b = new Float64Array(6);
const _c = new Float64Array(6), _dd = new Float64Array(6), _t = new Float64Array(6);

function rk4(h) {
  deriv(y, _a);
  for (let i = 0; i < 6; i++) _t[i] = y[i] + 0.5 * h * _a[i];
  deriv(_t, _b);
  for (let i = 0; i < 6; i++) _t[i] = y[i] + 0.5 * h * _b[i];
  deriv(_t, _c);
  for (let i = 0; i < 6; i++) _t[i] = y[i] + h * _c[i];
  deriv(_t, _dd);
  for (let i = 0; i < 6; i++) {
    y[i] += (h / 6) * (_a[i] + 2 * _b[i] + 2 * _c[i] + _dd[i]);
    if (y[i] < 0) y[i] = 0;
  }
}

// Cutting the link has to send people home as well. Commuting piles people into
// environment 2 (it ends up with 1.5 units against environment 1's 0.5), and an
// overcrowded environment 2 sustains the disease entirely on its own -- indeed
// b_2 * Sbar > gamma is precisely the condition that makes the rescue work in the first
// place. So "D = 0" has to mean two genuinely separate towns, each with its own unit of
// population, which is the baseline the whole analysis is scaled to.
function sendEveryoneHome() {
  for (let i = 0; i < 2; i++) {
    const n = y[i] + y[2 + i] + y[4 + i];
    if (n > 1e-12) { y[i] /= n; y[2 + i] /= n; y[4 + i] /= n; }
  }
}

// ---------------------------------------------------------------- layout
const PAD = 38, BOT = 104;
const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
const overlay = document.getElementById('overlay');
const headline = document.getElementById('headline');
const verdictEl = document.getElementById('verdict');
const legendEl = document.getElementById('legend');

let width, height, R, innerR, rimR, qInnerR;
const centres = [{ x: 0, y: 0 }, { x: 0, y: 0 }];
const FAR_ANGLE = [Math.PI, 0];        // direction away from the other neighbourhood
const Q_ARC = 1.25;                    // half-width of the protected crescent, radians

function panelBottom() {
  let b = 0;
  for (const id of ['ui', 'controls']) {
    const el = document.getElementById(id);
    if (el) b = Math.max(b, el.getBoundingClientRect().bottom);
  }
  return b || 160;
}

function layout() {
  const top = Math.min(panelBottom() + 20, height * 0.45);
  R = Math.max(70, Math.min((width - 4 * PAD) / 4, (height - top - BOT) / 2));
  innerR = R - 14;
  rimR = R - 8;
  qInnerR = rimR - 0.30 * R;
  const cy = top + (height - top - BOT) / 2;
  centres[0] = { x: width / 2 - (R + PAD), y: cy };
  centres[1] = { x: width / 2 + (R + PAD), y: cy };
}

function resize() {
  const dpr = window.devicePixelRatio || 1;
  width = window.innerWidth;
  height = window.innerHeight;
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);
  canvas.style.width = width + 'px';
  canvas.style.height = height + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  layout();
  layoutLabels();
  for (const ag of people) { ag.qth = undefined; retarget(ag); }
}
window.addEventListener('resize', resize);

// ---------------------------------------------------------------- the dots
let people = [];
const counts = [0, 0, 0, 0, 0, 0];               // bucket = patch * 3 + state

function retarget(ag) {
  const c = centres[ag.patch];
  if (ag.state === Q) {
    // the protected cannot commute: cluster them tight against the far edge of their own
    // circle, on the side away from the other neighbourhood
    if (ag.qth === undefined) {
      ag.qth = FAR_ANGLE[ag.patch] + (Math.random() - 0.5) * 2 * Q_ARC;
      ag.qr = rimR - (rimR - qInnerR) * Math.pow(Math.random(), 1.7);
    }
    const th = ag.qth + (Math.random() - 0.5) * 0.05;
    const r = ag.qr + (Math.random() - 0.5) * 4;
    ag.tx = c.x + r * Math.cos(th);
    ag.ty = c.y + r * Math.sin(th);
  } else {
    ag.qth = undefined;
    const r = innerR * Math.sqrt(Math.random());
    const th = Math.random() * 2 * Math.PI;
    ag.tx = c.x + r * Math.cos(th);
    ag.ty = c.y + r * Math.sin(th);
  }
  ag.wait = 0.5 + Math.random() * 1.5;
}

function makePeople() {
  people = [];
  for (let i = 0; i < NTOT; i++) {
    const patch = i < NTOT / 2 ? 0 : 1;
    const c = centres[patch];
    const r = innerR * Math.sqrt(Math.random());
    const th = Math.random() * 2 * Math.PI;
    const ag = { patch, state: S, x: c.x + r * Math.cos(th), y: c.y + r * Math.sin(th),
                 tx: 0, ty: 0, wait: 0, qth: undefined, qr: 0 };
    retarget(ag);
    people.push(ag);
  }
  recount();
}

function recount() {
  counts.fill(0);
  for (const ag of people) counts[ag.patch * 3 + ag.state]++;
}

// how many dots each bucket should hold: largest-remainder rounding of the densities,
// normalised so the six buckets always add up to exactly NTOT people
function targetCounts() {
  const dens = [y[0], y[4], y[2], y[1], y[5], y[3]];   // -> S1 I1 Q1 S2 I2 Q2
  let total = 0;
  for (const d of dens) total += d;
  const raw = dens.map((d) => (total > 0 ? d * NTOT / total : 0));
  const out = raw.map(Math.floor);
  const left = NTOT - out.reduce((p, q) => p + q, 0);
  const order = raw.map((r, i) => [r - Math.floor(r), i]).sort((p, q) => q[0] - p[0]);
  for (let k = 0; k < left; k++) out[order[k % 6][1]]++;
  return out;
}

// reassign dots until the buckets match. Changing compartment in place is preferred to
// crossing the gap, so a dot only travels when the equations really do move a person.
function reconcile() {
  const want = targetCounts();
  const diff = counts.map((c, i) => c - want[i]);
  const cost = (from, to) => {
    const samePatch = ((from / 3) | 0) === ((to / 3) | 0);
    return samePatch ? 0 : (from % 3 === to % 3 ? 1 : 2);
  };

  for (let level = 0; level <= 2; level++) {
    for (let from = 0; from < 6; from++) {
      while (diff[from] > 0) {
        let to = -1;
        for (let t = 0; t < 6; t++) if (diff[t] < 0 && cost(from, t) === level) { to = t; break; }
        if (to < 0) break;
        const ag = people.find((a) => a.patch * 3 + a.state === from);
        if (!ag) { diff[from] = 0; break; }
        ag.patch = (to / 3) | 0;
        ag.state = to % 3;
        ag.qth = undefined;                    // pick a fresh spot in its new home
        retarget(ag);
        counts[from]--; counts[to]++;
        diff[from]--; diff[to]++;
      }
    }
  }
}

const SEEK = 2.6;                      // 1/seconds: how fast a dot chases its target
function movePeople(dt) {
  const k = Math.min(1, SEEK * dt);
  for (const ag of people) {
    ag.wait -= dt;
    if (ag.wait <= 0) retarget(ag);
    ag.x += (ag.tx - ag.x) * k;
    ag.y += (ag.ty - ag.y) * k;
  }
}

// ---------------------------------------------------------------- labels
function makeLabel(cls) {
  const el = document.createElement('div');
  el.className = cls;
  overlay.appendChild(el);
  return el;
}
const LABELS = [0, 1].map(() => ({
  title: makeLabel('circle-title'),
  params: makeLabel('circle-params'),
  counts: makeLabel('circle-counts'),
}));

function layoutLabels() {
  for (let i = 0; i < 2; i++) {
    const c = centres[i], L = LABELS[i], p = ENVS[i];
    L.title.textContent = p.name;
    L.title.style.left = c.x + 'px';
    L.title.style.top = (c.y - R - 24) + 'px';
    L.params.textContent =
      `β=${p.beta}  γ=${GAMMA}  ω=${p.omega}  v=${p.v}    S* = ` +
      (p.Sstar === 1 ? '1' : p.Sstar.toFixed(2));
    L.params.style.left = c.x + 'px';
    L.params.style.top = (c.y + R + 14) + 'px';
    L.counts.style.left = c.x + 'px';
    L.counts.style.top = (c.y + R + 42) + 'px';
  }
}

let countTimer = 0;
function updateCounts() {
  for (let i = 0; i < 2; i++) {
    const b = i * 3;
    LABELS[i].counts.innerHTML =
      `<span style="color:${COLORS.S}">${counts[b + S]}</span> / ` +
      `<span style="color:${COLORS.I}">${counts[b + I]}</span> / ` +
      `<span style="color:${COLORS.Q}">${counts[b + Q]}</span>`;
  }
}

// ---------------------------------------------------------------- drawing
function draw() {
  ctx.clearRect(0, 0, width, height);

  ctx.strokeStyle = 'rgba(0,0,0,0.55)';
  ctx.lineWidth = 2;
  for (const c of centres) {
    ctx.beginPath();
    ctx.arc(c.x, c.y, R, 0, Math.PI * 2);
    ctx.stroke();
  }

  if (D > 0) {                          // the commuting channel, thickening with D
    ctx.save();
    ctx.setLineDash([6, 7]);
    ctx.strokeStyle = 'rgba(0,0,0,0.35)';
    ctx.lineWidth = Math.min(6, 1 + 2.2 * Math.log10(1 + D / 0.05));
    ctx.beginPath();
    ctx.moveTo(centres[0].x + R, centres[0].y);
    ctx.lineTo(centres[1].x - R, centres[1].y);
    ctx.stroke();
    ctx.restore();
  }

  for (const [state, color] of [[S, COLORS.S], [Q, COLORS.Q], [I, COLORS.I]]) {
    ctx.fillStyle = color;
    ctx.beginPath();
    for (const ag of people) {
      if (ag.state !== state) continue;
      ctx.moveTo(ag.x + DOT_R, ag.y);
      ctx.arc(ag.x, ag.y, DOT_R, 0, Math.PI * 2);
    }
    ctx.fill();
  }
}

// ---------------------------------------------------------------- loop
let lastFrameTime = performance.now();

function animate() {
  const now = performance.now();
  let dt = (now - lastFrameTime) / 1000;
  if (dt > 0.1) dt = 0.1;
  lastFrameTime = now;

  const h = dt * TIME_SCALE / SUBSTEPS;
  for (let k = 0; k < SUBSTEPS; k++) rk4(h);

  reconcile();
  movePeople(dt);
  draw();

  countTimer -= dt;
  if (countTimer <= 0) { updateCounts(); countTimer = 0.15; }

  requestAnimationFrame(animate);
}

// ---------------------------------------------------------------- ui
function updateText() {
  headline.innerText = D === 0
    ? 'No commuting  (D = 0)'
    : `Commuting rate D = ${D < 1 ? D.toFixed(3) : D.toFixed(2)}`;

  let msg, color;
  if (D === 0) {
    msg = 'neighbourhoods isolated → the disease dies in both';
    color = COLORS.ink;
  } else if (D < D_C) {
    msg = `D < D_c = ${D_C.toFixed(2)} → the disease persists`;
    color = COLORS.I;
  } else {
    msg = `D > D_c = ${D_C.toFixed(2)} → too much mixing, the disease dies`;
    color = COLORS.ink;
  }
  verdictEl.innerText = msg;
  verdictEl.style.color = color;

  legendEl.innerHTML =
    `<b style="color:${COLORS.S}">■</b> susceptible &nbsp; ` +
    `<b style="color:${COLORS.I}">■</b> infected &nbsp; ` +
    `<b style="color:${COLORS.Q}">■</b> protected`;
}

function seedInfection() {
  for (let i = 0; i < 2; i++) {
    const take = Math.min(SEED, y[i]);
    y[i] -= take;
    y[4 + i] += take;
  }
}

function reset() {
  y = initialState();
  makePeople();
  reconcile();
}

function changeD(dir) {
  if (D === 0) {
    if (dir > 0) D = D_MIN;
  } else {
    const d = D * (dir > 0 ? D_FACTOR : 1 / D_FACTOR);
    if (d < D_MIN * 0.999) { D = 0; sendEveryoneHome(); }
    else D = Math.min(D_MAX, d);
  }
  updateText();
}

document.addEventListener('keydown', (e) => {
  if (e.code === 'Space') { e.preventDefault(); seedInfection(); }
  else if (e.key === 'r' || e.key === 'R') { e.preventDefault(); reset(); }
  else if (e.code === 'ArrowUp' || e.key === 'w' || e.key === 'W') { e.preventDefault(); changeD(+1); }
  else if (e.code === 'ArrowDown' || e.key === 's' || e.key === 'S') { e.preventDefault(); changeD(-1); }
});

for (const [id, fn] of [['btn-d-plus', () => changeD(+1)],
                        ['btn-d-minus', () => changeD(-1)],
                        ['btn-seed', seedInfection],
                        ['btn-reset', reset]]) {
  const el = document.getElementById(id);
  if (el) el.addEventListener('click', (e) => { e.preventDefault(); fn(); el.blur(); });
}

width = window.innerWidth;
height = window.innerHeight;
updateText();                    // size the panels before the circles are laid out
layout();
y = initialState();
makePeople();
resize();
reconcile();
updateCounts();
animate();
