# Poisson switching: an exact closed form

In the paper the predator–prey model had a *logistic* flow on the extinction boundary, so the
turning points of the cycle could only be found numerically, and the stochastic case needed a
simulated stationary distribution.

Here the boundary flow is **linear**, and that changes everything: the whole problem can be solved
with pen and paper, including the Poisson case. We end up with an explicit *rational* function
$\lambda(\alpha)$ and an explicit critical rate $\alpha_c$ — no transcendental equation at all.

---

## 1. The model and the boundary flow

The SIQS model, with a protected compartment $Q$ that people enter at rate $v$ and leave at rate
$\omega$:

$$\frac{dS}{dt} = \omega Q - vS - \beta S I + \gamma I, \qquad
  \frac{dQ}{dt} = vS - \omega Q, \qquad
  \frac{dI}{dt} = I(\beta S - \gamma).$$

The population is conserved, $S+Q+I=1$, so we may eliminate $Q = 1-S-I$ and work with two variables.

We want to know whether a *small* amount of disease can grow. So we sit on the disease-free boundary
$I = 0$, where $Q = 1-S$ and

$$\frac{dS}{dt} = \omega(1-S) - vS = \omega - (\omega+v)\,S .$$

Two shorthands, and these are the only ones we will need:

$$\boxed{\ \mu = \omega + v\ } \qquad \text{and} \qquad \boxed{\ S^{*} = \frac{\omega}{\mu}\ }$$

so that the boundary flow is just

$$\frac{dS}{dt} = -\mu\,(S - S^{*}).$$

$S^*$ is the susceptible fraction the environment settles at, and $\mu$ is the rate at which it gets
there. Each environment $E\in\{1,2\}$ has its own $\beta_E,\gamma_E,\mu_E,S^*_E$.

### The two facts we need about this flow

**Fact A (where you end up).** The equation $dS/dt = -\mu(S-S^*)$ is exponential relaxation. If a
dwell in environment $E$ starts at $S_{\rm in}$ and lasts a time $\tau$, it ends at

$$S_{\rm out} \;=\; S^{*}_E + \bigl(S_{\rm in} - S^{*}_E\bigr)\,e^{-\mu_E \tau}. \tag{A}$$

*Proof.* Write $u = S - S^*_E$. Then $du/dt = -\mu_E u$, so $u(\tau) = u(0)e^{-\mu_E\tau}$. Add
$S^*_E$ back. $\square$

**Fact B (the area underneath).** Over that same dwell,

$$\int_0^{\tau} S(t)\,dt \;=\; S^{*}_E\,\tau \;-\; \frac{S_{\rm out}-S_{\rm in}}{\mu_E}. \tag{B}$$

*Proof.* Rearrange the boundary equation to make $S$ the subject:
$S = S^*_E - \frac{1}{\mu_E}\frac{dS}{dt}$. Integrate both sides from $0$ to $\tau$. The first term
gives $S^*_E\tau$; the second gives $-\frac{1}{\mu_E}\bigl(S(\tau)-S(0)\bigr)$. $\square$

Fact B is the workhorse. It says: *to know the time-average of $S$ over a stretch, you only need the
endpoints.* It holds for any $\tau$ whatsoever — including a random one. That is why the stochastic
case is no harder than the deterministic one.

### How the disease responds

Since $dI/dt = I(\beta_E S - \gamma_E)$, we have $\dfrac{d\ln I}{dt} = \beta_E S(t) - \gamma_E$.
So the change in $\ln I$ over a dwell is

$$\Delta \ln I \;=\; \int_0^\tau\!\bigl(\beta_E S(t) - \gamma_E\bigr) dt
\;\overset{(B)}{=}\; \bigl(\beta_E S^*_E - \gamma_E\bigr)\tau \;-\; \frac{\beta_E}{\mu_E}\bigl(S_{\rm out}-S_{\rm in}\bigr). \tag{C}$$

The first piece is the boring one: it is the growth rate the disease would have if the environment
never changed, times the time spent. Define

$$a_E \;=\; \beta_E S^{*}_E - \gamma_E .$$

We *assume throughout that $a_1 < 0$ and $a_2 < 0$*: the disease dies in each fixed environment.
The second piece of (C) is the interesting one — it is the bonus (or penalty) the disease collects
because $S$ was still *moving* when the environment changed.

---

## 2. Periodic switching: the cycle in closed form

Let each dwell last exactly $\tau = 1/\alpha$, and label the turning points as in `criticalRate.md`:
environment 1 runs from $S_b$ up or down to $S_a$, then environment 2 runs from $S_a$ back to $S_b$.

By Fact A, writing $E_1 = e^{-\mu_1/\alpha}$ and $E_2 = e^{-\mu_2/\alpha}$,

$$S_a = S^*_1 + (S_b - S^*_1)E_1, \qquad S_b = S^*_2 + (S_a - S^*_2)E_2 .$$

Expand the brackets:

$$S_a = S^*_1(1-E_1) + E_1 S_b, \qquad S_b = S^*_2(1-E_2) + E_2 S_a .$$

Substitute the second into the first:

$$S_a = S^*_1(1-E_1) + E_1\Bigl[S^*_2(1-E_2) + E_2 S_a\Bigr]
\;\Longrightarrow\;
S_a\bigl(1 - E_1E_2\bigr) = S^*_1(1-E_1) + E_1S^*_2(1-E_2),$$

and by the same manipulation the other way round,

$$S_a = \frac{S^*_1(1-E_1) + E_1S^*_2(1-E_2)}{1-E_1E_2},
\qquad
S_b = \frac{S^*_2(1-E_2) + E_2S^*_1(1-E_1)}{1-E_1E_2}. \tag{1}$$

Now — and this is the step that makes everything collapse — we only ever need the **difference**.
Subtracting,

$$S_a - S_b = \frac{S^*_1(1-E_1) - S^*_2(1-E_2) + E_1S^*_2(1-E_2) - E_2S^*_1(1-E_1)}{1-E_1E_2}.$$

Group the $S^*_1$ terms and the $S^*_2$ terms:

$$S_a - S_b = \frac{S^*_1(1-E_1)(1-E_2) - S^*_2(1-E_2)(1-E_1)}{1-E_1E_2},$$

so the factor $(1-E_1)(1-E_2)$ is common and

$$\boxed{\;S_a - S_b = \bigl(S^{*}_1 - S^{*}_2\bigr)\,\frac{(1-E_1)(1-E_2)}{1-E_1E_2}\;} \tag{2}$$

A sanity check: as $\alpha\to0$ (very long dwells) $E_1,E_2\to0$ and $S_a-S_b\to S^*_1-S^*_2$ — the
system reaches each fixed point in turn, as it should. As $\alpha\to\infty$, $E_i\to1$ and the swing
shrinks to zero — the system barely moves.

### The Floquet exponent

Add up (C) over the two dwells, noting that environment 1 takes $S_b\to S_a$ and environment 2 takes
$S_a\to S_b$, and divide by the cycle length $2/\alpha$:

$$\lambda^{\rm det}(\alpha)
= \frac{\alpha}{2}\left[\frac{a_1+a_2}{\alpha} - \frac{\beta_1}{\mu_1}(S_a-S_b) - \frac{\beta_2}{\mu_2}(S_b-S_a)\right]
= \frac{a_1+a_2}{2} - \frac{\alpha}{2}\left(\frac{\beta_1}{\mu_1}-\frac{\beta_2}{\mu_2}\right)(S_a - S_b).$$

Substituting (2) gives the *fully explicit* deterministic answer:

$$\lambda^{\rm det}(\alpha) = \frac{a_1+a_2}{2}
\;-\;\frac{\alpha}{2}\left(\frac{\beta_1}{\mu_1}-\frac{\beta_2}{\mu_2}\right)
\bigl(S^*_1-S^*_2\bigr)\frac{(1-E_1)(1-E_2)}{1-E_1E_2}. \tag{3}$$

Setting $\lambda = 0$ still needs a numerical root because of the exponentials $E_i = e^{-\mu_i/\alpha}$.
That is the last numerical step in the whole problem — and the stochastic case will not even need it.

---

## 3. Poisson switching

Now let the environment change as a Poisson process of rate $\alpha$: each dwell time $\tau$ is an
independent exponential random variable with mean $1/\alpha$, i.e. density $\alpha e^{-\alpha\tau}$.

### 3.1 The one fact about exponentials we need

Fact A involves $e^{-\mu\tau}$. With $\tau$ random we need its average:

$$\mathbb{E}\bigl[e^{-\mu\tau}\bigr]
= \int_0^\infty e^{-\mu\tau}\,\alpha e^{-\alpha\tau}\,d\tau
= \alpha\int_0^\infty e^{-(\alpha+\mu)\tau}d\tau
= \frac{\alpha}{\alpha+\mu}.$$

So define, in exact parallel with $E_1, E_2$:

$$\varepsilon_1 = \frac{\alpha}{\alpha+\mu_1}, \qquad \varepsilon_2 = \frac{\alpha}{\alpha+\mu_2}.$$

### 3.2 Why the averages obey the same equations

Fact A is *linear* in $S_{\rm in}$:

$$S_a = S^*_1 + (S_b - S^*_1)\,e^{-\mu_1\tau_1}.$$

The dwell time $\tau_1$ is drawn fresh when environment 1 begins. The value $S_b$ at which that dwell
starts was produced by *earlier* dwells only. So $\tau_1$ and $S_b$ are **independent**, and the
average of the product is the product of the averages:

$$\mathbb{E}[S_a] = S^*_1 + \bigl(\mathbb{E}[S_b] - S^*_1\bigr)\,\mathbb{E}\bigl[e^{-\mu_1\tau_1}\bigr]
= S^*_1 + \bigl(\mathbb{E}[S_b] - S^*_1\bigr)\varepsilon_1 .$$

This is the *only* place independence is used, and it is why the averages close on themselves — we
never need the full distribution of $S_b$, unlike the logistic case in the paper, where
$\ln \mathcal X$ is not linear in the starting point.

After many cycles the process forgets its initial condition and settles into a stationary state, in
which the average value of $S$ at the start of every environment-1 dwell is the same number. Call the
stationary averages $\langle S_a\rangle$ and $\langle S_b\rangle$. They satisfy

$$\langle S_a\rangle = S^*_1 + (\langle S_b\rangle - S^*_1)\varepsilon_1, \qquad
  \langle S_b\rangle = S^*_2 + (\langle S_a\rangle - S^*_2)\varepsilon_2,$$

which are **exactly** the deterministic equations with $E_i \to \varepsilon_i$. Therefore we may
reuse (2) verbatim:

$$\langle S_a\rangle - \langle S_b\rangle
= \bigl(S^{*}_1-S^{*}_2\bigr)\frac{(1-\varepsilon_1)(1-\varepsilon_2)}{1-\varepsilon_1\varepsilon_2}. \tag{4}$$

### 3.3 The exponentials cancel

Here is where the linear boundary flow pays off completely. Substitute
$\varepsilon_i = \alpha/(\alpha+\mu_i)$ into (4):

$$1 - \varepsilon_i = 1 - \frac{\alpha}{\alpha+\mu_i} = \frac{\mu_i}{\alpha+\mu_i},$$

$$1 - \varepsilon_1\varepsilon_2 = 1 - \frac{\alpha^2}{(\alpha+\mu_1)(\alpha+\mu_2)}
= \frac{(\alpha+\mu_1)(\alpha+\mu_2)-\alpha^2}{(\alpha+\mu_1)(\alpha+\mu_2)}
= \frac{\alpha(\mu_1+\mu_2)+\mu_1\mu_2}{(\alpha+\mu_1)(\alpha+\mu_2)} .$$

Dividing, the awkward factors $(\alpha+\mu_1)(\alpha+\mu_2)$ cancel outright:

$$\frac{(1-\varepsilon_1)(1-\varepsilon_2)}{1-\varepsilon_1\varepsilon_2}
= \frac{\mu_1\mu_2}{\alpha(\mu_1+\mu_2)+\mu_1\mu_2}
= \frac{1}{1 + \alpha\left(\dfrac{1}{\mu_1}+\dfrac{1}{\mu_2}\right)},$$

where the last step is just dividing top and bottom by $\mu_1\mu_2$. So

$$\boxed{\;\langle S_a\rangle - \langle S_b\rangle
= \frac{S^{*}_1 - S^{*}_2}{\,1 + \alpha\left(\frac{1}{\mu_1}+\frac{1}{\mu_2}\right)}\;} \tag{5}$$

No exponentials survive. The average swing of the susceptible pool is a simple saturating function
of the switching rate.

### 3.4 The invasion exponent

Over one full cycle (one dwell $\tau_1$ in environment 1, then one dwell $\tau_2$ in environment 2),
equation (C) gives

$$\Delta\ln I = a_1\tau_1 + a_2\tau_2 - \frac{\beta_1}{\mu_1}\bigl(S_a - S_b\bigr) - \frac{\beta_2}{\mu_2}\bigl(S_b' - S_a\bigr),$$

where $S_b'$ is the value at the end of environment 2, i.e. at the start of the *next* cycle. Take
averages. Two things happen:

* $\mathbb{E}[\tau_1] = \mathbb{E}[\tau_2] = 1/\alpha$;
* in the stationary state $\mathbb{E}[S_b'] = \mathbb{E}[S_b] = \langle S_b\rangle$, because the
  distribution is the same at the start of every cycle. So
  $\mathbb{E}[S_b' - S_a] = -(\langle S_a\rangle - \langle S_b\rangle)$.

Hence

$$\mathbb{E}\bigl[\Delta\ln I\bigr] = \frac{a_1+a_2}{\alpha}
- \left(\frac{\beta_1}{\mu_1}-\frac{\beta_2}{\mu_2}\right)\bigl(\langle S_a\rangle - \langle S_b\rangle\bigr).$$

The growth rate per unit *time* is the average gain per cycle divided by the average duration of a
cycle, $\mathbb{E}[\tau_1+\tau_2] = 2/\alpha$. (This is just the law of large numbers: after $n$
cycles $\ln I$ has changed by about $n\,\mathbb{E}[\Delta\ln I]$ and a time of about
$n\,\mathbb{E}[\tau_1+\tau_2]$ has elapsed.) So

$$\lambda^{\rm stoch}(\alpha) = \frac{a_1+a_2}{2}
- \frac{\alpha}{2}\left(\frac{\beta_1}{\mu_1}-\frac{\beta_2}{\mu_2}\right)\bigl(\langle S_a\rangle-\langle S_b\rangle\bigr),$$

identical in form to the deterministic result — and now, using (5),

$$\boxed{\;
\lambda^{\rm stoch}(\alpha) = \frac{a_1+a_2}{2}
\;+\;\frac{\alpha}{2}\cdot
\frac{\left(\dfrac{\beta_1}{\mu_1}-\dfrac{\beta_2}{\mu_2}\right)\bigl(S^{*}_2 - S^{*}_1\bigr)}
{1+\alpha\left(\dfrac{1}{\mu_1}+\dfrac{1}{\mu_2}\right)}\;}
\tag{6}$$

(the sign flip in the last bracket is just absorbing the minus sign). **This is a rational function
of $\alpha$** — a Michaelis–Menten curve rising from $\frac{a_1+a_2}{2}<0$ at $\alpha=0$ to a plateau
at $\alpha\to\infty$.

### 3.5 The critical rate, in closed form

Because (6) is rational, $\lambda = 0$ is a *linear* equation for $\alpha$. Write the two positive
quantities

$$A = -(a_1+a_2) = \gamma_1+\gamma_2-\beta_1S^*_1-\beta_2S^*_2 \;>\;0
\qquad\text{(how badly the disease loses, on average)}$$

$$G = \left(\frac{\beta_1}{\mu_1}-\frac{\beta_2}{\mu_2}\right)\bigl(S^{*}_2-S^{*}_1\bigr)
\qquad\text{(the gain from switching)}$$

Multiplying (6) by $2\bigl[1+\alpha(\tfrac1{\mu_1}+\tfrac1{\mu_2})\bigr]$ and setting it to zero:

$$-A\left[1+\alpha\left(\frac{1}{\mu_1}+\frac{1}{\mu_2}\right)\right] + \alpha G = 0
\;\Longrightarrow\;
\alpha\left[G - A\left(\frac{1}{\mu_1}+\frac{1}{\mu_2}\right)\right] = A,$$

$$\boxed{\;\alpha_c^{\rm stoch} = \frac{A}{\;G - A\left(\dfrac{1}{\mu_1}+\dfrac{1}{\mu_2}\right)}\;}
\tag{7}$$

Two consistency checks fall out immediately.

* Letting $\alpha\to\infty$ in (6) gives the plateau
  $\lambda_\infty = -\tfrac{A}{2} + \tfrac{G}{2}\bigl(\tfrac1{\mu_1}+\tfrac1{\mu_2}\bigr)^{-1}$.
  This is positive exactly when the denominator of (7) is positive. So $\alpha_c$ is finite and
  positive precisely when infinitely fast switching rescues the disease, and it **diverges** as one
  approaches that boundary — the same divergence the paper had to establish numerically.
* $\lambda^{\rm stoch}$ is monotonically increasing in $\alpha$ (a Michaelis–Menten curve with
  positive $G$), so the root is unique.

### Worked example

$\beta_1=1,\ \gamma_1=1,\ \omega_1=9,\ v_1=1$ and $\beta_2=1.5,\ \gamma_2=1,\ \omega_2=1,\ v_2=1$:

$$\mu_1=10,\ S^*_1=0.9,\ a_1=-0.1;\qquad \mu_2=2,\ S^*_2=0.5,\ a_2=-0.25$$
$$A = 0.35,\qquad G = (0.1-0.75)(0.5-0.9) = 0.26, \qquad \tfrac1{\mu_1}+\tfrac1{\mu_2} = 0.6$$
$$\alpha_c^{\rm stoch} = \frac{0.35}{0.26 - 0.35\times0.6} = \frac{0.35}{0.05} = \mathbf{7} \ \text{exactly.}$$

The periodic threshold, from the numerical root of (3), is $\alpha_c^{\rm det} \approx 2.406$. Random
switching must be about three times faster.

---

## 4. Why random switching is always worse: a one-line proof

Compare (2) and (4): they are the same function

$$\Phi(x,y) = \frac{(1-x)(1-y)}{1-xy}$$

evaluated at $(E_1,E_2)$ for periodic switching and at $(\varepsilon_1,\varepsilon_2)$ for Poisson
switching. Two observations finish the argument.

**(i) $\Phi$ decreases in each argument.** Differentiate the $x$-dependent part:

$$\frac{\partial}{\partial x}\frac{1-x}{1-xy}
= \frac{-(1-xy) - (1-x)(-y)}{(1-xy)^2}
= \frac{-1+xy+y-xy}{(1-xy)^2}
= \frac{y-1}{(1-xy)^2} \;<\;0 \quad\text{for } 0<y<1 .$$

Multiplying by the positive factor $(1-y)$, $\partial\Phi/\partial x<0$; by symmetry
$\partial\Phi/\partial y<0$.

**(ii) Randomness raises both arguments.** The function $\tau \mapsto e^{-\mu\tau}$ is convex, so by
Jensen's inequality the average of the function exceeds the function of the average:

$$\varepsilon_i = \mathbb{E}\bigl[e^{-\mu_i\tau}\bigr] \;>\; e^{-\mu_i\mathbb{E}[\tau]} = e^{-\mu_i/\alpha} = E_i .$$

(For $\mu/\alpha=1$, for instance, $\tfrac12 > e^{-1}=0.368$.)

Putting (i) and (ii) together, $\Phi(\varepsilon_1,\varepsilon_2) < \Phi(E_1,E_2)$: the *average*
swing of the susceptible pool under Poisson switching is strictly smaller than the swing under
periodic switching of the same mean rate. Since the disease's bonus term is proportional to that
swing (with a positive coefficient $G$ whenever switching helps at all),

$$\lambda^{\rm stoch}(\alpha) < \lambda^{\rm det}(\alpha) \quad\text{for every } \alpha
\qquad\Longrightarrow\qquad \alpha_c^{\rm stoch} > \alpha_c^{\rm det}.$$

The intuition is the diminishing-returns argument from the paper, now made exact: a dwell twice as
long does *not* push $S$ twice as far, because $S$ saturates at $S^*$, but the disease keeps paying
$a_E$ per unit time regardless. Spreading the dwell times out therefore costs the disease.

So: **clockwork seasons sustain a disease that erratic ones cannot.**

---

## 5. Remark: you need two things to vary, and they must be anti-correlated

Write $\rho_E = \beta_E/\gamma_E$ and $S^*_E = \omega_E/\mu_E$. Extinction in each environment is
$\rho_E S^*_E<1$; fast-switching rescue is $\bar\rho\,\bar S^*>1$ with
$\bar\rho = \frac{\beta_1+\beta_2}{\gamma_1+\gamma_2}$ and
$\bar S^* = \frac{\omega_1+\omega_2}{\mu_1+\mu_2}$. Both are *mediants*
$\frac{p_1+p_2}{q_1+q_2}$, which always lie strictly between $\frac{p_1}{q_1}$ and $\frac{p_2}{q_2}$;
so if either $\rho$ or $S^*$ is the same in both environments, or if the two are ordered the *same*
way, then $\bar\rho\bar S^* \le \max_E \rho_E S^*_E < 1$ and no switching can help. Rescue requires
the more transmissible environment to be the better-protected one. (This is also why the plain
homogeneous SIR model of `../../SIRSQ_switch/readme.md` can never show the effect: it has only one
ratio, so there is nothing to anti-correlate.)
