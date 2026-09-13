# Commuting between two neighbourhoods: a window of persistence

Instead of one population whose environment changes in *time*, take two neighbourhoods that differ
permanently and let people commute between them. The switching is now spatial: an individual
experiences environment 1, then environment 2, because they travelled.

Two modelling choices define the problem:

* **Susceptible and infectious people commute, at the same rate $D$. The protected $Q$ do not** —
  someone isolating, or freshly vaccinated where they live, stays there.
* **The recovery rate $\gamma$ is the same in both neighbourhoods.** Only the transmission rate
  $\beta$ and the protection rates $\omega, v$ differ.

The result is a **window**:

| commuting rate | outcome |
|---|---|
| $D = 0$ | dies — by assumption, in each isolated neighbourhood |
| $0 < D < D_c$ | **survives** |
| $D > D_c$ | dies |
| $D \to \infty$ | dies — and this is a *theorem*, not a parameter choice (§3) |

with

$$\frac{1}{D_c} = \frac{1}{a_1} + \frac{1}{a_2}, \qquad a_i = \beta_i \bar S - \gamma,
\qquad \bar S = \frac{2}{\dfrac{1}{S^*_1}+\dfrac{1}{S^*_2}} .$$

This is the mirror image of switching in time, where persistence needs the switching to be *fast*.

---

## 1. The model

Neighbourhood $i=1,2$ has its own $\beta_i, \omega_i, v_i$ and the shared $\gamma$. As in the
well-mixed model write

$$\mu_i = \omega_i + v_i, \qquad S^*_i = \frac{\omega_i}{\mu_i},$$

so $S^*_i$ is the susceptible fraction neighbourhood $i$ settles at when it is alone and
disease-free. Letting $j$ denote the other neighbourhood,

$$\frac{dS_i}{dt} = \omega_i Q_i - v_i S_i - \beta_i S_i I_i + \gamma I_i + D\,(S_j - S_i)$$
$$\frac{dQ_i}{dt} = v_i S_i - \omega_i Q_i \qquad\qquad\text{(no diffusion term)}$$
$$\frac{dI_i}{dt} = \beta_i S_i I_i - \gamma I_i + D\,(I_j - I_i)$$

Densities are scaled so that each neighbourhood holds one unit of people when isolated. The grand
total $\sum_i(S_i+Q_i+I_i)=2$ is conserved, but the split between the two is not, because $S$ and
$I$ move and $Q$ does not.

**Assumption (dies alone).** In isolation, neighbourhood $i$ is disease-free at $S_i = S^*_i$, and
a small outbreak there grows at rate $\beta_iS^*_i-\gamma$. We require this to be negative:

$$\boxed{\;\beta_i S^*_i < \gamma \qquad (i=1,2)\;} \tag{1}$$

It is convenient to write this using each neighbourhood's **herd-immunity threshold**

$$S_{c,i} = \frac{\gamma}{\beta_i},$$

the susceptible fraction above which the disease would grow there. Condition (1) is then simply
$S^*_i < S_{c,i}$: each neighbourhood keeps its susceptible pool below its own threshold.

Label the neighbourhoods so that

$$S^*_1 > S^*_2 :$$

neighbourhood 1 is the **open** one (few people protected), neighbourhood 2 the **protected** one.

---

## 2. Commuting levels the susceptible pools — at *any* $D>0$

Set $I_1=I_2=0$ and look for the steady state. The $Q$ equation gives, in each neighbourhood
separately,

$$\omega_i Q_i = v_i S_i . \tag{2}$$

Substitute that into the $S$ equation. The terms $\omega_iQ_i - v_iS_i$ cancel **exactly**, and all
that survives is

$$\frac{dS_i}{dt} = D\,(S_j - S_i) = 0 \qquad\Longrightarrow\qquad \boxed{S_1 = S_2 \equiv \bar S}$$

This is the step to pause on. The only thing that changes $S_i$ locally is the exchange with
$Q_i$, and at steady state that exchange must balance *on its own*, in each neighbourhood. So the
diffusive flux has nothing left to balance against and must itself vanish — which forces the two
susceptible densities to be equal, however small $D$ is. There is no gradual approach to the mixed
state: at $D=0^+$ the levelling is already complete. (This is exactly what freezing $Q$ buys us,
and it is why the threshold below turns out to be an upper one.)

### What $\bar S$ is

From (2), $Q_i = (v_i/\omega_i)\bar S$, so the number of people living in neighbourhood $i$ is

$$S_i + Q_i = \bar S\left(1+\frac{v_i}{\omega_i}\right)
= \bar S\,\frac{\omega_i+v_i}{\omega_i} = \frac{\bar S}{S^*_i} .$$

The two must add to 2, so $\bar S/S^*_1 + \bar S/S^*_2 = 2$ and

$$\boxed{\;\bar S = \frac{2}{\dfrac{1}{S^*_1}+\dfrac{1}{S^*_2}}\;}
\qquad\text{the \emph{harmonic} mean of the two isolated pools.} \tag{3}$$

A harmonic mean lies strictly between its inputs, so $S^*_2 < \bar S < S^*_1$. Commuting therefore
**raises** the susceptible density in the protected neighbourhood and lowers it in the open one.
That is the whole source of the effect: neighbourhood 2 is topped up by commuters, and if the
top-up carries it past its own herd-immunity threshold $S_{c,2}$, the disease can grow there.

Note also what came along for the ride: the protected neighbourhood ends up holding *more people*
($\bar S/S^*_2$ of them, which is larger). $Q$ is a trap nobody can diffuse out of, so susceptibles
keep flowing in to refill it. That redistribution is what makes the average harmonic rather than
arithmetic, and §3 shows it is exactly what closes the door at large $D$.

---

## 3. Infinitely fast commuting always kills the disease

At large $D$ the infection is shared equally as well, so its growth rate is the plain average of
the two local rates (derived properly in §4). With $a_i = \beta_i\bar S-\gamma$,

$$\lambda(D\to\infty) = \frac{a_1+a_2}{2} = \frac{\beta_1+\beta_2}{2}\,\bar S - \gamma .$$

Substitute (3); the factors of 2 cancel and this is positive only if

$$\frac{\beta_1+\beta_2}{\dfrac{1}{S^*_1}+\dfrac{1}{S^*_2}} \;>\; \gamma . \tag{4}$$

Now look at that left-hand side. It has the form

$$\frac{p_1+p_2}{q_1+q_2}, \qquad p_i = \beta_i, \quad q_i = \frac{1}{S^*_i},$$

which is called a **mediant** of the two fractions $p_1/q_1 = \beta_1S^*_1$ and
$p_2/q_2 = \beta_2S^*_2$. A mediant always lies strictly between the two fractions it is built
from. (Quick proof: if $p_1/q_1 < p_2/q_2$ then $p_1q_2 < p_2q_1$; adding $p_1q_1$ to both sides
and dividing by $q_1(q_1+q_2)$ gives $p_1/q_1 < (p_1+p_2)/(q_1+q_2)$, and the other inequality
follows the same way.)

But assumption (1) says $\beta_1S^*_1 < \gamma$ **and** $\beta_2S^*_2 < \gamma$. A number squeezed
between two numbers that are both below $\gamma$ is itself below $\gamma$, so (4) can never hold:

$$\boxed{\;\lambda(D\to\infty) < 0 \quad\text{always.}\;}$$

**A disease that dies in each isolated neighbourhood also dies under infinitely fast commuting,
whatever you choose for $\beta_i, \omega_i, v_i$.** No parameter tuning escapes it.

This is the same mediant obstruction that rules the effect out in the plain homogeneous SIR model
(see §5 of `../../SIQS_switching/math/poissonSwitching.md`): with $\gamma$ uniform and $Q$ frozen,
everything collapses onto the single ratio $\beta_iS^*_i$, and there is nothing left to
anti-correlate. Switching in time escapes precisely because the two environments are averaged
*differently* there — arithmetically in $\omega$ and $\mu$ separately, giving
$\bar S_{\rm time} = (\omega_1+\omega_2)/(\mu_1+\mu_2)$, which is not the harmonic mean.

So: dies at $D=0$, dies at $D=\infty$. Everything interesting has to happen in between.

---

## 4. Finite $D$

Put a small $I_i$ into the disease-free state. Since $S_i = \bar S$ there, the linearised equations
are

$$\frac{dI_1}{dt} = a_1 I_1 + D(I_2-I_1), \qquad
  \frac{dI_2}{dt} = a_2 I_2 + D(I_1-I_2),
\qquad a_i = \beta_i\bar S-\gamma . \tag{5}$$

(The $S_i$ perturbations drop out at linear order: they enter $dI_i/dt$ only multiplied by $I_i$,
which is already small.)

The key point is that $a_i \neq \beta_iS^*_i-\gamma$: the disease is judged against the **mixed**
pool $\bar S$, not against each neighbourhood's own equilibrium. Since $\bar S < S^*_1$, we have
$a_1 < \beta_1S^*_1-\gamma < 0$ — the open neighbourhood can only get worse. But
$\bar S > S^*_2$, so $a_2$ can be positive. Explicitly,

$$a_2 > 0 \iff \bar S > S_{c,2},$$

and together with (1) this is the **interleaving condition** for the whole effect:

$$\boxed{\;S^*_2 \;<\; S_{c,2} \;<\; \bar S \;<\; S^*_1 \;<\; S_{c,1}\;} \tag{6}$$

The four marks alternate rather than nest. Reading it left to right: neighbourhood 2 keeps itself
below its own threshold, but only just, so the top-up from commuting pushes it over; meanwhile
neighbourhood 1 sits far below its threshold and stays there. Note that $S_{c,2} < S_{c,1}$ forces
$\beta_2 > \beta_1$ — **the protected neighbourhood must also be the more transmissible one**. It
is the same anti-correlation that the temporal problem needs, in spatial clothing.

### The growth rate

System (5) is $\dot{\vec I} = M \vec I$ with

$$M = \begin{pmatrix} a_1 - D & D \\ D & a_2 - D \end{pmatrix},$$

so $\operatorname{tr}M = a_1+a_2-2D$ and, using $(a_1-D)(a_2-D)-D^2$,

$$\det M = a_1a_2 - D(a_1+a_2).$$

The eigenvalues are $\tfrac12\bigl[\operatorname{tr}\pm\sqrt{\operatorname{tr}^2-4\det}\bigr]$, and
the discriminant simplifies beautifully:

$$\operatorname{tr}^2-4\det = (a_1+a_2-2D)^2 - 4a_1a_2 + 4D(a_1+a_2) = (a_1-a_2)^2 + 4D^2 ,$$

because expanding the square produces $-4D(a_1+a_2)$, which cancels the $+4D(a_1+a_2)$. Hence

$$\boxed{\;\lambda_{\max}(D) = \frac{a_1+a_2}{2} - D + \sqrt{\left(\frac{a_1-a_2}{2}\right)^2+D^2}\;}
\tag{7}$$

Read off the two ends:

* $D\to0^+$: the square root is $|a_1-a_2|/2$, so $\lambda_{\max}\to\max(a_1,a_2) = a_2 > 0$.
  The infection stays put in neighbourhood 2 and grows at that neighbourhood's own rate.
* $D\to\infty$: $\sqrt{c^2+D^2}\approx D+c^2/2D$, so $\lambda_{\max}\to\tfrac12(a_1+a_2) < 0$ by §3.

And in between $\lambda_{\max}$ falls **monotonically**, since

$$\frac{d}{dD}\left[-D+\sqrt{c^2+D^2}\right] = -1 + \frac{D}{\sqrt{c^2+D^2}} < 0 .$$

So commuting by the infectious is purely a suppressant: it drags the growth rate from the good
neighbourhood's value down towards the average, delocalising the infection out of the one place
where it is winning. Because the curve starts positive, ends negative, and is monotone, it crosses
zero exactly once.

### The critical rate

Set (7) to zero:

$$\sqrt{\left(\frac{a_1-a_2}{2}\right)^2+D^2} = D - \frac{a_1+a_2}{2} .$$

The right-hand side is positive because $a_1+a_2<0$, so we may square both sides:

$$\left(\frac{a_1-a_2}{2}\right)^2+D^2 = D^2 - D(a_1+a_2) + \left(\frac{a_1+a_2}{2}\right)^2 .$$

The $D^2$ cancel, and since $\bigl[(a_1-a_2)^2-(a_1+a_2)^2\bigr]/4 = -a_1a_2$, this is
$-a_1a_2 = -D(a_1+a_2)$, giving

$$\boxed{\;D_c = \frac{a_1a_2}{a_1+a_2}, \qquad\text{equivalently}\qquad
\frac{1}{D_c} = \frac{1}{a_1}+\frac{1}{a_2}\;} \tag{8}$$

the "parallel resistors" combination. ($D_c>0$ because $a_1a_2<0$ and $a_1+a_2<0$.) The window
widens without bound as $a_1+a_2\to0^-$, i.e. as the parameters approach the mediant bound of §3 —
the one thing they can never cross.

---

## 5. Why this is the opposite of switching in time

Both mechanisms rescue the disease by exposing it to an average of the two environments rather than
to either one. The difference is in how the average is delivered.

**In time**, the environments arrive one after the other. The disease has to survive the bad
environment long enough to reach the good one, so the switching has to be **fast**:
$\alpha > \alpha_c$. And the averaging is arithmetic, which leaves room for the effect.

**In space**, both environments exist at once. The disease does not need to wait for anything — it
just needs to *stay* in the neighbourhood that suits it. So commuting has to be **slow**:
$D < D_c$. Moving susceptibles around is what creates the good neighbourhood, and that happens at
$D=0^+$; moving infectious people around only smears the infection back out.

The price of freezing $Q$ is that the spatial average becomes harmonic, and §3 shows the harmonic
mean is exactly enough averaging to cancel the benefit. So the effect exists only at finite $D$,
squeezed between two extremes that both kill it.

---

## 6. Worked example

$$\beta = (1.05,\ 3.0), \qquad \gamma = 1, \qquad \omega = (9,\ 3), \qquad v = (1,\ 7)$$

* $S^*_1 = 9/10 = 0.9$ and $S^*_2 = 3/10 = 0.3$, so neighbourhood 2 is the protected one.
* Herd-immunity thresholds $S_{c,1} = 1/1.05 = 0.952$ and $S_{c,2} = 1/3 = 0.333$.
* Isolated growth rates $\beta_iS^*_i-\gamma = -0.055$ and $-0.100$: **the disease dies in both.**
* Harmonic mean $\bar S = 2/(1/0.9+1/0.3) = 0.45$.

Check the interleaving (6): $\;0.3 < 0.333 < 0.45 < 0.9 < 0.952$ ✓.

After mixing, $a_1 = 1.05(0.45)-1 = -0.5275$ and $a_2 = 3.0(0.45)-1 = +0.35$. So

$$\lambda(D\to0^+) = +0.35, \qquad \lambda(D\to\infty) = -0.0888,
\qquad D_c = \frac{(-0.5275)(0.35)}{-0.5275+0.35} = 1.04 .$$

The disease dies if nobody commutes, thrives at moderate commuting, and dies again once commuting
is fast. `../diffusion_threshold.py` checks Eq. (7) against direct integration of the full
nonlinear six-variable system over $D \in [0.01, 20]$; theory and simulation agree to
three or four decimals, apart from the smallest $D$, where levelling the susceptible pools takes a
time of order $1/D$ and the run is not much longer than that transient.
