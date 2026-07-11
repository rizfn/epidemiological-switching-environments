# Derivation of Critical Switching Rate for Epidemiological Rescue

Based on the provided epidemiological switching-rescue model and the Floquet theory approach presented in the paper *"A loser in both environments can survive by switching between them"*, this document derives the critical switching rate $\alpha_c$ required for the disease to persist in a dynamic environment, despite being driven to extinction in any static environment.

---

## 1. The Epidemiological Model

We consider a variation of the SIS model that includes a protected/quarantine state $Q$, consisting of individuals who are vaccinated or isolated and thus cannot be infected. 
* Individuals enter $Q$ from the susceptible pool $S$ with rate $v$.
* Individuals leave $Q$ and return to $S$ with rate $\omega$.
* Infection occurs at rate $\beta S I$, and recovery occurs at rate $\gamma I$.

The governing differential equations are:
$$\frac{dS}{dt} = \omega Q - v S - \beta S I + \gamma I$$
$$\frac{dQ}{dt} = v S - \omega Q$$
$$\frac{dI}{dt} = \beta S I - \gamma I$$

Assuming a conserved population $S + I + Q = 1$, we substitute $Q = 1 - S - I$ into the first equation to reduce the system to two dimensions:
$$\frac{dS}{dt} = \omega(1 - S - I) - v S - \beta S I + \gamma I$$
$$\frac{dI}{dt} = I(\beta S - \gamma)$$

To simplify notation in our derivation, we define the total rate of leaving the susceptible state (excluding infection) as $\mu = \omega + v$. Thus, the boundary flow for $S$ when $I=0$ becomes:
$$\frac{dS}{dt} = \omega - \mu S$$

## 2. Extinction in a Fixed Environment

In a static environment, the disease-free equilibrium occurs when $I^* = 0$. Substituting this into the $S$ equation yields the steady-state susceptible population:
$$S^* = \frac{\omega}{\omega + v} = \frac{\omega}{\mu}$$

For the disease to be eliminated (i.e., for $I=0$ to be a stable fixed point), the per-capita growth rate of the infection must be negative at this equilibrium:
$$\frac{1}{I}\frac{dI}{dt} \bigg|_{I=0} = \beta S^* - \gamma < 0$$
$$\beta \left( \frac{\omega}{\omega+v} \right) < \gamma \implies \frac{\gamma}{\beta} \left( \frac{v}{\omega} + 1 \right) > 1$$

We assume that environments 1 and 2 **both** satisfy this extinction condition, meaning the disease dies out if the environment remains statically in either state.

## 3. Fast Switching Limit (Time-Averaged Environments)

When the environment switches infinitely fast between Environment 1 (with parameters $\omega_1, v_1, \beta_1, \gamma_1$) and Environment 2 (with parameters $\omega_2, v_2, \beta_2, \gamma_2$) with an even duty cycle ($p = 0.5$), the system experiences time-averaged parameters:
$$\bar{\beta} = \frac{\beta_1 + \beta_2}{2}, \quad \bar{\gamma} = \frac{\gamma_1 + \gamma_2}{2}, \quad \bar{\omega} = \frac{\omega_1 + \omega_2}{2}, \quad \bar{v} = \frac{v_1 + v_2}{2}$$

The disease survives in the fast-switching limit if the condition is reversed for the time-averaged system:
$$\frac{\bar{\gamma}}{\bar{\beta}} \left( \frac{\bar{v}}{\bar{\omega}} + 1 \right) < 1$$

## 4. Deterministic Periodic Switching & Floquet Theory

To find the critical switching rate $\alpha_c$ that bridges the gap between slow switching (extinction) and fast switching (survival), we apply Floquet theory to the periodic boundary orbit $I=0$. 

Let the environment switch periodically with rate $\alpha$. Each environment is held for a dwell time $\tau = 1/\alpha$. A full cycle takes $T = 2/\alpha$.

### 4.1 Boundary Flow of Susceptibles
On the disease-free boundary ($I=0$), the dynamics of $S$ in environment $E \in \{1, 2\}$ are governed by:
$$\frac{dS}{dt} = \omega_E - \mu_E S$$

This is a linear ordinary differential equation, which can be solved exactly for an initial condition $S_{in}$:
$$\mathcal{S}_E(t; S_{in}) = S^*_E + (S_{in} - S^*_E)e^{-\mu_E t}$$
where $S^*_E = \omega_E / \mu_E$ is the fixed point of environment $E$.

### 4.2 Periodic Cycle Endpoints
After a transient period, the susceptible population settles into a periodic orbit. Let $S_a$ be the population at the end of Environment 1, and $S_b$ be the population at the end of Environment 2. These turning points satisfy:
$$S_a = \mathcal{S}_1(1/\alpha; S_b)$$
$$S_b = \mathcal{S}_2(1/\alpha; S_a)$$

By defining $E_1 = e^{-\mu_1/\alpha}$ and $E_2 = e^{-\mu_2/\alpha}$, we can solve for $S_a$ and $S_b$ algebraically:
$$S_a = \frac{S^*_1 (1 - E_1) + S^*_2 E_1 (1 - E_2)}{1 - E_1 E_2}$$
$$S_b = \frac{S^*_2 (1 - E_2) + S^*_1 E_2 (1 - E_1)}{1 - E_1 E_2}$$

### 4.3 Integration of the Boundary Orbit
To evaluate whether a small infection $I$ will grow, we must calculate the Floquet exponent $\lambda(\alpha)$, which represents the average per-capita growth rate of $I$ over one full cycle:
$$\lambda(\alpha) = \frac{\alpha}{2} \int_0^{2/\alpha} (\beta(t)S(t) - \gamma(t)) dt$$
$$\lambda(\alpha) = \frac{\alpha}{2} \left[ \int_0^{1/\alpha} (\beta_1 \mathcal{S}_1(t; S_b) - \gamma_1) dt + \int_0^{1/\alpha} (\beta_2 \mathcal{S}_2(t; S_a) - \gamma_2) dt \right]$$

We utilize an integral identity for the boundary flow. Since $\frac{dS}{dt} = \omega_E - \mu_E S$, we have $S = S^*_E - \frac{1}{\mu_E}\frac{dS}{dt}$. Integrating this from $t=0$ to $t=\tau$ yields:
$$\int_0^\tau S(t) dt = S^*_E \tau - \frac{S(\tau) - S(0)}{\mu_E}$$

Applying this to each environment dwell:
* **Environment 1:** (Starts at $S_b$, ends at $S_a$)
    $$\int_0^{1/\alpha} S_1(t) dt = S^*_1 \frac{1}{\alpha} - \frac{S_a - S_b}{\mu_1}$$
* **Environment 2:** (Starts at $S_a$, ends at $S_b$)
    $$\int_0^{1/\alpha} S_2(t) dt = S^*_2 \frac{1}{\alpha} - \frac{S_b - S_a}{\mu_2} = S^*_2 \frac{1}{\alpha} + \frac{S_a - S_b}{\mu_2}$$

### 4.4 The Floquet Exponent
Substituting the integrals into the expression for $\lambda(\alpha)$:
$$\lambda(\alpha) = \frac{\alpha}{2} \left[ \beta_1 \left( \frac{S^*_1}{\alpha} - \frac{S_a - S_b}{\mu_1} \right) - \frac{\gamma_1}{\alpha} + \beta_2 \left( \frac{S^*_2}{\alpha} + \frac{S_a - S_b}{\mu_2} \right) - \frac{\gamma_2}{\alpha} \right]$$

Grouping the terms, we obtain a beautifully closed-form exact expression analogous to Equation 7 from the reference paper:
$$\lambda(\alpha) = \frac{(\beta_1 S^*_1 - \gamma_1) + (\beta_2 S^*_2 - \gamma_2)}{2} - \frac{\alpha}{2} \left( \frac{\beta_1}{\mu_1} - \frac{\beta_2}{\mu_2} \right) (S_a - S_b)$$

The critical switching rate $\alpha_c^{det}$ is the value at which the disease neither grows nor decays ($\lambda = 0$). Thus, $\alpha_c^{det}$ is the solution to the transcendental equation:
$$\frac{(\beta_1 S^*_1 - \gamma_1) + (\beta_2 S^*_2 - \gamma_2)}{2} = \frac{\alpha_c^{det}}{2} \left( \frac{\beta_1}{\omega_1 + v_1} - \frac{\beta_2}{\omega_2 + v_2} \right) (S_a - S_b)$$

*Note: The first term on the left side is simply the average of the per-capita growth rates at the respective fixed points. Because both environments are lethal to the disease, this term is strictly negative. Therefore, for a solution to exist, the dynamic correction term on the right must perfectly balance this baseline decay.*

---

## 5. Stochastic (Poisson) Switching

If the environmental switching is modeled as a Poisson process where the environment changes randomly with an average rate $\alpha$, the dwell times $\tau_1, \tau_2$ are exponentially distributed random variables with mean $1/\alpha$. 

Following the methodology of the paper, the expected log-increment (or in our linear boundary case, the direct density increment) over a full cycle must be zero in the stationary state: $\mathbb{E}[S_{out} - S_{in}] = 0$. 

Defining the stationary mean increment over an environment-1 dwell as $\mathcal{A}(\alpha) = \mathbb{E}[S_a - S_b]$, the stochastic invasion exponent becomes:
$$\lambda^{stoch}(\alpha) = \frac{(\beta_1 S^*_1 - \gamma_1) + (\beta_2 S^*_2 - \gamma_2)}{2} - \frac{\alpha}{2} \left( \frac{\beta_1}{\mu_1} - \frac{\beta_2}{\mu_2} \right) \mathcal{A}(\alpha)$$

The critical stochastic switching rate $\alpha_c^{stoch}$ is the root of $\lambda^{stoch}(\alpha_c) = 0$. Because long spells in the unfavorable environment yield a linear penalty to survival, but bounded benefits in the favorable environment due to the boundary fixed point, we mathematically expect $\alpha_c^{stoch} > \alpha_c^{det}$, requiring faster switching to compensate for variance.