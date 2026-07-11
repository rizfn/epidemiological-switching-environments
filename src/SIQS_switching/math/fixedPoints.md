### Epidemiological switching-rescue

Similar to an SIS model with Susceptible $S$, infectious $I$.

Add a protected state $Q$, of people who vaccine/isolate and so cannot be infected. People enter $Q$ with rate $v$ and leave with rate $\omega$.

$$ \frac{dS}{dt} = \omega Q - v S - \beta S I + \gamma I $$
$$ \frac{dP}{dt} = v S - \omega Q $$
$$ \frac{dI}{dt} = \beta S I - \gamma I $$

As the population is conserved, $S+I+Q=1$. Substituting for $Q$ in the first equation.

$$ \frac{dS}{dt} = \omega(1 - S - I) - v S - \beta S I + \gamma I $$
$$ \frac{dI}{dt} = I(\beta S - \gamma) $$

Fixed points:

$I=0$, or $S=\gamma/\beta$

$$
\begin{align*}
0 &= \omega(1-\gamma/\beta-I) - v\gamma/\beta - \gamma I +\gamma I \\
\frac{v\gamma}{\beta} &= \omega(1-\gamma/\beta-I) \\
\frac{v\gamma}{\beta\omega} &= (1-\gamma/\beta-I) \\
I &= 1 - \frac{v\gamma}{\beta\omega} - \frac{\gamma}\beta
\end{align*}
$$

For the system to die out in the fixed environment, we need $I<0$, and so

$$ \frac{\gamma}{\beta} \left( \frac{v}\omega +1 \right) > 1$$


Suppose the population switches between Environment 1 (with $\omega_1, v_1, \beta_1, \gamma_1$) and Environment 2 (with $\omega_2, v_2, \beta_2, \gamma_2$) in a duty cycle of $p$. 

$$
\begin{align*}
    \frac{dI}{dt} &= p \Big[ I(\beta_1 S - \gamma_1) \Big] + (1-p) \Big[ I(\beta_2 S - \gamma_2) \Big] \\
    &= I \Big[ (p\beta_1 + (1-p)\beta_2)S - (p\gamma_1 + (1-p)\gamma_2) \Big]
\end{align*}
$$

By defining time-averaged parameters:

$$
\begin{align*}
\bar{\beta} &= p\beta_1 + (1-p)\beta_2 \\
\bar{\gamma} &= p\gamma_1 + (1-p)\gamma_2 \\
\bar{\omega} &= p\omega_1 + (1-p)\omega_2 \\
\bar{v} &= p v_1 + (1-p)v_2
\end{align*}
$$

We get

$$
\begin{equation*}
    \frac{dI}{dt} = I(\bar{\beta} S - \bar{\gamma})
\end{equation*}
$$

Doing the same for $S$:

$$
\begin{align*}
    \frac{dS}{dt} &= p \Big[ \omega_1(1 - S - I) - v_1 S - \beta_1 S I + \gamma_1 I \Big] \\
    &\quad + (1-p) \Big[ \omega_2(1 - S - I) - v_2 S - \beta_2 S I + \gamma_2 I \Big] \\
    \frac{dS}{dt} &= (1 - S - I)[p\omega_1 + (1-p)\omega_2] - S[p v_1 + (1-p)v_2] \\
    &\quad - SI[p\beta_1 + (1-p)\beta_2] + I[p\gamma_1 + (1-p)\gamma_2] \\
    \frac{dS}{dt} &= \bar{\omega}(1 - S - I) - \bar{v} S - \bar{\beta} S I + \bar{\gamma} I
\end{align*}
$$

As the structure is identical, the condition for switching to rescue the system implies:

$$ \frac{\bar\gamma}{\bar\beta} \left( \frac{\bar{v}}{\bar{\omega}} +1 \right) < 1$$


### Trial parameter set:

Even duty cycle, $p=0.5$.

$\gamma=v=1$, in both environments.

Death condition:

$$ \frac{1}{\beta} \left( \frac1\omega +1 \right) > 1$$
$$ \left( \frac1\omega +1 \right) > \beta$$

Survival condition:

$$ \left( \frac1{\bar\omega} +1 \right) < \bar\beta$$


Trying $\beta_1=1, \omega_1=9, \beta_2=1.5, \omega_2=1$:


Environment 1: 

$$1 + \frac1{9} > 1$$

Environment 2:

$$ 1 + 1 > 1.5$$

Switching:

$$ \frac15 + 1 < 1.25$$

So it works!