Similar to an SIS model with Susceptible $S$, infectious $I$.

Add a protected state $Q$, of people who vaccine/isolate and so cannot be infected. People enter $Q$ with rate $v$ and leave with rate $\omega$.

$$ \frac{dS}{dt} = \omega Q - v S - \beta S I + \gamma I $$
$$ \frac{dP}{dt} = v S - \omega Q $$
$$ \frac{dI}{dt} = \beta S I - \gamma I $$