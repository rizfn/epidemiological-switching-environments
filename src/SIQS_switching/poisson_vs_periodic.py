"""
Invasion exponent of the SIQS model under periodic and Poisson environmental switching.

Theory (closed forms from math/poissonSwitching.md) against direct RK4 simulation of the
full nonlinear model.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

plt.rcParams.update({
    "font.size": 20, "axes.labelsize": 24, "xtick.labelsize": 20, "ytick.labelsize": 20,
    "legend.fontsize": 18, "axes.linewidth": 1.4, "xtick.major.width": 1.4,
    "ytick.major.width": 1.4, "xtick.major.size": 7, "ytick.major.size": 7,
    "figure.autolayout": True,
})

# --- environment 1 (protective: large S*), environment 2 (transmissive: small S*)
BETA = (1.0, 1.5)
GAMMA = (1.0, 1.0)
OMEGA = (9.0, 1.0)
V = (1.0, 1.0)

MU = tuple(OMEGA[i] + V[i] for i in range(2))
SSTAR = tuple(OMEGA[i] / MU[i] for i in range(2))
A_E = tuple(BETA[i] * SSTAR[i] - GAMMA[i] for i in range(2))


def swing(alpha, stochastic):
    """<S_a> - <S_b>, Eq. (2) / (4)."""
    if stochastic:
        e = [alpha / (alpha + MU[i]) for i in range(2)]
    else:
        e = [np.exp(-MU[i] / alpha) for i in range(2)]
    return (SSTAR[0] - SSTAR[1]) * (1 - e[0]) * (1 - e[1]) / (1 - e[0] * e[1])


def lam_theory(alpha, stochastic):
    return (A_E[0] + A_E[1]) / 2 - (alpha / 2) * (
        BETA[0] / MU[0] - BETA[1] / MU[1]
    ) * swing(alpha, stochastic)


def lam_stoch_rational(alpha):
    """Eq. (6): the exponential-free form, as a check on the algebra."""
    G = (BETA[0] / MU[0] - BETA[1] / MU[1]) * (SSTAR[1] - SSTAR[0])
    return (A_E[0] + A_E[1]) / 2 + (alpha / 2) * G / (1 + alpha * (1 / MU[0] + 1 / MU[1]))


def alpha_c_closed_form():
    """Eq. (7)."""
    A = -(A_E[0] + A_E[1])
    G = (BETA[0] / MU[0] - BETA[1] / MU[1]) * (SSTAR[1] - SSTAR[0])
    return A / (G - A * (1 / MU[0] + 1 / MU[1]))


def simulate(alpha, stochastic, t_max=6000.0, dt=0.002, seed=0):
    """Growth rate of ln I from the full nonlinear model, linearised by rescaling I."""
    rng = np.random.default_rng(seed)
    S, I = SSTAR[0], 1e-12
    env = 0
    t_switch = rng.exponential(1 / alpha) if stochastic else 1 / alpha
    t, log_gain, n = 0.0, 0.0, int(t_max / dt)

    def deriv(S, I, e):
        dS = OMEGA[e] * (1 - S - I) - V[e] * S - BETA[e] * S * I + GAMMA[e] * I
        return dS, I * (BETA[e] * S - GAMMA[e])

    for _ in range(n):
        if t >= t_switch:
            env ^= 1
            t_switch = t + (rng.exponential(1 / alpha) if stochastic else 1 / alpha)
        k1 = deriv(S, I, env)
        k2 = deriv(S + dt / 2 * k1[0], I + dt / 2 * k1[1], env)
        k3 = deriv(S + dt / 2 * k2[0], I + dt / 2 * k2[1], env)
        k4 = deriv(S + dt * k3[0], I + dt * k3[1], env)
        S += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        I += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t += dt
        if I > 1e-9 or I < 1e-15:          # keep I infinitesimal, bank the growth
            log_gain += np.log(I / 1e-12)
            I = 1e-12
    log_gain += np.log(I / 1e-12)
    return log_gain / t_max


def main():
    alphas = np.logspace(np.log10(0.3), np.log10(60), 200)
    ac_det = brentq(lambda a: lam_theory(a, False), 0.5, 100)
    ac_sto = alpha_c_closed_form()
    print(f"alpha_c det   = {ac_det:.4f}")
    print(f"alpha_c stoch = {ac_sto:.4f} (closed form)")
    assert abs(lam_stoch_rational(3.7) - lam_theory(3.7, True)) < 1e-12

    a_sim = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 20.0, 50.0])
    sim_det = [simulate(a, False) for a in a_sim]
    sim_sto = [np.mean([simulate(a, True, seed=s) for s in range(6)]) for a in a_sim]

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.axhline(0, color="0.6", lw=1.2, zorder=0)
    ax.plot(alphas, [lam_theory(a, False) for a in alphas], color="#1f5c99", lw=3,
            label="periodic (theory)")
    ax.plot(alphas, [lam_stoch_rational(a) for a in alphas], color="#c73e1d", lw=3,
            label="Poisson (theory)")
    ax.plot(a_sim, sim_det, "o", ms=11, mfc="none", mew=2.5, color="#1f5c99",
            label="periodic (sim.)")
    ax.plot(a_sim, sim_sto, "s", ms=10, mfc="none", mew=2.5, color="#c73e1d",
            label="Poisson (sim.)")
    for xc, c in ((ac_det, "#1f5c99"), (ac_sto, "#c73e1d")):
        ax.axvline(xc, color=c, ls=":", lw=2)
    ax.annotate(rf"$\alpha_c^{{\rm det}}={ac_det:.2f}$", (ac_det * 0.92, -0.145),
                color="#1f5c99", ha="right", va="bottom", fontsize=19, rotation=90)
    ax.annotate(rf"$\alpha_c^{{\rm stoch}}={ac_sto:.2f}$", (ac_sto * 0.92, -0.145),
                color="#c73e1d", ha="right", va="bottom", fontsize=19, rotation=90)
    ax.set_xscale("log")
    ax.set_xlabel(r"switching rate $\alpha$")
    ax.set_ylabel(r"invasion exponent $\lambda$")
    ax.legend(frameon=False, loc="lower right")
    ax.set_xlim(alphas[0], alphas[-1])

    name = (f"lambda_vs_alpha_beta_{BETA[0]}-{BETA[1]}_gamma_{GAMMA[0]}-{GAMMA[1]}"
            f"_omega_{OMEGA[0]}-{OMEGA[1]}_v_{V[0]}-{V[1]}")
    for ext in ("svg", "pdf"):
        fig.savefig(f"plots/{name}.{ext}")
    print("saved", name)


if __name__ == "__main__":
    main()
