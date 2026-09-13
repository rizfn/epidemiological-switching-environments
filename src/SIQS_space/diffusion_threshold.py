"""
SIQS_space: two neighbourhoods coupled by commuting.

S and I commute at rate D; the protected Q do not. The recovery rate is the SAME in both
neighbourhoods. The disease dies in each neighbourhood on its own, and dies again under
infinitely fast commuting -- but survives in a window of intermediate D.

Theory (math/criticalDiffusion.md) against direct RK4 integration of the full nonlinear
six-variable system.
"""
import os

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 20, "axes.labelsize": 24, "xtick.labelsize": 20, "ytick.labelsize": 20,
    "legend.fontsize": 17, "axes.linewidth": 1.4, "xtick.major.width": 1.4,
    "ytick.major.width": 1.4, "xtick.major.size": 7, "ytick.major.size": 7,
})

# neighbourhood 1 is the open one (large S*), neighbourhood 2 the protected one (small S*)
BETA = (1.05, 3.0)
GAMMA = 1.0                                  # the same in both, by construction
OMEGA = (9.0, 3.0)
V = (1.0, 7.0)

SSTAR = tuple(OMEGA[i] / (OMEGA[i] + V[i]) for i in range(2))
SBAR = 2.0 / sum(1.0 / s for s in SSTAR)     # harmonic mean, Eq. (2)
A = tuple(BETA[i] * SBAR - GAMMA for i in range(2))       # Eq. (5)

D_SIM = np.array([0.01, 0.05, 0.2, 0.5, 1.04, 2.0, 5.0, 20.0])


def lam_theory(D):
    """Eq. (6)."""
    return 0.5 * (A[0] + A[1]) - D + np.sqrt((0.5 * (A[0] - A[1])) ** 2 + D ** 2)


def D_c():
    """Eq. (7)."""
    return A[0] * A[1] / (A[0] + A[1])


def simulate(D, t_max=6000.0, dt=0.002):
    """Growth rate of the total infected, from the full nonlinear six-variable system."""
    b1, b2 = BETA
    g = GAMMA
    (w1, w2), (v1, v2) = OMEGA, V
    S1, S2 = SSTAR
    y = (S1, 1 - S1, 1e-12, S2, 1 - S2, 1e-12)
    log_gain = 0.0

    def f(S1, Q1, I1, S2, Q2, I2):
        return (w1 * Q1 - v1 * S1 - b1 * S1 * I1 + g * I1 + D * (S2 - S1),
                v1 * S1 - w1 * Q1,
                b1 * S1 * I1 - g * I1 + D * (I2 - I1),
                w2 * Q2 - v2 * S2 - b2 * S2 * I2 + g * I2 + D * (S1 - S2),
                v2 * S2 - w2 * Q2,
                b2 * S2 * I2 - g * I2 + D * (I1 - I2))

    for _ in range(int(t_max / dt)):
        k1 = f(*y)
        k2 = f(*(y[i] + dt / 2 * k1[i] for i in range(6)))
        k3 = f(*(y[i] + dt / 2 * k2[i] for i in range(6)))
        k4 = f(*(y[i] + dt * k3[i] for i in range(6)))
        y = tuple(y[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6))
        tot = y[2] + y[5]
        if tot > 1e-9 or tot < 1e-15:        # keep I infinitesimal, bank the growth
            log_gain += np.log(tot / 1e-12)
            y = (y[0], y[1], y[2] / tot * 1e-12, y[3], y[4], y[5] / tot * 1e-12)
    return (log_gain + np.log((y[2] + y[5]) / 1e-12)) / t_max


def sim_points():
    """Integration is slow, so keep the results next to the script."""
    cache = "data/lambda_vs_D.csv"
    if os.path.exists(cache):
        return list(np.loadtxt(cache, delimiter=",")[:, 1])
    sim = [simulate(D) for D in D_SIM]
    os.makedirs("data", exist_ok=True)
    np.savetxt(cache, np.column_stack([D_SIM, sim]), delimiter=",")
    return sim


def main():
    Dc = D_c()
    Sc = tuple(GAMMA / b for b in BETA)      # herd-immunity thresholds
    print(f"S*_i          = {SSTAR}       (disease-free susceptible pools, isolated)")
    print(f"S_c,i = g/b_i = {Sc}       (herd-immunity thresholds)")
    print(f"Sbar (harmonic mean) = {SBAR:.4f}")
    print(f"isolated rates b_i S*_i - g = "
          f"{tuple(BETA[i]*SSTAR[i]-GAMMA for i in range(2))}  (both < 0: dies alone)")
    print(f"a_i = b_i Sbar - g   = {A}")
    print(f"lambda(D->0+) = {max(A):+.5f}   lambda(D->inf) = {0.5*(A[0]+A[1]):+.5f}")
    print(f"D_c = {Dc:.5f}")

    sim = sim_points()
    for D, s in zip(D_SIM, sim):
        print(f"   D={D:<6} theory={lam_theory(D): .6f}  sim={s: .6f}")

    Ds = np.logspace(-2.3, 1.6, 500)
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9.5, 9.0),
                                   gridspec_kw={"height_ratios": [1, 2.2]})

    # --- top: where the susceptible pools sit relative to the herd-immunity thresholds
    ax0.axhline(0, color="0.35", lw=1.6, zorder=1)
    # equilibria labelled below the line, herd-immunity thresholds above it
    below = [(SSTAR[1], r"$S^*_2$", "#1f5c99"), (SSTAR[0], r"$S^*_1$", "#1f5c99")]
    above = [(Sc[1], r"$S_{c,2}$", "#c73e1d"), (SBAR, r"$\bar S$", "black"),
             (Sc[0], r"$S_{c,1}$", "#c73e1d")]
    for pos, lab, col in below + above:
        ax0.plot([pos], [0], "|", ms=24, mew=3, color=col, zorder=2)
    for pos, lab, col in below:
        ax0.annotate(lab, (pos, -1.15), color=col, fontsize=21, ha="center")
    for pos, lab, col in above:
        ax0.annotate(lab, (pos, 0.3), color=col, fontsize=21, ha="center")
    ax0.annotate("", xy=(SBAR, 1.5), xytext=(SSTAR[1], 1.5),
                 arrowprops=dict(arrowstyle="->", lw=2.2, color="0.25"))
    ax0.annotate("commuting lifts neighbourhood 2\npast its own threshold",
                 (SSTAR[1], 1.9), fontsize=16, ha="left", color="0.25")
    ax0.set_xlim(0.22, 1.02)
    ax0.set_ylim(-1.8, 4.2)
    ax0.set_yticks([])
    ax0.set_xlabel("susceptible fraction")
    for s in ("left", "right", "top"):
        ax0.spines[s].set_visible(False)

    # --- bottom: the persistence window
    ax1.axhline(0, color="0.6", lw=1.2, zorder=0)
    ax1.axvspan(Ds[0], Dc, color="#c73e1d", alpha=0.09, zorder=0)
    ax1.semilogx(Ds, [lam_theory(D) for D in Ds], color="#1f5c99", lw=3,
                 label="theory, Eq. (6)")
    ax1.semilogx(D_SIM, sim, "o", ms=11, mfc="none", mew=2.5, color="#c73e1d",
                 label="simulation")
    ax1.axvline(Dc, color="k", ls=":", lw=2)
    ax1.annotate(rf"$D_c={Dc:.2f}$", (Dc * 1.25, 0.22), fontsize=20)
    ax1.annotate("disease persists", (0.02, -0.055), fontsize=18, color="#c73e1d")
    ax1.annotate(r"$\lambda\to\frac{a_1+a_2}{2}<0$", (6.0, -0.055), fontsize=18, color="0.3")
    ax1.set_xlabel(r"commuting rate $D$")
    ax1.set_ylabel(r"invasion exponent $\lambda_{\max}$")
    ax1.set_xlim(Ds[0], Ds[-1])
    ax1.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    name = (f"lambda_vs_D_beta_{BETA[0]}-{BETA[1]}_gamma_{GAMMA}"
            f"_omega_{OMEGA[0]}-{OMEGA[1]}_v_{V[0]}-{V[1]}")
    for ext in ("svg", "pdf"):
        fig.savefig(f"plots/{name}.{ext}")
    print("saved", name)


if __name__ == "__main__":
    main()
