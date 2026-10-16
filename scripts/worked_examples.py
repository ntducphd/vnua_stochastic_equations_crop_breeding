# -*- coding: utf-8 -*-
"""Worked examples quoted in the review. Prints each number with the sentence it supports and
writes the same text to worked_examples_output.txt in this folder.

1  Fixation from x0 = 0.5 (fixation table and the Wright-Fisher worked example, Ne = 100, s = 0.15).
   Exact chain (wf_exact.transition; moments from monte_carlo_check.exact_moments):
   P_fix, and the mean T1 and standard deviation of the number of generations to fixation among
   runs that fix. Diffusion: Kimura's u(x) = (1 - exp(-4 Ne s x)) / (1 - exp(-4 Ne s)) and the
   conditional mean time of Kimura and Ohta for the generator L = M d/dx + (V/2) d2/dx2 with
   M(x) = s x (1 - x) and V(x) = x (1 - x) / (2 Ne):
     T1(p) = int_p^1 psi u (1 - u) dx + [(1 - u(p)) / u(p)] int_0^p psi u^2 dx,
     psi(x) = 2 int_0^1 G(y) dy / [V(x) G(x)],  G(x) = exp(-4 Ne s x),
   integrated with scipy.integrate.quad. With S = 4 Ne s the integrands simplify to
     psi u (1 - u)                 = (4 Ne / S) (1 - e^(-S x)) (1 - e^(-S (1 - x))) / [x (1 - x) (1 - e^(-S))]
     [(1 - u(p)) / u(p)] psi u^2   = (4 Ne / S) [(1 - e^(-S (1 - p))) / (1 - e^(-S p))]
                                     e^(-S (p - x)) (1 - e^(-S x))^2 / [x (1 - x) (1 - e^(-S))],
   which are evaluated with expm1 and stay finite for any S. Computing 1 - u(p) by subtraction
   would lose the second term once e^(-S p) < 1e-16: at S = 80, 1 - u(0.5) = 4e-18, yet the term
   contributes 0.5 generation. The direct form, with 1 - u computed as in kimura_1mu, is
   evaluated as well and must agree. Checks: the neutral limit 4 Ne ln 2 (from s -> 0), and the
   strong-selection approximation (ln(4 Ne s) + gamma) / s, gamma = 0.5772 (Euler's constant).
2  Branching process. A copy leaves a Poisson(1 + s) number of copies (the large-population limit
   of the binomial sampling of the chain). The establishment probability pi solves
   1 - pi = exp(-(1 + s) pi); it is compared with 2s and with Kimura's u(1/(2 Ne)) for Ne = 100
   and 1,000. With k independent copies (edited events, introgressed carriers), the allele
   establishes with probability 1 - (1 - pi)^k.
3  Replicates needed by a simulation: R = (SD / SE)^2 runs that fix give a standard error SE of
   the mean conditional fixation time; R = p (1 - p) / SE^2 runs give SE for a probability p.
4  Monte Carlo check of the exact values (monte_carlo_check.run: 100,000 replicates per scenario,
   seed 20260926).

Usage: python -I -B scripts/worked_examples.py [replicates per scenario for part 4]
"""
import math
import os
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # python -I leaves the script folder off the path
import monte_carlo_check as mc  # noqa: E402
from wf_exact import ROWS, transition  # noqa: E402

GAMMA = 0.5772156649015329
EXAMPLE = ("Worked example (Ne = 100, s = 0.15)", 100, 0.15)
OUTFILE = os.path.join(HERE, "worked_examples_output.txt")
QUAD = dict(epsabs=0.0, epsrel=1e-11, limit=500)


# --------------------------------------------------------------------------- diffusion
def kimura_u(x, Ne, s):
    S = 4 * Ne * s
    return x if S == 0 else math.expm1(-S * x) / math.expm1(-S)


def kimura_1mu(x, Ne, s):
    """1 - u(x) = e^(-S x) (1 - e^(-S (1 - x))) / (1 - e^(-S)), S = 4 Ne s. Written this way because
    1 - u(x) computed by subtraction rounds to 0 once e^(-S x) falls below 1e-16 (S = 80, x = 0.5)."""
    S = 4 * Ne * s
    return 1 - x if S == 0 else math.exp(-S * x) * math.expm1(-S * (1 - x)) / math.expm1(-S)


def tbar1_diffusion(Ne, s, p=0.5):
    """Kimura-Ohta conditional mean time to fixation from p (generations), simplified integrands."""
    S = 4 * Ne * s
    if S == 0:
        return -4 * Ne * (1 - p) * math.log1p(-p) / p
    c = 4 * Ne / S / -math.expm1(-S)

    def a(x):      # psi u (1 - u)
        return c * -math.expm1(-S * x) * -math.expm1(-S * (1 - x)) / (x * (1 - x))

    k = -math.expm1(-S * (1 - p)) / -math.expm1(-S * p)

    def b(x):      # [(1 - u(p)) / u(p)] psi u^2, with e^(-S p) (e^(S x) - 1) = e^(-S (p - x)) (1 - e^(-S x))
        return c * k * math.exp(-S * (p - x)) * math.expm1(-S * x) ** 2 / (x * (1 - x))

    return quad(a, p, 1, **QUAD)[0] + quad(b, 0, p, **QUAD)[0]


def tbar1_direct(Ne, s, p=0.5):
    """The same integral written term by term (G, V, psi, u), as a check on the simplification."""
    S = 4 * Ne * s
    G = lambda x: math.exp(-S * x)                        # noqa: E731
    V = lambda x: x * (1 - x) / (2 * Ne)                  # noqa: E731
    IG = quad(G, 0, 1, **QUAD)[0]
    psi = lambda x: 2 * IG / (V(x) * G(x))                # noqa: E731
    u = lambda x: kimura_u(x, Ne, s)                      # noqa: E731
    v = lambda x: kimura_1mu(x, Ne, s)                    # noqa: E731
    first = quad(lambda x: psi(x) * u(x) * v(x), p, 1, **QUAD)[0]
    second = quad(lambda x: psi(x) * u(x) ** 2, 0, p, **QUAD)[0]
    return first + v(p) / u(p) * second


# --------------------------------------------------------------------------- branching process
def establishment(s):
    """Non-extinction probability of a Galton-Watson process with Poisson(1 + s) offspring."""
    lam = 1 + s
    return brentq(lambda q: -math.expm1(-lam * q) - q, s / (10 * lam), 1.0, xtol=1e-15, rtol=1e-14)


def copies_needed(pi, target):
    return math.ceil(math.log(1 - target) / math.log(1 - pi) - 1e-9)


# --------------------------------------------------------------------------- report
def main(reps):
    out = []
    say = out.append
    say("Worked examples for the review (worked_examples.py)")
    say("python %s, numpy %s, scipy %s" % (platform.python_version(), np.__version__, scipy.__version__))

    # 1 ---------------------------------------------------------------------------------------------
    say("")
    say("1. Fixation from x0 = 0.5: exact Wright-Fisher chain and diffusion (generations)")
    say("%-38s %4s %5s %6s | %8s %8s | %7s %7s %6s | %6s" % ("scenario", "Ne", "s", "4Nes", "Pfix_ex", "Pfix_dif",
                                                            "T1_ex", "T1_dif", "dif-ex", "SD_ex"))
    scen = [(n, Ne, s) for n, Ne, s in ROWS] + [EXAMPLE]
    res = {}
    for name, Ne, s in scen:
        u, t_ex, sd = mc.exact_moments(Ne, s)
        t_dif = tbar1_diffusion(Ne, s)
        if s:
            t_dir = tbar1_direct(Ne, s)
            assert abs(t_dir - t_dif) < 1e-6 * t_dif, (name, t_dir, t_dif)
        res[(Ne, s)] = (u, t_ex, sd, t_dif)
        say("%-38s %4d %5.2f %6.1f | %8.5f %8.5f | %7.2f %7.2f %+6.2f | %6.2f"
            % (name[:38], Ne, s, 4 * Ne * s, u, kimura_u(0.5, Ne, s), t_ex, t_dif, t_dif - t_ex, sd))
    gaps = [(abs(res[(Ne, s)][0] - kimura_u(0.5, Ne, s)), Ne, s) for _, Ne, s in scen]
    g, gNe, gs = max(gaps)
    say("Largest |P_fix diffusion - exact| over these rows: %.4f (Ne = %d, s = %.2f: %.3f against %.3f)."
        % (g, gNe, gs, kimura_u(0.5, gNe, gs), res[(gNe, gs)][0]))
    tgap = max((abs(res[(Ne, s)][3] - res[(Ne, s)][1]), Ne, s) for _, Ne, s in scen)
    say("Largest |T1 diffusion - exact|: %.2f generations (Ne = %d, s = %.2f)." % tgap)
    for _, Ne, s in scen:
        if s == 0:
            lim = tbar1_diffusion(Ne, 1e-9)
            say("Neutral limit, Ne = %d: 4 Ne ln 2 = %.2f; the general integral at s = 1e-9 gives %.2f; "
                "exact chain %.2f." % (Ne, 4 * Ne * math.log(2), lim, res[(Ne, s)][1]))
    u, t_ex, sd, t_dif = res[(100, 0.15)]
    asym = (math.log(4 * 100 * 0.15) + GAMMA) / 0.15
    say("Worked example (Ne = 100, s = 0.15, x0 = 0.5): P_fix = %.6f (exact chain); conditional mean time "
        "to fixation %.2f generations in the exact chain (SD %.2f), %.2f from the diffusion integral, and "
        "%.2f from the strong-selection approximation (ln(4 Ne s) + gamma)/s with 4 Ne s = 60."
        % (u, t_ex, sd, t_dif, asym))
    big = tbar1_diffusion(10_000, 0.15)
    say("Check of the approximation at Ne = 10,000 (4 Ne s = 6,000, s = 0.15): integral %.2f, "
        "(ln(4 Ne s) + gamma)/s = %.2f." % (big, (math.log(6000) + GAMMA) / 0.15))

    # 2 ---------------------------------------------------------------------------------------------
    say("")
    say("2. Branching process: establishment probability pi of one copy, offspring Poisson(1 + s)")
    say("%6s | %8s %6s %7s | %13s %14s" % ("s", "pi", "2s", "pi/2s", "u(1/2Ne) 100", "u(1/2Ne) 1000"))
    pis = {}
    for s in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
        pi = establishment(s)
        pis[s] = pi
        say("%6.2f | %8.5f %6.2f %7.3f | %13.5f %14.5f"
            % (s, pi, 2 * s, pi / (2 * s), kimura_u(1 / 200, 100, s), kimura_u(1 / 2000, 1000, s)))
    say("For s = 0.01, pi = %.4f against 2s = 0.02 and u(1/(2Ne)) = %.4f (Ne = 100) and %.4f (Ne = 1,000); "
        "2s overstates pi by %.0f%% at s = 0.1 and %.0f%% at s = 0.5."
        % (pis[0.01], kimura_u(1 / 200, 100, 0.01), kimura_u(1 / 2000, 1000, 0.01),
           100 * (0.2 / pis[0.1] - 1), 100 * (1.0 / pis[0.5] - 1)))
    pi = pis[0.05]
    say("k independent copies at s = 0.05 (pi = %.4f): establishment probability 1 - (1 - pi)^k" % pi)
    say("   " + "; ".join("k = %d: %.3f" % (k, 1 - (1 - pi) ** k) for k in (1, 5, 10, 20)))
    say("   copies needed for a 50%%, 90%% and 95%% chance that at least one establishes: %d, %d and %d"
        % tuple(copies_needed(pi, q) for q in (0.5, 0.9, 0.95)))

    # 3 ---------------------------------------------------------------------------------------------
    say("")
    say("3. Replicates needed by a simulation study")
    se_t = 0.5
    r_t = (sd / se_t) ** 2
    say("Worked example: the SD of the conditional fixation time is %.2f generations (exact chain), so a "
        "standard error of %.1f generation needs R = (SD / %.1f)^2 = %.1f, that is %d runs that fix; with "
        "P_fix = %.4f that is also the number of runs." % (sd, se_t, se_t, r_t, math.ceil(r_t), u))
    p, se_p = 0.73, 0.01
    r_p = round(p * (1 - p) / se_p ** 2, 6)
    say("A fixation probability of %.2f (small-family row) estimated with a standard error of %.2f needs "
        "R = p (1 - p) / %.2f^2 = %s runs." % (p, se_p, se_p, format(math.ceil(r_p), ",")))

    # 4 ---------------------------------------------------------------------------------------------
    say("")
    say("4. Monte Carlo check of the exact chain (monte_carlo_check.py)")
    out.extend(mc.report(mc.run(reps)))

    # 5 ---------------------------------------------------------------------------------------------
    say("")
    say("5. Further numbers quoted in Sections 3, 4 and 7")
    for ne, s_, x0 in ((10, 0.05, 0.1), (50, 0.05, 0.1)):
        P = transition(ne, s_)
        n = 2 * ne
        k = int(round(x0 * n))
        v = np.zeros(n + 1)
        v[k] = 1.0
        for _ in range(10):
            v = v @ P
        Q = P[1:-1, 1:-1]
        loss = np.linalg.solve(np.eye(n - 1) - Q, P[1:-1, 0])[k - 1]
        say("Ne = %d, s = %.2f, x0 = %.1f: allele lost within ten generations with probability %.4f, "
            "eventually %.4f" % (ne, s_, x0, v[0], loss))
    say("Strong-selection asymptote at Ne = 200, s = 0.10: (ln 4Nes + gamma)/s = %.1f generations"
        % ((math.log(4 * 200 * 0.10) + GAMMA) / 0.10))
    P = transition(50, 0.0)
    n = 100
    Q = P[1:-1, 1:-1]
    b1 = P[1:-1, -1]
    Nm = np.linalg.inv(np.eye(n - 1) - Q)
    i = n // 2 - 1
    u1 = (Nm @ b1)[i]
    m1 = (Nm @ Nm @ b1)[i] / u1
    m2 = ((2 * Nm @ Nm @ Nm - Nm @ Nm) @ b1)[i] / u1
    sd50 = math.sqrt(m2 - m1 ** 2)
    rfix = (sd50 / 0.5) ** 2
    say("Ne = 50, s = 0, x0 = 0.5: SD of the conditional fixation time %.2f; a standard error of 0.5 generation "
        "needs %s runs that fix, %s runs in all at P_fix = %.2f"
        % (sd50, format(math.ceil(rfix), ","), format(math.ceil(rfix / u1), ","), u1))
    r10 = 0.5 * (1 - math.exp(-2 * 0.10))
    say("Haldane map function: r = %.3f at 10 cM; a BC6 plant that carries the target still carries the donor "
        "allele at that locus with probability (1 - r)^6 = %.3f" % (r10, (1 - r10) ** 6))
    say("a1 at r = 0.0008 stays with sh2 through six backcrosses with probability (1 - r)^6 = %.3f"
        % ((1 - 0.0008) ** 6))

    text = "\n".join(out) + "\n"
    print(text, end="")
    with open(OUTFILE, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("written to %s" % OUTFILE)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else mc.REPS)
