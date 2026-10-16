# -*- coding: utf-8 -*-
"""Exact Wright-Fisher values: the fixation table, the Wright-Fisher worked example, Kimura's
fixation probability at the values quoted in the text, and the a1-sh2 screening sizes.

Model: diploid Wright-Fisher population with 2*Ne gene copies and genic selection. The favoured
allele has relative fitness 1+s in the gamete pool, so its post-selection frequency is
p(1+s)/(1+ps), and the next generation is a binomial draw of 2*Ne copies. The chain is solved
exactly from its transition matrix; no random numbers are used.

For a start of i copies:
  P_fix(i)        = probability of absorption at 2*Ne (fundamental-matrix solution)
  T_fix(i)        = mean generations to fixation among runs that fix, (N u)_i / u_i,
                    where N = (I - Q)^-1 is the fundamental matrix of the transient states
  distribution at generation t = row i of the t-step transition matrix
Kimura's diffusion value P_fix(x0) = (1 - exp(-4 Ne s x0)) / (1 - exp(-4 Ne s)) is printed beside
the exact one. The diffusion conditional fixation times, the standard deviations of the exact
times, the branching-process examples and the replicate numbers are in worked_examples.py; the
Monte Carlo check of these values is monte_carlo_check.py.

The a1-sh2 screening size uses q = (r/2)(1 - r/2), the frequency of purple shrunken kernels on a
selfed A1 Sh2 / a1 sh2 ear (both gametes carry sh2 and at least one carries A1), and
n = ln 0.05 / ln(1 - q) kernels for a 95% chance of at least one.

Usage: python -I -B scripts/wf_exact.py
Requires numpy and scipy. Runs in a few seconds on one core.
"""
import math

import numpy as np
from scipy.stats import binom


def transition(Ne, s):
    n = 2 * Ne
    p = np.arange(n + 1) / n
    w = p * (1.0 + s) / (1.0 + p * s) if s else p
    return binom.pmf(np.arange(n + 1)[None, :], n, w[:, None])     # rows: from i, columns: to j


def exact(Ne, s, x0=0.5):
    n = 2 * Ne
    P = transition(Ne, s)
    tr = np.arange(1, n)                                             # transient states 1..n-1
    Q = P[np.ix_(tr, tr)]
    R_fix = P[tr, n]
    Nmat = np.linalg.inv(np.eye(n - 1) - Q)
    u = Nmat @ R_fix                                                 # P(fix) from each transient state
    t_fix = (Nmat @ u) / u
    i0 = int(round(x0 * n))
    return u[i0 - 1], t_fix[i0 - 1], P


def kimura(Ne, s, x0):
    S = 4 * Ne * s
    return x0 if S == 0 else (1 - math.exp(-S * x0)) / (1 - math.exp(-S))


def dist_at(P, i0, t):
    v = np.zeros(P.shape[0])
    v[i0] = 1.0
    for _ in range(t):
        v = v @ P
    return v


def quantile(v, q):
    n = len(v) - 1
    c = np.cumsum(v)
    return np.searchsorted(c, q) / n


ROWS = [  # scenario, Ne, s (start x0 = 0.5)
    ("Population under selection (Ne = 200)", 200, 0.10),
    ("Population under selection (Ne = 50)", 50, 0.10),
    ("Small population under selection", 30, 0.05),
    ("Population improvement (yield QTL)", 100, 0.02),
    ("Base population (drift)", 50, 0.00),
    ("Conservation line", 20, 0.00),
    ("Small family (inbreeding)", 10, 0.05),
]

if __name__ == "__main__":
    print("Fixation table (x0 = 0.5)")
    print("%-40s %4s %5s %6s | %8s %8s %7s | %7s" % ("scenario", "Ne", "s", "4Nes", "Pfix_dif", "Pfix_ex",
                                                     "|diff|", "Tfix_ex"))
    for name, Ne, s in ROWS:
        u, t, _ = exact(Ne, s)
        kd = kimura(Ne, s, 0.5)
        print("%-40s %4d %5.2f %6.1f | %8.5f %8.5f %7.4f | %7.2f" % (name, Ne, s, 4 * Ne * s, kd, u, abs(kd - u), t))
        if s == 0:
            print("%40s diffusion 4 Ne ln 2 = %.1f" % ("", 4 * Ne * math.log(2)))

    print("\nWorked example: Ne = 100, s = 0.15, x0 = 0.5")
    Ne, s = 100, 0.15
    u, t, P = exact(Ne, s)
    print("  P_fix exact %.6f, Kimura %.6f" % (u, kimura(Ne, s, 0.5)))
    print("  conditional mean time to fixation (exact chain) %.2f generations" % t)
    print("  strong-selection approximation (ln(4 Ne s) + 0.5772)/s = %.1f" % ((math.log(4 * Ne * s) + 0.5772) / s))
    i0 = Ne
    for tt in (2, 4, 6, 10):
        v = dist_at(P, i0, tt)
        x = np.arange(2 * Ne + 1) / (2 * Ne)
        mean = float(v @ x)
        sd = float(np.sqrt(v @ (x - mean) ** 2))
        print("  t = %2d: mean %.3f, sd %.3f, central 95%% %.3f-%.3f" % (tt, mean, sd, quantile(v, 0.025),
                                                                          quantile(v, 0.975)))
    x = 0.5
    det = [x]
    for _ in range(10):
        x = x + s * x * (1 - x)
        det.append(x)
    print("  deterministic recursion x(t+1) = x + s x (1-x): x6 = %.3f, x10 = %.3f" % (det[6], det[10]))
    print("  neutral SD approximation at t = 6: sqrt(0.25 (1 - exp(-6/200))) = %.4f" % math.sqrt(0.25 * (1 - math.exp(-6 / 200))))

    print("\nKimura P_fix from x0 = 0.5 (text and drift-selection figure c)")
    for S in (0.5, 1, 2, 8):
        print("  4Nes = %-4g P_fix = %.3f" % (S, (1 - math.exp(-S / 2)) / (1 - math.exp(-S))))
    print("Worked-example figure caption: P_fix at s = 0.15, x0 = 0.05: Ne = 10 -> %.3f, Ne = 50 -> %.3f" %
          (kimura(10, 0.15, 0.05), kimura(50, 0.15, 0.05)))

    print("\na1-sh2 screening sizes, q = (r/2)(1 - r/2)")
    for r in (0.0007, 0.0008, 0.0009):
        q = (r / 2) * (1 - r / 2)
        print("  r = %.4f: q = 1/%.0f, kernels for 95%% = %.0f, expected among 1,843 = %.2f"
              % (r, 1 / q, math.log(0.05) / math.log(1 - q), 1843 * q))
