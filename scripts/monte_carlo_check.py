# -*- coding: utf-8 -*-
"""Monte Carlo check of the exact Wright-Fisher values (fixation table and the worked example).

Each scenario is simulated forward as a diploid Wright-Fisher population of 2 Ne gene copies with
genic selection (post-selection frequency p(1 + s)/(1 + ps), then a binomial draw of 2 Ne copies),
starting from x0 = 0.5, until every replicate has fixed or been lost. The estimates (fixation
probability, and the mean and standard deviation of the number of generations to fixation among
the runs that fix) are compared with the exact chain of wf_exact.py, which is the target.
A Monte Carlo mean is consistent with the exact value when it lies within about two standard
errors of it.

One generator (numpy default_rng, seed 20260926) is shared by the scenarios in the order listed,
so a run with the same number of replicates reproduces every number.

Usage: python -I -B scripts/monte_carlo_check.py [replicates per scenario, default 100000]
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # python -I leaves the script folder off the path
from wf_exact import ROWS, transition  # noqa: E402

SEED = 20260926
REPS = 100_000
X0 = 0.5
MAXGEN = 200_000
SCENARIOS = [(name, Ne, s) for name, Ne, s in ROWS] + [("Worked example (Ne = 100, s = 0.15)", 100, 0.15)]


def exact_moments(Ne, s, x0=X0):
    """Exact P_fix and the mean and SD of the time to fixation among runs that fix.

    Transient states 1..2Ne-1 with block Q, A = I - Q, fundamental matrix N = A^-1, and u = N R the
    fixation probabilities (R: one-step probability of reaching 2Ne). The chain conditioned on
    fixation has fundamental matrix D^-1 N D with D = diag(u), so
      E[T | fix] = (N u)_i / u_i  and  E[T^2 | fix] = 2 (N N u)_i / u_i - E[T | fix].
    T counts generations up to and including the one in which fixation occurs."""
    n = 2 * Ne
    P = transition(Ne, s)
    tr = np.arange(1, n)
    A = np.eye(n - 1) - P[np.ix_(tr, tr)]
    u = np.linalg.solve(A, P[tr, n])
    Nu = np.linalg.solve(A, u)
    NNu = np.linalg.solve(A, Nu)
    i = int(round(x0 * n)) - 1
    mean = Nu[i] / u[i]
    second = 2 * NNu[i] / u[i] - mean
    return u[i], mean, math.sqrt(second - mean ** 2)


def simulate(Ne, s, reps, rng, x0=X0):
    """Generations to absorption of `reps` independent populations; returns the times of the runs
    that fixed and the fraction that fixed."""
    two_n = 2 * Ne
    i = np.full(reps, int(round(x0 * two_n)), dtype=np.int64)
    t = np.zeros(reps, dtype=np.int64)
    alive = np.ones(reps, dtype=bool)
    gen = 0
    while alive.any() and gen < MAXGEN:
        gen += 1
        p = i[alive] / two_n
        w = p * (1.0 + s) / (1.0 + p * s) if s != 0 else p
        i[alive] = rng.binomial(two_n, w)
        t[alive] += 1
        alive = (i > 0) & (i < two_n)
    fixed = i >= two_n
    return t[fixed], fixed.mean()


def run(reps=REPS, seed=SEED):
    rng = np.random.default_rng(seed)
    out = []
    for name, Ne, s in SCENARIOS:
        times, pfix = simulate(Ne, s, reps, rng)
        u, mean, sd = exact_moments(Ne, s)
        m = times.mean()
        se = times.std(ddof=1) / math.sqrt(len(times))
        out.append(dict(name=name, Ne=Ne, s=s, reps=reps, n_fixed=len(times), pfix_mc=pfix,
                        pfix_se=math.sqrt(pfix * (1 - pfix) / reps), pfix_exact=u, t_mc=m, t_se=se,
                        t_exact=mean, sd_mc=times.std(ddof=1), sd_exact=sd))
    return out


def report(rows, seed=SEED):
    lines = ["%d replicates per scenario, x0 = 0.5, seed %d (numpy default_rng, shared in the order below)"
             % (rows[0]["reps"], seed),
             "%-38s %4s %5s | %8s %8s | %7s %5s %7s %6s %5s | %6s %6s"
             % ("scenario", "Ne", "s", "Pfix_MC", "Pfix_ex", "T_MC", "SE", "T_ex", "diff", "z", "SD_MC", "SD_ex")]
    for r in rows:
        d = r["t_mc"] - r["t_exact"]
        lines.append("%-38s %4d %5.2f | %8.5f %8.5f | %7.2f %5.2f %7.2f %+6.2f %+5.1f | %6.2f %6.2f"
                     % (r["name"][:38], r["Ne"], r["s"], r["pfix_mc"], r["pfix_exact"], r["t_mc"], r["t_se"],
                        r["t_exact"], d, d / r["t_se"], r["sd_mc"], r["sd_exact"]))
    worst = max(rows, key=lambda r: abs(r["t_mc"] - r["t_exact"]))
    zmax = max(abs(r["t_mc"] - r["t_exact"]) / r["t_se"] for r in rows)
    pmax = max(abs(r["pfix_mc"] - r["pfix_exact"]) for r in rows)
    lines.append("Largest |Monte Carlo - exact| conditional mean time: %.2f generations (%s), at most %.1f "
                 "standard errors in any scenario; largest |Monte Carlo - exact| P_fix: %.4f."
                 % (abs(worst["t_mc"] - worst["t_exact"]), worst["name"], zmax, pmax))
    return lines


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else REPS
    print("\n".join(report(run(reps))))
