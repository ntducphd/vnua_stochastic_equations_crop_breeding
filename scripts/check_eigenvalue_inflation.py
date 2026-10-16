# -*- coding: utf-8 -*-
"""Inflation of the largest eigenvalue of an estimated covariance matrix (the multivariate-response
row of the quantitative-genetics table).

For q uncorrelated traits with unit variance, the sample covariance matrix estimated from n
observations has eigenvalues that fill [(1 - sqrt(q/n))^2, (1 + sqrt(q/n))^2] as q and n grow
with q/n fixed (Marchenko-Pastur), and the largest eigenvalue tends to the upper edge
(1 + sqrt(q/n))^2 although every true eigenvalue is 1. The script simulates the null case directly
(200 replicate matrices per (q, n), seed 20260926) and prints the mean largest eigenvalue beside
that limit, to show how close breeding-sized q and n come to it.

Usage: python -I -B scripts/check_eigenvalue_inflation.py
"""
import numpy as np

rng = np.random.default_rng(20260926)
REPS = 200

print("%6s %6s %7s %14s %18s %8s" % ("q", "n", "q/n", "simulated max", "(1 + sqrt(q/n))^2", "ratio"))
for q, n in ((10, 100), (20, 100), (50, 200), (5, 50), (40, 400), (100, 1000)):
    g = q / n
    tops = []
    for _ in range(REPS):
        X = rng.standard_normal((n, q))
        S = X.T @ X / n
        tops.append(np.linalg.eigvalsh(S)[-1])
    sim = float(np.mean(tops))
    limit = (1 + np.sqrt(g)) ** 2
    print("%6d %6d %7.3f %14.3f %18.3f %8.3f" % (q, n, g, sim, limit, sim / limit))

print()
print("Every true eigenvalue is 1. The mean largest sample eigenvalue approaches (1 + sqrt(q/n))^2 as q")
print("and n grow with q/n fixed, and stays below it for the few traits of a breeding index (q = 5-10).")
