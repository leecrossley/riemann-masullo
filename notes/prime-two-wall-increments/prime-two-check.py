"""Numerical checks for the prime 2 test. Requires NumPy and mpmath."""

import argparse
from decimal import Decimal, localcontext
from functools import lru_cache
from math import ceil, comb, factorial, log

import mpmath as mp
import numpy as np


HALF_BOX = 1 / 32
STEP = 1 / 512
MOLLIFIER_RADIUS = 1 / 256
BOXES = 8
SUPPORT = HALF_BOX + BOXES * STEP + MOLLIFIER_RADIUS
CENTRE = log(2)
SHIFT = 1 / 8


@lru_cache
def nodes(order):
    return np.polynomial.legendre.leggauss(order)


def mollifier_rule(order):
    points, weights = nodes(order)
    weights = weights * np.exp(-1 / (1 - points**2))
    return MOLLIFIER_RADIUS * points, weights / weights.sum()


def box_sum_cdf(y):
    x = (y / STEP + BOXES) / 2
    left = np.clip(np.minimum(x, BOXES - x), 0, BOXES / 2)
    polynomial = np.zeros_like(left)
    for k in range(BOXES // 2):
        polynomial += (-1)**k * comb(BOXES, k) * np.maximum(left - k, 0)**BOXES
    lower = polynomial / factorial(BOXES)
    return np.where(x <= BOXES / 2, lower, 1 - lower)


def bump(y, order):
    v, mass = mollifier_rule(order)
    cdf = box_sum_cdf(HALF_BOX - np.abs(y)[..., None] + v)
    return cdf @ mass


def real_weight(u):
    return np.exp(u / 2) * (1 - 1 / (np.exp(3 * u) - np.exp(u)))


def probe_rule(order):
    points, weights = nodes(order)
    v = SUPPORT * points
    weights = SUPPORT * weights * bump(v, order)
    left, right = CENTRE + v, CENTRE + SHIFT + v
    c = np.dot(weights, real_weight(left)) / np.dot(weights, real_weight(right))
    return c, np.r_[left, right], np.r_[weights, -c * weights]


def transform_on_line(gamma, c, order):
    t = np.asarray(gamma, dtype=float)
    v, mass = mollifier_rule(order)
    eta = np.cos(t[..., None] * v) @ mass
    base = (2 * HALF_BOX * np.sinc(HALF_BOX * t / np.pi)
            * np.sinc(STEP * t / np.pi)**BOXES * eta)
    return base * (np.exp(1j * CENTRE * t) - c * np.exp(1j * (CENTRE + SHIFT) * t))


def phi(x):
    r = np.pi * x**2
    return (4 * r**2 - 6 * r) * np.exp(-r)


def global_theta(x):
    large = np.maximum(x, 1 / x)
    value = sum(phi(n * large) for n in range(1, 7))
    return np.where(x >= 1, value, value / x)


def finite_theta(x, excluded):
    terms = [(1, 1)]
    for p in excluded:
        terms += [(d * p, -sign) for d, sign in terms]
    return sum(sign * global_theta(d * x) for d, sign in terms)


def dyadic_theta(x):
    value = np.zeros_like(x)
    while np.min(x) <= 6:
        value += phi(x)
        x = 2 * x
    return value


def global_generator(u):
    x = np.exp(np.abs(u))
    return np.sqrt(x) * sum(phi(n * x) for n in range(1, 6))


def shell_cutoff(v):
    t = np.clip(8 * (4 - np.exp(np.abs(v))), 0, 1)
    left, right = np.zeros_like(t), np.zeros_like(t)
    left[t > 0] = np.exp(-1 / t[t > 0])
    right[t < 1] = np.exp(-1 / (1 - t[t < 1]))
    return left / (left + right)


def global_projection_trace(alpha, shifts, u, weights, order):
    """Trace on a finite span of translates of the global Gaussian theta generator."""
    limit = np.max(np.abs(shifts)) + 3
    walls = [log(4), log(4 - 1 / 8)]
    edges = np.unique(np.r_[np.linspace(-limit, limit, ceil(8 * limit) + 1), walls, -np.array(walls)])
    half = np.diff(edges) / 2
    points, mass = nodes(order)
    v = ((edges[:-1] + half)[:, None] + half[:, None] * points).ravel()
    dv = (half[:, None] * mass).ravel()
    metric = np.sqrt(dv) * np.exp(alpha * np.abs(v)) * (1 + v**2)
    translates = global_generator(v[:, None] - shifts)
    q, r = np.linalg.qr(metric[:, None] * translates)
    keep = np.abs(v) < log(4)
    cut_weight = metric[keep] * shell_cutoff(v[keep])
    images = np.column_stack([
        cut_weight * (global_generator(v[keep, None] - shift - u) @ weights)
        for shift in shifts
    ])
    return np.trace(np.linalg.solve(r, q[keep].T @ images))


def cutoff_coefficient(log_cutoff, gamma, u, weights, theta, panel_order=16):
    points, mass = nodes(panel_order)
    edges = np.linspace(-log_cutoff, log_cutoff, ceil(2 * log_cutoff / 0.4) + 1)
    total = 0j
    for lower, upper in zip(edges[:-1], edges[1:]):
        half = (upper - lower) / 2
        v = (lower + upper) / 2 + half * points
        x = np.exp(v[:, None] - u)
        values = (np.sqrt(x) * theta(x)) @ weights
        total += half * np.dot(mass, np.exp(1j * gamma * v) * values)
    return total


def tail_bound():
    with localcontext() as ctx:
        ctx.prec = 60
        H = Decimal(3) * 10**12
        r = Decimal(1) / 128
        return Decimal(512)**8 / H**7 * (H.ln() / (1 - r) + Decimal(2).ln() / (1 - r)**2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zeros", help="Optional LMFDB text table with index and positive ordinate.")
    args = parser.parse_args()
    order = 96
    c, u, weights = probe_rule(order)
    c_again = probe_rule(160)[0]
    assert abs(c - c_again) < 3e-13
    assert abs(bump(np.array([0.0, SHIFT]), order) @ np.array([1., -c]) - 1) < 1e-13
    assert abs(np.dot(weights, real_weight(u))) < 1e-14

    mp.mp.dps = 48
    zero = mp.zetazero(1)
    gamma = float(zero.imag)
    mf = transform_on_line(gamma, c, order)
    assert abs(mf - np.dot(weights, np.exp(1j * gamma * u))) < 2e-13
    local_factor = zero * (zero - 1) * mp.pi**(-zero / 2) * mp.gamma(zero / 2) / 2
    local = complex(local_factor / (1 - 2**(-zero))) * mf
    target = log(2) / np.sqrt(2)
    print(f"Second bump coefficient: {c:.14f}")
    print(f"Local geometric trace target: {target:.15f}")
    print(f"Dyadic absolute tail majorant: {tail_bound():.10E}")
    print(f"Mellin transform at the first zero: {mf:.12g}")
    print(f"Limiting dyadic Mellin coefficient: {local:.12g}")

    for alpha in [1 / 8, 1 / 4]:
        for count in [9, 17, 33]:
            shifts = np.linspace(-4, 4, count)
            trace = global_projection_trace(alpha, shifts, u, weights, 12)
            repeat = global_projection_trace(alpha, shifts, u, weights, 18)
            assert abs(trace - repeat) < 5e-9
            print(f"Finite projection at half-width {alpha:g}, {count} translates: {trace:.9f}")

    fixed_primes = lambda x: finite_theta(x, (3, 5, 7))
    for scale in [4, 64, 4096, 2**36]:
        finite = cutoff_coefficient(log(scale), gamma, u, weights, fixed_primes)
        dyadic = cutoff_coefficient(log(scale), gamma, u, weights, dyadic_theta)
        print(f"Cutoff {scale:>11}: finite {finite:.9g}, dyadic {dyadic:.9g}")
    assert abs(finite) < 3e-12
    assert abs(dyadic - local) < 2e-7

    if args.zeros:
        table = np.loadtxt(args.zeros)
        assert np.array_equal(table[:, 0], np.arange(1, len(table) + 1))
        assert np.all(np.diff(table[:, 1]) > 0) and abs(table[0, 1] - gamma) < 1e-12
        values = transform_on_line(table[:, 1], c, order)
        geometric = -2 * values.real.sum()
        print(f"Truncated geometric zero sum at {len(table)} positive ordinates: {geometric:.15f}")
        print(f"Difference from the local target: {geometric - target:.3e}")


if __name__ == "__main__":
    main()
