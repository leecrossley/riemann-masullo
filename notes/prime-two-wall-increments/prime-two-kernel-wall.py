"""Check the three-wall identity and the odd Euler logarithm for the original pulse."""

import argparse
import importlib.util
from pathlib import Path

import mpmath as mp
import numpy as np


def load_wall():
    path = Path(__file__).with_name("prime-two-wall.py")
    spec = importlib.util.spec_from_file_location("prime_two_wall", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def potential(y, nodes, weights):
    small = y < 35
    values = np.empty_like(y)
    values[small] = -0.5 * (np.expm1(-y[small, None] * nodes) @ (weights / nodes))
    values[~small] = 0.5 * (np.euler_gamma + np.log(y[~small]))
    return values


def finite_eta(points):
    two = np.exp(-np.log(2) * points)
    three = np.exp(-np.log(3) * points)
    return 1 - two + three, np.log(2) * two - np.log(3) * three


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=1536)
    parser.add_argument("--steps", type=int, default=32)
    parser.add_argument("--order", type=int, default=192)
    args = parser.parse_args()
    if min(args.height, args.steps, args.order) < 2 or args.steps % 2:
        parser.error("Use positive quadrature sizes of at least two and an even step count.")

    beta = 0.75
    heat = 64.0
    half_width = 1 / (2 * np.log(heat))
    side_walls = np.array([beta - half_width, beta, beta + half_width])
    wall = load_wall()
    pulse = wall.load_pulse()
    _, v, signed = pulse.probe_rule(args.order)
    p = 2 * signed * np.cosh((beta - 0.5) * v)
    b = p * v / np.sinh(half_width * v)
    c = p * v / np.tanh(half_width * v)
    pole = np.dot(p, np.sinh((1 - beta) * v) / v)

    quad_nodes, quad_weights = np.polynomial.legendre.leggauss(32)
    unit_nodes, unit_weights = (1 + quad_nodes) / 2, quad_weights / 2
    test_y = np.array([0, 1e-10, 1, 12, 34.999, 35, 64])
    reference = [0 if y == 0 else float((mp.euler + mp.log(y) + mp.e1(y)) / 2) for y in test_y]
    np.testing.assert_allclose(potential(test_y, unit_nodes, unit_weights), reference, atol=2e-14, rtol=0)

    nodes, weights = np.polynomial.legendre.leggauss(20)
    length = half_width * (1 + nodes) / 2
    area_s = np.r_[beta - length, beta + length]
    area_weights = half_width * np.tile(weights, 2) / 2
    profile = (p * np.sinh((half_width - abs(area_s[:, None] - beta)) * v)
               / np.sinh(half_width * v))
    filters = np.vstack([p, b, c, profile])

    n = max(128, int(np.ceil(args.height / 2)))
    logs = np.log(np.arange(1, n, dtype=float))
    n_weights = np.exp(-side_walls[:, None] * logs)
    bernoulli = [float(mp.bernoulli(2 * k) / mp.factorial(2 * k)) for k in range(1, 17)]
    heights = np.arange(args.height * args.steps + 1) / args.steps
    samples = np.empty((4, len(heights)))
    for first in range(0, len(heights), 192):
        batch = slice(first, first + 192)
        t = heights[batch]
        transformed = filters @ np.cos(v[:, None] * t)
        original, outer, centre = transformed[:3]
        pts = side_walls[:, None] + 1j * t
        etas, _ = finite_eta(pts)
        heat_logs = potential(heat * np.abs(etas) ** 2, unit_nodes, unit_weights)
        samples[0, batch] = (outer * (heat_logs[0] + heat_logs[2]) - 2 * centre * heat_logs[1]) / np.pi

        eta, deriv = finite_eta(area_s[:, None] + 1j * t)
        density = abs(deriv) ** 2 * np.exp(-heat * abs(eta) ** 2)
        samples[1, batch] = 2 * heat / np.pi * (area_weights @ (transformed[3:] * density))

        odd_logs = np.empty_like(pts.real)
        for i, w in enumerate(pts):
            zeta, _ = wall.zeta_and_derivative(w, logs, n_weights[i], n, bernoulli)
            odd_logs[i] = np.log(abs((1 - np.exp(-np.log(2) * w)) * zeta))
        samples[2, batch] = original * odd_logs[1] / np.pi
        samples[3, batch] = (outer * (odd_logs[0] + odd_logs[2]) - 2 * centre * odd_logs[1]) / np.pi

    print(f"beta={beta:g}  X={heat:g}  h={2 * half_width:.10f}  expected pole contribution={-pole:.12f}")
    for stride in (2, 1):
        three, area, odd, combined = np.trapz(samples[:, ::stride], dx=stride / args.steps, axis=1)
        print(f"T={args.height} dt={stride / args.steps:g}  K=3 walls={three:.10f}  area={area:.10f}"
              f"  difference={three-area:.2g}")
        print(f"  odd logarithm={odd:.12f}  pole difference={odd+pole:.2g}  three-wall odd logarithm={combined:.2g}")
    if abs(three - area) > 4e-6 or abs(odd + pole) > 4e-6 or abs(combined) > 4e-6:
        raise ArithmeticError("The quadrature has not resolved the identities to the stated finite precision.")


if __name__ == "__main__":
    main()
