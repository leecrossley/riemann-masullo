"""Check the endpoint wall and the heat profile at its first dyadic zero."""

import argparse
from functools import partial
from math import erf
from pathlib import Path
from runpy import run_path

import mpmath as mp
import numpy as np


def eta_on_wall(s, n, bernoulli, zeta_and_derivative):
    a = np.log(2.0)
    log_n = np.log(np.arange(1, n, dtype=float))
    regular = s != 1
    eta = np.empty_like(s)
    derivative = np.empty_like(s)
    eta[~regular] = a
    derivative[~regular] = a * float(mp.euler) - a * a / 2
    if np.any(regular):
        sr = s[regular]
        zeta, dzeta = zeta_and_derivative(sr, log_n, np.exp(-sr[0].real * log_n), n, bernoulli)
        exponent = a * (1 - sr)
        factor = -np.expm1(exponent)
        eta[regular] = factor * zeta
        derivative[regular] = a * np.exp(exponent) * zeta + factor * dzeta
    return eta, derivative


def endpoint_integrals(periods, steps, transform, evaluate, target):
    period = 2 * np.pi / np.log(2)
    spacing = period / steps
    heights = spacing * np.arange(periods * steps + steps // 2 + 1)
    n = max(96, int(np.ceil(heights[-1] / 2)))
    heat = np.array([16.0, 64.0, 256.0, 1024.0, 4096.0])
    integrand = np.empty((len(heights), len(heat)))
    for start in range(0, len(heights), 192):
        batch = slice(start, start + 192)
        s = 1 + 1j * heights[batch]
        eta, derivative = evaluate(s, n)
        norm2 = np.abs(eta[:, None])**2
        quotient = np.broadcast_to(heat, (len(s), len(heat))).copy()
        np.divide(-np.expm1(-norm2 * heat), norm2, out=quotient, where=norm2 > 0)
        integrand[batch] = (transform(s - 0.5) * derivative * eta.conj()).real[:, None] * quotient / np.pi
    print(f"Endpoint vertical integral, limiting value {target:.12f}, height {heights[-1]:.3f}")
    print("    X       coarser mesh         finer mesh")
    coarse = np.trapz(integrand[::2], dx=2 * spacing, axis=0)
    fine = np.trapz(integrand, dx=spacing, axis=0)
    for x, first, second in zip(heat, coarse, fine):
        print(f"{x:5g} {first:18.12f} {second:18.12f}")


def first_dyadic_zero(transform, evaluate):
    a = np.log(2)
    rho = 1 + 2j * np.pi / a
    slope = a * abs(complex(mp.zeta(rho)))
    weight = transform(np.array([rho - 0.5]))[0]
    horizontal, wx = np.polynomial.legendre.leggauss(32)
    vertical, wy = np.polynomial.legendre.leggauss(80)
    half_width, half_height = 0.75, 7
    expected = erf(half_width) * erf(half_height)
    print(f"First dyadic zero, scaled rectangle, limiting fraction {expected:.12f}")
    print("    X           real part     imaginary part")
    for x in (256, 4096, 65536):
        scale = 1 / (slope * np.sqrt(x))
        mass = 0j
        for offset, row_weight in zip(horizontal, wx):
            s = rho + scale * (half_width * offset + 1j * half_height * vertical)
            eta, derivative = evaluate(s, 64)
            density = transform(s - 0.5) * np.abs(derivative)**2 * np.exp(-x * np.abs(eta)**2)
            mass += row_weight * np.dot(wy, density)
        fraction = mass * (half_width * half_height / (np.pi * slope**2)) / weight
        print(f"{x:5d} {fraction.real:19.12f} {fraction.imag:18.12f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--periods", type=int, default=192)
    parser.add_argument("--steps-per-period", type=int, default=1024)
    args = parser.parse_args()
    if args.periods < 1 or args.steps_per_period < 4 or args.steps_per_period % 4:
        parser.error("Use a positive number of periods and a positive multiple of four steps per period.")
    wall = run_path(str(Path(__file__).with_name("prime-two-wall.py")))
    pulse = wall["load_pulse"]()
    c, u, weights = pulse.probe_rule(128)
    offsets, masses = pulse.mollifier_rule(64)
    transform = partial(wall["pulse_transform"], pulse=pulse, c=c, offsets=offsets, masses=masses)
    bernoulli = [float(mp.bernpoly(2 * k, 0) / mp.factorial(2 * k)) for k in range(1, 17)]
    evaluate = partial(eta_on_wall, bernoulli=bernoulli, zeta_and_derivative=wall["zeta_and_derivative"])
    local = np.log(2) / np.sqrt(2)
    pole = weights @ (2 * np.cosh(u / 2))
    print(f"Predicted endpoint tail limits   left {pole:.12f}   midpoint {pole / 2:.12f}   right {0.0:.12f}")
    endpoint_integrals(args.periods, args.steps_per_period, transform, evaluate, (pole - local) / 2)
    first_dyadic_zero(transform, evaluate)


if __name__ == "__main__":
    main()
