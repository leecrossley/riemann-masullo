"""Check the odd-prime phase average for the original pulse at finite Dirichlet cutoff."""

import argparse
import importlib.util
import math
from pathlib import Path

import mpmath as mp
import numpy as np


def coefficient_grid(size, heat):
    path = Path(__file__).with_name("prime-two-increments.py")
    spec = importlib.util.spec_from_file_location("prime_two_increments", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.PrimeTwoGrid(size, heat)


def averaged_square(size, heat):
    grid = coefficient_grid(size, heat)
    denominator = [2**max(-i, 0) * 3**max(-j, 0) for i, j in zip(grid.i.tolist(), grid.j.tolist())]
    retained = np.array([q > 8 for q in denominator])
    powers = grid.j[retained]
    matches = powers[:, None] == powers[None, :]
    phases = np.exp(2j * np.pi * np.outer(np.arange(size) / size, powers))
    nodes, weights = np.polynomial.legendre.leggauss(18)
    h = 1 / np.log(heat)
    total = matching = difference = 0.0
    error = 0.0
    for node, weight in zip(nodes, weights):
        s = 0.75 + 0.025 * node
        increment = (grid.coefficients(s + h / 2) - grid.coefficients(s - h / 2))[retained]
        pairs = increment[:, None] * np.conj(increment[None, :])
        mean = np.mean(np.abs(phases @ increment)**2)
        error = max(error, abs(pairs[matches].sum() - mean))
        mass = 0.025 * weight * (1 - node**2)**2 / h
        total += mass * abs(increment.sum())**2
        matching += mass * mean
        difference += mass * pairs[~matches].sum().real
    if error > 2e-12 or abs(total - matching - difference) > 2e-12:
        raise ArithmeticError("The phase average and the direct pair sum differ.")
    return total, matching, difference, error


def linear_means():
    mp.mp.dps = 45
    s = mp.mpf(3) / 4
    prime = bytearray([1]) * 10001
    prime[:2] = bytes(2)
    for p in range(2, math.isqrt(len(prime) - 1) + 1):
        if prime[p]:
            count = (len(prime) - 1 - p*p) // p + 1
            prime[p*p::p] = bytes(count)
    for y in (10, 100, 1000, 10000):
        sifted = mp.altzeta(s) * mp.fprod(1 - mp.power(p, -s) for p in range(3, y + 1) if prime[p])
        print(f"Y={y:5d}  K first {float(sifted): .10f}")
    dyadic = (1 - mp.power(2, 1-s)) / (1 - mp.power(2, -s))
    print(f"All odd phases first, then K {float(dyadic): .10f}")


def second_moment():
    s = 0.75
    phases = np.exp(2j * np.pi * np.arange(96) / 96)
    x, y = 3**(-s), 5**(-s)
    eta = float(mp.altzeta(s))
    twisted = eta * (1-x) * (1-y) / ((1-x*phases[:, None]) * (1-y*phases[None, :]))
    exact = eta**2 * (1-x)/(1+x) * (1-y)/(1+y)
    error = abs(np.mean(np.abs(twisted)**2) - exact)
    if error > 1e-14:
        raise ArithmeticError("The averaged second moment differs from the Euler product.")
    print(f"Odd primes 3 and 5, second moment {exact:.10f}, quadrature error {error:.2e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=int, default=192)
    args = parser.parse_args()
    if args.grid < 96 or args.grid % 4:
        parser.error("Use a multiple of four of at least 96.")
    print("K=4, q>8, h=1/log X, wall weight (1-((s-0.75)/0.025)^2)^2 on [0.725,0.775]")
    print(" X        full / h       matching / h       differing / h       phase error")
    for heat in (16, 32, 64):
        total, matching, difference, error = averaged_square(args.grid, heat)
        print(f"{heat:2d}  {total:16.9e}  {matching:16.9e}  {difference:16.9e}  {error:.2e}")
    linear_means()
    second_moment()


if __name__ == "__main__":
    main()
