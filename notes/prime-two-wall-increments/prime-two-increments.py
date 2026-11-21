"""Check the coefficient derivative and the squared wall increment at Dirichlet cutoff four."""

import argparse
import importlib.util
import math
from pathlib import Path

import numpy as np


def load_pulse():
    path = Path(__file__).with_name("prime-two-check.py")
    spec = importlib.util.spec_from_file_location("prime_two_pulse", path)
    pulse = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pulse)
    return pulse


class PrimeTwoGrid:
    def __init__(self, size, heat):
        self.heat = heat
        angles = 2 * np.pi * np.arange(size) / size
        self.two = np.exp(1j * angles)[:, None]
        self.three = np.exp(1j * angles)[None, :]
        pulse = load_pulse()
        c, _, _ = pulse.probe_rule(96)
        frequencies = np.arange(-size // 4, size // 4 + 1)
        i, j = np.meshgrid(frequencies, frequencies, indexing="ij")
        logs = i * np.log(2) + j * np.log(3)
        absolute = np.abs(logs)
        support = ((np.abs(absolute - pulse.CENTRE) <= pulse.SUPPORT)
                   | (np.abs(absolute - pulse.CENTRE - pulse.SHIFT) <= pulse.SUPPORT))
        retain = support & ((i <= -2) | (j < 0))
        self.i, self.j, self.logs = i[retain], j[retain], logs[retain]
        offsets = np.abs(self.logs) - pulse.CENTRE
        self.weights = pulse.bump(offsets, 96) - c * pulse.bump(offsets - pulse.SHIFT, 96)

    def coefficients(self, sigma, density=False):
        two = 2 ** (-sigma) * self.two
        three = 3 ** (-sigma) * self.three
        four = two**2
        eta = 1 - two + three - four
        derivative = np.log(2) * two - np.log(3) * three + np.log(4) * four
        y = self.heat * np.abs(eta)**2
        if density:
            values = 2 * self.heat * np.abs(derivative)**2 * np.exp(-y)
        else:
            regulariser = np.divide(-np.expm1(-y), eta, out=np.zeros_like(eta), where=eta != 0)
            values = derivative * regulariser
        coeff = np.fft.fft2(values, norm="forward")[self.i, self.j]
        return self.weights * np.exp((sigma - 0.5) * self.logs) * coeff


def check(size, heat):
    grid = PrimeTwoGrid(size, heat)
    nodes, weights = np.polynomial.legendre.leggauss(20)
    step = 0.04
    middle = 0.75
    endpoints = grid.coefficients(middle + step / 2) - grid.coefficients(middle - step / 2)
    integrated = (step / 2) * sum(
        weight * grid.coefficients(middle + step * node / 2, density=True)
        for node, weight in zip(nodes, weights)
    )
    derivative_error = np.max(np.abs(endpoints - integrated))
    total = diagonal = 0.0
    for node, weight in zip(nodes, weights):
        centre = middle + 0.10 * node
        increment = grid.coefficients(centre + step / 2) - grid.coefficients(centre - step / 2)
        wall_weight = 0.10 * weight * (1 - node**2)**2 / step
        total += wall_weight * abs(increment.sum())**2
        diagonal += wall_weight * np.vdot(increment, increment).real
    if derivative_error > 1e-10:
        raise ArithmeticError(f"Coefficient discrepancy {derivative_error:.3g}")
    return derivative_error, total, diagonal


def arithmetic_check(size):
    grid = PrimeTwoGrid(size, 0.25)
    centre, step = 0.75, 0.04
    ratios = [
        (2**max(i, 0) * 3**max(j, 0), 2**max(-i, 0) * 3**max(-j, 0))
        for i, j in zip(grid.i.tolist(), grid.j.tolist())
    ]
    sums = np.zeros(len(ratios))
    coefficients = {1: 1}
    for j in range(23):
        product = {}
        for n, coefficient in coefficients.items():
            for k, sign in ((1, 1), (2, -1), (3, 1), (4, -1)):
                product[k * n] = product.get(k * n, 0) + sign * coefficient
        coefficients = product
        scale = (-grid.heat)**j / (math.factorial(j) * (j + 1)**2)
        for index, (numerator, denominator) in enumerate(ratios):
            for m, at_m in coefficients.items():
                if m % denominator:
                    continue
                n = (m // denominator) * numerator
                at_n = coefficients.get(n, 0)
                if at_n:
                    sums[index] += (scale * at_n * at_m * math.log(n) * math.sinh(step * math.log(m))
                                    * m**(-2 * centre))
    arithmetic = 2 * grid.heat * grid.weights * np.exp(-grid.logs / 2) * sums
    endpoints = grid.coefficients(centre + step / 2) - grid.coefficients(centre - step / 2)
    error = np.max(np.abs(arithmetic - endpoints))
    if error > 1e-12:
        raise ArithmeticError(f"Arithmetic increment discrepancy {error:.3g}")
    return error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=int, default=96)
    args = parser.parse_args()
    if args.grid < 64 or args.grid % 4:
        parser.error("Use a multiple of four of at least 64.")
    print("K=4, q>2, h=0.04, wall weight (1-((s-0.75)/0.10)^2)^2 on [0.65,0.85]")
    print("X   coefficient error   full energy/h    equal ratios     distinct ratios")
    for heat in (1.0, 4.0, 12.0):
        error, total, diagonal = check(args.grid, heat)
        print(f"{heat:2.0f}  {error:16.3e}  {total:15.8e}  {diagonal:15.8e}  {total-diagonal:15.8e}")
    print(f"X=0.25, arithmetic convolution through j=22, maximum error {arithmetic_check(args.grid):.3e}")


if __name__ == "__main__":
    main()
