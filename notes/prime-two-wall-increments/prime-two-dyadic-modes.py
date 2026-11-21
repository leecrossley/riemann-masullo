"""Check the dyadic block, the twisted Dirichlet polynomial and the two wall terms."""

import argparse
import importlib.util
from pathlib import Path

import mpmath as mp
import numpy as np


def load_walls():
    path = Path(__file__).with_name("prime-two-wall.py")
    spec = importlib.util.spec_from_file_location("prime_two_wall", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def valuations(indices, prime):
    powers = np.zeros_like(indices)
    factor = prime
    while factor <= indices[-1]:
        powers[factor - 1::factor] += 1
        factor *= prime
    return powers


def check_blocks():
    n = np.arange(1, 4097)
    parity = np.where(n % 2, 1, -1)
    powers = valuations(n, 2)
    largest = 0.0
    for angle in (0, 0.03, 0.8, np.pi, -2.2):
        z = np.exp(1j * angle)
        c = (1 - z) / (2 - z)
        sums = np.cumsum(parity * np.exp(1j * angle * powers))
        for m in range(1, 13):
            expected = c * (2 ** m - z ** m)
            largest = max(largest, abs(sums[2 ** m - 1] - expected))
    if largest > 3e-10:
        raise ArithmeticError("The direct block sums disagree with the dyadic formula.")
    print(f"Dyadic blocks through 4096, largest absolute error {largest:.2g}")


def check_dirichlet():
    mp.mp.dps = 32
    primes = {3: np.exp(0.4j), 5: np.exp(-0.6j)}
    n = np.arange(1, 2 ** 18 + 1)
    log_n = np.log(n)
    fixed = np.where(n % 2, 1, -1).astype(complex)
    for p, u in primes.items():
        fixed *= u ** valuations(n, p)
    powers = valuations(n, 2)

    def multiplier(w):
        return mp.fprod((1 - p ** (-w)) / (1 - u * p ** (-w)) for p, u in primes.items())

    residuals = []
    for length in (2 ** 12, 2 ** 15, 2 ** 18):
        log_k = np.log(length)
        errors = []
        for angle in (0, 0.7, length ** (-0.25)):
            z = mp.exp(mp.j * angle)
            c = (1 - z) / (2 - z)
            for w in (0.68 + 1.75j, 0.82 + 6j):
                terms = (fixed[:length] * np.exp(1j * angle * powers[:length])
                         * np.exp(-w * log_n[:length]))
                actual, derivative = np.sum(terms), np.sum(-log_n[:length] * terms)

                def continuation(v):
                    q = 2 ** (-v)
                    return multiplier(v) * (1 - 2 * z * q) * (1 - q) * mp.zeta(v) / (1 - z * q)

                main = c * multiplier(1) * mp.exp((1 - w) * log_k) / (1 - w)
                predicted = main + continuation(mp.mpc(w))
                predicted_derivative = main * (-log_k + 1 / (1 - w)) + mp.diff(continuation, w)
                errors.extend((abs(actual - complex(predicted)),
                               abs(derivative - complex(predicted_derivative)) / log_k))
        residuals.append(max(errors))
        print(f"K={length:6d}  maximum remainder in E or E'/log K  {residuals[-1]:.5g}")
    if residuals[-1] >= 0.005 or residuals[-1] >= residuals[0] / 3:
        raise ArithmeticError("The finite Dirichlet values disagree with the expansion.")


def check_walls(height, steps):
    wall = load_walls()
    pulse = wall.load_pulse()
    coefficient, nodes, weights = pulse.probe_rule(128)
    offsets, masses = pulse.mollifier_rule(64)
    pole_constant = np.dot(weights, np.exp(-nodes / 2))
    plateau = np.dot(weights, np.exp(nodes / 2)) - np.sqrt(2) * np.log(2)
    print(f"C={pole_constant:.12f}  (F(1/2) - sqrt(2) log 2)^2={plateau ** 2:.12f}")

    heat, centre, length = 64.0, 0.75, 512
    spacing = 1 / np.log(heat)
    sigmas = np.array([centre + spacing / 2, centre - spacing / 2])
    n = np.arange(1, length + 1)
    log_n = np.log(n)
    weights_n = np.where(n % 2, 1, -1) * np.exp(-sigmas[:, None] * log_n)
    ts = np.arange(height * steps + 1) / steps
    data = np.empty((2, 2, len(ts)))
    for j in range(0, len(ts), 128):
        batch = slice(j, j + 128)
        times = ts[batch]
        waves = np.exp(-1j * times[:, None] * log_n)
        for index, sigma in enumerate(sigmas):
            w = sigma + 1j * times
            transform = wall.pulse_transform(w - 0.5, pulse, coefficient, offsets, masses)
            value = np.einsum("ij,j->i", waves, weights_n[index])
            deriv = np.einsum("ij,j->i", waves, -log_n * weights_n[index])
            regularised = np.divide(-np.expm1(-heat * abs(value) ** 2), value,
                                    out=np.zeros_like(value), where=value != 0)
            data[index, 0, batch] = (transform / (1 - w)).real / np.pi
            data[index, 1, batch] = (transform * deriv * regularised).real / np.pi

    for stride in (2, 1):
        values = np.trapz(data[:, :, ::stride], dx=stride / steps, axis=2)
        pole_error = max(abs(values[:, 0] - pole_constant))
        q_right, q_left = values[:, 1] - pole_constant
        positive = q_right ** 2
        mixed = q_left ** 2 - 2 * q_right * q_left
        difference = (values[0, 1] - values[1, 1]) ** 2
        np.testing.assert_allclose(positive + mixed, difference, rtol=1e-12, atol=1e-14)
        print(f"Whole finite walls, K=512 X=64 u=1 dt={stride / steps:g} T={height}")
        print(f"  C integral error {pole_error:.2g}  right square {positive:.9g}"
              f"  left and mixed {mixed:.9g}  sum {difference:.9g}")
    if pole_error > 2e-7:
        raise ArithmeticError("The wall integral disagrees with the pulse moment.")
    if positive <= 0 or mixed >= 0:
        raise ArithmeticError("The two sampled wall terms do not have the reported signs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=1536)
    parser.add_argument("--steps", type=int, default=32)
    args = parser.parse_args()
    if args.height < 1536 or args.steps < 32 or args.steps % 2:
        parser.error("Choose a height of at least 1536 and an even step count of at least 32.")
    check_blocks()
    check_dirichlet()
    check_walls(args.height, args.steps)


if __name__ == "__main__":
    main()
