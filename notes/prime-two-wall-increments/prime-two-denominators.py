"""Compare the first heat coefficient in reduced and unreduced denominator sums."""

import argparse
from math import ceil, floor
from pathlib import Path
from runpy import run_path

import mpmath as mp
import numpy as np


def totients(limit):
    values = np.arange(limit + 1)
    for p in range(2, limit + 1):
        if values[p] == p:
            values[p::p] -= values[p::p] // p
    return values


def dyadic_boundary(limit, phi, z, zeta, dzeta):
    odd = np.arange(1, limit + 1, 2)
    lengths = np.frexp(limit // odd)[1]
    weights = -zeta / 2 * phi[odd] * odd**(-z) * 2.0**((1 - z) * lengths)
    slopes = weights * (dzeta / zeta - np.log(odd) - lengths * np.log(2))
    return weights.sum(), slopes.sum()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sigma", type=float, default=0.75)
    parser.add_argument("--limit", type=int, default=4096)
    parser.add_argument("--grid", type=int, default=32769)
    args = parser.parse_args()
    if not 0.5 < args.sigma < 1 or args.limit < 2 or args.grid < 3:
        parser.error("Choose 1/2 < sigma < 1, limit >= 2 and grid >= 3.")

    pulse = run_path(str(Path(__file__).with_name("prime-two-check.py")))
    c, nodes, weights = pulse["probe_rule"](128)
    pole = weights @ (2 * np.cosh(nodes / 2))
    pole_slope = weights @ (2 * nodes * np.sinh(nodes / 2))
    centre, shift, radius = (pulse[name] for name in ("CENTRE", "SHIFT", "SUPPORT"))
    grid = np.linspace(0, radius, args.grid)
    profile = np.concatenate([pulse["bump"](part, 64) for part in np.array_split(grid, 16)])
    centres = (centre, centre + shift, -centre, -centre - shift)
    bands = [(np.exp(u - radius), np.exp(u + radius)) for u in centres]

    z = 2 * args.sigma
    zeta = float(mp.zeta(z))
    dzeta = float(mp.zeta(z, derivative=1))
    v = 2.0**(1 - z)
    eta = (1 - v) * zeta
    deta = (1 - v) * dzeta + np.log(2) * v * zeta
    normalisation = 3 * float(mp.zeta(2))
    b = (2 * eta - zeta) / (normalisation * (2 - z))
    db = ((2 * deta - dzeta) / normalisation + b) / (2 - z)
    phi = totients(args.limit)
    reduced = unreduced = density = slope = 0.0
    print(f"sigma={args.sigma:g}  P={pole:.12f}  H'(1/2)={pole_slope:.12f}  leading={pole * b:.12f}")
    print("    Q          reduced     boundary     corrected      two-term     unreduced")
    for q in range(1, args.limit + 1):
        n = np.concatenate([np.arange(ceil(lo * q), floor(hi * q) + 1) for lo, hi in bands])
        x = n / q
        u = np.abs(np.log(x))
        left = np.interp(np.abs(u - centre), grid, profile, right=0.0)
        right = np.interp(np.abs(u - centre - shift), grid, profile, right=0.0)
        h = (left - c * right) / np.sqrt(x)
        log_n = np.log(n)
        even_q = q % 2 == 0
        same_parity = n % 2 == q % 2
        scale = q**(-z)
        unreduced += scale * (h @ np.where(same_parity, -log_n, log_n))
        primitive = np.gcd(n, q) == 1
        grouped = np.where(same_parity, -zeta * log_n + dzeta, eta * log_n - deta)
        reduced += scale * (h[primitive] @ grouped[primitive])

        average = eta if even_q else (eta - zeta) / 2
        d_average = deta if even_q else (deta - dzeta) / 2
        density += scale * phi[q] * average
        slope += scale * phi[q] * (d_average - np.log(q) * average)
        if q == args.limit or (q >= 16 and q & (q - 1) == 0):
            chain_value, chain_slope = dyadic_boundary(q, phi, z, zeta, dzeta)
            np.testing.assert_allclose([density, slope], [chain_value, chain_slope],
                                       rtol=5e-12, atol=5e-12)
            boundary = pole_slope * density - pole * slope
            smooth = q**(2 - z) * (pole * b * np.log(q) + pole_slope * b - pole * db)
            print(f"{q:5d}  {reduced:15.9f} {boundary:12.9f} {reduced - boundary:13.9f} "
                  f"{reduced - smooth:13.9f} {unreduced:13.9f}")


if __name__ == "__main__":
    main()
