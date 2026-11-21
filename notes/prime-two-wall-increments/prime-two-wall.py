"""Evaluate a finite part of the wall integral for the original prime 2 pulse."""

import argparse
import importlib.util
from pathlib import Path

import mpmath as mp
import numpy as np


def load_pulse():
    path = Path(__file__).with_name("prime-two-check.py")
    spec = importlib.util.spec_from_file_location("prime_two_pulse", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pulse_transform(z, pulse, c, offsets, masses):
    mollifier = np.einsum("ij,j->i", np.exp(z[:, None] * offsets), masses)
    bump = 2 * np.sinh(pulse.HALF_BOX * z) / z
    bump *= (np.sinh(pulse.STEP * z) / (pulse.STEP * z)) ** pulse.BOXES * mollifier
    return 2 * bump * (np.cosh(pulse.CENTRE * z) - c * np.cosh((pulse.CENTRE + pulse.SHIFT) * z))


def zeta_and_derivative(s, log_n, n_weight, n, bernoulli):
    phases = np.exp(-1j * s.imag[:, None] * log_n)
    value = np.einsum("ij,j->i", phases, n_weight)
    derivative = np.einsum("ij,j->i", phases, -log_n * n_weight)

    # Scale each rising factorial by the corresponding power of N.
    tail = n / (s - 1) + 0.5
    d_tail = -n / (s - 1) ** 2
    product = s / n
    d_product = np.full(s.shape, 1 / n, dtype=complex)
    for k, coefficient in enumerate(bernoulli, start=1):
        tail += coefficient * product
        d_tail += coefficient * d_product
        for offset in (2 * k - 1, 2 * k):
            factor = (s + offset) / n
            d_product = d_product * factor + product / n
            product = product * factor

    n_to_minus_s = np.exp(-s * np.log(n))
    return value + n_to_minus_s * tail, derivative + n_to_minus_s * (d_tail - np.log(n) * tail)


def verify_theta_signs():
    iv = mp.iv
    iv.dps = 35
    panels = 4096
    step = iv.mpf(2) / panels
    sums = [iv.mpf(0), iv.mpf(0)]
    for j in range(panels + 1):
        u = j * step
        scale = iv.pi * iv.exp(2 * u)
        terms = iv.mpf(0)
        for n in range(1, 4):
            q = n * n * scale
            terms += (4 * q * q - 6 * q) * iv.exp(-q)
        weight = 1 if j in (0, panels) else 4 if j % 2 else 2
        weighted = weight * 2 * iv.exp(u / 2) * terms
        for index, t in enumerate((14, 15)):
            sums[index] += weighted * iv.cos(t * u)

    # The proof note bounds the Simpson error and the omitted theta terms together.
    error = iv.mpf([-3, 3]) / 10**6
    enclosures = [total * step / 3 + error for total in sums]
    for t, enclosure in zip((14, 15), enclosures):
        print(f"xi(1/2 + {t}i) in {enclosure}")
    if not (enclosures[0] > 0 and enclosures[1] < 0):
        raise ArithmeticError("The enclosures do not establish the sign change.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=4096)
    parser.add_argument("--steps-per-unit", type=int, default=16)
    parser.add_argument("--sigma", type=float, default=0.75)
    parser.add_argument("--theta-signs", action="store_true", help="Check the Gaussian theta signs and exit.")
    parser.add_argument("--pole-split", action="store_true", help="Integrate both parts of the regularised pole subtraction.")
    args = parser.parse_args()
    if args.theta_signs:
        verify_theta_signs()
        return
    if args.height < 8 or args.steps_per_unit < 2 or args.steps_per_unit % 2:
        parser.error("Use a height of at least eight and an even number of steps per unit.")
    if not np.isfinite(args.sigma) or args.sigma <= 0.5 or args.sigma == 1:
        parser.error("Choose a finite wall to the right of 1/2, away from the zeros at real part 1.")

    pulse = load_pulse()
    c, u, weights = pulse.probe_rule(128)
    v, mollifier_weights = pulse.mollifier_rule(64)
    a = np.log(2)
    local = a / np.sqrt(2)
    pole = np.dot(weights, 2 * np.cosh(u / 2))
    heights = np.arange(args.height * args.steps_per_unit + 1) / args.steps_per_unit
    integrand = np.empty_like(heights)
    heat_parameters = np.array([1.0, 8.0, 64.0, 512.0])
    if args.pole_split:
        split_integrands = np.empty((len(heights), len(heat_parameters), 3))

    n = max(128, int(np.ceil(args.height / 2)))
    log_n = np.log(np.arange(1, n, dtype=float))
    n_weight = np.exp(-args.sigma * log_n)
    bernoulli = [float(mp.bernpoly(2 * k, 0) / mp.factorial(2 * k)) for k in range(1, 17)]
    for start in range(0, len(heights), 192):
        batch = slice(start, start + 192)
        s = args.sigma + 1j * heights[batch]
        z = s - 0.5
        value, derivative = zeta_and_derivative(s, log_n, n_weight, n, bernoulli)
        dyadic = a / np.expm1(a * (s - 1))
        eta_log_derivative = derivative / value + dyadic
        h = pulse_transform(z, pulse, c, v, mollifier_weights)
        integrand[batch] = (h * eta_log_derivative).real / np.pi
        if args.pole_split:
            two = np.exp(a * (1 - s))
            eta = (1 - two) * value
            multiplier = 1 + (1 / a - 1) * two
            pole_part = dyadic + multiplier * (1 - two) * derivative
            remainder = (1 - multiplier * eta) * derivative / value
            heat = -np.expm1(-np.abs(eta[:, None])**2 * heat_parameters)
            terms = np.column_stack((eta_log_derivative, pole_part, remainder))
            weighted = (h[:, None] * terms).real / np.pi
            split_integrands[batch] = heat[:, :, None] * weighted[:, None, :]

    reference = f"P_f - 2L={pole - 2 * local:.14f}" if args.sigma < 1 else f"L={local:.14f}"
    print(f"sigma={args.sigma:g}  {reference}")
    for height in sorted({args.height // 4, args.height // 2, args.height}):
        stop = height * args.steps_per_unit + 1
        for stride in (2, 1):
            values = integrand[:stop:stride]
            integral = stride / args.steps_per_unit * (values[1:-1].sum() + (values[0] + values[-1]) / 2)
            print(f"T={height:5d}  dt={stride / args.steps_per_unit:g}  integral={integral:.14f}")

    if args.pole_split:
        totals = np.trapz(split_integrands, dx=1 / args.steps_per_unit, axis=0)
        print("    X             total         pole part         remainder")
        for heat_parameter, (full, pole_part, remainder) in zip(heat_parameters, totals):
            print(f"{heat_parameter:5g} {full:17.12f} {pole_part:17.12f} {remainder:17.12f}")
        np.testing.assert_allclose(totals[:, 0], totals[:, 1] + totals[:, 2], rtol=5e-11, atol=5e-13)

    mp.mp.dps = 40
    sample = args.sigma + 1j * np.array([0, 14.125, args.height / 2, args.height])
    value, derivative = zeta_and_derivative(sample, log_n, n_weight, n, bernoulli)
    errors = [
        abs(estimate - complex(mp.zeta(complex(s), derivative=1) / mp.zeta(complex(s))))
        for s, estimate in zip(sample, derivative / value)
    ]
    print(f"Largest sampled zeta log-derivative error={max(errors):.2g}")


if __name__ == "__main__":
    main()
