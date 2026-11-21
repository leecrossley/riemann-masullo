"""Check the dyadic ordinate, the first ratios and the tail bound for the prime 2 pulse."""

from math import factorial

from mpmath import iv


def check_dyadic_ordinate():
    cutoff = 32
    panels = 8
    ordinate = 2 * iv.pi / iv.log(2)
    logs = [iv.log(n) for n in range(1, cutoff)]
    phases = [iv.cos(ordinate * log_n) for log_n in logs]
    remainder = iv.sqrt(1 + ordinate**2) * iv.sqrt(4 + ordinate**2)
    remainder /= 18 * iv.mpf(cutoff) ** (iv.mpf(3) / 2)
    if not remainder < iv.mpf(1) / 32:
        raise ArithmeticError("The Euler–Maclaurin remainder exceeds its allowance.")

    intervals = []
    for j in range(panels):
        sigma = iv.mpf([panels + j, panels + j + 1]) / (2 * panels)
        partial = sum(iv.exp(-sigma * log_n) * phase for log_n, phase in zip(logs, phases))
        pole = cutoff * (sigma - 1) / ((sigma - 1) ** 2 + ordinate**2)
        tail = iv.exp(-sigma * iv.log(cutoff)) * (pole + iv.mpf(1) / 2 + sigma / (12 * cutoff))
        intervals.append(partial + tail + iv.mpf([-1, 1]) / 32)

    lower = min(value.a for value in intervals)
    upper = max(value.b for value in intervals)
    if not lower > iv.mpf(1027) / 1000 or not upper < iv.mpf(1875) / 1000:
        raise ArithmeticError("The real part falls outside the stated enclosure.")
    print("1.027 < Re zeta(sigma + 2 pi i / log 2) < 1.875 throughout [1/2, 1].")


def check_exterior_bound():
    height = 3 * 10**12
    q = iv.mpf(1) / 128
    box_bound = iv.mpf(512) ** 8 / height**7
    box_bound *= iv.log(height) / (1 - q) + iv.log(2) / (1 - q) ** 2
    if not box_bound < iv.mpf(65) / 10**66 or not 2 * iv.exp(iv.mpf(4) / 3) < 8:
        raise ArithmeticError("The constants for the box or smooth density exceed their bounds.")

    order = 24_000
    exponent = 20_000
    numerator = 52 * 10 ** (exponent - 65) * 5120**order * factorial(order) ** 2
    if numerator >= height**order:
        raise ArithmeticError("The factorial bound exceeds the target.")
    print("The exterior residue sum has absolute value below 10^-20000.")


def check_first_ratios():
    centre = iv.log(2)
    centres = (centre, centre + iv.mpf(1) / 8)
    radius = iv.mpf(13) / 256
    selected = []
    for m in range(1, 4):
        for n in range(1, 5 * m // 2 + 1):
            if n == 2 * m or m == 2 * n:
                continue
            u = abs(iv.log(iv.mpf(n) / m))
            for bump, location in enumerate(centres, start=1):
                distance = abs(u - location)
                if distance < radius:
                    selected.append((n, m, bump))
                elif not distance > radius:
                    raise ArithmeticError("A ratio touches an unresolved support endpoint.")
    if selected != [(7, 3, 2)]:
        raise ArithmeticError("The first non-dyadic coefficient has unexpected terms.")
    print("Through denominator 3 the only other ratio is 7/3, in the second bump.")


if __name__ == "__main__":
    iv.dps = 35
    check_dyadic_ordinate()
    check_first_ratios()
    check_exterior_bound()
