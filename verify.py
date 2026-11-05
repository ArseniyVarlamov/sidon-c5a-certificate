#!/usr/bin/env python3
"""Exact, standalone verifier; does not import or trust the constructor.

All acceptance decisions use integers/Fraction. Decimal is output formatting only.
Every covering inequality is checked, including the last shift L*m.
"""
import argparse
import hashlib
import json
import time
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify(c):
    require(c["format"] == "sidon-centered-uniform-v1", "unknown format")
    n, m, depth, d, repair = (c[x] for x in ("n", "m", "L", "D", "repair"))
    require(all(type(x) is int for x in (n, m, depth, d, repair)), "noninteger parameter")
    require(n >= 1 and m == 2*n and depth >= 1 and d >= 1 and repair >= 0, "bad parameter")
    size = depth*m
    v = c["renewal_numerators"]
    require(len(v) == (size+1)//2 and all(type(x) is int for x in v), "bad multiplier data")
    # v may be arbitrary signed integers. No recurrence, convergence, positivity
    # or optimality assertion about v is needed for the certificate theorem.
    cover = [0] * (size+1)
    weight_sum = Fraction(0)
    kernel_energy = Fraction(0)
    energy_numerator = 0
    common = 2*n**4
    for ell in range(1, n+1):
        weight = Fraction(4*ell**3, n**4) if ell < n else Fraction(2*n-1, n*n)
        require(weight > 0, "negative weight")
        weight_sum += weight
        # p has value 1/(2*ell) on [n-ell,n+ell), zero elsewhere.
        left, right = n-ell, n+ell
        require(0 <= left < right <= m and left == m-right, "kernel symmetry/support")
        require(Fraction(right-left, 2*ell) == 1, "kernel mass")
        kernel_energy += m*weight/Fraction(2*ell)
        coeff = weight/Fraction(2*ell)*common
        require(coeff.denominator == 1, "common covering denominator")
        coeff = coeff.numerator

        # Reconstruct boundaries by adding pulses, independently of the
        # constructor's prefix-index formula. Integer floor is mathematical floor.
        changes = [0] * (size+1)
        for k, value in enumerate(v):
            first, end = 2*k+left, 2*k+right
            if first < size:
                changes[first] += value
                changes[min(end, size)] -= value
        total = 0
        boundary = []
        square_sum = 0
        for j in range(size):
            total += changes[j]
            value = total//ell + repair
            boundary.append(value)
            square_sum += value*value
        weight_num = weight*n**4
        require(weight_num.denominator == 1, "common energy denominator")
        energy_numerator += weight_num.numerator * square_sum

        # Rolling sums check all shifted supports; all tail entries are exactly d.
        padded = boundary + [d]*m
        running = sum(padded[left:right])
        for q in range(size+1):
            cover[q] += coeff*running
            if q < size:
                running += padded[q+right] - padded[q+left]

    require(weight_sum == 1, "mixing weights do not sum to 1")
    target = common*d
    minimum = min(cover)
    require(minimum >= target, "covering inequality failed at shift " + str(cover.index(minimum)))
    require(cover[-1] == target, "tail covering must equal one")
    a = kernel_energy
    b = 1 + 2*(Fraction(energy_numerator, n**4*d*d*m) - depth)
    require(a > 0 and b > 0, "nonpositive energy or boundary coefficient")
    gamma2 = a*b
    bound = Fraction(9428095, 10**7)
    passed = gamma2 < bound*bound
    with localcontext() as ctx:
        ctx.prec = 45
        decimal_gamma = str((Decimal(gamma2.numerator)/Decimal(gamma2.denominator)).sqrt())
    # A rigorous 15-decimal upper bound: smallest integer u with u^2 >= gamma^2*10^30.
    from math import isqrt
    scale = 10**15
    scaled_floor = gamma2.numerator*scale*scale//gamma2.denominator
    ceiling = isqrt(scaled_floor)
    if ceiling*ceiling*gamma2.denominator < gamma2.numerator*scale*scale:
        ceiling += 1
    upper = Fraction(ceiling, scale)
    require(gamma2 <= upper*upper, "square-root enclosure")
    return {
        "certificate_feasible": True,
        "verifies_upstream_0.9428095": passed,
        "m": m, "L": depth, "kernels": n,
        "covering_inequalities_checked": size+1,
        "covering_inequalities_passed": sum(x >= target for x in cover),
        "minimum_cover": str(Fraction(minimum, target)),
        "exact_cover_equalities": sum(x == target for x in cover),
        "a": str(a), "b": str(b), "gamma_squared": str(gamma2),
        "gamma_decimal_informational": decimal_gamma,
        "gamma_upper_rational": str(upper),
        "gap_below_upstream_squared": str(bound*bound-gamma2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("certificate", nargs="?", type=Path, default=Path(__file__).with_name("certificate.json"))
    ap.add_argument("--feasibility-only", action="store_true", help="allow feasible certificates above upstream")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    start = time.perf_counter()
    raw = args.certificate.read_bytes()
    result = verify(json.loads(raw))
    result["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    result["verification_seconds"] = time.perf_counter()-start
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    require(args.feasibility_only or result["verifies_upstream_0.9428095"], "feasible but not below upstream bound")


if __name__ == "__main__":
    main()
