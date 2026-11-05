#!/usr/bin/env python3
"""Construct an untrusted compact rational Hou--Zhao symmetric certificate.

Only the separate verify.py determines acceptance. Python standard library only.
The renewal recurrence guides construction; it is NOT a verifier assumption.
"""
import argparse
import json
import time
from pathlib import Path


def build(n, depth, denominator):
    m, size = 2 * n, 2 * n * depth
    # Floor-rounded renewal multiplier, stored in fixed-point integer units.
    v = [denominator * n * n // (2 * n - 1)]
    vp = [0, v[0]]
    for k in range(1, (size + 1) // 2):
        total = 2 * (vp[k] - vp[max(0, k - n + 1)])
        if k >= n:
            total += v[k - n]
        v.append(total // (2 * n - 1))
        vp.append(vp[-1] + v[-1])

    # Common denominator for lambda_l * p_l[i] is 2*n**4.
    cd = 2 * n**4
    cover = [0] * (size + 1)
    in_box = [0] * (size + 1)
    for ell in range(1, n + 1):
        coefficient = 4 * ell**2 if ell < n else n * (2 * n - 1)
        wp = [0]
        for j in range(size):
            lower = max(0, (j - n - ell + 2) // 2)
            upper = max(0, (j - n + ell) // 2 + 1)
            w = (vp[upper] - vp[lower]) // ell
            wp.append(wp[-1] + w)
        for q in range(size + 1):
            left, right = q + n - ell, q + n + ell
            count = max(0, min(size, right) - left)
            total = wp[min(size, right)] - wp[min(size, left)]
            total += (2 * ell - count) * denominator
            cover[q] += coefficient * total
            in_box[q] += coefficient * count
    repair = 0
    worst = 0
    for q, total in enumerate(cover):
        deficit = cd * denominator - total
        if deficit > 0:
            if in_box[q] == 0:
                raise ArithmeticError("unrepairable tail")
            required = (deficit + in_box[q] - 1) // in_box[q]
            if required > repair:
                repair, worst = required, q
    return {
        "format": "sidon-centered-uniform-v1",
        "n": n,
        "m": m,
        "L": depth,
        "D": denominator,
        "repair": repair,
        "renewal_numerators": v,
        "construction_diagnostics": {
            "repair_witness_shift": worst,
            "minimum_unrepaired_cover_numerator": min(cover),
            "cover_denominator": cd * denominator,
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=512)
    ap.add_argument("--depth", type=int, default=12)
    ap.add_argument("--denominator", type=int, default=10**18)
    ap.add_argument("--output", type=Path, default=Path(__file__).with_name("certificate.json"))
    args = ap.parse_args()
    if args.n < 1 or args.depth < 1 or args.denominator < 1:
        ap.error("positive parameters required")
    start = time.perf_counter()
    certificate = build(args.n, args.depth, args.denominator)
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "n": args.n, "m": 2*args.n,
                      "L": args.depth, "repair": certificate["repair"],
                      "construction_seconds": time.perf_counter() - start}))


if __name__ == "__main__":
    main()
