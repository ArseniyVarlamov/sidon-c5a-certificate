#!/usr/bin/env python3
"""Small direct Fraction cross-check and a corruption check on the real artifact."""
import json
from fractions import Fraction as Q
from pathlib import Path

from build_certificate import build
from verify import require, verify


def direct(c):
    n, m, size, d = c["n"], c["m"], c["L"]*c["m"], c["D"]
    kernels = []
    for ell in range(1, n+1):
        p = [Q(1, 2*ell) if n-ell <= i < n+ell else Q(0) for i in range(m)]
        weight = Q(4*ell**3, n**4) if ell < n else Q(2*n-1, n*n)
        w = []
        for j in range(size):
            # Literal convolution, with no prefix sums or event arrays.
            numerator = sum(value for k, value in enumerate(c["renewal_numerators"])
                            if 0 <= j-2*k < m and p[j-2*k] != 0)
            w.append(Q(numerator//ell + c["repair"], d))
        w += [Q(1)]*m
        kernels.append((weight, p, w))
    covers = [sum(weight*sum(p[i]*w[q+i] for i in range(m))
                  for weight, p, w in kernels) for q in range(size+1)]
    a = m*sum(weight*sum(x*x for x in p) for weight, p, w in kernels)
    b = 1 + 2*(sum(weight*sum(x*x for x in w[:size])
                    for weight, p, w in kernels)/m-c["L"])
    return covers, a, b, kernels


def main():
    checked = 0
    lifted = 0
    for n, depth in ((1, 2), (2, 2), (3, 3), (5, 2)):
        c = build(n, depth, 10**5)
        covers, a, b, kernels = direct(c)
        fast = verify(c)
        require(str(a) == fast["a"] and str(b) == fast["b"], "energy disagreement")
        require(min(covers) == Q(fast["minimum_cover"]), "cover disagreement")
        require(sum(x >= 1 for x in covers) == fast["covering_inequalities_passed"], "coverage disagreement")
        checked += len(covers)

        # Check the finite lifting formula at every x, including both endpoints,
        # for h=1,2,3 and the minimal allowed N. This is an indexing audit,
        # not a substitute for the universal proof in PROOF.md.
        for h in (1, 2, 3):
            m, L = c["m"], c["L"]
            H = m*h
            N = 2*L*H
            expanded = []
            for weight, p, w in kernels:
                K = [p[s//h]/h for s in range(H)]
                boundary = [w[j//h] for j in range(L*H)]
                weights = boundary + [Q(1)]*(N+H-1-2*L*H) + boundary[::-1]
                require(len(weights) == N+H-1, "expanded support length")
                expanded.append((weight, K, weights))
            total_norm = sum(weight*sum(z*z for z in weights) for weight, K, weights in expanded)
            require(total_norm == N+b*H-1, "lifted boundary norm")
            for x in range(N):
                actual = sum(weight*sum(K[s]*weights[x+s] for s in range(H))
                             for weight, K, weights in expanded)
                require(actual >= 1, "lifted covering failed")
                lifted += 1

    path = Path(__file__).with_name("certificate.json")
    corrupted = json.loads(path.read_text())
    corrupted["renewal_numerators"][0] = 0
    try:
        verify(corrupted)
    except ValueError as error:
        require("covering inequality failed" in str(error), "unexpected corruption failure")
        rejection = str(error)
    else:
        raise ValueError("corrupted certificate was accepted")
    print(json.dumps({"audit_passed": True,
                      "direct_rational_cover_checks": checked,
                      "lifted_endpoint_cover_checks": lifted,
                      "real_certificate_corruption_rejected": rejection}, indent=2))


if __name__ == "__main__":
    main()
