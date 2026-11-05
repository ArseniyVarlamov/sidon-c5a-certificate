# Certificate for C_5a ≤ 0.9428095

This repository holds a certificate and a checker for the upper bound 0.9428095 on the
Sidon set constant C_5a from Terence Tao's table of optimization constants
(https://teorth.github.io/optimizationproblems/constants/5a.html).

The bound was announced by John Akwei (report v6, 2026) in the vector-valued smoothing
framework of Hou and Zhao (arXiv:2607.01169). His certificate files are not public yet,
so I built an independent certificate in the same framework.

Result: all 12289 covering inequalities hold, and

    C_5a ≤ γ = sqrt(ab) ≤ 0.942809491151629 < 0.9428095.

## Files

- `certificate.json` — 6144 integers that define the boundary weights, the integer
  repair term and the parameters (n = 512 kernels, m = 1024, L = 12, D = 10^18).
- `verify.py` — checks the certificate with exact integer and rational arithmetic.
  No floating point is used for any decision. Standard library only.
- `build_certificate.py` — rebuilds the certificate from scratch.
- `reproduce.py` — rebuilds the certificate, checks it, and runs the audit below.
- `audit_checker.py` — compares the checker with direct rational sums on small cases
  and checks that a corrupted certificate is rejected.
- `PROOF.md` — why the finite checks imply the bound: a restatement of Hou–Zhao,
  Theorem 2.1, in the symmetric case.

## Run

    python3 verify.py       # about 2 seconds
    python3 reproduce.py    # about 6 seconds

Python 3.9 or newer.

## What is new and what is not

The kernel family (centred uniform kernels with mixing weights 4ℓ³/n⁴) is taken from
Akwei's report, Appendix D. The boundary weights, the depth L = 12 and the repair term are
built here. The slightly smaller γ is not claimed as a new record: the limit 2√2/3 of this
method is already described by Akwei and others.

Prepared with help from AI tools (Codex, Claude).

## License

MIT
