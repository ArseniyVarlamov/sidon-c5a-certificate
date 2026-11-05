# Independent verification of C_5a <= 0.9428095

The machine-checked part is a finite rational certificate. This file supplies the
ordinary mathematical proof that its inequalities imply the stated Sidon bound.
It is not a Lean formalization. In particular, checking the certificate does not
purport to machine-prove Cauchy--Schwarz or an asymptotic expansion.

## 1. Precisely what the certificate contains

Put n = 512, m = 2n = 1024, L = 12, M = Lm, and D = 10^18.
For 1 <= ell <= n, define the probability vector

\[
 p^{(\ell)}_i = \begin{cases}
 (2\ell)^{-1},& n-\ell\le i<n+\ell,\\
 0,&\text{otherwise},
 \end{cases}\qquad 0\le i<m,
\]

and mixing weights

\[
 \lambda_\ell=\begin{cases}
 4\ell^3/n^4,&\ell<n,\\
 (2n-1)/n^2,&\ell=n.
 \end{cases}
\]

Each kernel is symmetric under i -> m-1-i and nonnegative. Their masses and the
sum of the mixing weights are exactly 1; the verifier checks these conditions.
This kernel family is from Akwei V6, Appendix D, not a new family.

The compact file `certificate.json` contains integers V_k for
0 <= k < M/2, and the integer repair r = 43922. It defines the boundary weights by

\[
 W^{(\ell)}_j=
 \left\lfloor {1\over\ell}
   \sum_{\substack{0\le k<M/2\\n-\ell\le j-2k<n+\ell}} V_k\right\rfloor+r,
 \qquad w^{(\ell)}_j=W^{(\ell)}_j/D\quad(0\le j<M).
\]

Set w_j = 1 for j >= M, separately for every kernel. The floor is part of the
exact definition, not a floating-point rounding convention.

For reproducibility the constructor starts with

\[
 V_0=\left\lfloor {Dn^2\over 2n-1}\right\rfloor,
 \quad
 V_k=\left\lfloor{
 2\sum_{j=\max(0,k-n+1)}^{k-1}V_j+
 \boldsymbol{1}_{k\ge n}V_{k-n}
 \over 2n-1}\right\rfloor.
\]

The recurrence is inspired by the discrete renewal multiplier on page 25 of
Akwei V6. Placing this multiplier at the even sites gives the convolution above
for the centered uniform kernels. **Neither the recurrence nor a theorem about
its convergence is assumed by the verifier.** It treats the V_k as arbitrary
integers and checks the resulting finite boundary vectors directly.

## 2. The exact finite checks

For every integer q with 0 <= q <= M, check

\[
 \sum_{\ell=1}^n\lambda_\ell\sum_{i=0}^{m-1}
 p_i^{(\ell)}w_{q+i}^{(\ell)}\ge1.\tag{1}
\]

All **12289 of 12289** inequalities pass. The minimum is exactly 1, attained
only at q = M. Every later shift has all boundary entries equal to 1 and hence
also has covering sum exactly 1.

The check is over integers: the common denominator in (1) is 2n^4 D, with
numerator coefficient 4 ell^2 for ell < n and n(2n-1) for ell = n.
The constructor uses prefix sums; the independent verifier reconstructs each
boundary by adding pulses and checks the constraints with rolling sums.

Define

\[
 a=m\sum_\ell\lambda_\ell\sum_i(p_i^{(\ell)})^2,
 \quad
 b=1+2\left({1\over m}\sum_\ell\lambda_\ell
                     \sum_{j=0}^{M-1}(w_j^{(\ell)})^2-L\right).
\]

The exact results are

\[
 a={349525\over262144},\qquad
 b={1172814639933515173834491879581015223161385476153
 \over1759218604441600000000000000000000000000000000000}>0,
\]

\[
 ab={16397121480910475645380030968422173835019330342095093
 \over18446744073709551616000000000000000000000000000000000}.
\]

Integer cross-multiplication verifies

\[
 ab\le\left({942809491151629\over10^{15}}\right)^2
   <\left({9428095\over10^7}\right)^2.\tag{2}
\]

The displayed decimal 0.9428094911516289306... is informational only. Acceptance
uses the rational squared comparison, not a computed floating-point square root.

## 3. From a finite certificate to every sufficiently large Sidon set

This is a self-contained presentation of the symmetric smoothing argument of
[Hou--Zhao, Theorem 2.1, v3](https://arxiv.org/html/2607.01169v3#S2).
It uses no renewal identity, method-barrier theorem, numerical optimum, or
two-sided extension.

Let A be a Sidon subset of {0,...,N-1}, with k elements. Thus each nonzero ordered
difference occurs at most once. Fix an integer h >= 1 and H = mh, with N >= 2LH.
For every ell, define the probability kernel on the integers

\[
 K^{(\ell)}_s=p^{(\ell)}_{\lfloor s/h\rfloor}/h
 \quad(0\le s<H),\qquad K^{(\ell)}_s=0\quad\text{otherwise},
\]

and u_ell = 1_A * K_ell. The kernels remain symmetric under s -> H-1-s.
For the nonnegative autocorrelation C_ell(d) = sum_s K_ell(s)K_ell(s+d),
sum_d C_ell(d) = 1 and C_ell(-d) = C_ell(d). Therefore the Sidon condition gives

\[
 \|u_\ell\|_2^2
 =kC_\ell(0)+\sum_{x,y\in A,\ x\ne y}C_\ell(x-y)
 \le kC_\ell(0)+\sum_{d\ne0}C_\ell(d)
 =1+(k-1)C_\ell(0).
\]

Since C_ell(0) = h^(-1) sum_i p_i^2, the mixed energy obeys

\[
 E:=\sum_\ell\lambda_\ell\|u_\ell\|_2^2
 \le1+{a(k-1)\over H}.\tag{3}
\]

The support of u lies in J = {0,...,N+H-2}. On J assign Q_ell(j) the weight
w_floor(j/h) in the first LH sites, the reflected weight
w_floor((N+H-2-j)/h) in the last LH sites, and 1 in between.
These regions are disjoint. Each interval {x,...,x+H-1} with 0 <= x < N meets at
most one boundary region.

For x = qh+t near the left endpoint, with 0 <= t < h, its mixed covering sum is
exactly

\[
 \sum_\ell\lambda_\ell\sum_{s=0}^{H-1}K_\ell(s)Q_\ell(x+s)
 =\left(1-{t\over h}\right)
  \sum_\ell\lambda_\ell\sum_i p_i^{(\ell)}w_{q+i}^{(\ell)}
 +{t\over h}\sum_\ell\lambda_\ell\sum_i p_i^{(\ell)}w_{q+1+i}^{(\ell)}
 \ge1.\tag{4}
\]

The last boundary block uses the checked shift q+1 = M. At the right endpoint,
write x = N-1-(qh+t) and replace s by H-1-s. Symmetry reduces it to (4).
In the interior the sum is exactly 1.

Sum (4) over x in A and apply weighted Cauchy--Schwarz to obtain

\[
 k^2\le\left(\sum_\ell\lambda_\ell\|Q_\ell\|_2^2\right)E.
\]

Counting sites gives the exact identity

\[
 \sum_\ell\lambda_\ell\|Q_\ell\|_2^2
 =N+H-1-2LH+2h\sum_\ell\lambda_\ell\sum_{j=0}^{M-1}w_j^2
 =N+bH-1.
\]

Together with (3), this proves, for every admissible integer h,

\[
 \boxed{k^2\le(N+bH-1)\left(1+{a(k-1)\over H}\right).}\tag{5}
\]

For sufficiently large N choose a positive multiple H of m within m of
sqrt(a/b) N^(3/4). Then N >= 2LH eventually, because m,L,a,b are fixed and b>0.
The distinct positive differences give k(k-1)/2 <= N-1, hence k = O(sqrt(N)).
Substitution into (5) first yields k^2 <= N+O(N^(3/4)), and consequently
k <= sqrt(N)+O(N^(1/4)). Substituting this sharper estimate back into (5) yields

\[
 k^2\le N+2\sqrt{ab}\,N^{3/4}+O(N^{1/2}),
 \quad
 k\le\sqrt N+\sqrt{ab}\,N^{1/4}+O(1).
\]

Translating {0,...,N-1} to {1,...,N} preserves the Sidon property. Thus (2) gives

\[
 \boxed{C_{5a}\le\sqrt{ab}\le0.942809491151629<0.9428095.}
\]

This independently verifies the starred 0.9428095 bound. The tiny stronger
numerical value is not claimed as a new research record: Akwei already describes
convergence to 2 sqrt(2)/3, and reports an earlier claim at that limit by Chen.

## 4. Scope of the machine check

`verify.py` checks all finite hypotheses of the theorem and the rational squared
comparison. `audit_checker.py` additionally compares its result with literal
Fraction sums for 54 small covering constraints, checks the lifted endpoint
covering and norm formulas at 600 positions, and rejects a corrupted full-size
certificate. These checks audit implementation and indexing. The universal
argument in section 3 is an ordinary proof, not a finite enumeration or a
claim of full formal verification.
