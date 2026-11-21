# Squared wall increments

For the original prime 2 test in §12.3, the squared change in the large-denominator tail, integrated against a continuous nonnegative weight, is $O(h)$ when $h=1/\log X$ and the Dirichlet cutoff $K$ tends to infinity first. This follows from a bound on the total variation of the heat measure. Computing the coefficient of $h$ gives the term that would have to vanish for an $o(h)$ result. The argument uses neither RH nor Alpöge’s estimate. [Original test](prime-two-kernel.md#the-test) and [item 18 on issue #2469](https://github.com/mccorvie/riemann-masullo/issues/2469#issuecomment-7067444368).

## The tail

Fix $I=[d,e]\subset(1/2,1)$ and $0\le\omega\in C_c((d,e))$, extended by zero. Let $f$ be the smooth radial pulse with $f(\log 2)=1$, put $g(u)=f(u)+f(-u)$ and write $H(z)=\int_{\mathbb R}g(u)e^{zu}\,du$. For a positive rational $r=n/q$ in lowest terms, let $a_{K,X}(s,r)$ be its weighted Fourier coefficient on the wall $s$. At finite $K,X$ these are the coefficients of the regularised alternating Dirichlet polynomial, multiplied by $r^{s-1/2}g(\log r)$. [Coefficient definition and derivative](prime-two-increments.md#one-coefficient).

Set $S_d(X)=(\log X)^{1/(1-d)}(\log\log X)^{d/(1-d)}$ and $R(X)=\lfloor\exp(c_I S_d(X)/2)\rfloor$ for large $X$, with $c_I>0$ from the prime-phase estimate below. Write $T_{K,X}(s)=\sum_{q(r)>R(X)}a_{K,X}(s,r)$ and $\Delta_h T(s)=T(s+h/2)-T(s-h/2)$. The quantity in the comment is

$$
E_X=\lim_{K\to\infty}\int_{\mathbb R}\omega(s)\,|\Delta_h T_{K,X}(s)|^2\,ds,
\qquad h=\frac{1}{\log X}.
$$

At each fixed $X$, first take $K\to\infty$ in the full wall integral and in the finite block $q(r)\le R(X)$. Partial summation of the alternating series gives convergence on compact sets and polynomial bounds in height, uniform in $K$ on $I$. The rapid vertical decay of $H$ allows passage through the full integral. The finite block converges uniformly on $I$ by convergence in the prime phases. The cutoff is held fixed between the two walls. [Vertical integral](prime-two-denominators.md#every-wall) and [prime phases](prime-two-increments.md#prime-phases).

## The bound

Put $\eta(w)=(1-2^{1-w})\zeta(w)$. After the first limit, Cauchy–Green identifies the increment of the full wall integral with the mass between the walls of the complex measure

$$
d\mu_X(w)=\frac X\pi H(w-1/2)|\eta'(w)|^2 e^{-X|\eta(w)|^2}\,dA(w).
$$

Restrict this measure to $\Re w\in I$ and write $M_X=|\mu_X|(\Re w\in I)$ for its total variation. The $M_X$ are bounded uniformly for $X\ge1$. On $|\eta|>1/2$ the factor $Xe^{-X|\eta|^2}$ is bounded, and $\eta'$ grows at most polynomially in height. On $|\eta|\le1/2$, change variables to $z=\eta(w)$. Jensen’s formula bounds the number of preimages of any $|z|\le1/2$ in a unit height band at $j$ by $O_I(\log(2+|j|))$. The Gaussian in $z$ has mass one. Multiplying this count by the supremum of $|H|$ on the band gives a summable bound, including a uniform bound tending to zero for the tails in height. [Jensen discs and area estimate](prime-two-denominators.md#every-wall).

Let $A_h(s)=\{w:|\Re w-s|<h/2\}$. For $h$ small the walls lie in $I$ on the support of $\omega$. Write $d\mu_X=\phi\,d|\mu_X|$. Expanding $|\mu_X(A_h(s))|^2$ gives the factor $\phi(w)\overline{\phi(z)}$ for each pair in $A_h(s)^2$. For fixed $w,z$ the allowed centres form the intersection of two open intervals of length $h$. Fubini applies since $\omega$ is integrable, $|\phi|\le1$ almost everywhere and $M_X$ is finite. The integral of $\omega$ over the intersection is bounded by $h\|\omega\|_\infty$, giving

$$
\int_{\mathbb R}\omega(s)|\mu_X(A_h(s))|^2\,ds
\ \le\ h\|\omega\|_\infty M_X^2 .
$$

For the removed block $q(r)\le R(X)$, the coefficient derivative has $\ell^2$ norm $O_I(e^{-c_I S_d(X)})$ after the first limit. There are $O(R(X)^2)$ such ratios in the pulse. Integrating the derivative between the walls and applying Cauchy–Schwarz bounds their summed increment uniformly by $O_I(h R(X)e^{-c_I S_d(X)})=o(h)$. Its $L^2(\omega(s)\,ds)$ norm is consequently $o(h)$. The full increment has norm $O(\sqrt h)$ by the displayed inequality, so removing the block changes its squared norm by $o(h)$. [Coefficient estimate](prime-two-increments.md#prime-phases) and [denominator cutoff](prime-two-increments.md#distinct-ratios).

For the $O(h)$ conclusion alone, take any constant $B$ that bounds the removed increment divided by $h$ for all sufficiently large $X$ on $\{\omega>0\}$. Using $|a+b|^2\le2|a|^2+2|b|^2$ gives $E_X\le 2h\|\omega\|_\infty M_X^2+2B^2h^2\int\omega(s)\,ds$. Lean checks the change of integration order, this inequality for a complex density and the resulting $E_X=O(h)$ implication. The preceding area and coefficient estimates provide its two uniform bounds for the original test. [Lean proof](PrimeTwoWallEnergy.lean) and [build instructions](README.md#lean-checks).

## The coefficient

Divide the intersection integral by $h$. Its kernel at $x=\Re w$, $y=\Re z$ is bounded by $\|\omega\|_\infty$ and differs from $\omega(x)(1-|x-y|/h)_+$ by at most the modulus of continuity of $\omega$ at $h$. The uniform variation estimate allows the height tails to be discarded in this double integral with an error tending to zero uniformly in $X$.

On a bounded rectangle with boundary clear of zeros, the heat density tends to zero away from the finitely many zeros of $\eta$. The factor $1-2^{1-w}$ has no zeros in $I$. Near a zeta zero $\rho$ of multiplicity $m_\rho$, use a holomorphic coordinate $u$ with $\eta=u^{m_\rho}$. Rescaling $u$ by $X^{1/(2m_\rho)}$ gives a fixed integrable radial density of total mass $m_\rho$. Its width $X^{-1/(2m_\rho)}$ is $o(h)$ for each of these zeros. Thus a pair of zeros with the same real part $\beta$ contributes the kernel value $\omega(\beta)$ in the limit, and a pair with different real parts contributes zero. Restoring the uniform height tails and grouping the absolutely summable weights by real part gives

$$
E_X=h\sum_\beta\omega(\beta)|J_\beta|^2+o(h),
\qquad J_\beta=\sum_{\Re\rho=\beta}m_\rho H(\rho-\tfrac12).
$$

Conjugation makes each $J_\beta$ real. A continuous nonnegative weight can be positive at any specified interior $\beta$, and every term in the coefficient is nonnegative. Hence $E_X=o(h)$ for every such weight if and only if $J_\beta=0$ at every interior real part. The two dyadic bands cancel their common leading term and leave this same coefficient. The three-wall identity removes the pole at 1 exactly. I could not show the remaining line sums vanish. [Diagonal limit](prime-two-increments.md#two-scales), [two bands](prime-two-dyadic-modes.md#opposite-limits) and [three-wall calculation](prime-two-kernel-wall.md#original-order).
