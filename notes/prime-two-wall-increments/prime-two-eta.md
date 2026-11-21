# The inner wall

For the original prime 2 pulse I obtain an absolute error below $10^{-20000}$ in the proposed inner-wall value. Exact cancellation of the exterior residues remains unproved. Masullo’s §12.3 names the comparison between the global and local kernels. The pulse and the integral below come from applying that comparison to the rational field and trivial character at a fixed strip width. [Masullo, p. 514](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514) and [kernel comparison](prime-two-kernel-comparison.md).

Write $a=\log 2$, $L=a/\sqrt2$, $g(u)=f(u)+f(-u)$, $`F(z)=\int f(u)e^{zu}\,du`$ and $H(z)=F(z)+F(-z)$. The pulse has $f(a)=1$ and support in $(\log(3/2),\log(5/2))$. Put $P_f=H(1/2)$. Take $1/2<\sigma<1$ with $|\sigma-\Re\rho|\ge c_\sigma(1+|\Im\rho|)^{-2}$ for every nontrivial zero of $\zeta$ and some $c_\sigma>0$. Almost every wall qualifies. For a function $Q$ on the wall write

```math
I_\sigma(Q)=\frac{1}{2\pi}\int_\mathbb R H(\sigma-1/2+it)Q(\sigma+it)\,dt.
```

Wall separation gives polynomial growth of $\zeta'/\zeta$ there and $H$ decays faster than every inverse power. If $R_\sigma=\sum_{\Re\rho>\sigma}m_\rho H(\rho-1/2)$, including both signs of the ordinate and counting multiplicities, the contour shift gives $I_\sigma(\eta'/\eta)=P_f-2L-R_\sigma$. Reflection $\rho\mapsto1-\rho$ identifies $R_\sigma$ with the sum of $m_\rho F(\rho-1/2)$ on the two exterior half-planes. The trace comparison asks whether this particular sum is zero. [Contour and signs](prime-two-kernel-comparison.md#the-wall) and [original pulse](prime-two-kernel.md#the-test).

## Matching the pole

Write $A(s)=1-2^{1-s}$, $\eta=A\zeta$ and $D=A'/A=a/(2^{s-1}-1)$. Set $C(s)=1+(a^{-1}-1)2^{1-s}$ and $v=1-C\eta$. This choice gives $v(1)=0$ since $\eta(1)=a$. The constant Dirichlet coefficient of $v$ is also zero. For $N\ge1$ put

$$
M_N=AC\sum_{j=0}^{N-1}v^j,\qquad Q_N=D+\zeta'M_N .
$$

Here $M_N$ is entire, $M_N(1)=0$, $M_N'(1)=1$ and its constant Dirichlet coefficient is one. The geometric sum gives both identities

$$
1-\zeta M_N=v^N,\qquad \frac{\eta'}{\eta}-Q_N=\frac{\zeta'}{\zeta}v^N .
$$

These functions do approximate the logarithmic derivative on $\Re s=3$. There $|\eta-1|\le\sum_{n\ge2}n^{-3}\le1/4$, $|C-1|<1/4$ and $|v|\le9/16$. On that outer wall Fourier inversion gives $I_3(n^{-s})=g(\log n)/\sqrt n$. The original pulse selects only $n=2$, for which the coefficient of $\zeta'M_N$ is $-a$. Moving its integral inwards crosses just the pole at $1$, with residue $-M_N'(1)=-1$. Each fixed $N$ gives polynomial growth on the intervening strip, so the horizontal integrals tend to zero. Thus $I_\sigma(\zeta'M_N)=P_f-L$. The inner geometric expansion $D=-a\sum_{j\ge0}2^{j(s-1)}$ gives $I_\sigma(D)=-L$. Altogether

```math
I_\sigma(Q_N)=P_f-2L,\qquad I_\sigma\!\left(\frac{\zeta'}{\zeta}v^N\right)=-R_\sigma
\quad(N\ge1).
```

The second identity also follows by shifting the error itself. On the outer wall both $\zeta'/\zeta$ and $v^N$ have zero constant coefficient, so their product has zero coefficient at $2$. At $s=1$ the factor $v^N$ removes the pole. At every exterior zeta zero $\rho$ it has value one, so the residue of the weighted error is exactly $m_\rho H(\rho-1/2)$ for every $N$, including multiple zeros. This gives the minus sign when moving inwards. [Alternating zeta identity](https://dlmf.nist.gov/25.2.E3).

The approximation also fails locally on every inner wall. At $t_0=2\pi/a$ the factor $A(\sigma+it_0)$ is negative real and $C(\sigma+it_0)>0$. Euler–Maclaurin summation gives $\Re\zeta(\sigma+it_0)>1.027$ throughout $1/2\le\sigma\le1$. To check this, put $s=\sigma+it_0$ and $K=32$. To $\sum_{1\le n<K}n^{-s}$ add $K^{1-s}/(s-1)+K^{-s}/2+sK^{-s-1}/12$. The periodic second Bernoulli polynomial has absolute value at most $1/6$, giving an error at most $|s(s+1)|K^{-\sigma-1}/(12(\sigma+1))<1/32$. Here $K^{-it_0}=1$. Outward interval evaluation on eight equal subintervals gives $1.027<\Re\zeta<1.875$ after allowing for this error. [Euler–Maclaurin formula](https://dlmf.nist.gov/2.10.E1) and [interval calculation](prime-two-eta-check.py).

Consequently $\Re v(\sigma+it_0)>1$ for each $\sigma<1$ in this range. By continuity there is an interval of ordinates with $|v|\ge1+\epsilon_\sigma$. The analytic factor $H(s-1/2)\zeta'(s)/\zeta(s)$ has only isolated zeros there. On a smaller closed interval of positive length its modulus has a positive minimum. The integral of the absolute weighted error over that interval grows at least as $c(1+\epsilon_\sigma)^N$. Dominated convergence on the inner wall fails even on this bounded interval.

## A convergent expansion

For $X>0$ take $`Q_X^*=(\eta'/\eta)(1-e^{-X|\eta|^2})`$. Then $`I_\sigma(Q_X^*)\to I_\sigma(\eta'/\eta)`$ by dominated convergence. Writing $\eta_K(s)=\sum_{n=1}^K(-1)^{n-1}n^{-s}$ gives an expression for the actual integral in terms of the alternating coefficients and the pulse.

```math
I_\sigma(\eta'/\eta)=\lim_{X\to\infty}\lim_{K\to\infty}
\sum_{j\ge0}\frac{(-1)^jX^{j+1}}{(j+1)!}
I_\sigma\!\left(\eta_K'\eta_K^j\,\overline{\eta_K}^{\,j+1}\right).
```

For fixed $X,K$, Taylor expansion is uniform in $t$ since $\eta_K$ is bounded. Its sum is $\eta_K'(1-e^{-X|\eta_K|^2})/\eta_K$, with the quotient defined by continuity at zero. The quotient multiplying $\eta_K'$ has modulus at most $\sqrt X$. Partial summation of the alternating coefficients gives $|\eta_K'(\sigma+it)|\le c_\sigma(1+|t|)$ uniformly in $K$, and $\eta_K,\eta_K'$ converge pointwise. This justifies the limit in $K$ before the limit in $X$.

The iterated regularised integral has a value on every fixed inner wall. A uniform estimate for its area integral gives the same exterior sum on a wall containing no zero and adds half the residue of each zero on the wall itself to the sum being subtracted. [Every wall](prime-two-denominators.md#every-wall).

## Other ratios

Each finite Dirichlet polynomial in the expansion pairs with $H$ by Fourier inversion, using $g$ at logarithms of rational ratios. The right derivative of $I_\sigma(Q_X^*)$ at $X=0$ is $J(\sigma)=I_\sigma(\eta'\bar\eta)$, by domination by $|H\eta'\bar\eta|$. Its term indexed by $n,m$ in the square truncation is

```math
-(-1)^{n+m}\log(n)\,n^{-1/2}m^{1/2-2\sigma}g(\log(n/m)).
```

Summing all terms with $n=2m$ or $m=2n$ gives

$$
D_*(\sigma)=\frac{a}{\sqrt2}\eta(2\sigma)
-\left(\frac1{\sqrt2}+2^{1/2-2\sigma}\right)\eta'(2\sigma).
$$

Both sums on these ratios converge absolutely because $2\sigma>1$. The combined contribution of all other ratios is nonzero except possibly at a discrete set of inner walls. To prove this, extend $J$ holomorphically to $\Re w>0$ by integrating $H(w-1/2+it)\eta'(w+it)\eta(w-it)$ against $dt/(2\pi)$. The alternating series and partial summation give polynomial bounds uniform on compact subsets, and the Mellin factor decays rapidly there. The expression for $D_*$ also extends holomorphically.

For real $w>1$, the double sum converges absolutely. Grouping by denominator gives $J(w)-D_*(w)=\sum_{m\ge1} b_m m^{-2w}$, where $b_m$ sums $`-(-1)^{n+m}\sqrt{m/n}\log(n)\,g(\log(n/m))`$ over $n/m\ne2,1/2$. The support gives $b_m=O(m\log(2m))$. At $m=1,2$ the pulse selects only the dyadic ratios. At $m=3$ its only other ratio is $7/3$. Indeed $5/3$ is below the first bump and $\log(7/6)-1/8$ is strictly inside the second bump, of radius $13/256$. In terms of the original coefficient $c>0$ and bump $B$ this gives

```math
b_1=b_2=0,\qquad b_3=c\sqrt{3/7}\,\log(7)\,B(\log(7/6)-1/8)>0 .
```

Hence $3^{2w}(J(w)-D_*(w))\to b_3$ as $w\to+\infty$. A nonzero holomorphic function has isolated zeros, proving the assertion on the inner walls. At all orders the coefficient of each fixed ratio in the pulse apart from $2$ and $1/2$ tends to zero, including the whole $7/3$ contribution. Summing these other ratios before taking the limits gives $P_f-R_\sigma>1/600$. The equality $R_\sigma=0$ remains undetermined. [Summing the ratios](prime-two-ratios.md).

## Residue bound

The original pulse contains the probability density on $|u|<\delta=1/256$ proportional to $\phi(u/\delta)$, where $\phi(x)=e^{-1/(1-x^2)}$ on $|x|<1$ and is zero elsewhere. Write $y=1-x^2$. On the complex disc of radius $y/8$ about $x$, $|(1-z^2)-y|\le17y/64$, hence $\Re(1/(1-z^2))\ge 2/(5y)$. Cauchy's estimate gives

$$
|\phi^{(n)}(x)|\le n!(8/y)^n e^{-2/(5y)}\le20^n(n!)^2 .
$$

The middle expression tends to zero at both endpoints for fixed $`n`$. The last inequality follows by maximising $u^n e^{-u}$ and using $(n/e)^n\le n!$. The normalising integral of $\phi$ exceeds $e^{-4/3}$ by integration over $[-1/2,1/2]$. After rescaling, the density has $`n`$th derivative of $L^1$ norm at most $2e^{4/3}5120^n(n!)^2$. Integrating by parts $`n`$ times in its bilateral Laplace transform at $z=\rho-1/2$ improves the elementary bound $e^{\delta|\Re z|}$ by the factor

```math
2e^{4/3}\,\frac{5120^n(n!)^2}{|\Im\rho|^n}.
```

Platt and Trudgian's finite verification puts every exterior zero above $Y=3\times10^{12}$ in absolute ordinate. Bellotti and Wong's count gives $N(T)\le T\log T$ for $T\ge Y$. The eight box factors in the same pulse already give a sum of absolute weights below $6.5\times10^{-65}$, counting both signs and their multiplicities. Multiplying that estimate by the factor above, with $n=24000$ and $2e^{4/3}<8$, gives

$$
\sum_{\Re\rho>\sigma}m_\rho |H(\rho-1/2)|
\ \le\ \sum_{|\Re\rho-1/2|>\sigma-1/2}m_\rho|F(\rho-1/2)|
\ <\ 10^{-20000}.
$$

The last step is the integer inequality $`52\,10^{19935}5120^{24000}(24000!)^2<(3\times10^{12})^{24000}`$. The Lean file checks this inequality and the algebra of the reciprocal error. The contour shift, the derivative bound and the use of the published zero count are analytic arguments in this note. [Lean calculation](PrimeTwoEta.lean), [reproducible checks](prime-two-eta-check.py), [box estimate](prime-two-kernel.md#an-explicit-bound), [Platt and Trudgian](https://arxiv.org/abs/2004.09765) and [Bellotti and Wong](https://arxiv.org/abs/2412.15470).
