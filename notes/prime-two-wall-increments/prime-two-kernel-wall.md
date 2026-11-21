# One vertical line

For the original prime 2 pulse I get an integral on three walls which gives $J_\beta$ as the walls coalesce. Applied to logarithmic moduli, the test gives zero for every finite odd Euler product. I could not justify the same evaluation for the global kernel. A calculation on the actual zero modes gives the obstruction to a separate commutator proof. Equation (2) takes the original Dirichlet length to infinity at fixed heat. [Original pulse](prime-two-kernel.md#the-test), [two-band identity](prime-two-dyadic-modes.md#opposite-limits) and [Masullo, §12.3](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514).

## Three walls

Fix $1/2<\beta<1$, set $h=1/\log X$ and $d=h/2$, and take $X$ large enough that $[\beta-d,\beta+d]\subset(1/2,1)$. Let $g(v)=f(v)+f(-v)$, $H(z)=\int g(v)e^{zv}dv$ and $p(v)=g(v)\cosh((\beta-1/2)v)$. Use $\widehat p(t)=\int p(v)e^{itv}dv$. This transform is $\Re H(\beta-1/2+it)$, and conjugation pairs all ordinates on the line. For a zero $\rho$ write $m_\rho$ for its multiplicity. Our target is $J_\beta=\sum_{\Re\rho=\beta}m_\rho H(\rho-1/2)$.

For $|x-\beta|\le d$ put

```math
V(x,t)=\int p(v)\frac{\sinh((d-|x-\beta|)v)}{\sinh(dv)}e^{itv}\,dv
```

and put $V=0$ outside. The pulse support is compact and stays away from $v=0$. In either open half of the strip, $\partial_x^2 V=-\partial_t^2 V$. On the centre wall $V=\widehat p$, while $V$ vanishes on the two outside walls. Set $b_d(v)=p(v)v/\sinh(dv)$ and $c_d(v)=p(v)v\coth(dv)$. The jumps of the first horizontal derivative give the following functional on a real potential $U$.

```math
D_{\beta,d}(U)=\frac1{2\pi}\int_{\mathbb R}
\left[\widehat b_d(t)\bigl(U(\beta+d,t)+U(\beta-d,t)\bigr)-2\widehat c_d(t)U(\beta,t)\right]dt
=\frac1{2\pi}\iint V\,\Delta U .
\tag{1}
```

This is integration by parts for $\Delta V$. It gives zero for a harmonic potential whose value and first derivatives have at most polynomial growth on the strip. In particular, if $U=\log|E(x+it)|$ for a holomorphic function, a zero $\rho$ contributes $m_\rho V(\Re\rho,\Im\rho)$. This follows from $\Delta\log|w-\rho|=2\pi\delta_\rho$. A pole has the opposite sign. For the completed zeta function the formula applies on every wall, including one through a zero. One justification is its genus-one product. On a unit height interval near $T$ each nearby logarithmic singularity has integrable negative part bounded by 2, since $|x+it-\rho|\ge|t-\Im\rho|$. Factors with $|\rho|$ up to $2(|T|+3)$ therefore cost at most a polynomial in $|T|$ in local $L^1$, by the zero count. For the other factors the subtracted logarithm is $O((1+|T|)^2/|\rho|^2)$. These estimates allow (1) term by term against the rapidly decaying test. For an individual logarithm one may smooth in $x$ first and pass to its locally integrable trace on each wall. [Hadamard product and zero count](https://terrytao.wordpress.com/2021/02/12/246b-notes-4-the-riemann-zeta-function-and-the-prime-number-theorem/).

## Original order

Take $E=\eta_K(w)=\sum_{n\le K}(-1)^{n-1}n^{-w}$ and use the potential

```math
U_{K,X}(w)=\frac12\int_0^{X|\eta_K(w)|^2}\frac{1-e^{-y}}{y}\,dy .
```

Differentiating in the radial variable gives $\Delta U_{K,X}=2X|\eta'_K|^2e^{-X|\eta_K|^2}$, including at a zero by smooth extension. At each fixed $X$ the limit in $K$ passes through (1). Indeed $\eta_K,\eta'_K$ converge locally uniformly to $\eta,\eta'$ and have uniform polynomial bounds on the closed substrip by partial summation of the alternating series. Thus the right side of (1) uses exactly the Gaussian area density of the original trace after its first limit. [Original coefficient identity](prime-two-increments.md#one-coefficient).

As $X$ grows, $V$ and its $t$ derivatives decrease faster than every inverse power of $|t|$ uniformly in $d$ and $x$. For $0\le a\le1$ the multiplier $\sinh(adv)/\sinh(dv)$ and its $v$ derivatives are uniformly bounded on the fixed pulse support, by the power series of $\sinh y/y$ at zero. Its horizontal derivative costs $O(1/d)$ on either side of the centre. At bounded height the Gaussian measure places mass tending to $m_\rho$ in the disc of radius $h^2$ about $\rho$. The variation of $V$ across the disc is $O(h^2/d)=O(h)$. Every zero on the centre line gives its exact limiting weight. Every other zero at bounded height eventually lies outside the strip. Jensen's bound for the preimages of $\eta$ and the rapid decay of $V$ give uniformly negligible tails as the height tends to infinity. The factor $1-2^{1-w}$ of $\eta$ has no zero on this strip. This proves

$$
\lim_{X\to\infty}\lim_{K\to\infty}D_{\beta,1/(2\log X)}(U_{K,X})=J_\beta .
\tag{2}
$$

The cutoffs in (2) move with the original heat. The disc and height estimates use the unshifted kernel throughout. [Gaussian localisation and uniform tails](prime-two-pair-sieve.md#growing-sets).

## The odd primes

On the rational radial global theta generator $e^{u/2}\sum_{n\ge1}\Phi(ne^u)$, take $\Phi(x)=(x\partial_x)(x\partial_x+1)e^{-\pi x^2}$. At Mellin argument $w-1/2$ its transform is $M_\Phi(w)\zeta(w)=\xi(w)$ with $M_\Phi(w)=\tfrac12w(w-1)\pi^{-w/2}\Gamma(w/2)$. Keeping only powers of two gives the genuine dyadic generator with transform $M_\Phi(w)/(1-2^{-w})$. Its logarithmic modulus is harmonic on the test strip and (1) gives zero exactly. For the global generator (1) leaves precisely the odd factor $Z_o(w)=(1-2^{-w})\zeta(w)$. This compares two specified radial inputs. The full local Fourier trace has its own cutoff computation. [Global generator](prime-two-kernel.md#strip-widths) and [local trace](prime-two-full-trace.md#local-cutoffs).

At every fixed $d$ the exact result for the global generator is

```math
D_{\beta,d}(\log|\xi|)=D_{\beta,d}(\log|Z_o|)
=\sum_{|\Re\rho-\beta|\lt d}m_\rho V(\Re\rho,\Im\rho).
\tag{3}
```

A zero outside the two outer walls contributes exactly zero already at this spacing. In particular, the pole at 1 creates no term in (3). A finite odd Euler product gives zero even in each of the three separate integrals. The real logarithm of each factor has only cosine frequencies $j\log p$, $j\ge1$, all above $\log(5/2)$, whereas $b_d,c_d$ are supported at $|v|<\log(5/2)$. The infinite Euler series establishes the same frequency statement on a wall to the right of 1.

Here is a direct check of what goes wrong with using that series on a single interior wall. Write $B_p(s)=(2\pi)^{-1}\int\widehat p(t)\log|Z_o(s+it)|dt$ with the even test $p$ fixed. For $1/2<s<1$ one has

```math
B_p(s)=-\int p(v)\frac{\sinh((1-s)v)}{v}dv
+\sum_{\Re\rho>s}m_\rho\int p(v)\frac{\sinh((\Re\rho-s)v)}{v}e^{i\Im\rho\,v}dv .
\tag{4}
```

Integrate Green's identity on $s<x<2$ using the harmonic test $\int p(v)\sinh((x-s)v)e^{itv}dv/v$. Its value at $s$ is zero and its derivative there is $\widehat p$. On the wall at 2 both terms vanish by the absolutely convergent odd Euler series and the support of $p$, also after multiplying $p$ by $\cosh((2-s)v)$ or $\sinh((2-s)v)/v$. The sources in the strip are the pole at 1 and the zeros to its left. The product estimate used for (1) justifies Green's identity by testing the individual factors. Equivalently subtract the smooth gamma and dyadic logarithms from $\log|\xi|$, retaining the pole at 1. A zero exactly on $s$ has weight zero in (4).

Take $s=\beta=3/4$. The pole integral before its minus sign is $P=\int_0^\infty f(v)j(v)dv$ with $j(v)=\sinh(v/2)/v$. Put $l=1/8$. For the normalising weight $r$ of the pulse, $r(y+l)>e^{l/2}r(y)$ on its first bump because $r(y)e^{-y/2}$ is positive and strictly increasing there. Thus its translate coefficient $c$ is less than $e^{-l/2}$. Also $j(v)=\int_0^{1/2}\cosh(tv)dt$ has $j'(v)/j(v)<1/4$ for $0<v<\log(5/2)<1$, using $\tanh(tv)\le tv$. Since $j\ge1/2$ and the first bump has mass $1/16$,

$$
P>\frac{1-e^{-1/32}}{32}>\frac1{1056}.
$$

The absolute sum of the zero terms in (4) is less than $10^{-63}$ here. Write $\sinh((\Re\rho-s)v)/v=\int_0^{\Re\rho-s}\cosh(tv)dt$. Multiplying by $\cosh(v/4)$ gives an average of two real shifts of the pulse transform, both of absolute size at most $1/2$. The eight-box estimate is uniform in these shifts, and every relevant zero lies above the finite verified height. Thus $B_p(3/4)<-1/1100$, although every finite odd Euler product gave zero for this very integral. In the three-wall combination the pole term cancels exactly by harmonicity near the centre. The zero terms in (4) give (3). In fact, keep $p$ fixed at the chosen $\beta$ and take the two one-sided derivatives in $s$ in (4). A term with $\Re\rho=\beta$ has left derivative $-m_\rho\widehat p(\Im\rho)$ and right derivative zero. All other terms have matching derivatives there. Rapid Fourier decay is uniform for the $\cosh((\Re\rho-s)v)$ derivative multipliers on a compact substrip, so absolute summability permits the one-sided limits. Hence $B_p'(\beta+)-B_p'(\beta-)=J_\beta$. The missing statement is differentiability at every such wall after the smooth pole term is removed. [Pulse and box bound](prime-two-kernel.md#an-explicit-bound), [Platt and Trudgian](https://arxiv.org/abs/2004.09765) and [Bellotti and Wong](https://arxiv.org/abs/2412.15470).

## The commutator

The actual global space also gives a finite-rank test of the proposed cancellation. Fix a genuine zero $\rho$ of multiplicity $m$, let $z=\rho-1/2$ and take a width $\alpha>|\Re z|$ below $1/2$. Write $W_\alpha$ for the closed radial theta space in $X_\alpha=L^2(e^{2\alpha|u|}(1+u^2)^2du)$. The Mellin evaluations through order $m-1$ at $z$ vanish on $W_\alpha$ by unfolding every global theta input and continuity.

One can isolate this block without a basis for the other zeros. Write $G$ for the global Gaussian theta generator above and apply $q\mapsto e^{-zu}\int_u^\infty e^{zt}q(t)dt$ to it $m$ times. At each step the Mellin transform divides by $w-z$ in the centred variable $w$. The vanishing moment at that step gives the equal expression $-e^{-zu}\int_{-\infty}^u e^{zt}q(t)dt$, so both tails still decay faster than every exponential. The resulting vector $v$ has nonzero Mellin transform at $z$, and $(-\partial_u-z)^m v=G$. Adjoin its translates to $W_\alpha$ and take the closure $W'$. Modulo $W_\alpha$ the translate $v(u+t)$ satisfies $(\partial_t+z)^m[v(u+t)]=0$, so these translates span at most $m$ dimensions in the quotient. The first $m$ Mellin evaluations of the translates are $e^{-zt}$ times polynomials in $t$ of distinct degrees $0$ to $m-1$, with nonzero leading terms. They give the matching lower bound. Thus the difference of the two orthogonal projections has rank $m$. Its trace after convolution by $f$ is $m\int f(t)e^{zt}dt$, since the induced convolution on this quotient has that integral on each diagonal entry. Adding the reflected zero gives $mH(z)$. This finite change of the global space measures one block of the original projection. [Radial spaces and conventions](prime-two-full-trace.md#global-projection).

The same finite-rank comparison can be made at the exact width $\alpha=|\Re z|>0$. The rank there is $\min(m,2)$. The evaluations of orders zero and one are continuous on the wall because $\int_1^\infty u^{2j}(1+u^2)^{-2}du$ is finite for $j=0,1$. They give the lower bound up to the multiplicity. For the upper bound, $e^{-\alpha|t|}(1+|t|)^{-2}S_t v$ tends weakly to zero as $|t|\to\infty$. The normalised operators have norm at most one by the translation bound, and the pairings of compactly supported vectors vanish eventually. Density gives the weak limit and it persists on the quotient. On the side where $|e^{-zt}|=e^{\alpha|t|}$, a nonzero coefficient of degree at least two in the quotient polynomial would violate this limit. Hence its degree is at most one. For $\alpha<|\Re z|$ the same growth bound rules out every polynomial coefficient. Thus the traces of this particular projection difference immediately below, exactly at and immediately above the width are $0$, $\min(m,2)\int f(t)e^{zt}dt$ and $m\int f(t)e^{zt}dt$ respectively. No assertion about the span of all the other boundary modes is needed.

On the dual block use $e_j(u)=u^j e^{zu}/j!$, $0\le j<m$. For $S_t e(u)=e(u+t)$ set $T_k=\int k(t)S_tdt$ and $M e(u)=u e(u)$. The original reflected pulse vanishes near zero, so $k(t)=g(t)/t$ is smooth and compactly supported. On the full polynomial orbit $[T_k,M]=T_g$. Writing $C(z)=\int k(t)e^{zt}dt$ gives $C'(z)=H(z)$ and $T_k e_j=\sum_{l=0}^j C^{(l)}(z)e_{j-l}/l!$. Multiplication by $u$ sends $e_{m-1}$ to $m e_m$, outside the retained block. If $M_m$ drops that last term and $T_{k,m}$ is the restriction, the difference $`T_g|_m-[T_{k,m},M_m]`$ has just the last column, with entries given by $m\sum_{l=1}^m C^{(l)}(z)e_{m-l}/l!$. Its trace is $m C'(z)=mH(z)$. This is a finite matrix computation on an invariant quotient and does not require commuting an unbounded operator through an infinite trace.

Conjugate zeros can cancel only through the real parts of the pulse weights. Reflection across the centre gives $H(-z)=H(z)$, so it gives equal weights. Equations (2) and (3) leave the full signed sum of these blocks on the centre wall. The original two bands would give $o(h)$ for every continuous nonnegative interior weight if these line sums vanished for every $1/2<\beta<1$. I could not establish that cancellation. Lean checks the differential equation and endpoint values of the test profile, the compressed matrix commutator and the trace of its last column. The Green identities, Gaussian limit and quotient argument are analytic. At heights 1536 and 2048 both grids give $B_p(3/4)\simeq-0.00372920$. The three-wall and area integrals of the original potential at $K=3$, $X=64$ both give $0.15992422$, agreeing to $1.4\times10^{-12}$ on the longer grid. These are finite quadrature checks. [Lean source](PrimeTwoKernelWall.lean) and [numerical check](prime-two-kernel-wall.py).
