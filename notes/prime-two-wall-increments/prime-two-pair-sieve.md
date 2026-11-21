# Large ratio pairs

At $h=1/\log X$, a signed part of the original prime 2 sum over distinct ratios with both denominators above $R(X)$ is $o(h)$. It contains every pair for which the prime exponents differ at some odd prime $p\le(\log X)^2$. A separate absolute estimate removes pairs for which the numerator or denominator of the reduced quotient $r/r'$ is at most $R(X)^2$. I could not prove $o(h)$ for the remaining pairs. The cutoff in the Dirichlet polynomial tends to infinity first throughout these two estimates. [Original coefficients and cutoff](prime-two-increments.md#distinct-ratios) and [Masullo, §12.3](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514).

## Matching exponents

Fix $1/2<d<e<1$ and $0\le\omega\in C_c((d,e))$. Put $h=1/\log X$ and take $X$ large enough that the two walls over the support of $\omega$ stay in $I=[d,e]$. With the original pulse $g$, write $W_s(r)=r^{s-1/2}g(\log r)$, $`H(z)=\int g(v)e^{zv}\,dv`$ and $b_r(s)=a_{K,X}(s+h/2,r)-a_{K,X}(s-h/2,r)$. All sums over $r$ in the energy below have reduced denominator $q(r)>R(X)$. Here $R(X)=\lfloor e^{c_I S_d(X)/2}\rfloor$, with $S_d(X)=(\log X)^{1/(1-d)}(\log\log X)^{d/(1-d)}$ and the same $c_I>0$ as in the coefficient estimate. [Squared increments](prime-two-increments.md#prime-phases).

Choose a finite set $P$ of odd primes. Give $p\in P$ a unit phase $u_p$, put $\chi_u(r)=\prod_{p\in P}u_p^{v_p(r)}$ and integrate the phases against Haar probability measure. The identity

```math
\int_{\mathbb T^P}\left|\sum_r b_r(s)\chi_u(r)\right|^2du
=\sum_{\substack{r,r'\\ v_p(r)=v_p(r')\ \forall p\in P}} b_r(s)\overline{b_{r'}(s)}
\tag{1}
```

holds before either limit by character orthogonality and absolute convergence at finite $K,X$. Subtracting (1) from the squared sum at $u_p=1$ gives exactly the signed sum of pairs which differ at some $p\in P$. All these pairs are distinct. Lean checks the finite Fourier identity on a product of unit circles. [Lean proof](PrimeTwoPairSieve.lean).

## Odd twists

For fixed $X,P$, let $K\to\infty$. The alternating Dirichlet sum twisted at $P$ and its complex derivative tend locally uniformly, uniformly in the phases, to $f_u$ and $f'_u$, where

$$
f_u(z)=\eta(z)M_u(z),\qquad
M_u(z)=\prod_{p\in P}\frac{1-p^{-z}}{1-u_pp^{-z}},\qquad
\eta(z)=(1-2^{1-z})\zeta(z).
$$

For $\Re z>0$ each factor of $M_u$ has constant coefficient one and coefficient $u_p^j-u_p^{j-1}$ at $p^j$ for $j\ge1$. Its series converges absolutely there. In the convolution with the alternating series the multiplier is supported on odd integers, which preserve the sign $(-1)^{n-1}$ of the other factor. At $p^k$ the multiplier coefficients through $p^k$ sum to $u_p^k$. This recovers the original twisted truncation with the cutoff on the product index. Uniform absolute convergence of the multiplier on a slightly larger half-plane, together with partial summation of the alternating series, gives the local convergence and polynomial bounds in $|\Im z|$ for the truncated series and its derivative, with constants allowed to depend on $P$. These bounds justify the limit under the wall integrals at fixed $X,P$.

For the sum over all ratios, the twisted increment after this limit is the mass between the two walls for the complex measure

```math
d\mu_{X,u}(z)=\frac X\pi H(z-1/2)|f'_u(z)|^2e^{-X|f_u(z)|^2}\,dA(z).
\tag{2}
```

Indeed differentiating each original weighted coefficient gives twice $W_s(r)$ times the Fourier coefficient at $r$ of $X|\dot E_K|^2e^{-X|E_K|^2}$. Since $W_s(r)=(2\pi)^{-1}\int H(s-1/2+it)e^{-it\log r}dt$, the sum evaluates the density along $z_p=u_p p^{-it}$, with $u_p=1$ off $P$. Integration across the two walls gives (2). For $q(r)\le R(X)$ the absolute sum of the increments after the first limit is $O_I(hR(X)e^{-c_I S_d(X)})$ by Cauchy–Schwarz and the original coefficient bound, uniformly in the phases. Its $L^2(\omega)$ norm divided by $\sqrt h$ tends to zero. [Coefficient derivative](prime-two-increments.md#one-coefficient).

## Growing sets

Take $1/2<d_0<d$. The following argument allows $P=P_X$ to grow provided

$$
\sum_{p\in P_X}p^{-d_0}=o(\log X).
\tag{3}
$$

Uniformly in the phases and $\Re z\ge d_0$, the bounds $|M_u(z)|,|M_u(z)|^{-1}\le\exp(C_{d_0}\sum_{p\in P_X}p^{-d_0})=X^{o(1)}$ follow by multiplying the elementary upper and lower bounds for each Euler factor. The factors have no zeros or poles on this half-plane. On $I$ Cauchy's estimate bounds $|M'_u|$ by $X^{o(1)}$ as well.

To bound the high ordinates in (2), divide the area at $|f_u|=X^{-1/4}$. Above this value the exponential is at most $e^{-\sqrt X}$. On $I$, $|f'_u|\ll_I(1+|\Im z|)X^{o(1)}$ by partial summation and the derivative bound. The absolute integral on this part tends to zero by rapid decay of $H$. Below the cutoff use $w=f_u(z)$ as area variable, counting preimages with multiplicity. They are zeros of $\eta(z)-w/M_u(z)$. For $|w|\le X^{-1/4}$ the perturbation is bounded by $X^{-1/4+o(1)}$ uniformly on $\Re z\ge d_0$.

For each height $T$, centre a disc at $4+iT$, choose an inner radius between $4-d$ and $4-d_0$, and use outer radius $4-d_0$. The inner disc covers $I$ across a height interval of fixed positive length. At its centre $|\eta(4+iT)|\ge1-\sum_{n\ge2}n^{-4}>3/4$, and on the outer disc $|\eta|\ll_I 1+|T|$ by partial summation (with the constant also depending on $d_0$). Jensen's formula thus bounds the number of preimages in each band by $O_{I,d_0}(\log(2+|T|))$, uniformly in $u$ and large $X$. The Gaussian in $w$ has integral one. Summing this bound times the supremum of $|H|$ in each band gives bounded total variation and tails tending uniformly to zero as the height bound increases. [Area argument](prime-two-denominators.md#every-wall).

On a fixed rectangle choose the edges away from zeros and the vertical edges outside the wall intervals. Surround each zero $\rho$ of $\eta$ there by a disc of radius $h^2$. On its circumference $|\eta|\gg_\rho h^{2m_\rho}$, with $m_\rho$ the multiplicity. This exceeds $X^{-1/4+o(1)}$ for every fixed $m_\rho$. Rouché's theorem gives exactly $m_\rho$ preimages of every $|w|\le X^{-1/4}$ in this disc. There are no others on the rectangle, by the local lower bound for $|\eta|$ near each zero and its positive minimum away from fixed zero neighbourhoods. Integrating the small values in $w$ uses Gaussian mass $1-e^{-\sqrt X}$. Thus the mass near $\rho$ tends to $m_\rho H(\rho-1/2)$ uniformly in the phases, in discs of horizontal diameter $o(h)$.

For two positions $x,x'$ on the real axis, the kernel of the energy divided by $h$ is $h^{-1}\int_{|s-x|<h/2,\ |s-x'|<h/2}\omega(s)ds$. It is bounded by $`\|\omega\|_\infty`$. For two zero discs with different real parts it eventually vanishes, and for a common real part $`\beta`$ it tends uniformly on the pair of discs to $`\omega(\beta)`$. The uniform variation estimate lets us restore the large heights after taking this limit on the rectangle. In particular, with $`J_\beta=\sum_{\Re\rho=\beta}m_\rho H(\rho-1/2)`$,

$$
\sup_{u\in\mathbb T^{P_X}}\left|
\lim_{K\to\infty}\frac1h\int\omega(s)\left|\sum_{q(r)>R(X)} b_r(s)\chi_u(r)\right|^2ds
-\sum_\beta\omega(\beta)|J_\beta|^2\right|\longrightarrow0 .
\tag{4}
$$

For each fixed $X$ the passages in $K$ are uniform in $u$, including the removal of the finite set $q\le R(X)$. Both the untwisted square and its phase average in (1) therefore have the same normalised limit in (4). Their difference is $o(h)$ before normalisation. Choosing every odd $p\le(\log X)^2$ meets (3), since even summing over all integers costs $O((\log X)^{2(1-d_0)})$. More generally $p\le(\log X)^A$ works for each fixed $0<A<1/(1-d)$ by choosing $d_0$ sufficiently close to $d$. [Original wall energy](prime-two-increments.md#two-scales).

## Relative quotients

There is an absolute estimate which can be applied inside the remaining signed sum. At each fixed quotient $u=r/r'$, Cauchy–Schwarz bounds the sum of $|b_r b_{r'}|$ by $\sum_r|b_r|^2$, since multiplication by $u$ is injective. After taking $K\to\infty$ the bound is $O_I(h^2 e^{-2c_I S_d(X)})$, uniformly in the wall centre. The quotient of two ratios in the pulse lies in $(4/25,25/4)$. Thus the number of possible reduced quotients $a/b$ with $\min(a,b)\le R(X)^2$ is $O(R(X)^4)$. Their absolute contribution, integrated against $\omega$, is $O_{I,\omega}(h^2)$. The same bound holds for any subfamily by positivity of the absolute sum. One can use a limit superior in $K$ for this bound. For a finite set of quotient shifts the sums of products also converge by the $\ell^2$ convergence of the increment vector at fixed $X$. [Coefficient bound](prime-two-increments.md#distinct-ratios) and [finite quotient bound in Lean](PrimeTwoPairSieve.lean).

Consequently the entire possible coefficient of $h$ is carried by pairs with $q(r),q(r')>R(X)$ whose reduced quotient $a/b$ has $\min(a,b)>R(X)^2$ and no odd prime factor of $ab$ at most $(\log X)^2$. Their combined normalised limit is still $\sum_\beta\omega(\beta)|J_\beta|^2$. Grouping this residual by $v_2(a)-v_2(b)$ confines the possible term of order $h$ to two disjoint bands with opposite nonzero limits. [Two dyadic bands](prime-two-dyadic-modes.md).

## The remaining pairs

The residual class already contains positive subseries at the first nonzero Taylor order of the square. Fix a large $X$ and hold its $h,R,Y=(\log X)^2$ fixed while expanding in an auxiliary heat parameter $\lambda$ at zero. Choose $0<\varepsilon<1$ with $g(\log r)=1$ on $(2,2+\varepsilon)$. Take denominators that are odd primes $p>\max(R,Y)$, and numerators $n$ with $n/p$ in this interval and $n\equiv1\pmod{2\prod_{3\le l\le Y,\ l\ {\rm prime}}l}$. Each ratio is reduced. Every multiple of $(n,p)$ has positive sign in the coefficient of $\lambda$ in $b_{n/p}$ because both $n,p$ are odd. Its first multiple alone contributes $\gg h(\log p)^2p^{-2e}$ for $n\le K$, by the exact factor $2(\log n)\sinh(h\log p)p^{-2s}/\sqrt{n/p}$. [Wall averaging](prime-two-increments.md#wall-averaging).

Each half of the ratio interval supplies $\gg_Y p$ such numerators for each sufficiently large prime $p\le K/(2+\varepsilon)$. Pair the two halves. No resulting quotient has an odd prime factor at most $Y$. The product of their two positive sums diverges as $K\to\infty$, since $\sum_{p}p^{1-2e}(\log p)^2$ diverges. Removing $\min(a,b)\le R^2$ costs at most $O_I(R^4 h^2)$ pointwise at this Taylor order. Indeed the coefficient derivative at $\lambda=0$ in the wall derivative is $2W_s\widehat{|\dot E_K|^2}$, and $`\|\dot E_K\|_4`$ is bounded uniformly in $K,s$ on $I$. The same quotient bound applies. A nonzero $\omega\ge0$ preserves the divergence after integration. Thus even in the residual class absolute summation of the $\lambda^2$ coefficient loses cancellations. This calculation leaves the signed, fully resummed heat expression at $\lambda=X$ open.

## Heat factor

Increasing the prime cutoff further changes the heat factor even at an interior real point. Fix $1/2<\sigma<1$ and average all odd phases up to $Y$ after the Dirichlet limit. Orthogonality in one Euler factor gives the exact second moment

$$
\mathbb E_u |f_u(\sigma)|^2
=|\eta(\sigma)|^2\prod_{3\le p\le Y}\frac{(1-p^{-\sigma})^2}{1-p^{-2\sigma}}
=|\eta(\sigma)|^2\prod_{3\le p\le Y}\frac{1-p^{-\sigma}}{1+p^{-\sigma}} .
$$

Put $Y=(\log X)^A$. If $A<1/(1-\sigma)$ then $|M_u(\sigma)|^{-1}\le\exp(O_\sigma(Y^{1-\sigma}))=X^{o(1)}$ for every phase choice. Since $\eta(\sigma)>0$ by the alternating series, the mean of $1-e^{-X|f_u(\sigma)|^2}$ tends to one. If $A>1/(1-\sigma)$ it tends to zero. Indeed it is bounded by $X\mathbb E|f_u(\sigma)|^2$, and the logarithm of the product above is at most $-2\sum_{3\le p\le Y}p^{-\sigma}\le-cY^{1-\sigma}/\log Y$ for some $c>0$ and large $Y$. The prime lower bound follows from $4^m/(2m+1)\le {2m\choose m}\le(2m)^{\pi(2m)}$, since each prime power in the central binomial coefficient is at most $2m$. For example at $\sigma=3/4$ and $Y=(\log X)^5$, the original heat factor tends to one at phases equal to one, and its odd-phase average tends to zero. This computes the factor at one point. The target energy involves the derivative and the whole wall integral.

There is also a direct check on taking the odd primes all at once. At finite $K$, averaging through every odd prime $p\le K$ leaves only quotients which are powers of two. The support of the pulse bounds the number of those quotients, so their absolute contribution vanishes after division by $h$. Using this order to evaluate (4) would require justifying the change of limits. Already the linear odd-phase mean of the original alternating polynomial has different orders at the interior point $s=3/4$. For a fixed odd-prime cutoff $Y$, the limit in $K$ is $\eta(s)\prod_{3\le p\le Y}(1-p^{-s})$, which tends to zero as $Y\to\infty$ because $\sum_p1/p$ diverges. At each finite $K$, first averaging all odd primes leaves only $1-\sum_{1\le j,\ 2^j\le K}2^{-js}$. Its limit is $(1-2^{1-s})/(1-2^{-s})=-0.4667212022\ldots$. This gives an order mismatch for the linear mean at a zero-free point. The question about the squared tail still requires control of the signed pairs above.

At $K=4$, $X=64$, $h=1/\log64$ and the illustrative cutoff $q>8$, averaging the prime 3 phase on a 288-point circle agreed with explicitly selecting the matching exponents of the sampled Fourier coefficients to within $3\times10^{-18}$ at every sampled wall centre. With the weight in the script the normalised contribution from differing exponents was $-0.000231297709$. This finite calculation checks the sign and the projection in (1). [Numerical checks](prime-two-pair-sieve.py).
