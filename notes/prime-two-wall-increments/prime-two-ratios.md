# Summing the ratios

For the original prime 2 pulse, the full contribution at the fixed ratio $7/3$ tends to zero. This includes every pair $(7k,3k)$ and every order of the exponential regularisation. Every other fixed ratio in the pulse apart from $2$ and $1/2$ has the same cancellation. In the order that gives the actual inner-wall integral, those other ratios together contribute $P_f-R_\sigma>1/600$. The known pole term $P_f$ accounts for this difference to within $10^{-20000}$. Its exact equality to the difference is still open in this calculation. [Original pulse](prime-two-kernel.md#the-test) and [inner-wall integral](prime-two-eta.md).

Write $a=\log 2$, $L=a/\sqrt2$, $g(u)=f(u)+f(-u)$ and $`H(z)=\int g(u)e^{zu}\,du`$. Then $P_f=H(1/2)$. The fixed-ratio calculation holds for every $1/2<\sigma<1$. For the integral $I_\sigma$ take a wall with $|\sigma-\Re\rho|\ge c_\sigma(1+|\Im\rho|)^{-2}$ for every nontrivial zeta zero. Almost every inner wall qualifies. In the formula below $R_\sigma=\sum_{\Re\rho>\sigma}m_\rho H(\rho-1/2)$ counts both signs of the ordinate.

## One ratio

Set $\eta_K(s)=\sum_{n\le K}(-1)^{n-1}n^{-s}$, summing over positive integers. For $X>0$ put $\Phi_X(w)=(1-e^{-X|w|^2})/w$ with $\Phi_X(0)=0$. The previous calculation gives the actual integral as $\lim_{X\to\infty}\lim_{K\to\infty}I_\sigma(\eta_K'\Phi_X(\eta_K))$.

To collect all terms at a ratio, give each prime $p$ a unit complex phase $z_p$. Write $z_n=\prod_p z_p^{v_p(n)}$ and $z_{n/m}=z_n\overline{z_m}$. Replace $n^{-\sigma-it}$ by $n^{-\sigma}z_n$ in the two Dirichlet polynomials, calling the results $E_K$ and $\dot E_K$. The dot includes the derivative factor $-\log n$. On the product of the unit circles use normalised Haar measure $dz$ and define

```math
c_{K,X}(r)=\int \dot E_K\Phi_X(E_K)\overline{z_r}\,dz,\qquad
W_\sigma(r)=r^{\sigma-1/2}g(\log r).
```

Orthogonality collects exactly the terms with numerator divided by denominator equal to $r$. At finite $K,X$ the exponential series converges in the norm given by the sum of the absolute Fourier coefficients, since $E_K$ is a finite polynomial. Substituting $z_p=p^{-it}$ and using Fourier inversion therefore gives

$$
I_\sigma(\eta_K'\Phi_X(\eta_K))=\sum_{r\in\mathbb Q_{>0}}W_\sigma(r)c_{K,X}(r).
$$

At fixed $X$ the Taylor orders can also be summed directly under the vertical integral of $\eta'\Phi_X(\eta)$. The pulse bound with $n=\lfloor\sqrt{|t|/(4\cdot5120)}\rfloor$ and $n!\le n^n$ gives $|H(\sigma-1/2+it)|\le C e^{-c_0\sqrt{|t|}}$ for large $|t|$ and some $c_0>0$. The classical convexity bound gives $|\eta(\sigma+it)|^2=O_\sigma(|t|^{1-\sigma}\log^2|t|)=o(\sqrt{|t|})$. Since $\eta'$ has polynomial growth, $|XH\eta'\bar\eta|e^{|X||\eta|^2}$ is an integrable majorant for the sum of absolute Taylor integrands, uniformly on bounded sets of $X$. The issue with the order of limits arises when $X$ tends to infinity. [Pulse derivative bound](prime-two-eta.md#residue-bound) and [Titchmarsh, Theorem 11](https://brittlebooks.library.illinois.edu/brittlebooks_closed/Books2011-10/titced0001zetof/titced0001zetof.pdf).

For example, the individual pair $(7,3)$ has unweighted coefficient $-(j+1)\log 7$ in $`\eta_K'\eta_K^j\overline{\eta_K}^{\,j+1}`$ when $K\ge7$. Summing over $j$ gives $-(\log 7)Xe^{-X}$ before multiplication by $21^{-\sigma}W_\sigma(7/3)$. The prime-phase integral also sums every multiple of this pair. [Lean coefficient calculation](PrimeTwoRatios.lean) and [Hedenmalm, Lindqvist and Seip on the character space of Dirichlet series](https://arxiv.org/abs/math/9512211).

## All multiples

Let $Z_P=\prod_{p\le P}(1-p^{-\sigma}z_p)^{-1}$ and $Y_P=1/Z_P$. Their Fourier coefficients are $n^{-\sigma}$ and $\mu(n)n^{-\sigma}$ on the integers with prime factors at most $P$, and zero otherwise. Orthogonality and $\sum n^{-2\sigma}<\infty$ give limits $Z,Y$ in $L^2(dz)$. Products converge in $L^1$, so $ZY=1$ almost everywhere. The same argument with an extra factor $\log n$ gives an $L^2$ limit $\dot Z$ for $\dot Z_P$. Taking the limit of the finite Euler product derivative in $L^1$ gives

$$
\dot Z=-ZV,\qquad V=\sum_{p,j\ge1}(\log p)p^{-j\sigma}z_p^j,
$$

where $p$ runs over primes and $V$ converges in $L^2$. Put $v=2^{1-\sigma}>1$ and $A=1-vz_2$. Multiplication by $A$ gives the coefficient $(-1)^{n-1}n^{-\sigma}$ at each $n$. Thus $E_K,\dot E_K$ converge in $L^2$ to $E=AZ$ and $\dot E=avz_2Z+A\dot Z$. Since $|A|\ge v-1$, the quotient exists almost everywhere and has the expansion

$$
T=\frac{\dot E}{E}=-a\sum_{j\ge0}v^{-j}z_2^{-j}-V .
$$

For fixed $X$, the bounded continuous $\Phi_X$ takes $E_K\to E$ in probability to $\Phi_X(E_K)\to\Phi_X(E)$ in $L^2$. Multiplying by the derivatives gives convergence in $L^1$. Sending $X$ to infinity then gives $T(1-e^{-X|E|^2})\to T$ in $L^2$. It follows that $\lim_{X\to\infty}\lim_{K\to\infty}c_{K,X}(r)$ is the Fourier coefficient of $T$ at $r$. The expansion of $T$ contains only $1$, positive integer prime powers and inverse powers of $2$. The pulse selects precisely two of them.

| Ratio $r$ | Limiting coefficient | Contribution with $W_\sigma(r)$ |
| --- | --- | --- |
| $2$ | $`-a\,2^{-\sigma}`$ | $-L$ |
| $1/2$ | $`-a\,2^{\sigma-1}`$ | $-L$ |

The coefficient at $7/3$ is zero after this limit. On an outer wall $b>1$, expanding $avz_2/(1-vz_2)$ with $v=2^{1-b}<1$ instead gives a total coefficient $a2^{-b}$ at $2$ and zero at $1/2$, recovering the outer integral $L$. Haar measure extracts the coefficients at finite $K,X$. Only the finite polynomials and their smooth regularisations are evaluated at the particular phases $z_p=p^{-it}$ in this step.

## The total

Taking the fixed-ratio limits and then adding their contributions gives $-2L$. The actual order gives a different result on every wall satisfying the stated separation. Subtracting just the two rows of the table from the exact finite sum before taking its limits gives

```math
\lim_{X\to\infty}\lim_{K\to\infty}
\sum_{r\notin\{2,1/2\}}W_\sigma(r)c_{K,X}(r)
=I_\sigma(\eta'/\eta)+2L=P_f-R_\sigma .
```

This is strictly positive for the chosen pulse. Its normalisation gives $`P_f=\int w(u)f(u)\,du`$ with $w(u)=e^{-u/2}/(1-e^{-2u})$. On the pulse $w>1/2$, and $w(u+1/8)<e^{-1/16}w(u)$. The bump $B$ has integral $1/16$ and the right-hand translate has coefficient $-c$ with $0<c<1$. Hence $P_f>(1-e^{-1/16})/32>1/544$. The exterior sum is real by conjugation. The previous absolute bound $|R_\sigma|<10^{-20000}$ now gives $P_f-R_\sigma>1/600$. [Pole normalisation](prime-two-full-trace.md) and [absolute bound](prime-two-eta.md#residue-bound).

Even allowing the retained set of ratios to grow polynomially in $X$ misses this contribution when $K\to\infty$ first. For each positive integer $k$, $|E|^{-1}$ has finite $`2k`$th moment. Indeed the moment of each factor of $Y_P$ is $\sum_{j=0}^k\binom{k}{j}^2p^{-2j\sigma}\le1+\binom{2k}{k}p^{-2\sigma}$. The product is bounded uniformly in $P$, so an almost everywhere convergent subsequence and Fatou's lemma give the assertion for $Y/A=1/E$. For $c_X(r)=\lim_{K\to\infty}c_{K,X}(r)$, Cauchy–Schwarz gives the uniform bound

```math
\sup_r |c_X(r)-\widehat T(r)|
\le \|T\|_2
\left(\frac{k}{2eX}\right)^{k/2}
\left(\int |E|^{-2k}\,dz\right)^{1/2}.
```

Here $u^k e^{-u}\le(k/e)^k$ was applied to $u=2X|E|^2$. Since $W_\sigma$ is bounded, any choice of $O(X^d)$ other ratios contributes a sum tending to zero, by taking $k>2d$. The positive total comes from summing all ratios before removing the regularisation.

For a cutoff by reduced denominator, the same conclusion holds up to $\exp(C(\log X)^{1/(1-\sigma)})$ for each fixed $C>0$. Splitting each coefficient at finite $K,X$ gives a part of the tail with limit exactly $P_f$ and a second part with limit $-R_\sigma$. Counting the coprime numerators gives an explicit divergent boundary term in the first Taylor coefficient. [Large denominators](prime-two-denominators.md).

## Pole subtraction

Set $G(s)=(s-1)\zeta(s)/s$ with the removable value at $1$. Its logarithmic derivative adds $1/(s-1)-1/s$ to $\zeta'/\zeta$. The difference of the outer and inner integrals of that added term is exactly $H(1/2)=P_f$ by the residue at $1$. The zeros of $G$ on $\sigma<\Re s<3$ are just the exterior zeta zeros. Set $U_X=G'\Phi_X(G)$. Differentiating gives $\partial_{\bar s}U_X=X|G'|^2e^{-X|G|^2}$, including at the zeros by smooth extension. Cauchy–Green on the strip gives

```math
I_3(U_X)-I_\sigma(U_X)=
\frac X\pi\int_{\sigma<\Re s<3}
H(s-1/2)|G'(s)|^2 e^{-X|G(s)|^2}\,dA(s).
```

For fixed $X$ the horizontal integrals vanish by the rapid decay of $H$, polynomial growth of $G'$ and $|\Phi_X|\le\sqrt X$. On both vertical boundaries dominated convergence gives $I_3(G'/G)-I_\sigma(G'/G)=R_\sigma$ in the limit. The sign also follows by moving the right boundary inwards across the zeros of $G$.

The contribution near one such zero can be integrated exactly. At a zero $\rho$ of multiplicity $m$ choose a holomorphic coordinate $z$ with $z(\rho)=0$ and $G=z^m$. Take the neighbourhood $|z|<h^{1/m}$ inside this coordinate chart and the strip. Under this substitution $|G'|^2dA(s)=m^2|z|^{2m-2}dA(z)$. The circular mean of $H(s(z)-1/2)$ is $H(\rho-1/2)$. The resulting radial integral is $m(1-e^{-Xh^2})H(\rho-1/2)$ for every $X>0$. Thus a possible exterior zero contributes its full weighted multiplicity as $X\to\infty$.

Conjugate zeros contribute twice the real part of their weight. For every fixed real $x$ the function $\Re H(x+it)$ takes both signs at arbitrarily large $t$. It is the Fourier transform of $g(u)\cosh(xu)$, which vanishes near $u=0$, so all its even moments in $t$ vanish. If it were eventually of one sign, an interval of that sign beyond the last possible change would dominate the bounded initial part in sufficiently high even moments. Such an interval exists because the transform is real analytic and nonzero. This rules out a sign argument based on the pulse alone. The signed sum over the actual exterior zeros remains undetermined.
