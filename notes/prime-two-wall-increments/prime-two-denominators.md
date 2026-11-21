# Large denominators

For the original prime 2 pulse, an explicit part of the large-denominator sum gives exactly the pole weight $P_f$. The split takes place in each coefficient before either limit. On a wall containing no zero the other part gives $-R_\sigma$, the negative of the signed exterior residue sum. Its multiplier takes the value one at each possible exterior zero, and I could not show that these weights sum to zero. A zero on the wall itself contributes half its weight under this regularisation. The split holds above a reduced denominator cutoff as large as $\exp(C(\log X)^{1/(1-\sigma)})$ for any fixed $C>0$, where $X$ is the exponential regularisation parameter. [Masullo, §12.3](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514), [original pulse](prime-two-kernel.md#the-test) and [contour calculation](prime-two-eta.md).

Take $1/2<\sigma<1$, $g(u)=f(u)+f(-u)$, $`H(w)=\int_{\mathbb R}g(u)e^{wu}\,du`$, $P_f=H(1/2)$ and $L=(\log 2)/\sqrt2$. For any wall $b$ write $`I_b(V)=\int_{\mathbb R}H(b-1/2+it)V(b+it)\,dt/(2\pi)`$. Write $R_\sigma=\sum_{\Re\rho>\sigma}m_\rho H(\rho-1/2)$ over nontrivial zeta zeros, counting both signs of the ordinate. For the original regularised sum on every fixed inner wall, put $R_\sigma^*=R_\sigma+\tfrac12\sum_{\Re\rho=\sigma}m_\rho H(\rho-1/2)$. The two agree on every wall containing no zero, including the polynomially separated walls used for the global kernel comparison.

## Moving cutoff

Write each positive rational $r=n/q$ in lowest terms and put $W_\sigma(r)=r^{\sigma-1/2}g(\log r)$. On the product of the prime unit circles with normalised Haar measure, let $\chi_n=\prod_p u_p^{v_p(n)}$ and $\chi_{n/q}=\chi_n\overline{\chi_q}$. At finite $K$ the polynomial $E_K=\sum_{n\le K}(-1)^{n-1}n^{-\sigma}\chi_n$ has derivative $\dot E_K$ with factor $-\log n$. Let $c_{K,X}(r)$ be the coefficient at $\chi_r$ of $\dot E_K(1-e^{-X|E_K|^2})/E_K$, extended continuously at zero. The sum $\sum_r W_\sigma(r)c_{K,X}(r)$ is absolutely convergent at finite $`K,X`$ and is the vertical integral for the truncated alternating zeta function. [Derivation and normalisations](prime-two-ratios.md#one-ratio).

For $X>1$ put $Q_C(X)=\lfloor\exp(C(\log X)^{1/(1-\sigma)})\rfloor$. Then

```math
\lim_{X\to\infty}\lim_{K\to\infty}
\sum_{q>Q_C(X)}\ \sum_{\substack{n\ge1\\(n,q)=1}}
 W_\sigma(n/q)c_{K,X}(n/q)=P_f-R_\sigma^* .
```

To prove this, write $E=(1-2^{1-\sigma}u_2)Z$ for the $L^2$ limit of $E_K$, where $Z$ is the random Euler product and $Y=1/Z$ almost everywhere. Put $T=\dot E/E\in L^2$. For fixed $X$, the coefficients $c_{K,X}(r)$ tend to the coefficients $c_X(r)$ of $T(1-e^{-X|E|^2})$. The only frequencies of $T$ in the pulse are $2$ and $1/2$, each with weighted contribution $-L$. [Euler product and coefficients](prime-two-ratios.md#all-multiples) and [Hedenmalm, Lindqvist and Seip](https://arxiv.org/abs/math/9512211).

For any integer $k\ge1$, the $`2k`$th moment of one factor of $Y$ is at most $\exp(\min(2kx,k^2x^2))$, where $x=p^{-\sigma}$. Indeed $|1-xu_p|^{2k}\le(1+x)^{2k}\le e^{2kx}$, and its mean is $\sum_{j=0}^k\binom{k}{j}^2 x^{2j}\le \sum_{j\ge0}(k^2x^2)^j/j!$. Multiply the bounds for the finite products and pass to an almost everywhere converging subsequence by Fatou's lemma. With $M=k^{1/\sigma}$ this gives

$$
\log\int |Y|^{2k}\le 2k\sum_{p\le M}p^{-\sigma}+k^2\sum_{p>M}p^{-2\sigma}
=o_\sigma(k^{1/\sigma}).
$$

Here the elementary fact $\pi(t)=o(t)$ suffices. Apart from finitely many exceptions the primes avoid every prime divisor of a fixed primorial, so their upper density is at most $\prod_{p\le y}(1-1/p)$. Its reciprocal is at least $\sum_{n\le y}1/n$. Partial summation now gives $\sum_{p\le M}p^{-\sigma}=o(M^{1-\sigma})$ and $\sum_{p>M}p^{-2\sigma}=o(M^{1-2\sigma})$.

Since $|1-2^{1-\sigma}u_2|\ge 2^{1-\sigma}-1>0$, the same upper bound on the logarithm of the moment holds for $1/E$, with an additional $O_\sigma(k)$. Cauchy–Schwarz and $y^ke^{-y}\le(k/e)^k$, with $y=2X|E|^2$, give

```math
\sup_r |c_X(r)-\widehat T(r)|
\le \|T\|_2\left(\frac{k}{2eX}\right)^{k/2}
\left(\int |E|^{-2k}\right)^{1/2}.
```

Set $t=\log X$ and $m=1/(1-\sigma)$. Take $k=\lfloor A t^{m-1}\rfloor$ for fixed $A>0$. The logarithm of the right side is at most $-(A/2)t^m+o_{\sigma,A}(t^m)$. The pulse has bounded weight and $n/q<5/2$ on its support, giving $O(Q^2)$ ratios of reduced denominator at most $Q$. Choosing $A>4C$ makes the total error on the retained ratios tend to zero. For each $X$ this set is finite, so its limit in $K$ can be taken after summing it. Subtract it from the full sum at finite $`K,X`$ and use the full limit $P_f-2L-R_\sigma^*$ proved below. This gives the displayed tail identity with $K$ tending to infinity first. The argument covers any integer cutoff $2\le Q(X)\le Q_C(X)$ for all sufficiently large $X$. In particular every fixed $Q\ge2$ gives the same answer and a subsequent limit in $Q$ preserves it. Taking $Q$ to infinity at finite $`K,X`$ gives zero by absolute convergence.

## Every wall

At fixed $X$, taking $K$ to infinity in the full sum gives $I_\sigma(V_X)$ with $V_X=(\eta'/\eta)(1-e^{-X|\eta|^2})$ extended smoothly at zeros. The alternating series gives the required polynomial bounds for dominated convergence. Since $\partial_{\bar s}V_X=X|\eta'|^2e^{-X|\eta|^2}$, Cauchy–Green on $\sigma<\Re s<3$ gives

```math
I_3(V_X)-I_\sigma(V_X)=\frac X\pi\int_{\sigma<\Re s<3}
 H(s-1/2)|\eta'(s)|^2e^{-X|\eta(s)|^2}\,dA(s).
```

The horizontal integrals vanish for fixed $X$ by rapid decay of $H$ and $|(1-e^{-X|w|^2})/w|\le\sqrt X$. For the limit in $X$ it suffices to bound the absolute area integrals at large heights uniformly in $X\ge1$. On $|\eta|>1/2$ use $Xe^{-X|\eta|^2}\le4/e$ and the polynomial growth of $\eta'$. On $|\eta|\le1/2$ use the change of variables $w=\eta(s)$ counting preimages with multiplicity. In each horizontal band $j\le\Im s\le j+1$ of the strip, the number of preimages of any $|w|\le1/2$ is $O(\log(2+|j|))$.

For this bound apply Jensen's formula to $\eta(s)-w$ about $4+i(j+1/2)$ with inner radius $18/5$ and outer radius $15/4$. The smaller disc contains the band of the strip, and the larger lies in $\Re s\ge1/4$, where partial summation gives $|\eta(s)|\ll1+|s|$. At the centre $|\eta(s)-w|>3/8$ since $\sum_{n\ge2}n^{-4}<1/8$. In the change of variables the Gaussian $(X/\pi)e^{-X|w|^2}$ has total mass one. The absolute integral on the small-value part of each band is therefore bounded by a constant times $\log(2+|j|)\sup |H(s-1/2)|$ over that band. This is summable and its tails tend uniformly to zero.

On a bounded part of the strip, the density tends uniformly to zero on compact sets away from the zeros. At a zero $\rho$ of multiplicity $m_\rho$, take a holomorphic local coordinate $u$ with $\eta=u^{m_\rho}$. Then $|\eta'|^2dA(s)=m_\rho^2 |u|^{2m_\rho-2}dA(u)$, so the limiting mass is $m_\rho$, weighted by $H(\rho-1/2)$. At a zero on $\Re s=\sigma$ the boundary in the $u$ coordinate is a smooth curve with a tangent through zero. Its interior takes an angular fraction tending to $1/2$ on shrinking circles. The radial Gaussian therefore contributes half the multiplicity. These statements follow also by rescaling $u$ by $X^{-1/(2m_\rho)}$ in the local integral. [Local integral](prime-two-ratios.md#pole-subtraction).

The outer integral tends to $I_3(\eta'/\eta)=L$. The zeros of the factor $1-2^{1-s}$ on $\Re s=1$ give total weight $3L-P_f$ after omitting the zero at $s=1$ cancelled by the pole of $\zeta$. The other zeros give $`R_\sigma^*`$. Subtraction proves $`\lim_{X\to\infty}I_\sigma(V_X)=P_f-2L-R_\sigma^*`$ for every $1/2<\sigma<1$. This is a value of the prescribed regularised sum even on a wall passing through a zero. [Dyadic weights](prime-two-kernel-comparison.md#the-wall).

Absolute summability of the exterior weights also gives the jump in the large-denominator limit as the wall passes a real part $\beta>1/2$ inside the strip. Its right limit minus its left limit is exactly $\sum_{\Re\rho=\beta}m_\rho H(\rho-1/2)$, even if other real parts accumulate at $\beta$. This follows by taking the one-sided limits of each indicator in $`R_\sigma^*`$ under the sum.

On compact interior intervals the retained coefficients converge uniformly. For the full tail, a squared increment of the wall position recovers the sum of squared jumps in the original order. At the upper endpoint the same test gives $P_f/2$, halfway between its two one-sided limits. [Moving walls](prime-two-continuity.md).

## Pole subtraction

Set $a=\log2$, $A(s)=1-2^{1-s}$, $D=A'/A$ and $C_2(s)=1+(a^{-1}-1)2^{1-s}$. The constant Dirichlet coefficient of $C_2$ is one, and $aC_2(1)=1$ fixes the residue at the pole. Write $w=\eta_K(s)$ on the vertical line or $w=E_K$ on the prime circles, with $\dot w$ its Dirichlet derivative in either case. Put $U=1-\exp(-X|w|^2)$ and $\Phi=U/w$, extended by zero at $w=0$. Split the original regularised expression into

```math
M_{K,X}=U\big[D+C_2(\dot w-Dw)\big],\qquad
N_{K,X}=(1-C_2w)\big[\dot w\,\Phi-DU\big].
```

Since $w\Phi=U$ including at $w=0$, $`M_{K,X}+N_{K,X}=\dot w\,\Phi`$ identically. On the prime circles use $2^{1-s}=2^{1-\sigma}u_2$. The inner expansion of $D$ is an absolutely summable geometric series in $u_2^{-1}$. The two terms thus have absolutely summable Fourier coefficients at finite $`K,X`$, and their coefficients add separately at each original rational ratio. [Finite algebra in Lean](PrimeTwoRatios.lean).

On the vertical line, $M_{K,X}$ tends first to $(1-e^{-X|\eta|^2})(D+C_2A\zeta')$ because $\eta=A\zeta$. In the limit $X\to\infty$ its integral is $I_\sigma(D+C_2A\zeta')$. Both passages follow by domination using the bounded $U,D,C_2$, the polynomial bounds for $\eta_K,\eta_K'$ and the rapid decay of $H$. On an outer wall the integral of $C_2A\zeta'$ is $-L$ from the coefficient at $2$. Shifting this term inwards crosses only $s=1$, where $C_2 A\zeta'$ has residue $-aC_2(1)=-1$. The shift therefore adds $P_f$. The inner expansion of $D$ gives $-L$, so the full integral of $M$ tends to $P_f-2L$. [Contour signs](prime-two-eta.md#matching-the-pole).

On the prime circles the limit in $K$ is $(1-e^{-X|E|^2})M_\infty$, where $M_\infty=D+C_2 A\dot Z$ belongs to $L^2$. Indeed $D E=\dot A Z$, and the unweighted bracket converges in $L^2$. The heat factors are bounded and converge in probability, giving $L^1$ convergence of the products. In $C_2 A\dot Z$ only positive integer frequencies occur. Its coefficient at $2$ is $-a2^{-\sigma}$, since $\dot Z$ has zero constant coefficient and $C_2 A$ has constant coefficient one. Together with $D$ at $1/2$, the retained frequencies again give $-2L$. The moving cutoff estimate applies with $`\|M_\infty\|_2`$ in place of $`\|T\|_2`$. Subtracting the retained sum at finite $`K,X`$ gives the two limits below, where $`\mathcal S_{>Q}(V)`$ sums $`W_\sigma(r)\widehat V(r)`$ over the ratios of reduced denominator above $Q$.

$$
\lim_{X\to\infty}\lim_{K\to\infty}\mathcal S_{>Q_C(X)}(M_{K,X})=P_f,\qquad
\lim_{X\to\infty}\lim_{K\to\infty}\mathcal S_{>Q_C(X)}(N_{K,X})=-R_\sigma^* .
$$

On a polynomially separated wall the second limit also has a direct contour check. Its limiting integrand before weighting is $(1-C_2\eta)\zeta'/\zeta$. On the outer wall both factors have zero constant Dirichlet coefficient, so the first possible integer frequency is $4$, outside the pulse. At $s=1$ the multiplier vanishes because $\eta(1)=a$. At an exterior zero $\rho$ it is exactly one, giving the original residue $m_\rho H(\rho-1/2)$. Moving the contour inwards gives the minus sign. The retained coefficients have total limit zero because $M_\infty$ has the same two pulse coefficients as $T$. Replacing $C_2$ by any Dirichlet polynomial with the same two normalisations still gives multiplier one at every exterior zero. Taking a wall strictly between $1/2$ and $\sigma$ in the earlier absolute residue estimate gives $|R_\sigma^*|<10^{-20000}$, including the half weights. [Absolute residue bound](prime-two-eta.md#residue-bound).

At $\sigma=3/4$ numerical integration over $|t|\le4096$ gives $I_\sigma(N)$ on that finite interval as $0.186622950450$ at $X=1$, $-0.000511682834$ at $X=64$ and $0.000000041506$ at $X=512$. These are values after taking $K$ to infinity at fixed $X$. The integrals of $M$ and $N$ added to the original regularised integral in each numerical check. [Wall calculation with the pole split](prime-two-wall.py).

## First coefficient

Put $z=2\sigma$, $h(x)=x^{-1/2}g(\log x)$ on $x>0$ and extend $h$ by zero to $x\le0$. Its first two logarithmic moments are $`\int h(x)\,dx=P_f`$ and $`\int h(x)\log x\,dx=P_1=H'(1/2)`$. Write $B_Q(z)$ for the right derivative at $X=0$ of the contribution from $q\le Q$, after the limit in $K$ at fixed $X$. For a fixed ratio the derivative uses $\dot E\overline E$, since $`|\dot E\,\overline E|`$ is integrable and $(1-e^{-u})/u\le1$ for $u>0$. The right derivative of the full vertical integral is $J(\sigma)=I_\sigma(\eta'\overline\eta)$. The alternating series and its derivative have polynomial bounds on this wall, and $H$ decays rapidly. Thus $J(\sigma)-B_Q(z)$ is the first coefficient of the large-denominator contribution at fixed $Q$. [Vertical derivative](prime-two-eta.md#a-convergent-expansion).

At the ratio $n/q$ sum the pairs $(\ell n,\ell q)$ over all $\ell\ge1$. For $z>1$ the series with weight $\ell^{-z}(1+\log\ell)$ converges. The contribution to $B_Q$ is

$$
q^{-z}h(n/q)
\begin{cases}
 -\zeta(z)\log n+\zeta'(z),& n,q\text{ both odd},\\
 \eta(z)\log n-\eta'(z),& n,q\text{ of opposite parity}.
\end{cases}
$$

For two odd entries the alternating signs multiply to one for every $\ell$. For entries of opposite parity their product is $(-1)^\ell$. These are the only cases for a coprime pair.

## Coprime numerators

For any smooth compactly supported $h$ on the positive line, $\sum_{(n,q)=1}h(n/q)=\varphi(q)\int h+O_h(1)$ uniformly in $q$. Insert $\sum_{d\mid(n,q)}\mu(d)$ and put $M=q/d$. Poisson summation gives $\sum_j h(j/M)=M\int h+O_{h,N}(M^{-N})$ for any $N>1$. The errors summed over $M\mid q$ are bounded by a constant independent of $q$. For odd $q$ the even numerators give half the main term by applying the same argument to $h(2x)$. The odd numerators give the other half. All these statements also hold with $h(x)\log x$ in place of $h(x)$. [Divisor identities](https://dlmf.nist.gov/27.5.E4).

For even $q$ set $A_q(z)=\eta(z)$, and for odd $q$ set $A_q(z)=(\eta(z)-\zeta(z))/2=-2^{-z}\zeta(z)$. Define

$$
D_Q(z)=\sum_{q\le Q}\varphi(q)q^{-z}A_q(z).
$$

Splitting $\log n=\log q+\log(n/q)$ in the pair formula gives $`B_Q=P_1D_Q-P_fD'_Q+\sum_{q\le Q}e_q`$. The prime on $D_Q$ differentiates with respect to $z$. The finite pair formula defines $B_Q$ also for complex $z$ with $\Re z>1$. The error for each denominator satisfies $|e_q(z)|\ll q^{-\Re z}(1+\log q)$ locally uniformly there, by the two coprime counts. In particular $\sum_q e_q$ is holomorphic on this half-plane.

## Dyadic boundary

Write $q=2^j d$ with $d$ odd and $v=2^{1-z}$. After factoring out $\zeta(z)\varphi(d)d^{-z}/2$, the odd denominator gives $-v$ and each even denominator gives $(1-v)v^j$ for $j\ge1$, since $\varphi(2^j d)=2^{j-1}\varphi(d)$. With $J_d=\lfloor\log_2(Q/d)\rfloor$ the sum through the cutoff is $-v^{J_d+1}$. The infinite sum along each such sequence is zero for $\Re z>1$. In particular

```math
D_Q(z)=-\frac{\zeta(z)}2\sum_{\substack{d\le Q\\ d\ {\rm odd}}}
\varphi(d)d^{-z}v^{J_d+1}.
```

For $\Re z>2$ the full sum over $q$ and its derivative converge absolutely, so the infinite chains give $\sum_{q\ge1}\varphi(q)q^{-z}A_q(z)=0$ there. In the same half-plane the original pair sum converges absolutely and is $J(z/2)$. Thus $\sum_q e_q=J(z/2)$. On $\Re z>1$ the integral defining $J(z/2)$ uses $\eta'(z/2+it)\eta(z/2-it)$ and $H(z/2-1/2+it)$. It is holomorphic by partial summation of the alternating series and rapid decay of the pulse transform, locally uniformly in $z$. The identity theorem extends the equality to this half-plane. For real $1<z<2$ it follows that

$$
B_Q(z)=J(z/2)+P_1D_Q(z)-P_fD'_Q(z)+O_{z,f}\big(Q^{1-z}\log(2Q)\big).
$$

The finite chain identity and its derivative with $dv/dz=-(\log 2)v$ are checked in Lean. [Lean source](PrimeTwoRatios.lean).

## Boundary size

Define $b(z)=(1-2^{2-z})\zeta(z)/(3\zeta(2)(2-z))$, taking the removable value at $z=2$. For real $1<z<2$ the previous formula has the explicit expansion

$$
B_Q(z)=J(z/2)+Q^{2-z}\big[P_f b(z)\log Q+P_1 b(z)-P_f b'(z)\big]
+O_{z,f}\big(Q^{1-z}\log^2(2Q)\big).
$$

Indeed $\sum_{q\le y,\ q\ {\rm odd}}\varphi(q)=y^2/(3\zeta(2))+O(y\log(2y))$. Insert $\varphi(q)=\sum_{d\mid q}\mu(d)q/d$ and use $\sum_{k\le t,\ k\ {\rm odd}} k=t^2/4+O(t)$ for $t\ge1$ and $\sum_{d\ {\rm odd}}\mu(d)/d^2=4/(3\zeta(2))$. The same calculation over all integers gives $y^2/(2\zeta(2))+O(y\log(2y))$. Hence with $C(z)=(2-z)b(z)$ we have $\sum_{q\le y} A_q(z)\varphi(q)=C(z)y^2/2+E(y,z)$, where $E$ and its $z$ derivative are $O(y\log(2y))$ on compact subsets of $\Re z>1$.

For $\Re z>2$ use the zero sum over all denominators and apply partial summation to the tail. This gives

```math
D_Q(z)=b(z)Q^{2-z}+Q^{-z}E(Q,z)-z\int_Q^\infty E(y,z)y^{-z-1}\,dy .
```

The error integral converges locally uniformly on $\Re z>1$, giving the same formula throughout that half-plane. On differentiating, the last two terms contribute $O_{z,f}(Q^{1-z}\log^2(2Q))$ for real $z>1$. Substitution gives the expansion. Since $P_f>0$ and $b(z)<0$ for $1<z<2$, $B_Q$ tends to $-\infty$ and the complementary first coefficient tends to $+\infty$. Thus the pole weight occurs in the first coefficient with a divergent multiplier. Summing the heat factor at fixed $`K,X`$ still leaves the question whether its finite limiting contribution is exactly $P_f$.

At $\sigma=3/4$, $P_f b(3/2)\simeq-0.006038988168$. The table uses the original pulse. The last column sums the pairs by their unreduced denominator $m\le Q$, including every numerator in the support. Its inner alternating sum is $O_{f,N}(m^{-N}(1+\log m))$ for each $N$, by Poisson summation on the two parities. It therefore gives an independent series for $J$ by absolute convergence of the grouped rows and continuation from $z>2$. The displayed numbers are numerical checks of the first coefficient. [Reproduce the table](prime-two-denominators.py).

| $Q$ | $B_Q$ | $B_Q-(P_1D_Q-P_fD'_Q)$ | Unreduced sum |
| ---: | ---: | ---: | ---: |
| 512 | $-0.571259$ | $0.467076$ | $0.451295$ |
| 2048 | $-2.001050$ | $0.454138$ | $0.451293$ |
| 8192 | $-5.217237$ | $0.451331$ | $0.451293$ |
