# Squared increments

Two contributions to the squared wall increment vanish after division by $h$ for the original prime 2 pulse. One comes from equal rational ratios, including arbitrarily large denominators. The other comes from distinct ratios with at least one denominator below the growing cutoff derived from (3). I could not establish cancellation for pairs of distinct ratios with both denominators beyond that cutoff. The argument uses the finite Dirichlet sums, takes the limit in their length first, removes the heat parameter second and sends the wall separation to zero last. [Masullo, §12.3](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514) and [original pulse](prime-two-kernel.md#the-test).

## One coefficient

Let $I=[d,e]$ with $1/2<d<e<1$. Take $0\le\omega\in C_c((d,e))$, extend it by zero and use $h>0$ small enough that both walls stay in $I$ on its support. Write $\Delta_h A(s)=A(s+h/2)-A(s-h/2)$. For fixed $Q\ge2$, $F_{K,X}$ sums $a_{K,X}(s,r)$ over positive rationals $r$ with reduced denominator $q(r)>Q$. On the product of the prime unit circles with Haar probability measure put

$$
E_K=\sum_{n\le K}(-1)^{n-1}n^{-s}z_n,\qquad
a_{K,X}(s,r)=W_s(r)\widehat{\dot E_K\Phi_X(E_K)}(r),\qquad
W_s(r)=r^{s-1/2}g(\log r).
$$

Here $s$ is the real wall position, $z_n=\prod_p z_p^{v_p(n)}$, the dot differentiates in $s$, $g(u)=f(u)+f(-u)$ and $\Phi_X(w)=(1-e^{-X|w|^2})/w$, extended by zero. The hat extracts the coefficient of $z_r$. The Fourier series at finite $K,X$ is absolutely convergent by its exponential expansion, also locally after differentiation. The value of the original sum uses all phases equal to one after the coefficients have been multiplied by $W_s$. [Rational coefficients](prime-two-ratios.md#one-ratio) and [Hedenmalm, Lindqvist and Seip](https://arxiv.org/abs/math/9512211).

Rotate the prime phases by $z_p\mapsto z_p p^{-it}$. Then $E_K$ is holomorphic in $s+it$ and $i\partial_t z_r(t)=(\log r)z_r(t)$. The derivative of $W_s$ supplies exactly this logarithm. Applying $\partial_s+i\partial_t$ to $\dot E_K\Phi_X(E_K)$ gives the identity

$$
\partial_s a_{K,X}(s,r)=2W_s(r)\widehat{B_{K,X}(s)}(r),\qquad
B_{K,X}=X|\dot E_K|^2e^{-X|E_K|^2}.
\tag{1}
$$

Indeed $\partial_{\bar w}\Phi_X(w)=Xe^{-X|w|^2}$, including at zero. Equation (1) is coefficient by coefficient and holds before taking either limit.

## Prime phases

At each fixed $X$, $B_{K,X}\to B_X=X|\dot E|^2e^{-X|E|^2}$ in $L^2$ of the prime phases, uniformly in $s\in I$. To see the required moment convergence, raise a finite Dirichlet tail or its derivative to an integer power $j$ and use orthogonality. The squared norm of that power is at most the tail over $n>K$ of $d_j(n)^2(1+\log n)^{2j}n^{-2d}$, where $d_j(n)$ counts ordered factorizations into $j$ positive integers. The majorant is summable by the Euler product for $\sum d_j(n)^2 n^{-2d'}$ at any $1/2<d'<d$. Thus $E_K$ and $\dot E_K$ converge in each fixed even moment uniformly on $I$. Use the $L^4$ bounds for the difference of the squared derivatives, the $L^8$ bound for $\dot E$, and the Lipschitz bound $O(\sqrt X)$ for $w\mapsto e^{-X|w|^2}$ to obtain the asserted $L^2$ convergence.

There is a quantitative bound for $B_X$. For $X$ sufficiently large set $t=\log X$ and $S_d(X)=t^{1/(1-d)}(\log t)^{d/(1-d)}$. There is $c_I>0$ such that

$$
\sup_{s\in I}\|B_X(s)\|_2\le e^{-c_I S_d(X)}.
\tag{2}
$$

In the Euler product representation $E=(1-2^{1-s}z_2)\prod_p(1-p^{-s}z_p)^{-1}$ almost everywhere. The dyadic factor has modulus at least $2^{1-e}-1$. The mean of $|1-p^{-s}z_p|^{2k}$ is at most $\exp(\min(2kp^{-d},k^2p^{-2d}))$. Multiplication in finite products and Fatou's lemma bound the inverse moment of $E$. Splitting the primes at $k^{1/d}$ and applying $\pi(y)\ll y/\log y$ by partial summation gives, uniformly on $I$ for integers $k\ge2$,

$$
\log\int |E|^{-2k}\,dz\le C_I\left(k+\frac{k^{1/d}}{\log(2k)}\right).
$$

The prime bound here follows from $\prod_{m<p\le2m}p\le {2m\choose m}\le4^m$ and summation over dyadic intervals. Hölder's inequality and $y^ke^{-y}\le(k/e)^k$ give

$$
\|B_X\|_2\le X\|\dot E\|_8^2
\left[\left(\frac{k}{4eX}\right)^k\int |E|^{-2k}\,dz\right]^{1/4}.
$$

Choose $k=\lfloor\alpha(t\log t)^{d/(1-d)}\rfloor$. In the logarithm of the bracket the negative term $-kt$ is asymptotic to $-\alpha S_d(X)$. The inverse moment costs $O_I(\alpha^{1/d}S_d(X))$ for fixed $\alpha$, and $k\log k=o(S_d(X))$. Taking $\alpha>0$ sufficiently small proves (2). [Euler product and moment bound](prime-two-denominators.md#moving-cutoff).

## Distinct ratios

Integrate (1) across the wall interval. Parseval's identity and the triangle inequality in $\ell^2$ give

$$
\lim_{K\to\infty}\left(\sum_{r:q(r)>Q}|\Delta_h a_{K,X}(s,r)|^2\right)^{1/2}
\le 2h\sup_{v\in I,r}|W_v(r)|\sup_{v\in I}\|B_X(v)\|_2 .
\tag{3}
$$

The convergence on the left is uniform for the allowed centres by (1) and uniform $L^2$ convergence at fixed $X$. In particular, after integration against $\omega$ and division by $h$, all equal-ratio terms together are $O_{I,\omega}(h e^{-2c_I S_d(X)})$. This includes the entire infinite denominator tail.

At most $O(R^2)$ ratios in the pulse have reduced denominator at most $R$, since $r<5/2$ on its support. Cauchy–Schwarz in (3) therefore bounds the absolute increment summed over $Q<q(r)\le R$ by $O_I(Rh e^{-c_I S_d(X)})$ after the first limit. Choose $R(X)=\lfloor\exp(c_I S_d(X)/2)\rfloor$. For each fixed $h$ this block has vanishing increment as $X\to\infty$. The full tail $F_X=\lim_K F_{K,X}$ is bounded uniformly on $I$ by the vertical area estimate and the finitely many retained coefficients. Expanding the difference of the two squares shows that the combined contribution of all pairs containing this block vanishes after integration in the prescribed order. This also holds for smooth or bounded measurable compactly supported weights at this stage. [Uniform area estimate](prime-two-denominators.md#every-wall).

More precisely, write $J_b=\sum_{\Re\rho=b}m_\rho H(\rho-1/2)$ for the signed weight of the zeta zeros at real part $b$, with $H(z)=\int g(u)e^{zu}\,du$. The remaining arithmetic expression satisfies

$$
\lim_{h\downarrow0}\lim_{X\to\infty}\lim_{K\to\infty}
\frac1h\int\omega(s)
\sum_{\substack{r\ne r'\\q(r),q(r')>R(X)}}\Delta_h a_{K,X}(s,r)\,
\overline{\Delta_h a_{K,X}(s,r')}\,ds
=\sum_b\omega(b)|J_b|^2 .
\tag{4}
$$

All rational sums in (4) are taken at finite $K,X$. The inner limit for the cross terms follows by subtracting the equal-ratio sum from the square of the tail. On the zero side the increment is $\sum_bJ_b\mathbf1_{|s-b|<h/2}$ for almost every $s$. Two distinct centres eventually have disjoint intervals. For equal centres the normalised overlap tends to $\omega(b)$, and absolute summability justifies the limit over pairs even when the centres accumulate. Lean checks the normalised weighted interval integrals and this pair limit for any summable family of real jumps at distinct locations. It also checks that the limit vanishes precisely when every jump at a point where the weight is positive vanishes. The grouped $J_b$ are real by conjugation. [Ordered energy identity](prime-two-continuity.md#squared-increments) and [Lean proof](PrimeTwoIncrements.lean).

If some $J_b$ is nonzero, a smooth nonnegative weight positive at $b$ gives a strictly positive coefficient in (4). If every $J_b$ vanishes, $F=\lim_X F_X$ is constant on the open strip and every increment vanishes. The same criterion thus decides the assertion for all bounded measurable nonnegative weights of compact support, since they include the smooth tests.

At $h=1/\log X$ one can remove the signed contribution of all distinct pairs whose prime exponents differ at some odd prime up to $(\log X)^2$. The proof averages the original coefficients at those primes. [Large ratio pairs](prime-two-pair-sieve.md).

## Wall averaging

The wall integral itself can be carried out in the arithmetic expansion. To see which pairs it suppresses, let $d_{K,j+1}(n)$ be the signed coefficient of $n^{-s}z_n$ in $E_K^{j+1}$. Since $\dot E_KE_K^j=(E_K^{j+1})'/(j+1)$, expanding (1) at finite $K,X$ and integrating in $s$ gives

$$
\Delta_h a_{K,X}(s,r)=\frac{2X g(\log r)}{\sqrt r}
\sum_{j\ge0}\frac{(-X)^j}{j!(j+1)^2}
\sum_{n/m=r}d_{K,j+1}(n)d_{K,j+1}(m)
(\log n)\sinh(h\log m)\,m^{-2s}.
$$

Indeed $r^{s-1/2}(nm)^{-s}=r^{-1/2}m^{-2s}$, and the integral of $m^{-2v}\log m$ across the walls is $m^{-2s}\sinh(h\log m)$. This holds also at $m=1$, when both sides vanish. At finite $K$ the exponential series and its derivatives converge in the sum of absolute Fourier coefficients uniformly for $s\in I$. Multiplying by the conjugate expansion replaces the two factors depending on the wall centre by

$$
\int\omega(s)(mm')^{-2s}\,ds .
$$

This is positive for nonzero $\omega\ge0$. It depends only on the product of the two unreduced denominators, including when the two reduced ratios differ. There is no oscillation in $s$ from their separation.

In fact a positive subseries already diverges as $K\to\infty$ at the $X^2$ Taylor coefficient of the squared increment. Choose $\varepsilon>0$ so that $g(\log r)=1$ for $2<r<2+\varepsilon$. Take odd primes $p>Q$ and odd $n$ in that interval times $p$, choosing $\varepsilon<1$. The ratio $n/p$ is reduced, and every multiple of $(n,p)$ has positive sign in the $j=0$ term because $n+p$ is even. For $n\le K$ the single pair $(n,p)$ contributes at least a constant times $h(\log p)^2p^{-2e}$ to the coefficient of $X$, uniformly for the allowed centres. Each half of the ratio interval supplies $\gg p$ such numerators for large $p$. Pair a ratio from the lower half with one from the upper half, so all the pairs are distinct-ratio terms with positive product. Each half has lower bound proportional to $h\sum_{Q<p\le K/(2+\varepsilon)}p^{1-2e}(\log p)^2$, which diverges since $e<1$ and $\sum_p1/p$ diverges. Integrating the product against a nonzero $\omega$ preserves divergence. Thus taking absolute values in this Taylor coefficient loses cancellations even after grouping all multiples of each ratio. Equations (1) to (4) keep the heat exponential intact throughout their limits.

## Point evaluation

One might try to apply positivity of $B_X$ to bound the sum at the original phases by its $L^2$ norm. There is an explicit obstruction to this estimate even with the original pulse. Choose $M$ pairwise coprime odd integers $n_j=1+(M^2+j)2M!$, $1\le j\le M$. A common prime divisor would divide $2M!(i-j)$ and hence $2M!$, whereas $n_j\equiv1\pmod{2M!}$. The characters $z_{n_j}$ are independent uniform phases and are independent of $z_2$. The nonnegative polynomial

$$
D_M=|1+z_2|^2\left|\frac1M\sum_{j=1}^M z_{n_j}\right|^2
$$

has squared $L^2$ norm $6(2M^2-M)/M^4$, by matching the character indices in the fourth moment. The ratios $n_i/n_j$ tend to one uniformly as $M\to\infty$. The pulse equals one near both $\log 2$ and $-\log 2$ and vanishes near zero. For $i\ne j$ it therefore selects the frequencies $2n_i/n_j$ and $n_i/(2n_j)$, each with coefficient $1/M^2$ and reduced denominator $n_j$ or $2n_j$. At the original phases their weighted sum above any fixed cutoff tends to $2^{s-1/2}+2^{1/2-s}$. This example concerns the proposed norm inequality. For the actual $B_X$, any further cancellation must use more of the coefficients than their positivity and size in (2).

## Two scales

The interval calculation also gives a way to test (4) while $h$ and $X$ move together. Keep the limit in $K$ first. For the actual $F_X$, the choice $h_X=1/\log X$ gives

$$
\lim_{X\to\infty}\frac1{h_X}\int\omega|\Delta_{h_X}F_X|^2
=\sum_b\omega(b)|J_b|^2.
\tag{5}
$$

For $h_X=o(X^{-1/2})$ the same expression tends to zero unconditionally. Thus a calculation on this finer scale can return zero even before the cancellation in (4) has been decided.

To prove the two limits, use (1) for the finite set $q\le Q$. Its increment after the first limit is $O_{I,Q}(h\|B_X\|_{2,I})$, where the norm denotes the supremum in (2). For the full vertical integral the increment is the mass between the walls of $d\mu_X(z)=(X/\pi)H(z-1/2)|\eta'(z)|^2e^{-X|\eta(z)|^2}\,dA(z)$, with $\eta(z)=(1-2^{1-z})\zeta(z)$. Project this measure onto $x=\Re z$ on $I$. The kernel of the squared increment after division by $h$ is $h^{-1}\int_{|s-x|<h/2,\ |s-y|<h/2}\omega(s)\,ds$. It differs from $\omega(x)(1-|x-y|/h)_+$ by at most the modulus of continuity of $\omega$ at $h$. The total variations on $I$ are uniformly bounded and their tails at large heights tend uniformly to zero by the Jensen preimage count. Consequently one can discard large heights in the double integral at a cost tending to zero independently of $h,X$. This bound also makes the normalised $L^2(\omega)$ norm of the full increment uniformly bounded. The retained finite set contributes $o(1)$ in that norm. [Area bound](prime-two-denominators.md#every-wall).

On a bounded rectangle take the horizontal edges away from zeros and the vertical edges outside the support of every wall interval, moving those vertical edges slightly to avoid zeros. Away from finitely many zero discs the density tends uniformly to zero. At a zero $\rho$ of multiplicity $m$, in the variable $u=d_\rho X^{1/(2m)}(z-\rho)$, the density converges in $L^1(dA(u))$ to $m H(\rho-1/2)$ times $(m/\pi)|u|^{2m-2}e^{-|u|^{2m}}$, where $d_\rho=|\eta^{(m)}(\rho)/m!|^{1/m}$. This follows by Taylor expansion with an integrable majorant on the rescaled disc. For $h_X=1/\log X$ each of these horizontal widths is $o(h_X)$. The overlap tends to $\omega(b)$ for a pair at the same real part $b$ and to zero for distinct real parts. This proves (5) after restoring the height tail. For $h_X=o(X^{-1/2})$, every width is much larger than $h_X$. For two centres with the same real part, divide the horizontal difference by the larger width. Along a subsequence the ratio of the smaller width to it converges in $[0,1]$. The limiting difference has a continuous distribution because the two profile coordinates are independent and the wider one has a density. Its probability of lying in a window of width $o(1)$ tends to zero. Pairs at distinct centres still give zero. [Local profiles](prime-two-continuity.md#thin-strips).

One can also use an explicit finite Dirichlet cutoff in (5). Uniformly on $I$, for $K\ge2$ and $X\ge1$,

$$
\sup_I|F_{K,X}-F_X|\ll_{I,Q} K^{1/2-d}\big(\sqrt X\log(2K)+X\big).
$$

For $q\le Q$ use $\|E_K-E\|_2\ll_I K^{1/2-d}$, $\|\dot E_K-\dot E\|_2\ll_I K^{1/2-d}\log(2K)$, $|\Phi_X|\le\sqrt X$ and its Lipschitz bound $2X$ in the coefficient integrals. For the full vertical integral, partial summation in the alternating series bounds the two respective tails by $O_I((1+|t|)K^{-d})$ and $O_I((1+|t|)K^{-d}\log(2K))$. Integrate against the rapidly decreasing $H$. The difference of the normalised increment norms is therefore at most $2\|\omega\|_1^{1/2}h^{-1/2}\sup_I|F_{K,X}-F_X|$. In particular $X_n=e^n$, $h_n=1/n$ and $K_n=\lceil e^{n^2}\rceil$ give exactly the ordered coefficient in (5) as a single limit of the original $F_{K_n,X_n}$. Its value still depends on the sums $J_b$. [Finite regulariser bounds](prime-two-continuity.md#interior-walls).

## Finite check

With the actual alternating polynomial through $K=4$, the two prime circles give a numerical check of (1). On a $144\times144$ grid at $X=12$, integrating the right side from $0.73$ to $0.77$ agreed with the endpoint difference in each sampled pulse coefficient to within $2.2\times10^{-17}$. At $h=0.04$, using $\omega(s)=(1-((s-0.75)/0.10)^2)^2$ on $[0.65,0.85]$ and zero elsewhere, the normalised energy was $0.0008937331$. Equal ratios gave $0.0004641141$ and distinct ratios gave $0.0004296190$. At $X=1/4$ direct multiplication of the signed Dirichlet polynomials through Taylor order $j=22$ agreed with the endpoint increments to within $8\times10^{-19}$ on that grid. The vanishing required in (4) remains unproved. [Numerical check](prime-two-increments.py).
