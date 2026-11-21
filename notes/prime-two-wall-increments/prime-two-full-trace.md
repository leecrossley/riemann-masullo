# Prime 2 limits

## Global projection

Put $u=\log|x|$, $a=1/2$, and use Sobolev parameter four. The smooth test $f$ is supported in $\log(3/2)<u<\log(5/2)$, has $f(\log2)=1$, and satisfies $`\int (e^{u/2}-e^{-u/2}/(e^{2u}-1))f(u)\,du=0`$. These conditions isolate the prime 2 and cancel the real-place contribution against the two poles. Write $L=\log(2)/\sqrt2$ and $`M_f(\rho)=\int f(u)e^{(\rho-1/2)u}\,du`$. The earlier note constructs this $f$ and proves the trace formula used here for almost every $\alpha$ in the stated range. These widths avoid zeros on the strip walls and have the polynomial separation used in the density proof. [Test, density proof and bound](prime-two-kernel.md).

On $`X_\alpha=L^2(\mathbb R,e^{2\alpha|u|}(1+u^2)^2\,du)`$, let $W_\alpha$ be the closed global theta space and $\Pi_\alpha$ its orthogonal projection. Let $`A_f h(u)=\int f(t)h(u-t)\,dt`$ and let $M_{\kappa_\Lambda}$ multiply by the smooth norm cutoff. Normalise the compact norm-one quotient to have Haar mass one. The norm-dependent test annihilates the other angular characters, so the trace on this space accounts for the full global trace. [Masullo, Definition 10.4.1.5](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=457).

For a finite set $P$ of excluded odd primes, put $S_t h(u)=h(u+t)$ and $D_p=I-p^{-1/2}S_{\log p}$. Removing multiples of $p$ from a theta sum applies $D_p$ to the half-weighted function. For example, $e^{u/2}\sum_{n\ge1}\Phi(pne^u)=p^{-1/2}S_{\log p}(e^{u/2}\sum_{n\ge1}\Phi(ne^u))$. The relevant operator bound is

```math
\|(p^{-1/2}S_{\log p})^n\|\ \le\ p^{-n(1/2-\alpha)}(1+n\log p)^2 .
```

Indeed $\sqrt{1+(x+t)^2}\le(1+|t|)\sqrt{1+x^2}$, and the exponential part of the norm contributes at most $e^{\alpha|t|}$. The bound is summable in $n$ because $\alpha<1/2$. Hence $D_p^{-1}=\sum_{n\ge0}p^{-n/2}S_{n\log p}$ in operator norm. Translations preserve $W_\alpha$, so closure shows that both $D_p$ and its inverse preserve it. Thus $(\prod_{p\in P}D_p)W_\alpha=W_\alpha$ for every finite $P$. At the widths above the translates of the Gaussian theta generator in the earlier note already span $W_\alpha$ densely, so using just the finite-prime versions of that generator gives the same projection.

Let $\Pi_{\alpha,P}$ project onto that closed space, and define $`T(P,\Lambda)=\mathrm{Tr}(M_{\kappa_\Lambda}A_f\Pi_{\alpha,P})`$. For each finite pair, $T(P,\Lambda)=T(\varnothing,\Lambda)$. The earlier trace-class estimate on $A_f(I-\Pi_\alpha)$ and the ordinary convolution diagonal $f(0)=0$ therefore give

$$
\lim_{P}\lim_{\Lambda\to\infty}T(P,\Lambda)
=\lim_{\Lambda\to\infty}\lim_P T(P,\Lambda)
=-\sum_{|\Re\rho-1/2|<\alpha}m_\rho M_f(\rho)
=L+R_\alpha,\qquad
R_\alpha=\sum_{|\Re\rho-1/2|\ge\alpha}m_\rho M_f(\rho).
$$

Here $P$ exhausts the odd primes, $\rho$ ranges over distinct nontrivial zeta zeros and $m_\rho$ counts multiplicity. For this $f$, $|R_\alpha|<6.5\times10^{-65}$ follows from the certified height $3\times10^{12}$, the test’s eighth-power decay and the unconditional zero count. [Trace-class argument](strip-resolvent.md), [Platt and Trudgian](https://arxiv.org/abs/2004.09765) and [Bellotti and Wong](https://arxiv.org/abs/2412.15470).

## Limiting generator

One can also evaluate the full projection obtained by *first* taking the pointwise limit of the generators. This changes the space on which the trace is taken. With $`\Phi(x)=(x\,d/dx)(x\,d/dx+1)e^{-\pi x^2}`$, the remaining generator is $g_{(2)}(u)=e^{u/2}\sum_{j\ge0}\Phi(2^je^u)$. Its translates span all of $X_\alpha$ densely for every $0<\alpha<1/2$.

To prove this, set $H_0(u)=e^{u/2}e^{-\pi e^{2u}}$. The function $D_2 g_{(2)}=e^{u/2}\Phi(e^u)$ equals $(\partial_u^2-1/4)H_0$. The dyadic series converges in $X_\alpha$ by the same summable bound used for $D_p^{-1}$. If a bounded linear functional $\ell$ annihilates all translates of $g_{(2)}$, it annihilates all translates of $D_2g_{(2)}$. Then $q(t)=\ell(S_t H_0)$ satisfies $q''=q/4$ because $H_0$ and its first two derivatives belong to $X_\alpha$ and translation is strongly continuous. The translation bound gives $|q(t)|\le C e^{\alpha|t|}(1+|t|)^2$. Since $\alpha<1/2$, both coefficients in $q(t)=c_+e^{t/2}+c_-e^{-t/2}$ vanish.

Represent $\ell$ by integration against $`b(u)\,du`$, with $b$ square integrable for the inverse weight. Putting $x=e^{2u}$ in $\ell(S_t H_0)=0$ shows that the Laplace transform of $b((\log x)/2)x^{-3/4}$ vanishes at every positive argument $\pi e^{2t}$. Cauchy–Schwarz gives absolute convergence at each argument. Uniqueness of the Laplace transform now gives $b=0$ almost everywhere. Equivalently, weight this density by $e^{-x}$ to make a finite measure and push it to $[0,1]$ by $x\mapsto e^{-x}$. All its polynomial moments vanish, so the measure vanishes. This proves density.

The projection obtained from this limiting generator is the identity, with cutoff trace $f(0)\int\kappa_\Lambda=0$. At every finite $P$ the projection onto the closed global space is $\Pi_\alpha$. Forming the projection after the pointwise limit therefore changes the limiting full trace by $L+R_\alpha$. The local arithmetic kernel in §12.4 has a further additive Fourier cutoff.

## Lifted kernel

Masullo’s (11.3.1.4) cuts the second variable of a lifted kernel. We can also compare the two orders for this operation on the integrable global trace diagonal. Take the rational fundamental domain $D=\mathbb R_{>0}\times\prod_p\mathbb Z_p^\times$ for $\mathbb Q^\times$ acting on the ideles. Let $b\ge0$ be a fixed lift weight with $\sum_{r\in\mathbb Q^\times}b(rx)=1$, as in §1.3.1. At the orbital diagonal $(x,\xi x)$ its weight is

$$
c_b(x,\xi)=\sum_{r\in\mathbb Q^\times}b(rx)b(r\xi x).
$$

For each $x$, $\sum_\xi c_b(x,\xi)=1$ by summing nonnegative terms. For $x\in D$, the requirement that $\xi x$ be a unit at $p$ amounts to $v_p(\xi)=0$. The retained weight $w_P(x)=\sum_{\xi:\ v_p(\xi)=0\ (p\in P)}c_b(x,\xi)$ decreases to the sum over $\xi=\pm2^j$, $j\in\mathbb Z$, and $0\le w_P\le1$. [Masullo, lifted kernel and orbital diagonal](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=24).

Write $d_\alpha(u)$ for the geometric trace density. It is integrable. To see this directly, conjugate $-A_f(I-\Pi_\alpha)$ to unweighted $L^2$ and take a singular-value expansion $\sum s_n a_n\otimes\overline{b_n}$. Then $d_\alpha=\sum s_n a_n\overline{b_n}$ in $L^1$, since the singular values are summable and both vectors have norm one. This gives $`\mathrm{Tr}(M_h[-A_f(I-\Pi_\alpha)])=\int h\,d_\alpha`$ for bounded multipliers $h$. The same expansion gives the radially averaged kernel diagonal because local interval averages of the two $\ell^2$-valued functions $(\sqrt{s_n}a_n)_n$ and $(\sqrt{s_n}b_n)_n$ converge to their values almost everywhere by Lebesgue differentiation. The conjugating weights are continuous and cancel on the diagonal. The ordinary convolution has zero diagonal. Choose its compact lift $1_D(z)f(\log|z|)$, which also vanishes at each rational $z$ because $|z|=1$.

With $`d^\times x=du\,dk`$ on $D$, where $dk$ is Haar probability on the finite units, the retained orbital mass is $`\int_D\kappa_\Lambda(u)d_\alpha(u)w_P(x)\,d^\times x`$. In either order, and along any joint exhaustion, dominated convergence gives $`\int_D d_\alpha(u)w_\infty(x)\,d^\times x`$. With $b=1_D$ only $\xi=1$ contributes, so $w_P=1$ at every stage and both orders give $L+R_\alpha$. With $b=(1_D+1_{3D})/2$ the weights are $1/2$ at $\xi=1$ and $1/4$ each at $3$ and $1/3$. Once 3 is excluded, the retained part and the removed part each have half the full trace. Adding them restores the same total. These two choices show exactly what a normalised lift can establish about individual shells. The examples and the separate defect in the printed lift were identified in issue 1358. [Issue 1358, rational lifts](https://github.com/mccorvie/riemann-masullo/issues/1358#issuecomment-7054088460).

Taken literally, (11.3.2.2) averages the indicator of the norm section $\tau(\mathbb R_{>0})$. For fixed $x$ and $\xi$, $\xi gx$ lying in that section forces $g=\xi^{-1}\tau(|x|)x^{-1}$. A singleton has Haar measure zero in the infinite profinite unit group. That printed weight therefore vanishes. The calculation above uses a normalised lift and compares its entire diagonal. Identifying its retained piece with the local Fourier calculation is a further step. [Masullo, (11.3.2.2)](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=495) and [issue 1358](https://github.com/mccorvie/riemann-masullo/issues/1358).

## Local cutoffs

The footnote on p. 514 asks what happens if the physical norm window has a lower endpoint. Use $1/\Lambda\le |x|_S\le\Lambda$ with $\Lambda\ge1$ and the additive Fourier cutoff $|\omega|_S\le B$. Denote the resulting complete unfolded orbital integral by $U(\Lambda,B)$, divided by the norm-one volume as on p. 514. All $S$-units $\xi=\pm2^j$ enter this calculation. The change of variables $w=x\omega$ gives the exact identity

$$
U(\Lambda,B)=J(\Lambda B)-J(B/\Lambda).
$$

Here $J(t)$ is the one-sided integral with logarithmic weight $`[\log(t/|w|_S)]_+`$. The two terms give the length of the allowed interval of $\log|x|_S$ at fixed $w$. Summing over all the orbits gives the endpoint values

```math
J(0+)=0,\qquad J(+\infty)=L+A_\infty,\qquad
A_\infty=\int_{\mathbb R}\frac{e^{-u/2}}{1-e^{-2u}}f(u)\,du .
```

For this test $A_\infty$ equals the pole contribution $`P_f=\int (e^{u/2}+e^{-u/2})f(u)\,du`$ by construction. Its positivity follows because $f$ is a bump minus a smaller translate to the right and the real-place weight is positive and strictly decreasing on its support. Thus $P_f\simeq0.0137703202700$ and $L+P_f\simeq0.5038993920043$. It follows that

$$
\boxed{\lim_{B\to\infty}\lim_{\Lambda\to\infty}U(\Lambda,B)=L+P_f,\qquad
       \lim_{\Lambda\to\infty}\lim_{B\to\infty}U(\Lambda,B)=0.}
$$

In the first order the two arguments of $J$ go to infinity and zero. In the second order both go to infinity at fixed $\Lambda$, so their limits cancel. Subtracting the same $P_f$ gives $L$ and $-P_f$. In the printed one-sided integral only the term $J(\Lambda B)$ occurs. Keeping the lower endpoint introduces the extra term, whose size depends on $B/\Lambda$. Along a joint limit with $B/\Lambda\to0$ one recovers $L+P_f$, while $B/\Lambda\to\infty$ gives zero. For a fixed positive limiting ratio $c$, the answer is $L+P_f-J(c)$ because the orbital series below converges locally uniformly on positive scales. This distinction also holds for cutoffs smoothed by averaging the endpoints over bounded log shifts, by linearity and the uniform bounds below. The sharp multiplier and the dependence on the ratio were previously analysed in the project’s chapter audit. [Masullo, footnote 2](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514), [equations (12.4.1.22–24)](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=518) and [issue 2469, item 14](https://github.com/mccorvie/riemann-masullo/issues/2469#issuecomment-7057757662).
