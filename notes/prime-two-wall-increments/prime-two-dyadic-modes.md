# Two dyadic bands

For the original prime 2 test, take the Dirichlet length $`K`$ to infinity at each fixed $X$, then let $X$ grow with $h=1/\log X$. Among the remaining pairs there are two disjoint ranges of the dyadic exponent, one for each wall. Their signed totals tend to $+A_\omega$ and $-A_\omega$, with $A_\omega>\tfrac12\int\omega$ for every nonzero admissible weight. The pairs outside both ranges contribute $O_{I,\omega}(h^2)$ after the first limit. After centring the two band sums at the same heat, their terms of order $h$ have opposite mixed terms and equal half-squares of the residue weights. I could not prove that the sum of the squares vanishes. [Original test](prime-two-kernel.md#the-test), [pair cutoffs](prime-two-pair-sieve.md#relative-quotients) and [Masullo, §12.3](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514).

## The pairs

Fix $I=[d,e]\Subset(1/2,1)$, $0\le\omega\in C_c((d,e))$ and a centre $s$ in its support. Put $\sigma_\pm=s\pm h/2$, $\varepsilon=h^2$, $Y=(\log X)^2$, $`P=\{p:3\le p\le Y,\ p\text{ prime}\}`$ and take $X$ large enough that both walls lie in $I$, $2\varepsilon<h$ and $\varepsilon<1-e$. Use the original cutoff $R=R(X)=\lfloor\exp(c_I S_d(X)/2)\rfloor$ and $S_d(X)=(\log X)^{1/(1-d)}(\log\log X)^{d/(1-d)}$. Write $b_r=\Delta_h a_{K,X}(s,r)$ for the original rational coefficient increment and $q(r)$ for the reduced denominator. In each ordered pair require $q(r),q(r')>R$. Write $r/r'=a/b$ in lowest terms and impose $a,b>R^2$ and $p\nmid ab$ for every $p\in P$. Sum the products $b_r\overline{b_{r'}}$ over these pairs. The sums converge absolutely at finite $K,X$. All cutoffs below are applied before the first limit. [Coefficient and pulse conventions](prime-two-increments.md#one-coefficient).

The integer $k=v_2(r)-v_2(r')=v_2(a)-v_2(b)$ labels each pair. Define two sets of integers, using the wall at the indicated sign,

```math
\mathcal B_\pm(K,s)=\{k:K^{1-\sigma_\pm-\varepsilon}\lt |k|\lt K^{1-\sigma_\pm+\varepsilon}\}.
\tag{1}
```

The exponent of 2 itself is of this size. The set $\mathcal B_+$ contains the smaller values of $|k|$. Write $S_+(K,X),S_-(K,X),S_0(K,X)$ for the signed sums integrated against $`\omega(s)\,ds`$, selecting the respective bands or their complement at each $s$.

## A single wall

Give 2 the phase $z=e^{i\theta}$ and the odd primes in $P$ independent unit phases $u_p$. Put $\chi_u(n)=\prod_{p\in P}u_p^{v_p(n)}$ and

$$
E_K(w)=\sum_{n\le K}(-1)^{n-1}z^{v_2(n)}\chi_u(n)n^{-w},\qquad
M_u(w)=\prod_{p\in P}\frac{1-p^{-w}}{1-u_pp^{-w}} .
$$

For fixed $P$ and any $\alpha>0$ the sum of the twisted coefficients through $N$ is $c(z)M_u(1)N+O_{P,\alpha}(N^\alpha)$ uniformly on the phase circles, where $c(z)=(1-z)/(2-z)$. Without the odd twists, counting the integers with $v_2(n)=j$ introduces an error $O(1)$ in each of $O(\log(2N))$ classes. The mean is $1/2-\sum_{j\ge1}z^j/2^{j+1}=c(z)$. For the odd twists convolve with the coefficients $m_u(\ell)$ of $M_u$. They are supported on odd integers and $\sum_\ell |m_u(\ell)|\ell^{-\alpha}$ is bounded uniformly for fixed $P$. This absorbs the logarithmic count error (use an exponent less than $\alpha$) and the tail of the linear main term. On a complete block $N=2^m$, $m\ge1$, the calculation with all odd phases equal to one gives the exact signed sum $c(z)(2^m-z^m)$. [Lean proof of the block identity](PrimeTwoDyadicModes.lean).

On $w=\sigma+it$, $\sigma\in I$, partial summation gives

$$
E_K(w)=D_K(w)+F_{z,u}(w)+O_{P,I,\alpha}((1+|t|)K^{-\sigma+\alpha}),\qquad
D_K(w)=c(z)M_u(1)\frac{K^{1-w}}{1-w},
\tag{2}
$$

where $0<\alpha<d$ and $F_{z,u}(w)=M_u(w)(1-2z2^{-w})(1-2^{-w})\zeta(w)/(1-z2^{-w})$. To obtain the constant term, first compute the Euler factor for $\Re w>1$, then continue the integral of the error in the coefficient count to $\Re w>\alpha$. Differentiating that integral costs a factor $O_I(\log(2K))$ in the remainder. The bounds $`|F_{z,u}|+|F'_{z,u}|\ll_{P,I}1+|t|`$ follow from Euler summation for $\zeta$ and its derivative. The constants here may depend on $X$ through $P$.

Use $g(v)=f(v)+f(-v)$ for the original pulse and $H(w)=\int g(v)e^{wv}dv$. Write $I_\sigma(B)=(2\pi)^{-1}\int H(\sigma-1/2+it)B(\sigma+it)dt$. For $\Phi_X(y)=(1-e^{-X|y|^2})/y$, extended by zero, put

$$
V_\sigma=I_\sigma(E'_K\Phi_X(E_K)),\quad
C=I_\sigma((1-w)^{-1})=\int_0^\infty f(v)e^{-v/2}dv,\quad Q_\sigma=V_\sigma-C .
\tag{3}
$$

The expression for $C$ follows by writing $(1-w)^{-1}=\int_0^\infty e^{-(1-w)v}dv$ and applying Fourier inversion. Also $I_\sigma(1)=g(0)=0$. These identities remove both terms in $D'_K/D_K=-\log K+(1-w)^{-1}$.

We need a uniform bound on $V_\sigma$ even when $`E_K`$ is small. On each unit height band of $[d,8]+i\mathbb R$, the number of preimages of any $y\in\mathbb C$ under $`E_K`$ is $O_I(\log(2K)+\log(2+|y|))$, counting multiplicities. Indeed $|E_K(8+it)-E_K(9+it)|\ge2^{-9}-\sum_{n\ge3}n^{-8}>0$ for $K\ge2$, so one centre stays a fixed distance from $y$. Use Jensen's formula at that centre with outer radius equal to the distance to a fixed line $\Re w=d_0$, $0<d_0<d$. A slightly smaller disc covers $[d,8]$ across a height interval of fixed positive length for either centre. On the outer circle $|E_K-y|\ll K^{1-d_0}+|y|$. Finitely many intervals cover a unit band. Integrating this count against $(X/\pi)e^{-X|y|^2}dA(y)$ bounds the Gaussian area of $`E_K`$ there by $O_I(\log(2K))$ for $X\ge1$. Cauchy–Green relates $V_\sigma$ to the trace at 8 by integration of $`H(w-1/2)(X/\pi)|E'_K|^2e^{-X|E_K|^2}`$. At each fixed $`K`$ the horizontal integrals vanish by rapid decay of $`H`$. The trace at 8 is bounded uniformly since $`E_K`$ stays away from zero there. Summing the rapid decay of $`H`$ over the bands gives $`\sup_{\sigma,\theta,u}|Q_\sigma|\ll_I\log(2K)`$. [Wall area identity](prime-two-pair-sieve.md#odd-twists).

For each fixed $X$ and $0<\delta<1-e$, the normalised angular integral satisfies

```math
\sup_u\|Q_\sigma\|_{L^1(d\theta/2\pi)}
\ll_{I,X,\delta}K^{-(1-\sigma-\delta)}\log^2(2K).
\tag{4}
```

On the arc $|1-z|\le K^{-(1-\sigma-\delta)}$ use the preceding supremum times the arc length. Outside it, for $|t|\le K^{\delta/8}$, equation (2) gives $|D_K|\gg_{I,P}|1-z|K^{1-\sigma}/(1+|t|)$, since $|M_u(1)|$ is bounded below uniformly at fixed $P$. The rest of (2) is $O_{I,P}(1+|t|)$, so $|E_K|\ge |D_K|/2$ for large $`K`$. Subtracting $D'_K/D_K$ from $`E'_K/E_K`$ costs $`O_{I,P}((1+|t|)^2\log(2K)/(|1-z|K^{1-\sigma}))`$. The Gaussian error decays faster than any power of $`K`$ on this part. Above the height cutoff use $|\Phi_X|\le\sqrt X$, the finite bound $|E'_K|\ll_I K^{1-d}\log(2K)$ and rapid decay of $`H`$. After integrating in $t$, the angular integral of $1/|1-z|$ outside the arc is $O(\log(2K))$. This proves (4) uniformly in the wall.

The other estimate uses complex $\theta$. The wall function $Q_\sigma$ extends holomorphically to $|\Im\theta|\le cK^{-(1-\sigma)}$ with size $O_{I,X}(K^{C_1})$ for some fixed power $C_1$. In the integrand continue $`E_K`$ and its conjugate by using $e^{ij\theta}$ and $e^{-ij\theta}$ respectively, and conjugate all the other coefficients for the second factor. Each differs from its value on the real angle by $O_I(1)$. Indeed $\sum_{n\le K}v_2(n)e^{c_0v_2(n)}n^{-\sigma}\ll_I K^{1-\sigma}$ for $0<c_0<\log 2$, by summing over $n=2^jm$ and using $\sum_j j(e^{c_0}/2)^j<\infty$. If the two extensions are $E$ and $`E^*`$, completing a square gives $`\Re(EE^*)\ge-O_I(1)`$. The continuation of $\Phi_X$ is $`X E^*\int_0^1e^{-tXEE^*}dt`$. It and the continued derivative have at most polynomial growth in $`K`$, uniformly in height, so the integral against $`H`$ has the asserted bound.

## The square

Let $L$ be the increment summed over $q(r)\le R$ with factors $z^{v_2(r)}\chi_u(r)$. Its dyadic frequencies are bounded by $O(\log(2R))$ since the pulse has $2/5<r<5/2$. At each fixed $X$, $L$ is uniformly bounded in $`K`$ by convergence of the finite coefficient vector. The identity $V_{\sigma_+}-V_{\sigma_-}=\sum_r b_r z^{v_2(r)}\chi_u(r)$ follows by Fourier inversion. Average $|V_{\sigma_+}-V_{\sigma_-}-L|^2$ over the odd phases. This retains exactly the pairs with $p\nmid ab$ for $p\in P$. Finally subtract the polynomial $P_K$ from the retained pairs with $\min(a,b)\le R^2$ and both denominators above $R$. It too has degree $O(\log(2R))$, since the ratio of two pulse arguments is between $4/25$ and $25/4$. The resulting real function $G_K(\theta)$ has $`k`$th Fourier coefficient equal to the sum of the requested products at that $`k`$. The unshifted pair sum is $G_K(0)$.

With $Q_\pm=Q_{\sigma_\pm}$, split $G_K=T_K+A_K+B_K$ as follows, writing $\mathbb E$ for the odd-phase average.

$$
\begin{aligned}
T_K&=\mathbb E|L|^2-P_K,\\
A_K&=\mathbb E\bigl[|Q_+|^2-2\Re(Q_+\overline L)\bigr],\\
B_K&=\mathbb E\bigl[|Q_-|^2-2\Re(Q_+\overline{Q_-})+2\Re(Q_-\overline L)\bigr].
\end{aligned}
\tag{5}
$$

At fixed $X$ the polynomial $T_K$ converges in its coefficients uniformly in $s$. The limiting sum of the absolute values of these coefficients is $O_I(h^2)$. For $P_K$ use the bound for the finite set of quotient shifts and convergence of the increment vector in $\ell^2$. For $L$ the limiting sum of absolute coefficients before the phase average is $O_I(hR e^{-c_I S_d(X)})$. [Quotient estimate](prime-two-pair-sieve.md#relative-quotients) and [increment norm](prime-two-increments.md#distinct-ratios).

For $A_K$ take $\sigma=\sigma_+$ and for $B_K$ take $\sigma=\sigma_-$. Each product in its expression contains a factor $Q_\sigma$ or its conjugate. Equation (4), the supremum bounds and $\delta=\varepsilon/2$ bound the $L^1$ norm by $O_{I,X}(K^{-(1-\sigma-\varepsilon/2)}\log^3(2K))$. Every Fourier coefficient satisfies the same bound. Their absolute sum over $|k|\le K^{1-\sigma-\varepsilon}$ therefore tends to zero. In the complex strip of width $cK^{-(1-\sigma)}$ both factors are bounded by a fixed power of $`K`$, using $Q(\theta)\overline{Q(\bar\theta)}$ for the continued square. For the mixed term in $B_K$ the strip for $Q_-$ fits inside the one for $Q_+$. Shifting the angular integral bounds each coefficient by $O_{I,X}(K^{C_2}\exp(-c|k|K^{-(1-\sigma)}))$. Summing over $|k|\ge K^{1-\sigma+\varepsilon}$ gives a quantity tending to zero. All bounds in this paragraph are uniform in the allowed $s$ at fixed $X$.

Thus the absolute sum of the grouped coefficients of $A_K$ outside $\mathcal B_+$ tends to zero, and the same holds for $B_K$ outside $\mathcal B_-$. The two bands are disjoint. Eventually the bounded frequencies of $T_K$ miss both bands. Fourier inversion at $\theta=0$ now gives

```math
\begin{aligned}
\lim_K S_+(K,X)&=\int\omega(s)\lim_K A_K(s,0)\,ds,\\ \lim_K S_-(K,X)&=\int\omega(s)\lim_K B_K(s,0)\,ds,\\ \left|\lim_K S_0(K,X)\right|&\ll_{I,\omega} h^2.
\end{aligned}
\tag{6}
```

The limits on the right exist uniformly in $s$ for fixed $X$, since at $z=1$ the alternating polynomial converges to $\eta M_u$ on the walls, with polynomial bounds in height uniform in $K,u$. Equation (6) is for pairs in the original unshifted sum. Its proof takes the absolute value after grouping by $`k`$, which permits cancellation within each group.

## Opposite limits

At $z=1$ let $R_s=\sum_{\Re\rho>s}m_\rho H(\rho-1/2)$, summing over zeta zeros. For almost every $s$ in the interior of $I$, as $X\to\infty$ the two wall values $\lim_K V_{s\pm h/2}(0,u)$ tend to the same $P_f-2\log(2)/\sqrt2-R_s$ uniformly over $u\in\mathbb T^{P_X}$, where $P_f=H(1/2)$.

To see the uniformity, start at $\Re w=3$. There $\eta M_u$ stays a fixed distance from zero and its logarithmic derivative is bounded uniformly in $P,u$. Every term of $`M'_u/M_u`$ has frequency $`\log(p^j)`$ with $`p\ge3`$ and its integral against $`H`$ vanishes by the pulse support. The regularised trace there therefore tends uniformly to the ordinary trace of $`\eta'/\eta`$. Move to either of the two inner walls with the Gaussian area identity. On $`\Re w\ge d_0`$ for fixed $`1/2<d_0<d`$, the bounds $`|M_u|,|M_u|^{-1}\le X^{o(1)}`$ hold because $`\sum_{p\le Y}p^{-d_0}=o(\log X)`$. On bounded rectangles the small Gaussian values of $\eta M_u$ are confined to discs of radius $h^2$ about the zeros of $\eta$, with their full multiplicities by Rouché's theorem. On $[d,3]$ the high ordinates have uniformly small total variation by the Jensen estimate for $\eta-w/M_u$ at $|w|\le X^{-1/4}$ and the exponential bound for larger $|w|$. The same proof on this wider interval uses the centre $4+it$. These estimates give a uniform bound for the traces. For $s$ different from the real parts of the zeros, both walls eventually include the same discs on each bounded rectangle. At $u=1$ the contour formula identifies the limit, which is independent of the phases. [Uniform Gaussian area estimate](prime-two-pair-sieve.md#growing-sets) and [unshifted wall value](prime-two-eta.md).

Write $F(t)=\int_0^\infty f(v)e^{tv}dv$, so $P_f-C=F(1/2)$. The low term $L$ tends to zero uniformly in the iterated limit. In (5) the two $Q$'s have the same limit for almost every $s$. Dominated convergence in (6) gives

$$
\lim_X\lim_K S_+= A_\omega,\qquad
\lim_X\lim_K S_-=-A_\omega,\qquad
A_\omega=\int\omega(s)\left|F(1/2)-\sqrt2\log2-R_s\right|^2ds .
\tag{7}
$$

In particular the original signed bands each have a nonzero limit for nonzero $\omega$. The positive bump has mass $1/16$, the subtracted translate is nonnegative and the whole pulse is supported below $\log(5/2)$ on the positive axis. Hence $F(1/2)\le\sqrt{5/2}/16<1/8$. The unconditional eight-box estimate gives $|R_s|<6.5\times10^{-65}$ by applying the decay of this pulse above the height in Platt and Trudgian's finite verification. These inequalities give $A_\omega>\tfrac12\int\omega$ for nonzero $\omega$. More precisely, replacing $R_s$ by zero in (7) changes $A_\omega$ by at most $2\times10^{-64}\int\omega$, since $|F(1/2)|\le 2\sqrt{5/2}/16$. Quadrature of the pulse gives $(F(1/2)-\sqrt2\log 2)^2\simeq0.949274032664$. [Pulse bound](prime-two-kernel.md#an-explicit-bound), [Platt and Trudgian](https://arxiv.org/abs/2004.09765) and [Bellotti and Wong](https://arxiv.org/abs/2412.15470).

The coefficient at the scale asked for follows by combining (6) with the full residual identity. With $J_\beta=\sum_{\Re\rho=\beta}m_\rho H(\rho-1/2)$ it is

$$
\lim_{X\to\infty}\frac1h\lim_{K\to\infty}(S_++S_-)=
\sum_\beta\omega(\beta)|J_\beta|^2 .
\tag{8}
$$

Equation (7) proves cancellation of the terms of order one. On a wall interval containing a zero mass, the two values used in (5) can differ even as the interval shrinks. The centres of those intervals occupy a length of order $h$. Equation (8) measures exactly their possible contribution. Proving $o(h)$ for every $\omega$ still requires proving $J_\beta=0$ at each interior real part for this pulse. [Original residual identity](prime-two-pair-sieve.md#growing-sets).

## At order $h$

The first correction can also be assigned to each band. To use the same reference for both walls at fixed heat, put $Q_{X,u}(v)=\lim_K V_v(0,u)-C$, $U_X(v)=\mathbb E_u|Q_{X,u}(v)|^2$ and $M_X=\int\omega(v)U_X(v)dv$. Write $q(v)=F(1/2)-\sqrt2\log2-R_v$ and $q(\beta^-)=F(1/2)-\sqrt2\log2-\sum_{\Re\rho\ge\beta}m_\rho H(\rho-1/2)$. For every original weight,

```math
\lim_{X\to\infty}\frac{\lim_{K\to\infty}S_\pm(K,X)\mp M_X}{h}
=\frac12\sum_\beta\omega(\beta)|J_\beta|^2
\ \pm\sum_\beta\omega(\beta)\Re\!\left(\overline{q(\beta^-)}J_\beta\right).
\tag{9}
```

Here $M_X\to A_\omega$. Using the same $M_X$ on both sides of the split avoids requiring a rate for this convergence. It cancels exactly in the original sum.

For the plus band, (5) and (6) give $\lim_K S_+=\int\omega(s)U_X(s+h/2)ds+o(h)$. The terms involving $L$ are $o(h)$ uniformly after the first limit, by its exponentially small coefficient bound and the uniform bound for the wall trace. For $\omega\in C_c^1((d,e))$, translate $\omega$ by $-h/2$ in this integral and subtract $M_X$. Division by $h$ and dominated convergence give $-\tfrac12\int\omega'(s)|q(s)|^2ds$. On any compact interval in the strip the functions $U_X$ have uniformly bounded variation. Indeed the variation of $Q_{X,u}$ is bounded by the total variation of the Gaussian area measure between the extreme walls, uniformly over $u$ and all large $X$. Squaring in modulus costs at most $2\sup|Q_{X,u}|$, and averaging preserves the bound. Hence the integral of $|U_X(s+h/2)-U_X(s)|/h$ on a slightly smaller interval is uniformly bounded. Uniform approximation by smooth weights extends the limit to every compactly supported continuous $\omega$, giving $`\tfrac12\int\omega\,d|q|^2`$.

Since the residues are absolutely summable, $q$ has jump $J_\beta$ at $\beta$, and $d|q|^2$ has mass $2\Re(\overline{q(\beta^-)}J_\beta)+|J_\beta|^2$ there. This follows first for finitely many real parts and then in bounded variation, since deleting an absolutely summable tail changes $q$ by a quantity small in both supremum norm and variation. This gives the plus sign in (9). Subtract the plus-band limit from (8) to get the minus sign. The mixed terms in (9) cancel, while the two half-squares add. Their sum can vanish for every nonnegative weight only if each $J_\beta$ vanishes. A direct integral on three walls of the radial theta kernel gives the same $J_\beta$ and removes the pole at 1 exactly. [Kernel calculation](prime-two-kernel-wall.md).

Lean checks the dyadic block formula, the algebra of (5), the error bound for extracting one band of an absolutely summable series and the half-jump identity used in (9). The estimates (2) to (4) and the limits (6) to (9) are analytic arguments in this note. At $K=512$, $X=64$, $s=3/4$ and all phases one, a numerical check on the full walls before sieving gives $0.86765314$ for $|Q_+|^2$ and $-0.86583845$ for $|Q_-|^2-2\Re(Q_+\overline{Q_-})$. It also recovers $C=0.007819126722$ by vertical integration on both walls. The proof of (7) applies the actual pair cutoffs. [Lean source](PrimeTwoDyadicModes.lean) and [numerical check](prime-two-dyadic-modes.py).
