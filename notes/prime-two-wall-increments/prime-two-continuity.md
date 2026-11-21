# Moving walls

For the original prime 2 test, the large-denominator tail is continuous at an inner wall $\beta$ exactly when the sum of the weights of the zeros on that wall is zero. I could prove continuity away from the real parts of zeros and a quantitative estimate at almost every wall. I could not prove cancellation at every real part in $1/2<\beta<1$. The same finite sums give an exact example at the upper endpoint where uniform convergence fails. [Masullo, §12.3](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514) and [original denominator sum](prime-two-denominators.md).

## Interior walls

Fix an integer $Q\ge2$. Write $F_{K,X}(\sigma)$ for the original sum over reduced denominators $q>Q$ with Dirichlet cutoff $K$ and heat parameter $X$. Write $F_X=\lim_{K\to\infty}F_{K,X}$ and $F=\lim_{X\to\infty}F_X$. Here $a=\log2$, $L=a/\sqrt2$, $H(z)=\int_{\mathbb R}(f(u)+f(-u))e^{zu}\,du$ for the original pulse, and $P_f=H(1/2)>0$. For a zeta zero $\rho$ of multiplicity $m_\rho$ put $w_\rho=m_\rho H(\rho-1/2)$ and $J_b=\sum_{\Re\rho=b}w_\rho$. Conjugation makes each $J_b$ real. The regularised trace gives

$$
F(\sigma)=P_f-\sum_{b>\sigma}J_b-\tfrac12 J_\sigma,\qquad
F(\beta+)-F(\beta-)=J_\beta .
$$

The absolute sum of the weights on any closed substrip to the right of $1/2$ is finite. Grouping by real part and taking one-sided limits therefore works even if the real parts accumulate. Continuity throughout the open strip would imply $J_\beta=0$ at each real part and hence $F\equiv P_f$. Lean checks both one-sided limits, the continuity criterion for an absolutely summable family of complex weights, and the regrouping when each sum is zero. [All-wall trace](prime-two-denominators.md#every-wall) and [Lean proof](PrimeTwoContinuity.lean).

On any compact interval $I=[d,e]\subset(1/2,1)$ the contribution from $q\le Q$ converges uniformly to $-2L$. In the prime-circle representation the $L^2$ norms of the logarithmic derivative $T_\sigma$ are uniformly bounded there. For each fixed $k$ the inverse moments $\int |E_\sigma|^{-2k}$ are bounded by the Euler product estimate at $d$ times $(2^{1-e}-1)^{-2k}$. The coefficient estimate for the omitted heat factor is thus uniform on $I$. The earlier proof with $k$ growing also works on $I$ for a common cutoff $2\le Q(X)\le\exp(C(\log X)^{1/(1-d)})$. In particular for fixed $Q$,

$$
F_X(\sigma)=I_\sigma(V_X)+2L+o_I(1),\qquad
V_X=(\eta'/\eta)(1-e^{-X|\eta|^2}),\qquad \eta=(1-2^{1-s})\zeta .
$$

Here $I_\sigma(V)=\int H(\sigma-1/2+it)V(\sigma+it)\,dt/(2\pi)$. At fixed $X$ the limit in $K$ is uniform on $I$. On the prime circles the $L^2$ tails of $E_K,\dot E_K$ are bounded using $\sum_{n>K}(1+\log n)^2n^{-2d}$. On the vertical lines partial summation gives local uniform convergence of $\eta_K,\eta_K'$ and a common polynomial bound in the height. The regularising function $\Phi_X(w)=(1-e^{-X|w|^2})/w$ extends by zero, is bounded by $\sqrt X$ and is $2X$-Lipschitz. Its two complex partial derivatives have modulus at most $X$, using $1-(1+y)e^{-y}\le y$ for $y\ge0$. [Prime circles](prime-two-ratios.md#all-multiples) and [moving cutoff](prime-two-denominators.md#moving-cutoff).

At an inner wall $\beta$ with $|\beta-\Re\rho|\ge c(1+|\Im\rho|)^{-2}$, a zero contributing to $F(\beta+h)-F(\beta)$ must have $|\Im\rho|\ge (c/|h|)^{1/2}-1$. The pulse bound $|H(\rho-1/2)|\ll e^{-c_0\sqrt{|\Im\rho|}}$ and the zero count $N(T)\ll T\log(2T)$ then give $|F(\beta+h)-F(\beta)|\ll_{\beta,f} e^{-c_\beta |h|^{-1/4}}$ for small $|h|$ and some $c_\beta>0$. The intervals of radius $(1+|\Im\rho|)^{-2}$ about the real parts have finite total length. Almost every $\beta$ lies in only finitely many of them and is distinct from all their centres, so a positive separation constant exists there. [Pulse bound](prime-two-eta.md#residue-bound).

## Thin strips

The Cauchy–Green formula makes $I_{\beta+h}(V_X)-I_{\beta-h}(V_X)$ the integral over $|\Re s-\beta|<h$ of $(X/\pi)H(s-1/2)|\eta'(s)|^2e^{-X|\eta(s)|^2}$. Its absolute integral at large heights tends to zero uniformly in $X\ge1$ on compact interior intervals, by the preimage count for $\eta$ in unit height bands. On bounded regions the density tends to zero away from the zeros. This proves that $F_X(\sigma_X)\to F(\beta)$ for every $\sigma_X\to\beta$ if the wall at $\beta$ contains no zero. It also gives an explicit limit at walls containing zeros. [Uniform area estimate](prime-two-denominators.md#every-wall).

For a zero $\rho$ of multiplicity $m$ on the wall set $d_\rho=|\eta^{(m)}(\rho)/m!|^{1/m}$. With $s=\rho+z/(d_\rho X^{1/(2m)})$, the limiting area density near $\rho$ becomes $w_\rho$ times the probability density $(m/\pi)|z|^{2m-2}e^{-|z|^{2m}}$. This follows from Taylor expansion of $\eta$ and $\eta'$. On a sufficiently small disc, $|\eta(s)|\gg|s-\rho|^m$ and $|\eta'(s)|\ll |s-\rho|^{m-1}$ supply an integrable bound after rescaling. Write $\Psi_m(y)$ for the probability that $|\Re z|<y$ under this density. In particular $\Psi_1(y)=\operatorname{erf}(y)$.

For a fixed integer $k\ge1$ and $\lambda>0$ set $h_X=\lambda X^{-1/(2k)}$. The actual large-denominator tails satisfy

$$
\lim_{X\to\infty}\big[F_X(\beta+h_X)-F_X(\beta-h_X)\big]
=\sum_{\substack{\Re\rho=\beta\\m_\rho<k}}w_\rho
+\sum_{\substack{\Re\rho=\beta\\m_\rho=k}}w_\rho\Psi_k(\lambda d_\rho).
$$

The scaled half-width at $\rho$ is $\lambda d_\rho X^{1/(2m_\rho)-1/(2k)}$, giving respectively the whole profile, a finite fraction or zero as $m_\rho<k$, $m_\rho=k$ or $m_\rho>k$. One can take the limit on a bounded region first and then remove the height bound by the uniform area estimate. For simple zeros the contribution at scale $X^{-1/2}$ is therefore $\sum_{\Re\rho=\beta,\ m_\rho=1}H(\rho-1/2)\operatorname{erf}(\lambda|\eta'(\rho)|)$.

In fact equicontinuity of $F_X$ at $\beta$ as $X\to\infty$ requires cancellation separately at each pair $(m_\rho,d_\rho)$ on that wall. To see this, start at $k=1$ and induct on $k$ in the displayed formula. With all lower multiplicity sums zero, equicontinuity gives $\sum_{m_\rho=k,\,\Re\rho=\beta}w_\rho\Psi_k(\lambda d_\rho)=0$ for every $\lambda>0$. This mixture determines the finite signed or complex measure of the scales. In the variable $\log\lambda$ it is the cumulative distribution of the convolution of $\sum w_\rho\delta_{-\log d_\rho}$ with the law of $\log|\Re z|$. The latter has characteristic function

$$
\mathbb E|\Re z|^{it}=\Gamma(1+it/(2k))
\frac{\Gamma((1+it)/2)}{\sqrt\pi\,\Gamma(1+it/2)} .
$$

In polar coordinates $|z|^{2k}$ has the exponential density and the angle is uniform. The gamma integral and the beta integral give the formula, which has no zero for real $t$. Uniqueness of Fourier transforms of finite measures now forces each grouped weight at a given $d_\rho$ to vanish. Conversely, if those grouped weights vanish on the wall for all $m,d$, take any $X_j\to\infty$, $\sigma_j\to\beta$. Along a subsequence the shifts have fixed sign and $|\sigma_j-\beta|X_j^{1/(2m)}$ has a limit in $[0,\infty]$ for every $m$. The fraction cut from each radial profile depends only on $m,d$ and those limits. The grouped contributions on the wall cancel, and the others tend to their usual whole or zero values. This gives $F_{X_j}(\sigma_j)\to F(\beta)$ by the subsequence criterion. Continuity of $F$ itself requires only $J_\beta=0$. [Beta integral](https://dlmf.nist.gov/5.12.E1) and [gamma function](https://dlmf.nist.gov/5.2).

## Squared increments

There is a test for continuity that allows different scales at the same real part to cancel. Take $0\le\omega\in C_c((1/2,1))$ and set $\Delta_h B(s)=B(s+h/2)-B(s-h/2)$. For fixed $Q$ the original finite sums in their prescribed order give

$$
\lim_{h\downarrow0}\lim_{X\to\infty}\lim_{K\to\infty}
 \frac1h\int\omega(s)|\Delta_h F_{K,X}(s)|^2\,ds
=\sum_b\omega(b)|J_b|^2 .
$$

For fixed small $h$, uniform convergence in $K$ at fixed $X$ gives the first passage. The vertical area estimate and the bounded retained coefficients bound $F_X$ uniformly on compact inner intervals, so dominated convergence gives the passage in $X$. For almost every $s$, $\Delta_h F(s)=\sum_b J_b\mathbf1_{|s-b|<h/2}$. Upon expanding the square, the two intervals centred at distinct $b,b'$ eventually have no overlap. A repeated centre gives $\omega(b)$ after division by $h$ and passage to the limit. Absolute double summability justifies the calculation.

Consequently the integral of $|\Delta_h F|^2$ against $\omega$ is $O(h)$, and continuity on the open strip is equivalent to $o(h)$ for every such $\omega$. Equivalently $F$ has local fractional Sobolev regularity $H^s$ for every $0<s<1/2$, and membership in $H^{1/2}_{\rm loc}$ would prove the requested cancellation. Indeed the difference quotient integral uses the weight $h^{-1-2s}$. At $s=1/2$ the coefficient of $\log(1/\epsilon)$ in its integral over $\epsilon<h<h_0$ is exactly $\sum_b\omega(b)|J_b|^2$. If all jumps vanish then $F=P_f$ throughout the open strip. Taking $h$ to zero at fixed $K,X$ gives zero because the finite sums are locally smooth. With the prescribed order the absolute residue bound only gives $\sum_b\omega(b)|J_b|^2\le\|\omega\|_\infty 10^{-40000}$.

For the squared increment one can remove all equal-ratio terms and all pairs with at least one denominator below a growing cutoff. The remaining pairs give the same sum of squared jumps in the prescribed order. [Distinct ratios](prime-two-increments.md#distinct-ratios).

## Upper endpoint

For the same test, the same cut $q>Q$ and the same order of limits, the values just below, at and just above the endpoint are

$$
F(1^-)=P_f,\qquad F(1)=P_f/2,\qquad F(1^+)=0 .
$$

The first follows by absolute summability and the classical nonvanishing of $\zeta$ on $\Re s=1$. For $\sigma>1$ the absolutely convergent series on the vertical line and on the prime circles both select only the prime $2$ with weight $L$, so the large-denominator tail is zero. On the endpoint line the full regularised vertical integral gives half weight to the zeros of $\eta$ at $1+2\pi i j/a$ for $j\ne0$. Their total weight is $3L-P_f$ by Poisson summation. Thus the vertical value is $L-(3L-P_f)/2=(P_f-L)/2$. [Classical nonvanishing](https://dlmf.nist.gov/25.10#i) and [dyadic weights](prime-two-kernel-comparison.md#the-wall).

To get the tail at the endpoint one also needs its retained rational coefficients. On the prime circles write $u=e^{i\theta}$ for the phase of $2$ and $Z_o$ for the Euler product over odd primes at $\sigma=1$. Then $E=(1-u)Z_o/(1-u/2)$ and $U=1-e^{-X|E|^2}$ is even in $\theta$ conditional on the odd phases. The regularised derivative is $(D-V)U$ almost everywhere, with $D=au/(1-u)$ and $V=\sum_{p,j\ge1}(\log p)p^{-j}\chi_{p^j}\in L^2$. The limit in $K$ at fixed $X$ still follows from the bounded continuous regulariser and $L^2$ convergence of the Dirichlet polynomials.

With normalised angular measure, $\int DU=-a\int U/2$ since $\Re D=-a/2$ and the imaginary part is odd. The bound $|DU|\le2a\sqrt X\,|Z_o|$ justifies conditioning at fixed $X$. At frequency $u^j$ the extra term $D(u)(u^{-j}-1)$ is a Laurent polynomial bounded by $a|j|$, whose constant coefficient is $a$ for $j\ge1$ and zero otherwise. Since $Z_o\ne0$ almost everywhere, dominated convergence gives limiting angular coefficient $a(\mathbf1_{j\ge1}-1/2)$, independent of the odd phases. All coefficients of this part with nontrivial odd character therefore vanish. The term $VU$ tends to $V$ in $L^1$. At ratio $2$ its subtraction leaves $a/2-a/2=0$. At $1/2$ the coefficient is $-a/2$, giving weighted contribution $-L/2$. No other retained ratio contributes. Subtracting $-L/2$ from the vertical value proves the endpoint formula.

At each fixed $X$, $F_X$ is continuous through $1$ by the bounded regulariser, including on the prime circles. Thus for every $0<\epsilon<1/2$,

$$
\liminf_{X\to\infty}\sup_{1-\epsilon<\sigma<1}|F_X(\sigma)-F(\sigma)|\ge P_f/2 .
$$

This rules out uniform convergence up to the upper endpoint for the original tail. The endpoint calculation leaves open the possibility of continuity everywhere strictly inside the strip.

For a numerical check I used the actual first dyadic zero $1+2\pi i/a$. At $X=65536$ a rectangle with scaled half-width $0.75$ and scaled half-height $7$ gave $0.71115559+0.00000224i$ after dividing by its pulse weight. The simple-zero profile gives $\operatorname{erf}(0.75)\operatorname{erf}(7)=0.71115563$ in the limit. The check also integrates the original regularised expression along $\Re s=1$. [Reproducible calculation](prime-two-continuity.py).
