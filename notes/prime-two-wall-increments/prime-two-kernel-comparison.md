# Comparing the kernels

For the original prime 2 test, removing the pole, real-place and dyadic terms leaves one integral of the logarithmic derivative of the alternating zeta function on the right wall of the strip. I could not evaluate it exactly on an inner wall. There is also a polar-free input on which the two kernels have provably different Mellin coefficients. Both calculations use the rational field, the trivial character and Sobolev parameter four. [Masullo, §12.3](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514).

## The wall

Write $a=\log 2$, $L=a/\sqrt2$, $`F(z)=\int f(u)e^{zu}\,du`$ and $H(z)=F(z)+F(-z)$. Here $f(u)=B(u-a)-cB(u-a-1/8)$ is the original pulse. Put $P_f=H(1/2)$. The complete norm-first local integral gives $L+P_f$ before pole subtraction. Take $\sigma=1/2+\alpha$ for $0<\alpha<1/2$ and write

```math
I_\sigma(Q)=\frac1{2\pi}\int_{\mathbb R} H(\alpha+it)Q(\sigma+it)\,dt.
```

For the global calculation take a width with the polynomial separation from the zeros used in the radial projection proof. This holds for almost every width. The Mellin transform of its Gaussian theta generator is $\xi(1/2+z)$. The trace on the complement counts the zeros inside the strip with multiplicity. Integrating $F(s-1/2)\xi'(s)/\xi(s)$ around the strip and reflecting the left wall with $\xi(s)=\xi(1-s)$ gives the global geometric trace $G_\alpha(f)=-I_\sigma(\xi'/\xi)$. The two vertical integrals add after reflection. For convergence, the genus-one logarithmic derivative is a sum of $1/(s-\rho)+1/\rho$ and a constant. Wall separation bounds the terms with $|\rho|\le2|s|$ by a polynomial in $|s|$, and the rest by $O(|s|\sum_{|\rho|>2|s|}|\rho|^{-2})$. Horizontal sides can be chosen at heights a polynomial distance from all nearby ordinates. The rapid decay of $F$ then removes these sides. [Original pulse and projection proof](prime-two-kernel.md) and [local Fourier integral](prime-two-full-trace.md).

Let $\eta(s)=(1-2^{1-s})\zeta(s)$ and $D_2(s)=a/(2^{s-1}-1)$. Factoring the completed zeta function gives

$$
\frac{\xi'}{\xi}(s)=\frac1s+\frac1{s-1}-\frac{\log\pi}{2}+\frac12\frac{\Gamma'}\Gamma(s/2)+\frac{\eta'}\eta(s)-D_2(s).
$$

The inverse transform under $I_\sigma$ pairs with the even function $g(u)=f(u)+f(-u)$, which vanishes near zero. On $1/2<\sigma<1$ the term $1/s$ gives $F(-1/2)$ from $u>0$, and $1/(s-1)$ gives $-F(-1/2)$ from $u<0$. All constant terms give zero. The partial fraction expansion of the digamma function gives

```math
I_\sigma\!\left(\tfrac12\Gamma'(s/2)/\Gamma(s/2)\right)
=-\int_0^\infty \frac{e^{-u/2}}{1-e^{-2u}} f(u)\,du=-P_f .
```

The nonconstant fractions are $-1/(s+2n)$ for $n\ge0$. With their usual constant subtractions, the partial sums are bounded by $C_\sigma(1+|s|)$ on the wall. This permits integration against $H$ before taking the limit, and the inverse transforms sum geometrically on the support of $f$. The last equality follows from the choice of the pulse coefficient $c$. [Digamma expansion](https://dlmf.nist.gov/5.7.E6) and [pulse normalisation](prime-two-kernel.md).

On the same wall $D_2(s)=-a\sum_{n\ge0}2^{n(s-1)}$, with absolute convergence. Fourier inversion gives $I_\sigma(2^{n(s-1)})=2^{-n/2}g(-na)$. Only $n=1$ survives, so $I_\sigma(D_2)=-L$. Combining the terms gives

$$
G_\alpha(f)-L=P_f-2L-I_\sigma(\eta'/\eta).
$$

This obtains the proposed scalar identity precisely when $I_\sigma(\eta'/\eta)=P_f-2L$. On a wall $b>1$, using the same definition with $\alpha=b-1/2$, the absolutely convergent Euler series gives $I_b(\zeta'/\zeta)=-L$ and the opposite geometric expansion gives $I_b(D_2)=2L$, hence $I_b(\eta'/\eta)=L$. In shifting to the inner wall, the zeros of $\eta$ at $1+2\pi ik/a$, $k\ne0$, give total residue $3L-P_f$. In fact Poisson summation gives $\sum_k H(1/2+2\pi ik/a)=a(\sqrt2+1/\sqrt2)=3L$, and the omitted $k=0$ term is $P_f$. Any zeta zeros crossed in the shift supply additional residues. The alternating Dirichlet series for $\eta$ converges on the inner wall, but that gives no evaluation of this remaining integral of its logarithmic derivative.

I evaluated the integral between heights $-4096$ and $4096$ using Euler–Maclaurin summation for $\zeta$ and its derivative, and the exact box factors in the transform of $B$. The target in the first two rows is $P_f-2L$. The last row checks the signs against the known outer-wall value $L$. The mesh comparisons in the script apply to the finite integral. The infinite tail still enters the trace equality. [Numerical check](prime-two-wall.py).

| Wall $\sigma$ | Mesh | Finite integral | Target |
| --- | ---: | ---: | ---: |
| $5/8$ | $1/64$ | $-0.96648782319859$ | $-0.96648782319854$ |
| $3/4$ | $1/32$ | $-0.96648782319860$ | $-0.96648782319854$ |
| $5/4$ | $1/32$ | $0.49012907173439$ | $0.49012907173427$ |

For the same pulse, finite geometric sums for the reciprocal give the proposed value at every order, while the error retains the full residue at each exterior zero. In a convergent expansion the contribution from $7/3$ vanishes after summing all orders at that fixed ratio. The total from all ratios other than $2$ and $1/2$ is $P_f-R_\sigma$ in the original order of limits. The smooth density in the pulse bounds the sum of absolute exterior weights by $10^{-20000}$. [Inner-wall calculation](prime-two-eta.md) and [sum over ratios](prime-two-ratios.md).

## One input

The following calculation uses the actual local Fourier cutoff in (12.4.1.8) and the global orthogonal projection. The convolution test is still $f$. Only the vector on which the operators act is chosen here. [Masullo, §§12.4.1, steps I and II](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=516).

Take a nonnegative $\psi\in C_c^\infty((1,e^{1/64}))$ with $`\int\psi(x)x^{-1/2}\,dx=1`$, set $`D=x\,d/dx`$ and $\phi=D(D+1)\psi$, and extend $\phi$ by zero on the negative real axis. Put $U=\mathbb Z_2^\times$ and $`b=1_U-2\,1_{2U}`$. The local input $v=\phi\otimes b$ is smooth, compactly supported away from both coordinate origins, and polar-free. Its half-weighted sum over the units $`\{\pm2^j\}`$ on the radial section is

$$
h(u)=e^{u/2}\big[\phi(e^u)-2\phi(2e^u)\big].
$$

For $s=1/2+z$ its transform is

```math
M_h(z)=(1-2^{1-s})s(s-1)\int_0^\infty \psi(x)x^{s-1}\,dx .
```

Thus the two pole evaluations vanish. The Gaussian calculation below gives a zero $\rho_1=1/2+i\gamma_1$ with $14<\gamma_1<15$. The real part of the integral in this display exceeds $0.97$ there since $0<\gamma_1\log x<15/64$ on the support. Also $|1-2^{1-\rho_1}|\ge\sqrt2-1>0.4$ and $|\rho_1(\rho_1-1)|>196$. The even nonnegative bump $B$ has mass $1/16$ and support of radius $13/256$, where $\cos(\gamma_1u)>0.7$ by $\cos y\ge1-y^2/2$. Its translate factor has modulus $|1-ce^{i\gamma_1/8}|\ge1$, since $c>0$ and $\pi/2<\gamma_1/8<\pi$. It follows that $|F(i\gamma_1)M_h(i\gamma_1)|>196(0.4)(0.97)(0.7)/16>3$ for this normalisation of $\psi$. [First zero](https://www.lmfdb.org/zeros/zeta/) and [certified verification](https://arxiv.org/abs/2004.09765).

Here are the limits on this vector. Normalise the additive measure of $\mathbb Z_2$ to one. The Fourier transform of $b$ is $`-1_{|w|_2=2}+\tfrac12 1_{|w|_2=4}`$. Consequently the local Fourier cutoff $`|w_\infty||w_2|_2\le T`$ acts on just two real Fourier cutoffs, at $`T/2`$ and $`T/4`$. Its inverse transform tends to $`v`$ after radial periodisation in $`X_\alpha=L^2(e^{2\alpha|u|}(1+u^2)^2\,du)`$. To check the norm, the real Fourier tail outside $[-R,R]$ is bounded by $C_N R^{-N}/(1+|x|)$ for $R\ge1$ by one integration by parts and the rapid decay of $\widehat\phi$. The square of this bound is integrable against $`e^{2\alpha|\log x|}(1+\log^2 x)^2\,dx`$ on $x>0$. For each inverse dyadic annulus its values at $2^j$ are bounded and vanish below some $j_0$. Periodisation contributes translations by $j\log2$ with factor $2^{-j/2}$, whose norms for $j\ge0$ sum by $\sum_j2^{-j(1/2-\alpha)}(1+j\log2)^2<\infty$. The bounds apply to both signs. Bounded smoothing of the endpoint over a fixed range of log shifts preserves these estimates.

Write $`A_fh(u)=\int f(t)h(u-t)\,dt`$. Periodisation commutes with the convolution in (12.4.1.7) after including the half-weight, and the physical norm cutoff commutes with periodisation because $`|\pm2^j|_S=1`$. At fixed $`T`$ the periodised vector is in $`X_\alpha`$, so removing the physical window recovers its convolution in norm. Sending $`T`$ to infinity next therefore gives $A_fh$, with Mellin value $F(i\gamma_1)M_h(i\gamma_1)\ne0$.

On the global side let $\Pi_\alpha$ project onto the closed theta space. Mellin evaluation at $i\gamma_1$ is continuous on $`X_\alpha`$. It vanishes on each global theta input by Mellin unfolding and the zeta factor, hence on the closure. Thus $A_f\Pi_\alpha h$ has Mellin value zero there. Removing the physical cutoff gives this same vector in norm. This mismatch holds for every $0<\alpha<1/2$ and only uses the one verified critical-line zero. It rules out identifying the two actions on every common radial input by taking the norm-first limit. A scalar trace comparison could still hold.

## Theta signs

The global Gaussian generator gives a direct check of the zero used above. By Poisson symmetry its Mellin transform on the central line is the real integral

```math
\xi(1/2+it)=2\int_0^\infty e^{u/2}\cos(tu)\sum_{n\ge1}(4q_n^2-6q_n)e^{-q_n}\,du,\qquad q_n=\pi n^2 e^{2u}.
```

I kept $n\le3$ and $0\le u\le2$. For $x\ge0$, $(4x^2+6x)e^{-x/2}<18$. The omitted $n\ge4$ on this interval cost at most $216\sum_{n\ge4}e^{-3n^2/2}<10^{-8}$. On $u=2+v$, $v\ge0$, one has $q_n/2\ge75n^2(1+2v)$, giving an integrated tail below $10^{-30}$.

For $t=14,15$ the truncated integrand is bounded by $1000$ on every complex disc of radius $1/16$ centred in $[0,2]$. There $\Re q_n\ge|q_n|/2$, $|e^{u/2}|<3$ and $|\cos(tu)|<3$. Cauchy’s estimate bounds its fourth derivative by $`4!\,1000\,16^4`$. Composite Simpson quadrature with $4096$ subintervals of length $1/2048$ then has error below $2\times10^{-6}$. Evaluating the nodes with interval arithmetic and widening by $3\times10^{-6}$ encloses the first value in $[0.000198,0.000205]$ and the second in $[-0.000709,-0.000702]$. Continuity gives the required zero between them. Run `python3 prime-two-wall.py --theta-signs` to reproduce the enclosures. [Calculation](prime-two-wall.py).
