# Prime 2 test

## The test

Write $u=\log|g|$ and $`M_f(\rho)=\int_{\mathbb R} f(u)e^{(\rho-1/2)u}\,du`$. Normalise the norm-one subgroup of the rational idele class group to have Haar mass one. A test depending only on the norm kills all its nontrivial characters under convolution, so only the radial part of the global kernel contributes. For support in $\log(3/2)<u<\log(5/2)$, the classical explicit formula is

```math
\sum_\rho m_\rho M_f(\rho)
=\int r(u)f(u)\,du-\frac{\log 2}{\sqrt2}f(\log2),
\qquad r(u)=e^{u/2}\left(1-\frac1{e^{3u}-e^u}\right).
```

Here and below a sum over distinct zeros uses their full multiplicities $m_\rho$. In the local integral at 2 the only contributing coset is $2^{-1}\mathbb Z_2^\times$. On it $|g|_2=|1-g|_2=2$. Its multiplicative measure equals that of $\mathbb Z_2^\times$, giving the coefficient $\log(2)/\sqrt2$. Every other finite prime gives zero. The smooth integral above includes the pole terms and the real place with their signs. [Classical explicit formula, exercise 22(iii)](https://terrytao.wordpress.com/2021/02/12/246b-notes-4-the-riemann-zeta-function-and-the-prime-number-theorem/#smoothed) and [Masullo’s local measures, pp. 519–520](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=519).

Here is the particular smooth test used for the bound and the numerical checks. Let $A=1/32$, $h=1/512$ and $\delta=1/256$. Convolve the indicator of $[-A,A]$ with eight copies of the uniform probability density on $[-h,h]$ and then with the even probability density proportional to $\exp[-1/(1-(u/\delta)^2)]$ on $|u|<\delta$, extended by zero. The resulting bump $B$ is smooth, even, equals one near zero and has mass $2A$. Its support has radius $d=A+8h+\delta=13/256$.

Put $a=\log2$, $b=a+1/8$, and take

```math
c=\frac{\int r(u)B(u-a)\,du}{\int r(u)B(u-b)\,du},
\qquad f(u)=B(u-a)-cB(u-b).
```

Both supports lie inside the prescribed interval and are disjoint. The function $r$ is positive and strictly increasing there, hence $0<c<1$. This gives $f(a)=1$ and $\int rf=0$ exactly. Numerically $c=0.87616498956040$. The full spectral sum is $-\log(2)/\sqrt2$.

## The global trace

On the log line the ambient Hilbert space is $`X_\alpha=L^2(W_\alpha(u)\,du)`$, where $W_\alpha=e^{2\alpha|u|}(1+u^2)^2$. Let $Q_\alpha$ project orthogonally onto the complement of the closed global theta space. Masullo’s *geometric* projection is $1-Q_\alpha$. Write $`A_f h(v)=\int f(u)h(v-u)\,du`$. For a smooth norm cutoff $\kappa_\Lambda$ as in §9.3.3, the actual truncated geometric operator is $M_{\kappa_\Lambda}A_f(1-Q_\alpha)$. Linearity gives an identity at each finite cutoff.

```math
G_{\alpha,\Lambda}
=\mathrm{Tr}(M_{\kappa_\Lambda} A_f)
 -\mathrm{Tr}(M_{\kappa_\Lambda}A_f Q_\alpha)
=-\mathrm{Tr}(M_{\kappa_\Lambda}A_f Q_\alpha).
```

The first operator has kernel $\kappa_\Lambda(v)f(v-w)$ relative to $dw$. It has trace $`f(0)\int\kappa_\Lambda(v)\,dv=0`$. One can take this trace using its compact smooth kernel on a containing box. Passing to the weighted space conjugates by bounded multiplication operators on the box and leaves this diagonal integral unchanged. Thus the entire width dependence at finite cutoff is in the spectral projection and its ambient space. This uses the operator definition and avoids assigning a diagonal to an arbitrary representative of a Hilbert–Schmidt kernel. [Masullo, definitions 10.4.1.1 and 10.4.1.5](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=455) and [definition 10.4.3.7](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=473).

For the widths justified below, the dual of the quotient is the closed span of the exponential modes for the zeros inside the strip, with their derivatives through multiplicity minus one. On that span the resolvent estimate gives a trace-class convolution. Both projection orders on the ambient line are trace class. To recall the argument, orthogonalise finite translation-invariant spans of the modes. The generator’s resolvent restricted to their closure has squared column norms bounded by $C/(1+\gamma^2)$ for a zero of ordinate $\gamma$. These bounds sum over the zeros. Two integrations by parts factor smooth convolution on the span through the square of this Hilbert–Schmidt resolvent. For the other projection order, split the output into a smooth compact cutoff near the kink of the weight at zero and its complement. The compact piece has a smooth compact kernel up to bounded weight factors. On the complement two adjoint derivatives of the convolution are bounded, so the projected piece factors through the square of the adjoint of the restricted resolvent. Weighted duality gives the corresponding assertion for $A_f Q_\alpha$ on $`X_\alpha`$. The proof and domain details are in the linked note. [Resolvent proof](strip-resolvent.md).

As $\kappa_\Lambda\to1$ strongly and $0\le\kappa_\Lambda\le1$, the error $(M_{\kappa_\Lambda}-1)A_f Q_\alpha$ tends to zero in trace norm for each such fixed $\alpha$, by finite-rank approximation. The diagonal entries of the quotient convolution in the orthogonalised flags are $M_f(\rho)$, once per power retained. Consequently, with $G_\alpha=\lim_\Lambda G_{\alpha,\Lambda}$,

$$
G_\alpha=-\sum_{|\Re\rho-1/2|<\alpha}m_\rho M_f(\rho)
=\frac{\log 2}{\sqrt2}+\underbrace{\sum_{|\Re\rho-1/2|\ge\alpha}m_\rho M_f(\rho)}_{R_\alpha}.
$$

This sign agrees with (12.5.1.4). At a width such as $\alpha=1$ every nontrivial zero is inside. At a specified exceptional width the split of the classical zero sum remains valid for the open strip. Relating it to the full projection at that width requires the corresponding density assertion and a separate count of admissible modes on the wall. Adjoining a generator with one zero divided out to its full multiplicity $m$ gives a projection difference of rank $\min(m,2)$ at the exact width for Sobolev parameter four. [Single-zero calculation](prime-two-kernel-wall.md#the-commutator). [Masullo, p. 520](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=520).

## Strip widths

Here is an argument for the required density on the rational component. Take a half-width satisfying $|\Re(\rho-1/2)\pm\alpha|\ge c_\alpha(1+|\Im\rho|)^{-2}$ at both walls for every zero, with some $c_\alpha>0$. Almost every $\alpha>0$ qualifies. The sum of the lengths of the excluded intervals of radii $(1+|\Im\rho|)^{-2}$ is finite by the ordinary zero count. Borel–Cantelli leaves only finitely many exceptions for almost every width, and these are absorbed in $c_\alpha$ after excluding the countably many exact walls. Every $\alpha>1/2$ qualifies by a uniform positive distance. The generic-width idea also occurs in yikerman’s earlier proposal. [Earlier proposal](https://github.com/mccorvie/riemann-masullo/issues/1169#issuecomment-7053478224) and [zero count](https://arxiv.org/abs/2412.15470).

Let $`D=x\,d/dx`$ and $\Phi(x)=D(D+1)e^{-\pi x^2}$. Its value at zero and its integral over the real line both vanish (for the integral, integrate by parts). The even test therefore meets both adelic moment conditions. Its radial global theta function, divided by two to sum over positive integers, is $g(u)=e^{u/2}\sum_{n\ge1}\Phi(ne^u)$. The Gaussian Poisson formula gives rapid decay of $`g`$ at both ends, faster than every exponential in $|u|$. Its bilateral Laplace transform is $E(z)=\xi(1/2+z)$, with $\xi(s)=s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)/2$. Conversely, the transform of every radial global theta input vanishes at each zero inside the strip to its multiplicity. To see this in the critical strip, average the finite factors over local units. Each modified Euler factor is the usual one times a Laurent polynomial in $p^{-s}$, since the sequence of values on the valuation classes is eventually constant towards zero and eventually zero in the other direction. The archimedean Mellin integral is holomorphic for $\Re s>0$. Unfolding for $\Re s>1$ and continuation give the claimed vanishing. Interior evaluation of each derivative is continuous on $`X_\alpha`$ by Cauchy–Schwarz, so the whole geometric closure has these zeros.

For the reverse inclusion take $`F(z)=\int h(u)e^{zu}\,du`$ for $h\in X_\alpha$ whose transform has the prescribed zeros, and write $d_\epsilon(z)=e^{-\epsilon z^4}$. By Plancherel, $`\|h\|_{X_\alpha}^2`$ is comparable to the sum of the squared $L^2$ norms of $F$ and its first two derivatives on the two walls. Both $F$ and its boundary values are bounded, by Cauchy–Schwarz with $(1+u^2)^{-2}$. The quotient $J_\epsilon=d_\epsilon F/E$ is holomorphic in the strip and continuous up to the walls on finite rectangles.

On the two walls one has $1/|E(z)|\le\exp(C(1+|\Im z|)^2)$. The same estimate holds on horizontal segments across the strip along a sequence of positive and negative heights tending to infinity. For the horizontal segments choose heights in $[n,n+1]$ avoiding every nearby zero ordinate by at least $c(1+n)^{-2}$, possible by $N(T)=O(T\log T)$. For both estimates use the genus-one Hadamard product for $\xi$. At radius $r\ge2$ the zeros with $|\rho|\le2r$ cost $O(r\log^2(2r))$ in the lower bound for the logarithm of the product, using the polynomial separation and $\sum_{|\rho|\le2r}|\rho|^{-1}=O(\log^2(2r))$. The tail costs $O(r^2\sum_{|\rho|>2r}|\rho|^{-2})$. Finite order and Cauchy’s estimates also bound the first two derivatives of $E$ on the walls by $\exp(C(1+|\Im z|)^2)$. [Hadamard factorisation and zero count](https://terrytao.wordpress.com/2021/02/12/246b-notes-4-the-riemann-zeta-function-and-the-prime-number-theorem/).

On those rectangles the quartic damping dominates these bounds. Applying the maximum principle to $(1-z^2)^m J_\epsilon$ for each integer $m\ge0$ gives decay faster than any inverse polynomial in height uniformly across the strip. On both walls $J_\epsilon$ and its first two weak derivatives are in $L^2$, using the product rule and the bounds for $1/E$. Invert its Fourier transform on the central line. Shifting to each wall, then applying Plancherel there to the boundary values and their derivatives, puts the inverse $`q_\epsilon`$ in $`X_\alpha`$. It follows that $d_\epsilon F$ is the transform of $`g*q_\epsilon`$. Convolution by $`g`$ is bounded on $`X_\alpha`$ by its two-sided rapid decay and the translation bound for the weight. Approximating $`q_\epsilon`$ by compact smooth functions puts $`g*q_\epsilon`$ in the closed span of translates of the genuine global generator $`g`$. Finally $d_\epsilon F\to F$ in the two boundary Sobolev norms by dominated convergence. The multipliers and their first two derivatives are uniformly bounded for $0<\epsilon\le1$, as is seen by scaling $z$ by $\epsilon^{1/4}$. This proves the reverse inclusion. Taking the annihilator now gives exactly the closed span of the zero evaluations and their derivatives in the dual.

## An explicit bound

For $z=\rho-1/2$ each normalised box in $B$ supplies the factor $\sinh(hz)/(hz)$. For $|\Re z|\le1/2$ use $|\sinh(hz)|\le e^{h|\Re z|}$ and $|z|\ge|\Im\rho|$. The indicator has mass $2A$, the mollifier has mass one, both translates stay inside the test interval, and $c<1$. Together these give

```math
|M_f(\rho)|\le (2A)\,2\sqrt{5/2}\,(512/|\Im\rho|)^8
<\tfrac14(512/|\Im\rho|)^8.
```

Platt and Trudgian proved that all zeros with $0<|\Im\rho|\le H=3\times10^{12}$ lie on the critical line. Thus every omitted zero for any positive width has height above $H$. Bellotti and Wong’s explicit count gives $N(T)\le T\log T$ for $T\ge H$, with multiplicity and positive ordinates only. Sum the preceding bound over $2^jH<|\Im\rho|\le2^{j+1}H$, counting both signs. With $q=2^{-7}$ the result is

$$
|R_\alpha|\le\frac{512^8}{H^7}
\left(\frac{\log H}{1-q}+\frac{\log2}{(1-q)^2}\right)<6.5\times10^{-65}.
$$

The same majorant applies to the difference of the zero sums at two positive widths, since that difference uses a subset of the same high zeros. This is an absolute bound for this particular smooth test. It leaves room for exact cancellation or for a nonzero sum. [Platt and Trudgian](https://arxiv.org/abs/2004.09765) and [Bellotti and Wong](https://arxiv.org/abs/2412.15470).
