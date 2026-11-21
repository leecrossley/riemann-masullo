# The strip resolvent

## The model

Let $`H=L^2(\mathbb R,w(u)\,du)`$, where $w(u)=e^{-2\alpha|u|}(1+u^2)^{-s/2}$, $\alpha,s\ge0$. Let $`E`$ be the closed span of $u^j e^{\lambda u}$ for a discrete set of exponents satisfying $|\Re\lambda|\le\alpha$. At each exponent take the initial segment $0\le j<m_\lambda$ and require these vectors to belong to $H$. In the open strip all powers belong to $H$. On a wall only the powers with $2j<s-1$ do. List each exponent $m_\lambda$ times as $\lambda_n$ and suppose

$$
\sum_n\frac{1}{1+(\Im\lambda_n)^2}<\infty.
$$

For zeta take $\lambda=\rho-1/2$ and use the true multiplicity in the interior. The classical zero count $N(T)=O(T\log T)$ implies the displayed bound by summing over dyadic height intervals. The same proof applies to any fixed character with this counting bound. Identifying $`E`$ with the full dual of the adelic quotient requires the density assertion in Theorem 10.2.2.1. If a zero lies on a wall, the trace counts the modes actually in $`E`$. Widths with no zeros on their walls avoid this distinction. [Zero count](https://arxiv.org/abs/2107.06506) and [Masullo, p. 446](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=446).

## The estimate

Let $U_t h(u)=h(u+t)$. Since $|(\log(1+u^2))'|\le1$, the function $`\log w`$ is Lipschitz with constant $2\omega=2\alpha+s/2$. A change of variable gives $`\|U_t\|\le e^{\omega|t|}`$. Translation is strongly continuous by density of smooth compactly supported functions. It preserves each span obtained by ordering the powers $j$ from lowest to highest at each exponent.

Write $A$ for the translation generator on $H$, $`R=(b-A)^{-1}=\int_0^\infty e^{-bt}U_t\,dt`$ for $b>\omega$, and take $b\ge\alpha+1$. The restriction $R_E$ is the resolvent on $`E`$. Gram-Schmidt orthogonalisation along the finite invariant spans gives an orthonormal basis $e_n$ of $`E`$. The successive one-dimensional quotients have translation eigenvalues $e^{t\lambda_n}$, hence the diagonal of $R_E$ in this basis is $(b-\lambda_n)^{-1}$. This step also applies to repeated zeros because translation preserves the ordering of the powers of $u$. [Engel and Nagel, Theorem II.1.10](https://www.math.uni-tuebingen.de/de/forschung/agfa/members/engel-nagel_one-parameter-semigroups.pdf/@@download/file/engel-nagel_one-parameter-semigroups.pdf).

For $y=R_E e_n$, differentiation of $`\|U_t y\|^2\le e^{2\omega t}\|y\|^2`$ at $t=0$ from the right gives $`\Re\langle y,Ay\rangle\le\omega\|y\|^2`$. Since $e_n=by-Ay$,

```math
(b-\omega)\|R_Ee_n\|^2
\le \Re\langle R_Ee_n,e_n\rangle
=\frac{b-\Re\lambda_n}{(b-\Re\lambda_n)^2+(\Im\lambda_n)^2}
\le\frac{b+\alpha}{1+(\Im\lambda_n)^2}.
```

Summation proves that $R_E$ is Hilbert-Schmidt. For $k\in C_c^\infty(\mathbb R)$ put $`C_k=\int k(t)U_t\,dt`$ and $`\widehat k(z)=\int k(t)e^{tz}\,dt`$. Twice integrating by parts gives $`C_k|_E=R_E^2 C_{(b+\partial_t)^2k}|_E`$. The last factor is bounded and $R_E^2$ is trace class. On the same finite spans the diagonal of $C_k$ is $\widehat k(\lambda_n)$, so $`\mathrm{Tr}_E C_k=\sum_n\widehat k(\lambda_n)`$ with absolute convergence. [Bell, trace ideals](https://jordanbell.info/LaTeX/mathematics/traceclass/).

## Projection order

Let $P$ be the orthogonal projection of $H$ onto $`E`$. Invariance gives $C_kP=\iota_E(C_k|_E)P$, which is trace class. For the reverse order, take $\chi\in C_c^\infty(\mathbb R)$ equal to one near zero. The operator $\chi C_k$ is trace class on $H$. On unweighted $L^2$ its kernel $\chi(u)k(v-u)$ is smooth with compact support in both variables, so its Fourier series on a containing box gives an absolutely summable rank-one expansion. Conjugation to the weighted space introduces bounded multipliers $w^{1/2}(u)$ and $w^{-1/2}(v)$ on that box.

Put $J=(1-\chi)C_k$. On its support $`A^*=-\partial_u-(\log w)'`$, and the first two derivatives of $`\log w`$ are bounded there. The identities $`\partial_u C_k=-C_{k'}`$ and $`\partial_u^2 C_k=C_{k''}`$ show that $`J(H)\subset D((A^*)^2)`$ and $`D=(b-A^*)^2J`$ is bounded. These identities follow first on smooth compactly supported inputs and extend in the graph norms by closure. Hence $`J=(R^*)^2D`$. Invariance gives $`PR^*=R_E^*P`$ as maps into $`E`$. It follows that $`PJ=(R_E^*)^2PD`$ is trace class, as is $PC_k$.

Extending the basis of $`E`$ to one of $H$ shows that both projection orders have trace $`\mathrm{Tr}_E C_k`$. Uniformly bounded multiplication cutoffs converging pointwise to one converge strongly to the identity. Left multiplication by them therefore converges in trace norm on both operators, by finite-rank approximation. This calculation on the weighted line uses the actual orthogonal projection and includes both orders relevant to §10.4. Masullo also asserts convergence under arbitrary enumeration of the individual biorthogonal rank-one operators. The proof here uses orthogonal finite spans. [Masullo, pp. 456-463](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=456).

## The exterior zeros

For zeta, write $A(f)$ for the arithmetic side of the classical explicit formula and $`M_f(\rho)=\int f(u)e^{(\rho-1/2)u}\,du`$ for a smooth compact test supported in $u>0$. Splitting the absolutely convergent full sum gives

$$
\sum_{|\Re\rho-1/2|<\alpha}m_\rho M_f(\rho)
=A(f)-\sum_{|\Re\rho-1/2|\ge\alpha}m_\rho M_f(\rho).
$$

This corrected formula is unconditional. For $\alpha\ge1/2$ the exterior is empty by the classical zero-free boundary of the critical strip. At a fixed narrower width the exterior term vanishes for every such test exactly when the strip contains every non-trivial zero. Smooth translates of a test whose transform is nonzero at a chosen zero isolate its multiplicity by a Laplace transform and residue. The Lean theorems `fixed_strip_iff` and `exterior_vanishes_iff` give these statements with the classical full explicit formula as an input. The resolvent argument works for any subcollection of the zero modes. An identification with the arithmetic expression must use further information about the actual adelic subspace. Section 12.3 of the posted draft contains only a heading and a footnote requesting the section. [Tao, exercise 22(iii)](https://terrytao.wordpress.com/2021/02/12/246b-notes-4-the-riemann-zeta-function-and-the-prime-number-theorem/#smoothed), [Masullo, p. 514](https://math216.wordpress.com/wp-content/uploads/2026/11/the-music-of-the-primes.pdf#page=514) and [issue 587](https://github.com/mccorvie/riemann-masullo/issues/587).
