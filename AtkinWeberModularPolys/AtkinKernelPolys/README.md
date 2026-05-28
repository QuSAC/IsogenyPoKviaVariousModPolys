# Kernel polynomial formulas for the Atkin modular polynomials
This directory contains the Atkin kernel polynomial formulas
```math
k_{\ell}(Y,j) \in h_{\ell}^{-1} \cdot ℤ[Y,j][x],
```
where $h_{\ell} = h_{\ell}(Y) \in ℤ[Y]$ denotes the denominator polynomial, in the files `kernel_poly_ℓ.sage`. Given a root of the Atkin modular polynomial $\Phi_\ell^A(Y,j_0)$, they can be used to compute the (kernel polynomial of) the isogeny corresponding to the given root.
Note that, currently, these formulas are only conjectured and experimentally tested for thousands of iterations, though a proof of their correctness might be derivable from the interpolation technique explained below.

## How to use the formulas
**Input**: A $j$-invariant $j_0 \in \overline{K}$ and a root $f \in \overline{K}$ of the $\ell^{th}$ Atkin modular polynomial $\Phi_\ell^A(Y,j)$ at $j_0$, i.e. $\Phi_\ell^A(f, j_0) = 0$. For the formula to work, one has to ensure the following:

1. $j_0 \neq 0, 1728$
2. $\text{char}(K) \neq 2,3, \ell$
3. $f$ does not correspond to two non-equivalent dual loops, i.e. $\Phi_\ell^A(Y, j_0) \in \overline{K}[Y]$ does not have a double root at $f$

**Output**: The kernel polynomial 
```math
k_{\ell}(f, j_0) \in \overline{K}[x]
```
giving the $\ell$-isogeny $\phi \colon E_0 \to E_1$ corresponding to the root $f$, defined on the standard $j_0$-model
```math
E_0: y^2 = x^3 -3j_0(j_0-1728) \cdot x - 2 j_0 (j_0 - 1728)^2.
```
In particular, one has
```math
\Phi_\ell^A(f, j) = (j - j_0) \cdot (j - j(E_1)),
```
and this (generically) determines the kernel polynomial uniquely.

## Comments on the testing script
We also provide a testing script `test_kernel_polynomials.sage`, which defines the method `test_Atkin_kernel_polynomial` to test our kernel polynomial formulas both over a finite field and over the algebraic numbers $$\overline{ℚ}$$; for this testing script, please note the following:
1. For efficiency the tests restrict to $j_0$ in finite field extensions of degree $\leq 8$, respectively to $j_0 \in \overline{ℚ}$ with minimal polynomial of degree $\leq 8$.
2. For slightly larger $\ell$ or $e$, the test over $\overline{ℚ}$ is very slow; we therefore recommend to perform tests over finite fields.
3. For small characteristics $0 < p < 4\ell^2$ the test might run into an infinite loop since it cannot find a simple root of the Atkin kernel polynomial; for example, for $p = 13$, $j_0 = 5 \in 𝔽_{13}$ and $\ell = 3$ we have $$\Phi_\ell^A(Y, j_0) = (Y-5)^2 \cdot (Y-9)^2$$.

## The interpolation technique used to derive the formulas

Unfortunately, it is not feasible to directly compute the kernel polynomial as a factor of the $\ell^{th}$ division polynomial on the (generic) standard model for $j_0$ defined over the function field $ℚ(Y)[j]/(\Phi_\ell^A(Y,j)) =: ℚ(Y)_\Phi$.

To circumvent this, the idea is as follows: The existence of a generic kernel polynomial formula 
```math
k_{\ell}(Y,j) \in ℚ(Y)_\Phi[x]
```
implies that, for any $y \in ℤ$, the polynomial 
```math
k_\ell(y,j) \in ℚ_{\Phi, y}[x]
```
over the number field 
```math
ℚ_{\Phi, y} :=  ℚ[j]/(\Phi_\ell^A(y,j))
```
gives a kernel polynomial formula for the two roots of $\Phi_\ell^A(y,j) \in ℚ[j]$, assuming that we can evaluate $k_\ell(Y,j)$ at $Y=y$ and that $\Phi_\ell^A(y,j) \in ℚ[j]$ is still irreducible.

However, we also obtain a kernel polynomial formula by factoring the $\ell$-th division polynomial of the standard model for $j_0 = \overline{j} := (j \bmod \Phi_\ell^A(y,j)) \in ℚ_{\Phi,y}$, and these two formulas should coincide since the kernel polynomial corresponding to $y$ is unique as we assume $\Phi_\ell^A(y,j)$ to be irreducible.

To avoid integers $y$ at which we cannot evaluate $k_\ell(Y,j)$, and to reduce to a polynomial interpolation for our coefficients, we introduce an "overestimated" denominator polynomial:
Given that the kernel polynomial should not be defined only if we are at a special $j$-invariant or if $y$ encodes two non-equivalent dual loops, a suitable such denominator is given by the resultant
```math
r_\ell(Y) = \text{res}_j \left(\Phi_\ell^A(Y,j), \tfrac{\partial}{\partial Y}\Phi_\ell^A(Y,j)\right) \in ℤ[Y],
```
which captures all $y \in \overline{ℚ}$ that are double roots of $\Phi_\ell^A(Y,j_0)$ for some $j_0 \in \overline{ℚ}$.

In summary, we followed the following interpolation approach for the kernel polynomial formula $k_\ell(Y,j)$:
1. Sample (enough) small integers $y \in ℤ$ such that $\Phi_\ell^A(y,j) \in ℚ[j]$ is irreducible (which also implies $r_\ell(y) \neq 0$).
2. Factor the $\ell^{th}$ division polynomial on $E_0 \colon y^2 = x^3 - 3 \overline{j} (\overline{j} - 1728) - 2 \overline{j} (\overline{j} - 1728)^2$ over the number field $ℚ_{\Phi, y}$, and let $k_\ell(y,j) \in ℚ_{\Phi, y}[x]$ denote the unique factor of degree $\lceil (\ell+1)/2 \rceil$ that we obtain.
3. Represent the coefficients of $k_\ell(y,j)$ uniquely with polynomials of degree $\leq 1$ in $j$.
4. Multiply these representatives by $r_\ell(y)$ and save the result in a list for later interpolation.
5. Use the saved interpolation data to compute, coefficientwise, the polynomial $r_\ell(Y) \cdot k_\ell(Y,j)$, and divide each coefficient by $r_\ell(Y)$ to obtain the desired kernel polynomial formula.

Writing $n_\ell(Y,j) = h_\ell(Y) \cdot k_\ell(Y,j) \in ℤ[Y,j][x]$ for the final numerator polynomial, the number of required interpolation points is
```math
\deg_Y n_\ell(Y,j)(x=0) + \deg_Y r_\ell(Y) - \deg_Y h_\ell(Y) + 1,
```
and we give these numbers explicitly in the table below; note that we did not know these numbers in advance, but simply stopped computing further interpolation points once we were certain that all coefficient polynomials stabilized (in their degree).
$\ell$ | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 41 | 47 | 59 | 71
:--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | 
Points | 7 | 10 | 21 | 36 | 78 | 105 | 171 | 210 | 300 | 465 | 528 | 903 | 1176 | 1830 | 2628
