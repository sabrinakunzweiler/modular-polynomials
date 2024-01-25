# Computing modular polynomials by deformation

This repository contains the code accompanying our paper

> S. Kunzweiler and D. Robert, Computing modular polynomials by deformation

In particular, we provide an implementation for computing the modular polynomial $\phi_\ell$
for any prime $\ell \equiv 3 \pmod{4}$.

## Examples

For instance, let $\ell = 11$, then the modular polynomial can be computed as follows.

```
sage: load('main.sage')
sage: phi11 = modular_polynomial(11)
```
The modular polynomial is represented as a list of coefficients. The monomial coefficient of $X^i \cdot Y^j$ 
is `phi11[i][j]` (where we assume $i \geq j$). 

For primes $p$ satisfying $2^n\cdot \ell \mid p+1$, where  $2^n - \ell = a^2 + 4b^2$, one can also compute 
the modular polynomial modulo $p$ directly. For example with the code below,
one computes the modular polynomial $\phi_{103} \in \mathbb{F}_{79103}[x,y]$.

```
sage: ell = 103
sage: [n,u,v] = [7,5,0]
sage: p = 79103
sage: phi103 = modular_polynomial_modp(p,ell,n,u,v)
```
