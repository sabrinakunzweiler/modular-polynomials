"""
Computing modular polynomials by deformation
==========================================================================

This is the code accompanying our paper

.. S. Kunzweiler and D. Robert, Computing modular polynomials by deformation

In particular, we provide an implementation for computing the modular polynomial phi_ell
for any prime ell = 3 (mod 4). 

 """

from elliptic_utilities import *
from kummer_isogenies.product_isogeny_kummer import *


load("lifting-isogenies.sage")
load("formulae_kummer_lift.sage")

# set number of threads to allow parallel computation
# otherwise its 1
os.environ['SAGE_NUM_THREADS'] = '10'

half = 1

def prime_bound(ell):
        """
        bounds for the logarithmic height of the coefficients of the modular polynomial
        """
        bound1 = 18*ell + 6*ell*log(ell)
        bound2 = 16*ell + 6*ell*log(ell) + 14*sqrt(ell)*log(ell)
        return min(bound1,bound2)


def suitable_n(ell):
        """
        find u,v,n such that: u^2+4*v^2 = 2^n-ell
        """
        assert ell % 4 == 3, "method is only implemented for ell = 3 (mod 4)"
        n = ceil(log(ell,2))
        found_n =False
        factors = []
        
        while not found_n:
            factors = factor(2^n-ell)
            if all([el[0]^el[1] % 4 != 3 for el in factors]):
                found_n = True
            else:
                n = n+1
        [u,v] = two_squares(2^n-ell)
        assert u^2 + v^2 == 2^n - ell, [n,u,v]
        #we want u^2 + 4*v^2 = 2^n - ell
        if v % 2 == 0:
            return[n, u, ZZ(v/2)]
        else:
            return [n, v, ZZ(u/2)]
 

def next_suitable_prime(ell, n, cofactor):
        while True:
            cofactor = cofactor + 1
            p = cofactor*ell*2^(n+1)-1
            if p.is_prime():
                return [p,cofactor]


def mod_poly_from_invariants(j_invariants, j00):
    '''
    input: list of all j_invariants ell-isogenous to j00 + eps
    output: modular polynomial 
    '''
    S.<Y> = j_invariants[0].parent()[]
    mod_poly_eps = prod([Y-jk for jk in j_invariants])
    #K = j00.base_ring() 
    K = j00.parent() #strange bug: for large primes, the computation over F_p doesn't work, but we need F_{p^2}
    R.<X,Y> = K[]
    mod_poly = 0
    for m in mod_poly_eps.monomials():
        d = m.degree()
        c = mod_poly_eps.monomial_coefficient(m)
        newc = sum([c.coefficients()[i]*(X-j00)^(c.exponents()[i]) for i in range(len(c.exponents()))])
        mod_poly = mod_poly + newc*Y**d

    return mod_poly


@parallel
def modular_polynomial_modp(p,ell,n,u,v):
    '''
    Input: a prime p of the form p = 2^(n+1)*ell*cofactor - 1, and integers u,v with 2^n-ell = u^2 + v^2
    Output: modular polynomial for degree ell over the finite field F_p
    '''
    Fp = FiniteField(p**2)
    global half
    half = Fp(1/2)
    E = EllipticCurve([0,Fp(6),0,1,0])
    [P_ell,Q_ell] = ell_torsion_basis(ell, p+1, E)
    [P_2n, Q_2n] = ell_torsion_basis(2**(n+1), p+1, E)
    j_invariants = [0]*(ell+1)
    gamma_P_2n = eval_gamma(u,v,P_2n,E)
    gamma_Q_2n = eval_gamma(u,v,Q_2n,E)
    for k in range(ell+1):
        if k == ell:
            Pk = Q_ell
        else:
            Pk = P_ell + k*Q_ell
        phik = E.isogeny(Pk, model = "montgomery")
        Ek = phik.codomain()
        Ak = montgomery_coefficient(Ek)
        ker = [gamma_P_2n, gamma_Q_2n, phik(P_2n),phik(Q_2n)]
        jk = lift_isogeny(E,Ek,ker, ell,n)
        j_invariants[k] = jk
    j00 = Fp(E.j_invariant())
    mod_poly = mod_poly_from_invariants(j_invariants, j00)
    return mod_poly


def modular_polynomial(ell,ncpu=8,verbose=False):
        """
        Input: integer ell = 3 (mod 4)
        Output: coefficients of the modular polynomial Phi_ell(X,Y) 
        where for j<=i, coefficients[i][j] is the monomial coefficient of X^iY^j.
        """
        [n,u,v] = suitable_n(ell)
        if verbose:
            print("parameters: ", [n,u,v])
        bound = prime_bound(ell)
        #print("bound on logarithmic height of the primes:", bound.n())
        cofactor = 0
        log_primes = 0
        suitable_primes = []
        M = 1
        coefficients = [[0 for j in range(i+1)] for i in range(ell+2)]

        while log_primes < log(2) + bound:
            [p, cofactor] = next_suitable_prime(ell,n,cofactor)
            #print("prime:", p)
            suitable_primes.append(p)
            log_primes = log_primes + log(p)

        while suitable_primes:
            next_primes = suitable_primes[:ncpu]
            suitable_primes = suitable_primes[ncpu:]
            mod_polys = list(modular_polynomial_modp([(p,ell,n,u,v) for p in next_primes]))
            M_small = product(next_primes)
            M_new = M*M_small
            for i in range(ell+2):
                for j in range(i+1):
                    #CRT of asubset of primes
                    cij = CRT([ZZ(mod_poly[1][i,j]) for mod_poly in mod_polys],[el[0][0][0] for el in mod_polys])
                    #update coefficients from previous computation
                    cij = CRT([cij, coefficients[i][j]],[M_small,M])   
                    if cij > M_new/2:
                        cij = cij -M_new
                    elif cij < -M_new/2:
                        cij = cij +M_new
                    coefficients[i][j] = cij
            M = M_new
            if verbose:
                print("progress: ", next_primes, "done")
        return coefficients