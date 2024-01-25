from sage.all import *

#########################################
######### Montgomery Curves #############
#########################################

def montgomery_coefficient(E):
    assert E.a1() == 0 and E.a3() == 0
    if E.a6() == 0 and E.a4() == 1:
        return E.a2()
    else:
        E_new = E.montgomery_model()
        print("transformed into Montgomery model")
        return E_new.a2()
        
def montgomery_double(P,E):
    [x1,y1,z1] = P
    A = montgomery_coefficient(E)
    assert z1 == 1, "we expect points of the form [x,y,1]"
    m = (3*x1**2 + A*(x1+x1) + 1)/(y1 + y1)
    x3 = m**2 - (A + x1 + x1)
    y3 = m * (x1 - x3) - y1
    return E([x3,y3,1])


#############################################
### Elliptic Curve: Torsion and Isogenies ###
#############################################

def ell_torsion_basis(ell, maxorder, E):
    """
    Input: integer ell, elliptic curve E of order maxorder (divisible by ell)
	Output: a basis (P_ell,Q_ell) for E[ell]
    """
    while True:
        P = E.random_element()
        P_ell = ZZ(maxorder/ell)*P
        if order(P_ell) == ell:
            while True:
                Q = E.random_element()
                Q_ell = ZZ(maxorder/ell)*Q
                if P_ell.weil_pairing(Q_ell,ell).multiplicative_order() == ell:
                    return [P_ell, Q_ell]


def compute_phi_end_i(E):
    """
    input: E: y^2 = x^3 + 6*x^2 + x over some finite field 
    output: phi: E -> E0: y^2 = x^3 + x
            end_i: non-trivial automorphism of E0
    """

    phi = E.isogeny(E([0,0]), model = "montgomery")
    E0 = phi.codomain()
    end_i_on_E0 = E0.automorphisms()[2]
    P = E0.random_point()
    assert end_i_on_E0(end_i_on_E0(P)) == -P, "wrong automorphism"
    return [phi, end_i_on_E0]


def eval_gamma(u,v,P,E):
    """
    input: integers u,v; P a point on E:y^2=x^3+6*x^2+x
    output: gamma(P), where gamma = u + v *iota (with iota=2i in End(E))
    """
    
    [phi,end_i_on_E0] = compute_phi_end_i(E)
    newP = phi.dual()(end_i_on_E0(phi(P)))
    return u*P + v*newP